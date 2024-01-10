import numpy as np
import matplotlib.pyplot as plt
import symlib
import scipy.stats as stats
import scipy.interpolate as interp

import palette
from palette import pc

CALC_RS = True

def err_bounds(x, y, bins):
    n_min = 20
    bs = stats.binned_statistic
    n, edges, _ = bs(x, y, "count", bins)
    lo, _, _ = bs(x, y, lambda z: np.quantile(z, 0.5-0.68/2), bins)
    md, _, _ = bs(x, y, lambda z: np.quantile(z, 0.5), bins)
    hi, _, _ = bs(x, y, lambda z: np.quantile(z, 0.5+0.68/2), bins)
    
    mid = 10**((np.log10(edges[1:]) + np.log10(edges[:-1]))/2)
    ok = n > n_min

    return mid[ok], lo[ok], md[ok], hi[ok]

def m23_vmax_ratio(m_ratio):
    """ m23_vmax_ratio converts an msub/mpeak ratio to the equivalent
    vmax/vpeak ratio according to the fit in Mansfield et al. 2023.
    """
    eta, mu = 0.184, 0.274
    return 2**mu * m_ratio**eta / (1 + m_ratio)**mu

def j_vdb_16_vmax_ratio(m_ratio):
    """ j_vdb_16_vmax_ratio converts an msub/mpeak ratio to the equivalent
    vmax/vpeak ratio according to the fit in Jiang & van den Bosch 2016.
    """
    eta, mu = 0.35, 0.45
    return 2**mu * m_ratio**eta / (1 + m_ratio)**mu


def g_vdb_19_vmax_ratio(m_ratio, c):
    p0 = 2.980
    p1 = 0.310
    p2 = -0.223
    p3 = -3.308
    p4 = -0.079
    q0 = 0.176
    q1 = -0.008
    q2 = 0.452
    mu = p0 + p1*c**p2*np.log10(m_ratio) + p3*c**p4 - 0.2
    eta = q0 + q1*c**q2*np.log10(m_ratio) - 0.05
    return 2**mu * m_ratio**eta / (1 + m_ratio)**mu

def f_nfw(x):
    return np.log(1+x) - x/(1+x)

def is_conv_vdb(param, scale, f_bound, r_half, rs0, c0):
    eps = param["eps"]/param["h100"] * scale
    """
    print(f_bound > 1.79/f_nfw(c0)*(eps/rs0)*(r_half/rs0))
    print("f_bound", f_bound)
    print(1.79/f_nfw(c0))
    print("eps", eps)
    print("r_half", r_half)
    print("rs0", rs0)
    print("c0", c0)
    print(f_bound / (1.79/f_nfw(c0)*(eps/rs0)*(r_half/rs0)))
    exit(0)
    """
    return f_bound > 1.79/f_nfw(c0)*(eps/rs0)*(r_half/rs0)


def numerical_limit(n_sub):
    A_lim, b_lim = 67.8, 0.257
    n_peak_lim = (n_sub/A_lim)**(1/b_lim)
    return n_sub/n_peak_lim

def main():
    palette.configure(True)

    suites = ["SymphonyMilkyWay", "SymphonyMilkyWayHR", "SymphonyLMC", 
              "SymphonyLCluster"]
    base_dir = "/sdf/home/p/phil1/ZoomIns"

    v_rs, v_c = [], []
    m_rs, m_c = [], []
    n_rs, n_c = [], []
    c_infall_c, c_infall_rs = [], []
    is_conv_c = []

    for suite in suites:
        n_hosts = symlib.n_hosts(suite)
        for i_host in range(n_hosts):
            print(suite, i_host)
            sim_dir = symlib.get_host_directory(base_dir, suite, i_host)
            h, hist = symlib.read_subhalos(sim_dir)
            c = symlib.read_cores(sim_dir)

            param = symlib.simulation_parameters(sim_dir)
            mp = param["mp"]/param["h100"]

            I = np.arange(len(h), dtype=int)
            snap = np.arange(h.shape[1], dtype=int)
            scale = symlib.scale_factors(sim_dir)
            a_infall = scale[hist["first_infall_snap"]]
            vmax_infall = h["vmax"][I,hist["first_infall_snap"]]
            m_infall = h["mvir"][I,hist["first_infall_snap"]]
            mhost_infall = h["mvir"][0,hist["first_infall_snap"]]
            c_infall = h["cvir"][I,hist["first_infall_snap"]]
            rs_infall = h["rvir"][I,hist["first_infall_snap"]]/c_infall

            for i in range(len(h)):
                if m_infall[i]/mhost_infall[i] > 0.1: continue

                ok_c = c["ok"][i]
                ok_rs = c["ok_rs"][i] & c["ok"][i]
            
                m_c.append(c["m_bound"][i,ok_c]/m_infall[i])
                m_rs.append(h["mvir"][i,ok_rs]/m_infall[i])

                v_c.append(c["vmax"][i,ok_c]/vmax_infall[i])
                v_rs.append(h["vmax"][i,ok_rs]/vmax_infall[i])
            
                c_infall_c.append(c_infall[i]*np.ones(np.sum(ok_c)))
                #c_infall_c.append(10*np.ones(np.sum(ok_c)))
                c_infall_rs.append(c_infall[i]*np.ones(np.sum(ok_rs)))
                n_c.append(c["m_bound"][i,ok_c]/mp)
                n_rs.append(h["mvir"][i,ok_rs]/mp)

                f_bound = c["m_bound"][i,ok_c]/m_infall[i]
                is_conv_1 = is_conv_vdb(
                    param, scale[ok_c], f_bound,
                    c["r50_bound"][i,ok_c], c_infall[i], rs_infall[i])
                is_conv_2 = f_bound > (0.32 * (m_infall[i]/mp/1e3)**-0.8)
                is_conv_c.append(is_conv_1 & is_conv_2)

    m_c = np.hstack(m_c)
    m_rs = np.hstack(m_rs)
    n_c = np.hstack(n_c)
    n_rs = np.hstack(n_rs)
    v_c = np.hstack(v_c)
    v_rs = np.hstack(v_rs)
    c_infall_c = np.hstack(c_infall_c)
    c_infall_rs = np.hstack(c_infall_rs)
    is_conv_c = np.hstack(is_conv_c)

    if CALC_RS:
        m_c, v_c, n_c, c_infall_c = m_rs, v_rs, n_rs, c_infall_rs

    v_g19_c = g_vdb_19_vmax_ratio(m_c, c_infall_c)
    v_g19_rs = g_vdb_19_vmax_ratio(m_rs, c_infall_rs)

    cut_exp = [1, 2, 3, 4][::-1]
    ls = [":", "-.", "--", "-"][::-1]

    gs_kw = {"height_ratios": [3, 1.5]}
    fig, ax = plt.subplots(nrows=2, sharex=True,
                           gridspec_kw=gs_kw)

    c_med = np.median(c_infall_c[n_c > 1e4])
    bins = 10**np.linspace(-4, 0, 40)
    #ax[0].plot(bins, j_vdb_16_vmax_ratio(bins), c=pc("r"), label=r"${\rm Idealised\ sims\ (P10)}$")
    ax[0].plot(bins, g_vdb_19_vmax_ratio(bins, c_med), c=pc("r"), label=r"${\rm Idealised\ sims\ (GvdB19)}$")
    ax[1].plot([1e-3, 1], [1, 1], c=pc("r"))

    for i in range(len(cut_exp)):
        ok_c = n_c > 10**cut_exp[i]
        ok_rs = n_rs > 10**cut_exp[i]
        m_g19_c, v_g19_c_lo, v_g19_c_md, v_g19_c_hi = err_bounds(
            m_c[ok_c], v_g19_c[ok_c], bins)
        m_g19_rs, v_g10_rs_lo, v_g19_rs_md, v_g19_rs_hi = err_bounds(
            m_rs[ok_rs], v_g19_rs[ok_rs], bins)
        M_c, v_c_lo, v_c_md, v_c_hi = err_bounds(
            m_c[ok_c], v_c[ok_c], bins)
        M_rs, v_rs_lo, v_rs_md, v_rs_hi = err_bounds(
            m_rs[ok_rs], v_rs[ok_rs], bins)
    
        ok_conv = ok_c & (m_c < 0.05)
        print(np.sum(ok_conv))

        c_med = np.median(c_infall_c[ok_c])
        ax[0].plot(M_c, v_c_md, ls[i], c=pc("b"), lw=1)
        ax[1].plot(M_c, v_c_md/g_vdb_19_vmax_ratio(M_c, c_med), ls[i],
                   c=pc("b"), lw=1)
        ok = M_c > numerical_limit(10**cut_exp[i])
        ax[0].plot(M_c[ok], v_c_md[ok], ls[i], c=pc("b"),
                   label=r"$n > 10^%d$" % cut_exp[i])
        ax[1].plot(M_c[ok], (v_c_md/g_vdb_19_vmax_ratio(M_c, c_med))[ok],
                   ls[i], c=pc("b"))
    
    ax[1].set_xlim(1e-3, 1)
    ax[0].set_ylim(0.05, 1)
    ax[1].set_ylim(0.5, 1.05)
    lo, hi = plt.xlim()
    ax[0].set_yscale("log")
    ax[1].set_xscale("log")
    ax[0].legend(loc="lower right", fontsize=16)
    ax[1].set_xlabel(r"$m/m_{\rm infall}$")
    ax[0].set_ylabel(r"$H=v_{\rm max}/v_{\rm max,infall}$")
    ax[1].set_ylabel(r"$H/H_{\rm GvdB19}$")

    if CALC_RS:
        plt.savefig("../plots/core_tracking/vacc_macc_rs.pdf")
    else:
        plt.savefig("../plots/core_tracking/vacc_macc.pdf")

    """
    plt.figure()

    eval_ratios = [0.3, 0.1, 0.03, 0.01]
    colors = [pc("r"), pc("o"), pc("b"), pc("p")]

    for j in range(len(eval_ratios)):
        lims = 10**np.linspace(1.5, 5, 80)
        ratio = np.zeros(len(lims))
        eval_ratio = eval_ratios[j]
        for i in range(len(lims)):
            bins = 10**np.linspace(-4, 0, 20)
            ok = n_c > lims[i]
            _, v_g19_lo, v_g19_md, v_g19_hi = err_bounds(
                m_c[ok], v_g19_c[ok], bins)
            M, v_lo, v_md, v_hi = err_bounds(
                m_c[ok], v_c[ok], bins)
            
            if len(M) == 0 or min(M) > eval_ratio:
                ratio[i] = np.nan
                continue
                
            f = interp.interp1d(np.log10(M), np.log10(v_md))
            f_g19 = interp.interp1d(np.log10(M), np.log10(v_g19_md))

            ratio[i] = 10**(f(np.log10(eval_ratio))-f_g19(np.log10(eval_ratio)))
        ok = ~np.isnan(ratio)
        plt.plot(lims[ok], ratio[ok], c=colors[j], label=r"$m_{\rm sub}/m_{\rm infall}=%.2f$" % eval_ratio)
    plt.xscale("log")
    lo, hi = plt.xlim()
    plt.xlim(lo, hi)
    plt.plot([lo, hi], [1, 1], "--", c="k", lw=1.5)
    plt.xlabel(r"$n_{\rm sub}$")
    plt.ylabel(r"$H/H_{\rm ideal}(>n_{\rm sub})$")
    plt.legend(loc="upper left")

    plt.savefig("../plots/core_tracking/nsub_vmax.pdf")
    """
    
if __name__ == "__main__": main()
