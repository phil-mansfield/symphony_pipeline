# Standarf library importat
import numpy as np
import scipy.interpolate as interpolate
import matplotlib.pyplot as plt

# Colossus for cosmology
from colossus.cosmology import cosmology
from colossus.halo import mass_so

# symlib: a library for interacting with the Symphony simulation
# See this link to install, don't worry about downloading yet
# http://web.stanford.edu/group/gfc/symphony/build/html/quickstart.html
import symlib

# This is my plotting library. Download it here
# https://github.com/phil-mansfield/palette
# and add it to your PYTHONPATH variable
import palette
from palette import pc

# when subhalos are larger than this fraction of their host's mass, they
# rapidly merge. We don't analyze them here.
max_ratio = 0.1
# This is where our Symphony data is stored
base_dir = "/sdf/home/p/phil1/ZoomIns"
# True if we want to combine the high resolution Symphony simulations and low
# reoslution Symphony simulations
COMBINE_HR = True

# determine whether we make a the "vmax" plto with orphan model comparions or
# the "mass" plot with just subhalo masses.
#ratio_type = "vmax"
ratio_type = "mass"

def plot_arrow_lines(ax, xc, yc, x_arrow, dlogy_arrow, color):
    """ Plots a curve with downward pointing arrows underneath it. ax is the 
    pyplot axis, xc and yc are the x and y values of the curve, x_arrow is
    the x values that the arrows will be drawn at, dlogy_arrow is the
    logarithmic distance ebtween the arrows and the curve, and color is a
    color string.
    """
    for i in range(len(x_arrow)):
        j = np.searchsorted(xc, x_arrow[i])
        tip_start = 10**(np.log10(yc[j]) + dlogy_arrow)
        tip_end = 10**(np.log10(tip_start) + dlogy_arrow*0.2)
        ax.plot([xc[j]], [tip_start], marker=r"$\downarrow$",
                color=color, markersize=20)
        

def kaplan_meier(dt, censored, t_eval, decreasing=False):
    """ kaplan_meier calculates Kaplan-Meier statsitics for a given set of
    survival times, dt. dt are the survival times, censored is a boolean
    array that's true if the observation ended before the data point died, and
    t_eval are the times that you want to evaluate at. By default, survival
    curves are calculated starting at lower t values (i.e. P_death(<t)), but
    you can reverse that by setting decreasing to False.
    
    Two arrays are returned. The first is the survival curve and the second is
    the 1-sigma errors on the survival curve estiamted via Greenwood's formula.
    """
    if decreasing:
        dt, t_eval = -dt, -t_eval
    dt_death = dt[~censored]
    dt = np.sort(dt)
    ti, di = np.unique(dt_death, return_counts=True)
    ni = len(dt) - np.searchsorted(dt, ti, side="left")

    hi = 1 - di/ni
    Si = np.cumprod(hi)

    # Greenwood forumla
    diff = ni-di
    diff[(ni-di) <= 0] = 1
    err = Si * np.sqrt(np.cumsum(di / (ni*diff)))
    
    f_S = interpolate.interp1d(ti, Si, fill_value=np.nan,
                               bounds_error=False)
    f_err = interpolate.interp1d(ti, err, fill_value=np.nan,
                                 bounds_error=False)

    return f_S(t_eval), f_err(t_eval)
    
def calc_mpeak_pre(h, hist):
    """ calc_mpeak_pre computes Mpeak prior to first infall for all subhaloes.
    Requires h, a symlib.SUBHALO_DTYPE array and hist, a symlib.HISTORY_DTYPE
    array.
    """
    mpeak_pre = np.zeros(len(h))
    for i in range(1, len(h)):
        mpeak_pre[i] = np.max(h["mvir"][i,:hist["first_infall_snap"][i]+1])
    mpeak_pre[0] = np.max(h["mvir"][0,:])
    return mpeak_pre

def find_quantile(x, cdf, p):
    """ find_quantile inverts a CDF to find the value of x where cdf reaches
    given quantile, p.
    """
    f = scipy.inter1d(np.log10(cdf/cdf[-1]), np.log10(x))
    return 10**f(np.log10(p))

def j_vdb_16_vmax_ratio(m_ratio):
    """ j_vdb_16_vmax_ratio converts an msub/mpeak ratio to the equivalent
    vmax/vpeak ratio according to the fit in Jiang & van den Bosch 2016.
    """
    eta, mu = 0.3, 0.4
    return 2*mu * m_ratio**eta / (1 + m_ratio)**eta

def survival_time():
    """ The left panel of surival curve plots
    """

    # Symphony simulation suites that we use.
    suite = "SymphonyMilkyWay"
    suite_lr = "SymphonyMilkyWayLR"
    suite_hr = "SymphonyMilkyWayHR"

    
    npeak = [] # array on n_{peak,pre-finall} for all subhalos
    m_c = [] # array of disruption masses with particle-tracking
    m_rs = [] # array of disruption masses with Rockstar
    censored_c = [] # True if subhalo survived to end os sim in pt
    censored_rs = [] # True if subhalo survived to end of sim in rockstar

    # Same as above for the low-res simulations
    npeak_lr = []
    t_lr = []
    m_lr = []
    censored_lr = []

    # Same as above for the high-res simulations
    npeak_hr = []
    t_hr = []
    m_hr = []
    censored_hr = []

    # Build up the arrays one suite at a time.
    for suite_i in [suite, suite_lr, suite_hr]:

        n_hosts = symlib.n_hosts(suite_i)

        for i_host in range(n_hosts):
            print(suite_i, i_host)
            sim_dir = symlib.get_host_directory(base_dir, suite_i, i_host)
        
            param = symlib.simulation_parameters(sim_dir)
            mp = param["mp"]/param["h100"]

            # Read in subhalo data.
            h, hist = symlib.read_subhalos(sim_dir)
            c = symlib.read_cores(sim_dir)
            mpeak_pre = calc_mpeak_pre(h, hist)

            # Calculate the last snapshot where we can see each subhalo.
            last_snap_rs = np.zeros(len(h), dtype=int)
            last_snap_c = np.zeros(len(h), dtype=int)
            snap = np.arange(h.shape[1], dtype=int)
            for i in range(1, len(h)):
                ok_rs, ok_c = c["ok_rs"][i], c["ok"][i]
                if np.sum(ok_c) == 0:
                    last_snap_c[i] = hist["first_infall_snap"][i]
                else:
                    last_snap_c[i] = np.max(snap[ok_c])
                if np.sum(ok_rs) == 0:
                    last_snap_rs[i] = hist["first_infall_snap"][i]
                else:
                    last_snap_rs[i] = np.max(snap[ok_rs])
            last_snap_c[0], last_snap_rs[0] = len(snap)-1, len(snap)-1
            
            idx = np.arange(len(c), dtype=int)
            if ratio_type == "mass":
                # mass ratios
                dm_c = c["m_bound"][idx,last_snap_c]/mpeak_pre
                dm_rs = h["mvir"][idx,last_snap_rs]/mpeak_pre
            elif ratio_type == "vmax":
                # vmax ratios, note that this is relative to vmax @ Mpeak, but
                # also mpeak is calculated prior to first infall
                vmpeak = np.zeros(len(h))
                snap = np.arange(h.shape[1], dtype=int)
                for i in range(1, len(h)):
                    ok = snap <= hist["first_infall_snap"][i]
                    i_peak = np.argmax(h["mvir"][i,ok])
                    vmpeak[i] = h["vmax"][i,i_peak]
                vmpeak[0] = 1
                dm_c = c["vmax"][idx,last_snap_c]/vmpeak
                dm_rs = h["vmax"][idx,last_snap_rs]/vmpeak
            else:
                assert(0)
    
            # subhalo needs to not be too big and needs to have been in its
            # host for at least one snapshot
            ratio_ok = ((hist["merger_ratio"] < max_ratio) &
                        (hist["first_infall_snap"] < last_snap_c) &
                        (hist["first_infall_snap"] < last_snap_rs))
            ratio_ok[0] = False

            # append arrays before hstacking them
            if suite_i == suite:
                censored_rs.append(c["ok_rs"][ratio_ok,-1])
                censored_c.append(c["ok"][ratio_ok,-1])
                m_rs.append(dm_rs[ratio_ok])
                m_c.append(dm_c[ratio_ok])
                npeak.append((mpeak_pre/mp)[ratio_ok])
            elif suite_i == suite_lr:
                censored_lr.append(c["ok"][ratio_ok,-1])
                m_lr.append(dm_c[ratio_ok])
                npeak_lr.append((mpeak_pre/mp)[ratio_ok])
            elif suite_i == suite_hr:
                censored_hr.append(c["ok"][ratio_ok,-1])
                m_hr.append(dm_c[ratio_ok])
                npeak_hr.append((mpeak_pre/mp)[ratio_ok])
                
    # combine appended arrays
    censored_rs = np.hstack(censored_rs)
    censored_c = np.hstack(censored_c)
    m_rs = np.hstack(m_rs)
    m_c = np.hstack(m_c)
    npeak = np.hstack(npeak)

    censored_lr = np.hstack(censored_lr)
    m_lr = np.hstack(m_lr)
    npeak_lr = np.hstack(npeak_lr)

    censored_hr = np.hstack(censored_hr)
    m_hr = np.hstack(m_hr)
    npeak_hr = np.hstack(npeak_hr)

    # masses that the survival curves will be evaluated at
    n_eval = 100
    m_eval = 10**np.linspace(-4, 1, n_eval)

    # what resolution bins do we want to plot?
    lims = [(10**4.5, 10**5)] # (lower limit npeak, upper limit npeak)
    labels = [r"$10^{4.5}<N_{\rm peak}<10^5$"]
    # colors to use for each bin
    colors = [pc("r")]

    fig1, ax_m = plt.subplots()

    # For each limit, calculate survival curves.
    for i in range(len(lims)):
        low, high = lims[i]
        ok = (npeak > low) & (npeak < high)

        S_m_c, err_m_c = kaplan_meier(m_c[ok], censored_c[ok], m_eval,
                                      decreasing=True)

        S_m_rs, err_m_rs = kaplan_meier(m_rs[ok], censored_rs[ok], m_eval,
                                        decreasing=True)
        
        ok_lr = (npeak_lr > low) & (npeak_lr < high)
        S_m_lr, err_m_lr = kaplan_meier(m_lr[ok_lr], censored_lr[ok_lr],
                                        m_eval, decreasing=True)

        ok_hr = (npeak_hr > low) & (npeak_hr < high)
        S_m_hr, err_m_hr = kaplan_meier(m_hr[ok_hr], censored_hr[ok_hr],
                                        m_eval, decreasing=True)


        ok = ~np.isnan(S_m_rs)
        ax_m.plot(m_eval[ok], S_m_rs[ok], "-", c=pc("b"),
                  label=r"${\rm Rockstar}$")
        ax_m.fill_between(m_eval[ok], (S_m_rs+err_m_rs)[ok],
                          (S_m_rs-err_m_rs)[ok], color=pc("b"), alpha=0.2)

        ok = ~np.isnan(S_m_c)
        ax_m.plot(m_eval[ok], S_m_c[ok], "-", c=pc("r"),
                  label=r"${\rm Particle{-}tracking}$")
        ax_m.fill_between(m_eval[ok], (S_m_c+err_m_c)[ok],
                          (S_m_c-err_m_c)[ok], color=pc("r"), alpha=0.2)

    # Everything from here on in this function is just plotting nonsense.

    ax_m.legend(loc="upper left", fontsize=17)
    ax_m.set_ylim(0, 1.25)
    if ratio_type == "mass":
        ax_m.set_xlabel(r"$m_{\rm sub,disrupt}/m_{\rm sub,peak}$")
        ax_m.set_ylabel(r"${\rm Pr}(<m_{\rm sub,disrupt}/m_{\rm sub,peak})$")
        ax_m.set_xscale("log")
        ax_m.text(1.2e-4, 0.87, r"$10^{4.5} < n_{\rm sub,peak} < 10^5$", fontsize=17)
        ax_m.set_xlim((1e-4, 1))
    elif ratio_type == "vmax":
        ax_m.set_xlabel(r"$V_{\rm max,disrupt}/V_{\rm max,Mpeak}$")
        ax_m.set_ylabel(r"${\rm Pr}(V_{\rm max,disrupt}/V_{\rm max,Mpeak})$")
        ax_m.text(0.032, 0.87, r"$10^{4.5} < n_{\rm sub,peak} < 10^5$", fontsize=17)
        ax_m.set_xlim((0.03, 1))
        ax_m.set_xscale("log")

    if ratio_type == "mass":
        fig1.savefig("../plots/core_tracking/survival.pdf")
    elif ratio_type == "vmax":
        fig1.savefig("../plots/core_tracking/survival_vmax.pdf")

def survival_quantile(x, y, y_err, p):
    """ survival_quantile takes in a survival curve, (x, y), the error on y,
    y_err, and a target quantile, p, and computes the x values (lower 1-sigma,
    median, upper 1-sigma) that it occurs at.
    """
    y_high = y + y_err
    y_low = y - y_err

    for i in range(1, len(x)):
        if np.isnan(y[i]) or np.isnan(y[i-1]): continue
        y_low[i] = max(y_low[i-1], y_low[i])
        y_high[i] = max(y_high[i-1], y_high[i])

    f_y = interpolate.interp1d(np.log10(y), np.log10(x))
    f_y_high = interpolate.interp1d(np.log10(y_high), np.log10(x))
    f_y_low = interpolate.interp1d(np.log10(y_low), np.log10(x))

    return (
        10**f_y_low(np.log10(p)),
        10**f_y(np.log10(p)),
        10**f_y_high(np.log10(p))
    )

def convergence_n():
    """ right panel of the survival curve plot. Very similar to above except in
    the places I comment.
    """

    suite = "SymphonyMilkyWay"
    suite_hr = "SymphonyMilkyWayHR"

    npeak = []
    m_c = []
    m_rs = []
    censored_c = []
    censored_rs = []

    npeak_hr = []
    m_c_hr = []
    m_rs_hr = []
    censored_c_hr = []
    censored_rs_hr = []

    # Collect all the data again
    for suite_i in [suite, suite_hr]:
        n_hosts = symlib.n_hosts(suite_i)
        for i_host in range(n_hosts):
            print(suite_i, i_host)
            sim_dir = symlib.get_host_directory(base_dir, suite_i, i_host)
        
            param = symlib.simulation_parameters(sim_dir)
            mp = param["mp"]/param["h100"]

            h, hist = symlib.read_subhalos(sim_dir)
            c = symlib.read_cores(sim_dir)
            mpeak_pre = hist["mpeak_pre"]

            last_snap_rs = np.zeros(len(h), dtype=int)
            last_snap_c = np.zeros(len(h), dtype=int)
            snap = np.arange(h.shape[1], dtype=int)
            for i in range(1, len(h)):
                ok_rs, ok_c = c["ok_rs"][i], c["ok"][i]
                if np.sum(ok_c) == 0:
                    last_snap_c[i] = hist["first_infall_snap"][i]
                else:
                    last_snap_c[i] = np.max(snap[ok_c])
                if np.sum(ok_rs) == 0:
                    last_snap_rs[i] = hist["first_infall_snap"][i]
                else:
                    last_snap_rs[i] = np.max(snap[ok_rs])
            last_snap_c[0], last_snap_rs[0] = len(snap)-1, len(snap)-1
            
            idx = np.arange(len(c), dtype=int)
            if ratio_type == "mass":
                dm_c = c["m_bound"][idx,last_snap_c]/mpeak_pre
                dm_rs = h["mvir"][idx,last_snap_rs]/mpeak_pre
            elif ratio_type == "vmax":
                vmpeak = np.zeros(len(h))
                snap = np.arange(h.shape[1], dtype=int)
                for i in range(1, len(h)):
                    ok = snap <= hist["first_infall_snap"][i]
                    i_peak = np.argmax(h["mvir"][i,ok])
                    vmpeak[i] = h["vmax"][i,i_peak]
                vmpeak[0] = 1
                dm_c = c["vmax"][idx,last_snap_c]/vmpeak
                dm_rs = h["vmax"][idx,last_snap_rs]/vmpeak
            else:
                assert(0)
    
            ratio_ok = ((hist["merger_ratio"] < max_ratio) &
                        (hist["first_infall_snap"] < last_snap_c) &
                        (hist["first_infall_snap"] < last_snap_rs))
            ratio_ok[0] = False

            if suite_i == suite:
                censored_rs.append(c["ok_rs"][ratio_ok,-1])
                censored_c.append(c["ok"][ratio_ok,-1])
                m_rs.append(dm_rs[ratio_ok])
                m_c.append(dm_c[ratio_ok])
                npeak.append((mpeak_pre/mp)[ratio_ok])
            elif suite_i == suite_hr:
                censored_rs_hr.append(c["ok_rs"][ratio_ok,-1])
                censored_c_hr.append(c["ok"][ratio_ok,-1])
                m_rs_hr.append(dm_rs[ratio_ok])
                m_c_hr.append(dm_c[ratio_ok])
                npeak_hr.append((mpeak_pre/mp)[ratio_ok])
                
    censored_rs = np.hstack(censored_rs)
    censored_c = np.hstack(censored_c)
    m_rs = np.hstack(m_rs)
    m_c = np.hstack(m_c)
    npeak = np.hstack(npeak)

    if COMBINE_HR:
        censored_rs = np.hstack([censored_rs, np.hstack(censored_rs_hr)])
        censored_c = np.hstack([censored_c, np.hstack(censored_c_hr)])
        m_rs = np.hstack([m_rs, np.hstack(m_rs_hr)])
        m_c = np.hstack([m_c, np.hstack(m_c_hr)])
        npeak = np.hstack([npeak, np.hstack(npeak_hr)])

    # now we construct a moving window of 200 subhaloes to analyze
    order = np.argsort(npeak)
    censored_rs, censored_c = censored_rs[order], censored_c[order]
    m_rs, m_c, npeak = m_rs[order], m_c[order], npeak[order]
    
    n_eval = 100
    m_eval = 10**np.linspace(-4, 1, n_eval)

    n_window = 200
    n_npeak = 50
    npeak_eval = 10**np.linspace(np.log10(300), np.log10(1e6), n_npeak)

    c_10, c_50, c_90 = np.zeros(n_npeak), np.zeros(n_npeak), np.zeros(n_npeak)
    rs_10, rs_50, rs_90 = np.zeros(n_npeak), np.zeros(n_npeak), np.zeros(n_npeak)

    fig1, ax_m = plt.subplots()

    for i in range(len(npeak_eval)):
        # for each central bin, we find the center
        j = np.searchsorted(npeak, npeak_eval[i])

        # if we're out of range, give up
        if j < n_window or j + n_window > len(npeak):
            c_10[i], c_50[i], c_90[i] = np.nan, np.nan, np.nan
            rs_10[i], rs_50[i], rs_90[i] = np.nan, np.nan, np.nan
            continue

        # extract m and censorship values within the window
        m_rs_i = m_rs[j-n_window: j+n_window]
        m_c_i = m_c[j-n_window: j+n_window]
        cen_rs_i = censored_rs[j-n_window: j+n_window]
        cen_c_i = censored_c[j-n_window: j+n_window]

        # calculate KM
        S_m_c, err_m_c = kaplan_meier(m_c_i, cen_c_i, m_eval, decreasing=True)
        S_m_rs, err_m_rs = kaplan_meier(m_rs_i, cen_rs_i, m_eval,
                                        decreasing=True)

        # I'm not currently using the 0.5 quantile values.
        _, rs_50[i], _ = survival_quantile(m_eval, S_m_rs, err_m_rs, 0.5)
        _, rs_90[i], _ = survival_quantile(m_eval, S_m_rs, err_m_rs, 0.9)
        _, c_50[i], _ = survival_quantile(m_eval, S_m_c, err_m_c, 0.5)
        _, c_90[i], _ = survival_quantile(m_eval, S_m_c, err_m_c, 0.9)
        

    # Everything from here on in this function is just plotting nonsense.

    x_arrow = [10**2.5, 10**3, 10**3.5, 10**4, 10**4.5]
    dlog_arrow = -np.log10(1.1)
    if ratio_type == "mass": dlog_arrow *= 2
    ok = (~np.isnan(c_90)) & (~np.isnan(c_90)) & (~np.isnan(c_90))
    plot_arrow_lines(ax_m, npeak_eval[ok], rs_90[ok],
                     x_arrow, dlog_arrow, pc("b"))
    ax_m.plot(npeak_eval[ok], rs_90[ok], "-", c=pc("b"),
              label=r"${\rm Rockstar}$")
    ax_m.plot(npeak_eval[ok], c_90[ok], "-", c=pc("r"),
              label=r"${\rm Particle{-}tracking}$")
    plot_arrow_lines(ax_m, npeak_eval[ok], c_90[ok],
                     x_arrow, dlog_arrow, pc("r"))

    ax_m.set_xscale("log")
    ax_m.set_yscale("log")
    ax_m.set_xlabel(r"$n_{\rm peak,sub}$")
    if ratio_type == "mass":
        ax_m.set_ylabel(r"$m_{\rm sub,90}/m_{\rm sub,peak}$")
    elif ratio_type == "vmax":
        ax_m.set_ylabel(r"$v_{\rm max,90}/v_{\rm max,Mpeak}$")

    xlo, xhi = np.min(npeak_eval[ok]), np.max(npeak_eval[ok])
    ylo, yhi = ax_m.get_ylim()
    ax_m.set_xlim(xlo, xhi)
    ax_m.set_ylim(ylo, 1)

    if ratio_type == "mass":
        npeak_ref = np.array([xlo, xhi])
        for log_n_disrupt in [1, 2, 3, 4, 5]:
            ax_m.plot(npeak_ref, 10**log_n_disrupt/npeak_ref, 
                      "--", lw=1, c=pc("a"))
            if log_n_disrupt > 1 and log_n_disrupt < 5:
                text = r"$n_{\rm sub,90} = 10^%d$" % log_n_disrupt if log_n_disrupt == 2 else r"$10^%d$" % log_n_disrupt
                if log_n_disrupt != 4:
                    ax_m.text(10**(log_n_disrupt+1)*0.6, 0.2, text,
                              color=pc("a"), fontsize=17,
                              horizontalalignment="left")
                else:
                    ax_m.text(3e4, 0.2, text,
                              color=pc("a"), fontsize=17,
                              horizontalalignment="left")

    elif ratio_type == "vmax":
        ax_m.fill_between([xlo, xhi], [0.466]*2, [0.544]*2, alpha=0.3,
                          color=pc("p"), label=r"${B19}$")
        m_low_s16 = 1/23.94 # reff/rvir < 0.025
        m_high_s16 = 1/8.6 # reff/rvir > 0.04
        v_low = j_vdb_16_vmax_ratio(m_low_s16)
        v_high = j_vdb_16_vmax_ratio(m_high_s16)
        ax_m.fill_between([xlo, xhi], [v_low]*2, [v_high]*2,
                          alpha=0.3, color=pc("o"),
                          label=r"${\rm S16\,+\,JvdB16}$")
        ax_m.set_ylim(0.1, 1)
        
    ax_m.legend(loc="lower left", frameon=True, fontsize=17)

    if ratio_type == "mass":
        fig1.savefig("../plots/core_tracking/n_converge.pdf")
    elif ratio_type == "vmax":
        fig1.savefig("../plots/core_tracking/n_converge_vmax.pdf")

def main():
    palette.configure(True)
    # The left panel of the survival curve plot
    survival_time()
    # The right panel of the survival curve plot
    convergence_n()

if __name__ == "__main__": main()
