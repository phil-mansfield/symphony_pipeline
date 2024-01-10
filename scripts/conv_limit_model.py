import numpy as np
import symlib
import matplotlib.pyplot as plt
import palette
from palette import pc
import scipy.stats as stats
import survival

base_dir = "/sdf/home/p/phil1/ZoomIns"
#suites = ["SymphonyMilkyWay"]
#suites = ["SymphonyLMC", "SymphonyMilkyWay", "SymphonyLCluster", 
#          "SymphonyGroup", "SymphonyMilkyWayHR"]
suites = ["SymphonyLCluster"]

def kaplan_meier_contours(n_peak, m_final, censored, bins,
                          p=(0.5-0.68/2, 0.5, 0.5+0.68/2),
                          n_min=100, m_eval=10**np.linspace(-4, 0, 200)):
    
    mid = np.sqrt(bins[1:]*bins[:-1])
    n = len(bins) -1
    ok = np.zeros(n, dtype=bool)
    out = [np.zeros(n) for _ in p]

    for i in range(n):
        low, high = bins[i], bins[i+1]
        idx = np.where((low < n_peak) & (n_peak <= high))[0]
        if len(idx) < n_min: continue
        
        ok[i] = True
        S, err = survival.kaplan_meier(m_final[idx], censored[idx], m_eval,
                                       decreasing=True)
        for j in range(len(p)):
            out[i] = m_eval[np.searchsorted(S, p[j])]
    
    return mid[ok], [x[ok] for x in out]

def mpeak_loss(m, ok):
    out = np.zeros(len(m))
    idx = np.where(ok)[0][::-1]
    mpeak = m[idx[0]]
    for i in idx:
        if m[i] > mpeak: mpeak = m[i]
        out[i] = mpeak
    return out

def get_lims(sim_dir):
    param = symlib.simulation_parameters(sim_dir)
    mp = param["mp"]/param["h100"]
    h, hist = symlib.read_subhalos(sim_dir)
    c = symlib.read_cores(sim_dir)
    
    lim_eps = np.zeros(len(h))
    lim_disc = np.zeros(len(h))
    lim_relax = np.zeros(len(h))
    ok = np.zeros(len(h), dtype=bool)
    npeak = np.zeros(len(h))

    for i in range(1, len(h)):
        h0 = hist[i]
        if  (h0["conv_snap_eps"] == -1 or h0["conv_snap_discrete"] == -1 or
             h0["conv_snap_relax"] == -1 or h0["merger_ratio"] > 0.1):
            continue

        m = mpeak_loss(c["m_bound"][i], c["ok"][i])/mp

        npeak[i] = hist[i]["mpeak_pre"]/mp
        #npeak[i] = m[hist["first_infall_snap"][i]]
        lim_eps[i] = m[h0["conv_snap_eps"]]
        n_infall = m[np.where(c["ok"][i])[0][0]]
        lim_disc[i] = m[h0["conv_snap_discrete"]] 
        lim_relax[i] = m[h0["conv_snap_relax"]]
        ok[i] = True

    return lim_eps, lim_disc, lim_relax, npeak, ok

def median_curve(x_bins, x, y):
    low, _, _ = stats.binned_statistic(
        x, y, lambda x: np.quantile(x, 0.5-0.68/2), bins=x_bins)
    med, _, _ = stats.binned_statistic(x, y, "median", bins=x_bins)
    high, _, _ = stats.binned_statistic(
        x, y, lambda x: np.quantile(x, 0.5+0.68/2), bins=x_bins)
    return low, med, high

def final_masses(c, hist, conv_limit_name):
    

def main():
    palette.configure(True)

    lims_eps = []
    lims_disc = []
    lims_relax = []
    npeaks = []

    for suite in suites:
        n_host = symlib.n_hosts(suite)
        for i_host in range(n_host):
            print(suite, i_host)
            sim_dir = symlib.get_host_directory(base_dir, suite, i_host)
            
            lim_eps, lim_disc, lim_relax, npeak, ok = get_lims(sim_dir)
            lims_eps.append(lim_eps[ok])
            lims_disc.append(lim_disc[ok])
            lims_relax.append(lim_relax[ok])
            npeaks.append(npeak[ok])
            

    lim_eps, lim_disc = np.hstack(lims_eps), np.hstack(lims_disc)
    lim_relax, npeak = np.hstack(lims_relax), np.hstack(npeaks)    
    lim_eps = lim_eps#/npeak
    lim_relax = lim_relax#/npeak
    lim_disc = lim_disc#/npeak

    #lim_eps = np.minimum(lim_eps, 1.0)
    #lim_disc = np.minimum(lim_disc, 1.0)
    #lim_relax = np.minimum(lim_relax, 1.0)

    fig, ax = plt.subplots()

    lim_max = np.maximum(lim_relax, np.maximum(lim_disc, lim_eps))

    #ax.plot(npeak, lim_eps, ".", c=pc("r"), alpha=0.1)
    #ax.plot(npeak, lim_disc, ".", c=pc("o"), alpha=0.1)
    #ax.plot(npeak, lim_relax, ".", c=pc("b"), alpha=0.1)
    #ax.plot(npeak, lim_max, ".", c=pc("k"), alpha=0.1)

    edges = 10**np.linspace(np.log10(300), np.log10(2e5), 21)
    mids = np.sqrt(edges[1:]*edges[:-1])

    eps_lo, eps_med, eps_hi = median_curve(edges, npeak, lim_eps)
    disc_lo, disc_med, disc_hi = median_curve(edges, npeak, lim_disc)
    relax_lo, relax_med, relax_hi = median_curve(edges, npeak, lim_relax)
    max_lo, max_med, max_hi = median_curve(edges, npeak, lim_max)

    n, _ = np.histogram(npeak, bins=edges)
    ok = n > 20

    ax.plot(mids[ok], eps_med[ok], pc("o"),
            label=r"${\rm force{-}softening\ limit}$")
    #ax.fill_between(mids, eps_lo, eps_hi, color=pc("r"), alpha=0.2)
    ax.plot(mids[ok], disc_med[ok], pc("b"),
            label=r"${\rm discreteness\ limit}$")
    #ax.fill_between(mids, disc_lo, disc_hi, color=pc("o"), alpha=0.2)
    ax.plot(mids[ok], relax_med[ok], pc("r"),
            label=r"${\rm relaxation\ limit}$")
    #ax.fill_between(mids, relax_lo, relax_hi, color=pc("b"), alpha=0.2)
    ax.plot(mids[ok], max_med[ok], pc("k"),
            label=r"${\rm maximum\ limit}$")
    ax.fill_between(mids[ok], max_lo[ok], max_hi[ok], color=pc("k"), alpha=0.2)

    b_med, a_med = np.polyfit(np.log10(mids[ok]), np.log10(max_med[ok]), 1)
    b_lo, a_lo = np.polyfit(np.log10(mids[ok]), np.log10(max_lo[ok]), 1)
    b_hi, a_hi = np.polyfit(np.log10(mids[ok]), np.log10(max_hi[ok]), 1)

    print(10**a_med, b_med)
    print(10**a_lo, b_lo)
    print(10**a_hi, b_hi)

    cut_expected = 10**(a_med + b_med*np.log10(npeak))
    log_disp = np.std(np.log10(lim_max) - np.log10(cut_expected))
    print(log_disp)

    cut_med = 10**(a_med + b_med*np.log10(mids))
    cut_hi = 10**(a_hi + b_hi*np.log10(mids))
    cut_lo = 10**(a_lo + b_lo*np.log10(mids))

    plt.plot(mids, cut_med, "--", c=pc("p"), lw=2, label=r"${\rm fit}$")
    plt.plot(mids, cut_lo, ":", c=pc("p"), lw=2)
    plt.plot(mids, cut_hi, ":", c=pc("p"), lw=2)

    ax.legend(loc="lower right", fontsize=17)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(100, None)

    ax.set_xlabel(r"$n_{\rm peak}$")
    ax.set_ylabel(r"$n_{\rm lim}$")

    fig.savefig("../plots/core_tracking/conv_limit_model.pdf")

if __name__ == "__main__": main()
