import numpy as np
import symlib
import matplotlib.pyplot as plt
import palette
from palette import pc
import scipy.stats as stats
import survival

base_dir = "/sdf/home/p/phil1/ZoomIns"
#suites = ["SymphonyMilkyWay"]
suites = ["SymphonyLMC", "SymphonyMilkyWay", "SymphonyLCluster",
          #"SymphonyGroup",
          "SymphonyMilkyWayHR"]
#suites = ["SymphonyMilkyWay"]
#suites = ["SymphonyLCluster"]

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

        m_eval_i, S = m_eval[~np.isnan(S)], S[~np.isnan(S)]

        for j in range(len(p)):
            k = np.searchsorted(S, p[j])
            if k < len(m_eval_i):
                out[j][i] = m_eval_i[k]
            else:
                out[j][i] = 1
    
    return out, ok

def mpeak_loss(m, ok):
    out = np.zeros(len(m))
    idx = np.where(ok)[0][::-1]
    mpeak = m[idx[0]]
    for i in idx:
        if m[i] > mpeak: mpeak = m[i]
        out[i] = mpeak
    return out

LIM_DTYPE = [
    ("n_eps", "f4"), ("n_disc", "f4"), ("n_relax", "f4"), ("n_max", "f4"),
    ("c_eps", "?"),  ("c_disc", "?"),  ("c_relax", "?"),  ("c_max", "?"),
    ("npeak", "f4"),
    ("ok", "?")
]

def get_lims(sim_dir):
    param = symlib.simulation_parameters(sim_dir)
    mp = param["mp"]/param["h100"]
    h, hist = symlib.read_subhalos(sim_dir)
    c = symlib.read_cores(sim_dir)
    
    lim = np.zeros(len(h), dtype=LIM_DTYPE)
    n_snap = h.shape[1]

    for i in range(1, len(h)):
        if  hist[i]["merger_ratio"] > 0.1: continue

        n = mpeak_loss(c["m_bound"][i], c["ok"][i])/mp

        lim["npeak"][i] = hist[i]["mpeak_pre"]/mp
        #lim["npeak"][i] = n[hist["first_infall_snap"][i]]
        lim["ok"][i] = True

        for name, m_name, c_name in zip(
                ["conv_snap_eps", "conv_snap_discrete", "conv_snap_relax"],
                ["n_eps", "n_disc", "n_relax"],
                ["c_eps", "c_disc", "c_relax"]
        ):
            if hist[i][name] != -1:
                final_snap = hist[i][name]
            elif c["ok"][i, -1]:
                final_snap = n_snap-1
            else:
                final_snap = hist[i]["disrupt_snap"]
            """
            if m_name == "n_disc":
                lim[i][m_name] = 0.32*(lim["npeak"][i]/1e3)**-0.8
            else:
                lim[i][m_name] = n[final_snap]/lim["npeak"][i]
            """
            lim[i][m_name] = n[final_snap]/lim["npeak"][i]
            lim[i][c_name] = final_snap != hist[i][name]

    lim["n_max"] = np.max(
        [lim["n_eps"], lim["n_disc"], lim["n_relax"]], axis=0)
    lim["c_max"] = lim["c_eps"] & lim["c_disc"] & lim["c_relax"]
    
    return lim

def median_curve(x_bins, x, y):
    low, _, _ = stats.binned_statistic(
        x, y, lambda x: np.quantile(x, 0.5-0.68/2), bins=x_bins)
    med, _, _ = stats.binned_statistic(x, y, "median", bins=x_bins)
    high, _, _ = stats.binned_statistic(
        x, y, lambda x: np.quantile(x, 0.5+0.68/2), bins=x_bins)
    return low, med, high
    

def main():
    palette.configure(True)

    lim = []

    for suite in suites:
        n_host = symlib.n_hosts(suite)
        for i_host in range(n_host):
            print(suite, i_host)
            sim_dir = symlib.get_host_directory(base_dir, suite, i_host)
            
            lim.append(get_lims(sim_dir))

    lim = np.hstack(lim)

    fig, ax = plt.subplots()

    edges = 10**np.linspace(np.log10(300), np.log10(2e5), 21)
    mids = np.sqrt(edges[1:]*edges[:-1])
    
    n_min = 100
    (max_lo, max_med, max_hi), ok = kaplan_meier_contours(
        lim["npeak"], lim["n_max"], lim["c_max"], edges, n_min=n_min)
    (_, eps_med, _), _ = kaplan_meier_contours(
        lim["npeak"], lim["n_eps"], lim["c_eps"], edges, n_min=n_min)
    (_, disc_med, _), _ = kaplan_meier_contours(
        lim["npeak"], lim["n_disc"], lim["c_disc"], edges, n_min=n_min)
    (_, relax_med, _), _ = kaplan_meier_contours(
        lim["npeak"], lim["n_relax"], lim["c_relax"], edges, n_min=n_min)

    _, disc_med_2, _ = median_curve(edges, lim["npeak"], lim["n_disc"])
    
    print(edges[:-1][ok])
    print(disc_med[ok])
    print(disc_med_2[ok])

    ax.plot(mids[ok], eps_med[ok]*mids[ok], pc("o"),
            label=r"${\rm force{-}softening\ limit}$")
    ax.plot(mids[ok], disc_med[ok]*mids[ok], pc("b"),
            label=r"${\rm discreteness\ limit}$")
    ax.plot(mids[ok], relax_med[ok]*mids[ok], pc("r"),
            label=r"${\rm relaxation\ limit}$")
    ax.plot(mids[ok], max_med[ok]*mids[ok], pc("k"),
            label=r"${\rm maximum\ limit}$")
    ax.fill_between(mids[ok], max_lo[ok]*mids[ok], max_hi[ok]*mids[ok],
                    color=pc("k"), alpha=0.2)

    c_med, b_med, a_med = np.polyfit(np.log10(mids[ok]),
                                     np.log10(max_med[ok]*mids[ok]), 2)
    c_lo, b_lo, a_lo = np.polyfit(np.log10(mids[ok]),
                                  np.log10(max_lo[ok]*mids[ok]), 2)
    c_hi, b_hi, a_hi = np.polyfit(np.log10(mids[ok]),
                                  np.log10(max_hi[ok]*mids[ok]), 2)

    print(a_med, b_med, c_med)
    print(a_lo, b_lo, c_lo)
    print(a_hi, b_hi, c_hi)

    cut_med = 10**(a_med + b_med*np.log10(mids) + c_med*np.log10(mids)**2)
    cut_hi = 10**(a_hi + b_hi*np.log10(mids) + c_hi*np.log10(mids)**2)
    cut_lo = 10**(a_lo + b_lo*np.log10(mids) + c_lo*np.log10(mids)**2)

    plt.plot(mids, cut_med, "--", c=pc("p"), lw=2, label=r"${\rm fit}$")
    plt.plot(mids, cut_lo, ":", c=pc("p"), lw=2)
    plt.plot(mids, cut_hi, ":", c=pc("p"), lw=2)

    ax.legend(loc="upper left", fontsize=17)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(1e2, 1e4)

    ax.set_xlabel(r"$n_{\rm peak}$")
    ax.set_ylabel(r"$n_{\rm lim}$")

    fig.savefig("../plots/core_tracking/conv_limit_model.pdf")

if __name__ == "__main__": main()
