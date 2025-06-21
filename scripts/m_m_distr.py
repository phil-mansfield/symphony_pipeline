import numpy as np
import matplotlib.pyplot as plt
import palette
from palette import pc
import symlib

palette.configure(False)
USE_NPEAK = False
ONE_HALO_TEST = False

def main():
    n, npeak, r = [], [], []

    r_bins = np.array([1/25, 1/10, 1/5, 1/2.5, 1])
    n_bins = np.array([3e2, 1e3, 3e3, 1e4, 3e4, 1e5])

    r_names = ["0.04", "0.1", "0.2", "0.4", "1"]
    n_names = [r"3\times10^2", r"1\times10^3", r"3\times10^3",
               r"1\times10^4", r"3\times10^4", r"1\times10^5"]
    
    print("r_bins", r_bins)
    print("n_bins", n_bins)
    
    base_dir = "/sdf/home/p/phil1/ZoomIns"
    suites = ["SymphonyMilkyWay", "SymphonyMilkyWayHR"]

    for suite in suites:
        n_host = symlib.n_hosts(suite)
        for i_host in range(n_host):
            print(suite, i_host)
            sim_dir= symlib.get_host_directory(base_dir, suite, i_host)
            param = symlib.simulation_parameters(sim_dir)
            mp = param["mp"]/param["h100"]
            
            rs, _ = symlib.read_rockstar(sim_dir)
            sf, hist = symlib.read_symfind(sim_dir)
            ok = (hist["merger_ratio"] < 0.15) & (sf["ok"][:,-1])

            n.append(sf["m"][ok,-1]/mp)
            npeak.append(hist[ok]["mpeak_pre"]/mp)
            r.append(np.sqrt(np.sum(sf["x"][ok,-1,:]**2, axis=1)) /
                     rs[0,-1]["rvir"])

            if ONE_HALO_TEST:
                break
        
    n = np.hstack(n)
    npeak = np.hstack(npeak)
    r = np.hstack(r)

    colors = [pc("r"), pc("o"), pc("g"), pc("b"), pc("p")]

    bins = 10**np.linspace(-3, 0, 200)

    figax_r = [plt.subplots() for _ in range(len(r_bins)-1)]
    figax_n = [plt.subplots() for _ in range(len(n_bins)-1)]

    n_tot = 0
    for i_r in range(len(r_bins)-1):
        for i_n in range(len(n_bins)-1):
            in_bin_r = (r < r_bins[i_r+1]) & (r > r_bins[i_r])
            if not USE_NPEAK:
                in_bin_n = (n < n_bins[i_n+1]) & (n > n_bins[i_n])
            else:
                in_bin_n = (npeak < n_bins[i_n+1]) & (npeak > n_bins[i_n])
            ok = in_bin_r & in_bin_n
            
            print("%4d" % np.sum(ok), end=" ")
            n_tot += np.sum(ok)
            if np.sum(ok) < 10: continue

            fig, ax = figax_r[i_r]
            c = colors[i_n]
            n_type = r"n_{\rm peak}" if USE_NPEAK else "n"
            ax.hist(n[ok]/npeak[ok],
                    bins=bins, cumulative=True, density=True,
                    histtype="step", lw=3, color=c,
                    label=r"$%s<%s<%s$" % (n_names[i_n], n_type, n_names[i_n+1]))

            fig, ax = figax_n[i_n]
            c = colors[i_r]
            ax.hist(n[ok]/npeak[ok],
                    bins=bins, cumulative=True, density=True,
                    histtype="step", lw=3, color=c,
                    label=r"$%s<r/R_{\rm vir}<%s$" % (r_names[i_r], r_names[i_r+1]))
        print()
    print(n_tot)

    out_dir = "../plots/binned_m_m"
    for i_r in range(len(figax_r)):
        fig, ax = figax_r[i_r]
        fig.suptitle(r"$%s<r/R_{\rm vir}<%s$" % (r_names[i_r], r_names[i_r+1]))
        ax.set_xlabel(r"$\mu \equiv m/m_{\rm infall}$")
        ax.set_ylabel(r"$P(<\mu)$")
        ax.set_xlim(1e-3, 1)
        ax.set_xscale("log")
        ax.legend(loc="upper left", fontsize=17)
        if USE_NPEAK:
            fig.savefig("%s/p_rbin_%d.png" % (out_dir, i_r))
        else:
            fig.savefig("%s/rbin_%d.png" % (out_dir, i_r))

    for i_n in range(len(figax_n)):
        fig, ax = figax_n[i_n]
        if not USE_NPEAK:
            fig.suptitle(r"$%s<n<%s$" % (n_names[i_n], n_names[i_n+1]))
        else:
            fig.suptitle(r"$%s<n_{\rm peak}<%s$" %
                         (n_names[i_n], n_names[i_n+1]))
        ax.set_xlabel(r"$\mu \equiv m/m_{\rm infall}$")
        ax.set_ylabel(r"$P(<\mu)$")
        ax.set_xlim(1e-3, 1)
        ax.set_xscale("log")
        ax.legend(loc="upper left", fontsize=17)
        if USE_NPEAK:
            fig.savefig("%s/p_nbin_%d.png" % (out_dir, i_n))
        else:
            fig.savefig("%s/nbin_%d.png" % (out_dir, i_n))
            
if __name__ == "__main__": main()
