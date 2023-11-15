import numpy as np
import symlib
import matplotlib.pyplot as plt
import palette
from palette import pc

base_dir = "/sdf/group/kipac/g/cosmo/ki21/phil1/simulations/ZoomIns"
halo_finder = "symfind"

def shmf():
    fig, ax = plt.subplots(1, 2, figsize=(16,8), sharey=True)
    ax_mvir, ax_mpeak = ax

    suites = ["SymphonyMilkyWayDiskDMO", "SymphonyMilkyWayDisk"]
    styles = ["--", "-"]
    r_cuts = [1.0, 1/3, 1/10]
    colors = [pc("r"), pc("o"), pc("b")]

    ax_mpeak.plot([], [], label=r"${\rm unevolved}$", c="k")
    ax_mpeak.plot([], [], label=r"$r < R_{\rm vir}$", c=pc("r"))
    ax_mpeak.plot([], [], label=r"$r < R_{\rm vir}/3$", c=pc("o"))
    ax_mpeak.plot([], [], label=r"$r < R_{\rm vir}/10$", c=pc("b"))
    ax_mpeak.plot([], [], label=r"${\rm DMO}$", ls="--", c=pc("a"))
    ax_mpeak.plot([], [], label=r"${\rm Disk}$", ls="-", c=pc("a"))
    ax_mpeak.legend(loc="upper right", fontsize=16)

    m_bins = 10**np.linspace(9, 12, 21)

    for i_suite, suite in enumerate(suites):
        n_host = symlib.n_hosts(suite)

        unev_n = np.zeros(len(m_bins)-1)
        mpeak_n = np.zeros((len(r_cuts), len(m_bins)-1))
        mvir_n = np.zeros((len(r_cuts), len(m_bins)-1))

        for i_host in range(n_host):
            sim_dir = symlib.get_host_directory(base_dir, suite, i_host)
            param = symlib.simulation_parameters(sim_dir)
            mp = param["mp"]/param["h100"]

            if halo_finder == "rockstar":
                h, hist = symlib.read_rockstar(sim_dir)
                host, h, hist = h[0], h[1:], hist[1:]
            elif halo_finder == "symfind":
                host = symlib.read_rockstar(sim_dir)[0][0]
                h, hist = symlib.read_symfind(sim_dir)
                h, hist = h[1:], hist[1:]
            r = np.sqrt(np.sum(h["x"]**2, axis=2))

            unev_n += np.histogram(hist["mpeak_pre"], bins=m_bins)[0]

            for i_cut, cut in enumerate(r_cuts):
                ok = h["ok"][:,-1] & (r[:,-1] < host["rvir"][-1]*cut)
                mpeak_n[i_cut] += np.histogram(
                    hist["mpeak_pre"][ok], bins=m_bins)[0]
                mvir_n[i_cut] += np.histogram(
                    h["m"][ok,-1], bins=m_bins)[0]

        unev_n = np.cumsum(unev_n[::-1])[::-1] / n_host
        mpeak_n = np.cumsum(mpeak_n[:,::-1],axis=1)[:,::-1] / n_host
        mvir_n = np.cumsum(mvir_n[:,::-1], axis=1)[:,::-1] / n_host
             
        # Plotting
        ls = styles[i_suite]
        
        m = m_bins[:-1]
        ax_mpeak.plot(m, unev_n, c="k", ls=ls)
        
        for i_cut in range(len(r_cuts)):
            ax_mpeak.plot(m, mpeak_n[i_cut], c=colors[i_cut], ls=ls)
            ax_mvir.plot(m, mvir_n[i_cut], c=colors[i_cut], ls=ls)

    ax_mpeak.set_yscale("log")
    ax_mpeak.set_xscale("log")
    ax_mpeak.set_xlim(1e9, 1e12)
    ax_mpeak.set_ylim(0.01, 100)
    ax_mpeak.set_xlabel(r"$m_{\rm peak}\ (M_\odot)$")
    ax_mpeak.set_ylabel(r"$N(>m_{\rm peak})$")

    ax_mvir.set_yscale("log")
    ax_mvir.set_xscale("log")
    ax_mvir.set_xlim(1e9, 1e12)
    ax_mvir.set_ylim(0.01, 100)
    ax_mvir.set_xlabel(r"$m_{\rm sub}\ (M_\odot)$")
    ax_mvir.set_ylabel(r"$N(>m_{\rm sub})$")
    
    fig.savefig("../plots/disk_shmf_test_%s.pdf" % halo_finder)

def match_sims():
    base_dir = "/sdf/group/kipac/g/cosmo/ki21/phil1/simulations/ZoomIns"
    suite_1, suite_2 = "SymphonyMilkyWayDisk", "SymphonyMilkyWayDiskDMO"
    suite_3 = "SymphonyMilkyWay"
    name_1, name_2 = r"${\rm Disk}$", r"${\rm DMO\ 4K}$"
    name_3 = r"${\rm DMO\ 8K}$"
    host_i = 0
    out_dir = "../plots/disk_tests/"

    sim_dir_1 = symlib.get_host_directory(base_dir, suite_1, host_i)
    sim_dir_2 = symlib.get_host_directory(base_dir, suite_2, host_i)
    sim_dir_3 = symlib.get_host_directory(base_dir, suite_3, host_i)
    
    h_1, hist_1 = symlib.read_rockstar(sim_dir_1)
    h_2, hist_2 = symlib.read_rockstar(sim_dir_2)
    h_3, hist_3 = symlib.read_rockstar(sim_dir_3)

    h_cmov_1, _ = symlib.read_rockstar(sim_dir_1, comoving=True)
    h_cmov_2, _ = symlib.read_rockstar(sim_dir_2, comoving=True)
    h_cmov_3, _ = symlib.read_rockstar(sim_dir_3, comoving=True)
    
    match_32, _ = symlib.match_subhalos(h_cmov_3, h_cmov_2, hist_3, hist_2)
    match_31, _ = symlib.match_subhalos(h_cmov_3, h_cmov_1, hist_3, hist_1)
    print(match_32[:20])
    print(match_31[:20])
    print(match_32[14])
    print(match_31[14])

def symfind_tests():
    i_host = 0
    i_subs = [10, 10]
    snap = 162
    suites = ["SymphonyMilkyWayDiskDMO", "SymphonyMilkyWayDisk"]
    names = [r"${\rm DMO{-}4K}$", r"${\rm Disk{-}8K}$"]
    suffix = ["4k", "4k_disk"]

    for i_suite in range(2):
        i_sub = i_subs[i_suite]
        suite = suites[i_suite]

        sim_dir = symlib.get_host_directory(base_dir, suite, i_host)
        param = symlib.simulation_parameters(sim_dir)
        mp = param["mp"]/param["h100"]

        rs, hist = symlib.read_rockstar(sim_dir)
        sf, _ = symlib.read_symfind(sim_dir)
        part = symlib.Particles(sim_dir)
    
        p = part.read(snap)[i_sub]

        print(rs["x"][i_sub,snap], "%.2f" % rs[i_sub,snap]["rvir"])
        print("%g" % rs[i_sub,snap]["m"], len(p))

        fig, ax = plt.subplots()
        r_max = rs[0,snap]["rvir"]*1.5

        symlib.plot_circle(ax, rs[0,snap]["x"][0], rs[0,snap]["x"][1],
                           rs[0,snap]["rvir"])
        ax.plot(p["x"][:,0], p["x"][:,1], ".", c="k", alpha=0.2)
        ax.set_xlim(-r_max, r_max)
        ax.set_ylim(-r_max, r_max)

        fig.savefig("../plots/disk_tests/particles_%s.png" % suffix[i_suite])


def main():
    palette.configure(False)
    shmf()
    #match_sims()
    #symfind_tests()
    
if __name__ == "__main__": main()
