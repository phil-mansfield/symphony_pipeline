import numpy as np
import matplotlib.pyplot as plt
import symlib
import palette
import matplotlib as mpl

from palette import pc
palette.configure(True)

def plot_mass_funcs():
    # Get file locations.
    base_dir = "/sdf/home/p/phil1/ZoomIns/"
    #suite = "SymphonyMilkyWay"
    #halo_name = "Halo416"
    suite = "SymphonyMilkyWayHR"
    halo_name = 0
    #suite = "TestSIDM"
    #halo_name = 0
    
    fig, ax = plt.subplots()
    
    # Set up histogram bins
    param = symlib.simulation_parameters(suite)
    m_min = param["mp"]/param["h100"]*300
    m_max = 1e12 
    n_bins = 200
    bins = 10**np.linspace(np.log10(m_min), np.log10(m_max), n_bins+1)

    # Read in halo catalogs
    sim_dir = symlib.get_host_directory(base_dir, suite, halo_name)
    rs, hist = symlib.read_rockstar(sim_dir)
    sf, _    = symlib.read_symfind(sim_dir)

    # Compute subhalo radii
    r_rs = np.sqrt(np.sum(rs["x"][:,-1]**2, axis=1))
    r_sf = np.sqrt(np.sum(sf["x"][:,-1]**2, axis=1))

    # Loop over R_max cutoffs
    mult = [0.25, 0.5, 1.0]
    colors = [pc("b"), pc("o"), pc("r")]
    for i in range(len(mult)):
        # Select halos in R_max cutoff
        r_max = rs["rvir"][0,-1]*mult[i]
        ok_rs = rs["ok"][:,-1] & (r_rs < r_max)
        ok_rs_2 = rs["ok"][:,-1] & (r_rs < r_max) & (sf["f_core_rs"][:,-1] > 0)
        ok_sf = sf["ok"][:,-1] & (r_sf < r_max)

        # Get cumulative mass functions with numpy array tricks
        n_rs, _ = np.histogram(hist["mpeak_pre"][ok_rs], bins=bins)
        n_rs_2, _ = np.histogram(hist["mpeak_pre"][ok_rs_2], bins=bins)
        n_sf, _ = np.histogram(hist["mpeak_pre"][ok_sf], bins=bins)
        
        N_rs = np.cumsum(n_rs[::-1])[::-1]
        N_rs_2 = np.cumsum(n_rs_2[::-1])[::-1]
        N_sf = np.cumsum(n_sf[::-1])[::-1]
    
        left_bins = bins[:-1]
    
        plt.plot(left_bins, N_rs, "--", c=colors[i])
        plt.plot(left_bins, N_rs_2, ":", c=colors[i], lw=2)
        plt.plot(left_bins, N_sf, "-", c=colors[i])

    # Plotting nonsense
    plt.plot([], [], pc("r"), label=r"$r < R_{\rm vir}$")
    plt.plot([], [], pc("o"), label=r"$r < R_{\rm vir}/2$")
    plt.plot([], [], pc("b"), label=r"$r < R_{\rm vir}/4$")
    plt.plot([], [], pc("k"), label=r"${\rm Symfind}$")
    plt.plot([], [], "--", c=pc("k"), label=r"${\rm Rockstar}$")
    plt.plot([], [], ":", c=pc("k"), label=r"${\rm Rockstar}\ (f_{\rm core} > 0)$", lw=2)
    plt.legend(loc="upper right", fontsize=17)
    
    ax.set_xlim(m_min, 1e12)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$M_{\rm sub,peak}$")
    ax.set_ylabel(r"$N(>M_{\rm sub,peak})$")
    fig.savefig("../plots/sidm/mass_func_%s.png" % suite)

def particle_counts():
    base_dir = "/sdf/home/p/phil1/ZoomIns/"
    suite = "TestSIDM"
    halo_name = 0

    sim_dirs = [symlib.get_host_directory(base_dir, suite, halo_name)]
    
    suite = "SymphonyMilkyWayHR"
    for i in range(symlib.n_hosts(suite)):
        sim_dirs.append(symlib.get_host_directory(base_dir, suite, i))
    

    for sim_dir in sim_dirs:
        print(sim_dir)
        part = symlib.Particles(sim_dir)
        p = part.read(235, mode="all")
        n = np.array([len(pi) for pi in p], dtype=int)
        print(n[:6])
        print(n[-10:])
        print()

def inspect_subhalos():
    plot_dir = "../plots/sidm/"
    base_dir = "/sdf/home/p/phil1/ZoomIns/"
    #suite = "SymphonyMilkyWayHR"
    suite = "TestSIDM"
    halo_name = 0

    sim_dir = symlib.get_host_directory(base_dir, suite, halo_name)
    param = symlib.simulation_parameters(sim_dir)
    mp = param["mp"]/param["h100"]
    
    rs, hist = symlib.read_rockstar(sim_dir)
    sf, _ = symlib.read_symfind(sim_dir)

    ok = rs["ok"][:,-1] & (sf["f_core_rs"][:,-1] == 0) & (~sf["ok"][:,-1])
    targets = np.where(ok)[0]

    print(len(targets))
    #print(targets)
    targets = np.hstack([targets[:3], targets[-3:]])
    print(targets)
    print("n_rs")
    print(hist["mpeak_pre"][targets]/mp)
    
    part = symlib.Particles(sim_dir)
    fig, ax = plt.subplots(figsize=(8, 8))
    for idx in targets:
        print(idx)
        ax.cla()
        
        mode = "current"
        p = part.read(235, mode=mode, halo=idx)
        cores = p[part.core_indices(235, mode=mode, halo=idx)]

        host, sub = rs[0], rs[idx]
        lim = 1.5*host["rvir"][-1]
        
        symlib.plot_circle(
            ax, host["x"][-1,0], host["x"][-1,1],
            host["rvir"][-1], c=pc("b"), ls=":"
        )
        symlib.plot_circle(
            ax, host["x"][-1,0], host["x"][-1,1],
            host["rvir"][-1], c=pc("b"), ls="-"
        )

        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)

        norm = mpl.colors.LogNorm(vmin=1, vmax=10000)
        kwargs = {"extent": [-lim, lim, -lim, lim],
                  "norm": norm, "cmap": "inferno", "gridsize": 200}
        ax.hexbin(p["x"][:,0], p["x"][:,1], **kwargs)
        ax.plot(cores["x"][:,0], cores["x"][:,1], ".", c=pc("b"))
        
        fig.savefig("%s/halo_%04d.png" % (plot_dir, idx))
        
def main():
    #plot_mass_funcs()
    #particle_counts()
    inspect_subhalos()
    
if __name__ == "__main__": main()
