import numpy as np
import matplotlib.pyplot as plt
import symlib

try:
    import palette
    palette.configure(False)
    from palette import pc
except:
    pc = lambda x: x
    pass

def main():
    base_dir = "/sdf/home/p/phil1/ZoomIns"
    sim_dir = symlib.get_host_directory(base_dir, "SymphonyMilkyWay", 0)
    rs, hist = symlib.read_rockstar(sim_dir)

    param = symlib.simulation_parameters(sim_dir)
    mp = param["mp"]/param["h100"]
    
    i_sub = 1
    snap = hist["first_infall_snap"][i_sub]
    
    min_r, max_r = 0.05*rs["rvir"][i_sub,snap], 10**0.5*rs["rvir"][i_sub,snap]
    
    fig, ax = plt.subplots(1, 3, figsize=(21, 7))
    
    bins = [30, 90, 150]
    colors = [pc("r"), pc("o"), pc("g")]
    part = symlib.Particles(sim_dir)
    p = part.read(snap, mode="current")[i_sub]
    p["x"] -= rs["x"][i_sub,snap]
    p["v"] -= rs["v"][i_sub,snap]
    
    for j in range(len(bins)):
        r_bins = 10**np.linspace(np.log10(min_r), np.log10(max_r), bins[j])
        r = np.sqrt(np.sum((p["x"])**2, axis=1))
        dm, _ = np.histogram(r, bins=r_bins, weights=mp*np.ones(r.shape))
        
        m = np.zeros(len(r_bins))
        m[0] = mp*np.sum(r < min_r)
        m[1:] = np.cumsum(dm) + m[0]
        
        ok = m > 0
        dlnr = np.log10(r_bins[1]) - np.log10(r_bins[0])
        dlnm_dlnr = np.gradient(np.log10(m[ok]), dlnr, edge_order=2)
        
        rho_enc = m / (r_bins**3 * 4*np.pi/3)

        ax[0].plot(r_bins[ok], m[ok], c=colors[j])
        ax[1].plot(m[ok], rho_enc[ok], c=colors[j])
        ax[2].plot(r_bins[ok], dlnm_dlnr, c=colors[j])

    ax[0].set_xscale("log")
    ax[0].set_yscale("log")
    ax[1].set_xscale("log")
    ax[1].set_yscale("log")
    ax[2].set_xscale("log")
    ax[0].set_xlabel(r"$r\ ({\rm kpc})$")
    ax[0].set_ylabel(r"$M(<r)\ M_\odot$")
    ax[1].set_xlabel(r"$M(<r)\ M_\odot$")
    ax[1].set_ylabel(r"$\rho(<r)\ M_\odot ({\rm kpc})^{-3}$")
    ax[2].set_xlabel(r"$r\ ({\rm kpc})$")
    ax[2].set_ylabel(r"$d\ln M/d\ln r$")
    
    fig.savefig("../plots/mass_profile_test.png")
        
if __name__ == "__main__": main()
