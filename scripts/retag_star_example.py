import numpy as np
import matplotlib.pyplot as plt
import symlib
import time

try:
    import palette
    palette.configure(True)
except:
    pass

base_dir = "/sdf/group/kipac/g/cosmo/ki21/phil1/simulations/ZoomIns"
gal_halo = symlib.GalaxyHaloModel(
    symlib.UniverseMachineMStarFit(),
    symlib.Jiang2019RHalf(),
    symlib.PlummerProfile(),
    symlib.Kirby2013Metallicity(),
    no_scatter=False
)

def compute_profile(mp, p, r_bins=10**np.linspace(0, np.log10(250), 40)):
    n_bins = len(r_bins) - 1
    mass = np.zeros(n_bins, dtype=np.float64)
    
    for i in range(1, len(p)):
        r = np.sqrt(np.sum(p[i]["x"]**2, axis=1))
        mass += np.histogram(r, bins=r_bins, weights=mp[i])[0]

    vol = 4*np.pi/3*(r_bins[1:]**3 - r_bins[:-1]**3)
    r_mid = np.sqrt(r_bins[1:]*r_bins[:-1])

    return r_mid, mass/vol

def main():
    sim_dir = symlib.get_host_directory(base_dir, "SymphonyMilkyWay", 0)
    
    part = symlib.Particles(sim_dir)
    p = part.read(235, mode="smooth")

    t0 = -time.time()
    mp, ranks, _, _, _ = symlib.tag_stars(sim_dir, gal_halo)
    t0 += time.time()

    r, rho = compute_profile(mp, p)
    plt.plot(r, rho, c="tab:blue", lw=4)

    state, t1, profs = None, 0, []
    resamples = 17
    for i in range(resamples):
        t1 -= time.time()
        mp, _, _, _, state = symlib.retag_stars(
            sim_dir, gal_halo, ranks, state=state)
        t1 += time.time()
    
        r, rho = compute_profile(mp, p)
        plt.plot(r, rho, c="tab:red", lw=1.5, alpha=0.5)

    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel(r"$r\,({\rm kpc})$")
    plt.ylabel(r"$\rho\,{M_\odot\,{\rm kpc}^{-3}}$")

    print("%.2f %.2f" % (t0, t1/resamples))

    plt.savefig("../plots/stellar_halo/resample_test.png")

if __name__ == "__main__": main()
