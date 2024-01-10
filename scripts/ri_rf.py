import numpy as np
import matplotlib.pyplot as plt
import palette
from palette import pc
import symlib
import matplotlib as mpl

base_dir = "/sdf/group/kipac/g/cosmo/ki21/phil1/simulations/ZoomIns"

palette.configure(False)

def main():
    sim_dir = symlib.get_host_directory(base_dir, "SymphonyCluster", 0)

    scale = symlib.scale_factors(sim_dir)
    snap_i = np.searchsorted(scale, 0.5)
    snap_f = len(scale) - 1
    print(snap_i, snap_f)

    part = symlib.Particles(sim_dir)
    p_i = part.read(snap_i, mode="all")[0]
    p_f = part.read(snap_f, mode="all")[0]

    ok = p_i["ok"] & p_f["ok"]
    p_i, p_f = p_i[ok], p_f[ok]

    r_i = np.sqrt(np.sum(p_i["x"]**2, axis=1))
    r_f = np.sqrt(np.sum(p_f["x"]**2, axis=1))

    norm = mpl.colors.LogNorm(vmin=1, vmax=10000)
    low, high = np.log10(1e0), np.log10(1e4)
    plt.plot([low, low], [high, high], "--", lw=2, c="w")
    kwargs = {"extent": [low, high, low, high],
              "norm": norm, "cmap": "inferno", "gridsize": 200}
    plt.hexbin(np.log10(r_i), np.log10(r_f), **kwargs)
    plt.xlabel(r"$\log_{10}\,r(a = 0.5)/1\,{\rm kpc}$")
    plt.ylabel(r"$\log_{10}\,r(a = 1.0)/1\,{\rm kpc}$")

    plt.savefig("../plots/ri_rf.png")

if __name__ == "__main__": main()
