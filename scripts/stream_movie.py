import matplotlib.pyplot as plt
import symlib
import matplotlib as mpl
import cache_stars
import numpy as np

try:
    import palette
    palette.configure(True)
except:
    pass

def plot_host(i_host):
    base_dir = "/sdf/home/p/phil1/ZoomIns"
    suite = "SymphonyMilkyWayHR"

    sim_dir = symlib.get_host_directory(base_dir, suite, i_host)
    param = symlib.simulation_parameters(sim_dir)
    mp = param["mp"]/param["h100"]

    m_min = mp*5e3

    rs, hist = symlib.read_rockstar(sim_dir)
    sf, hist = symlib.read_symfind(sim_dir)
    lim = rs["rvir"][0,-1]

    fig, ax = plt.subplots(1, 3, figsize=(24.5, 8), sharey=True)
    fig.subplots_adjust(wspace=0.05, hspace=0.05)

    for icol in range(3):
        if icol == 0: ax[icol].set_ylabel(r"$Y\ ({\rm kpc})$")
        ax[icol].set_xlabel(r"$X\ ({\rm kpc})$")

        symlib.plot_circle(
            ax[icol], rs["x"][0,-1,0], rs["x"][0,-1,1],
            rs["rvir"][0,-1], c="white"
        )
        ax[icol].set_xlim(-lim, lim)
        ax[icol].set_ylim(-lim, lim)
            

    part = symlib.Particles(sim_dir)
    p = part.read(235, mode="smooth")
    mp_star, Fe_H, m_star, r_half = cache_stars.read_stars(suite, i_host)

    mp_dm = [np.ones(len(x))*mp if x is not None else None for x in mp_star]

    resolved = hist["mpeak"] > m_min
    surviving = sf["ok"][:,-1]
    low_mass = hist["merger_ratio"] < 0.15
    stripped = (sf["m"][:,-1]/hist["mpeak"] < 1/15.0)
    sub_flags = [
        resolved | (~resolved),
        resolved & low_mass,
        resolved & low_mass & ((~surviving) | stripped)
     ]

    part = symlib.Particles(sim_dir)

    norm = mpl.colors.LogNorm(vmin=1, vmax=1e5)
    kwargs = {"extent": [-lim, lim, -lim, lim],
              "norm": norm, "cmap": "inferno", "gridsize": 300}

    for icol in range(3):
        x, y, m = [], [], []
        for i_sub in range(1, len(sf)):
            if not sub_flags[icol][i_sub]: continue

            ps = p[i_sub]
            x.append(ps["x"][:,0])
            y.append(ps["x"][:,1])
            m.append(mp_star[i_sub])

        x, y, m = np.hstack(x), np.hstack(y), np.hstack(m)
        print("%g" % np.sum(m))
        poly_coll = ax[icol].hexbin(x, y, m, **kwargs)
    print()

    #plt.colorbar(poly_coll)
    fig.savefig("../plots/stellar_streams/stellar_streams_%s_%d.png" %
                (suite, i_host))

def main():
    for i_host in range(symlib.n_hosts("SymphonyMilkyWayHR")):
        #if i_host > 4: continue
        plot_host(i_host)

if __name__ == "__main__": main()
