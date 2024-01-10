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

base_dir = "/sdf/home/p/phil1/ZoomIns"
suite = "SymphonyMilkyWayHR"
# Change this to the place where you want to save your plots.
plot_dir = "../plots/stellar_streams/"

def plot_host(i_host):
    fig, ax = plt.subplots(figsize=(10, 8))
    fig.subplots_adjust(wspace=0.05, hspace=0.05)

    # Directory where subhalos and particles are stored.
    sim_dir = symlib.get_host_directory(base_dir, suite, i_host)

    # Read in Rockstar subhalos to get virial radii.
    rs, hist = symlib.read_rockstar(sim_dir)
    lim = rs["rvir"][0,-1]            

    # Read in particles. You get p, a list of structured arrays representing
    # the "smoothly" accreted matter around each subhalo.
    part = symlib.Particles(sim_dir)
    p = part.read(235, mode="smooth")
    # Read in star data. Same format for the first two arrays, second two are
    # subhalo-by-subhalo. You can generate this directly for a specific
    # galaxy-halo model, but here I've cached UniverseMachine predictions
    # combined with several observed galaxy size and metallicity relations.
    mp_star, Fe_H, m_star, r_half = cache_stars.read_stars(suite, i_host)

    # Don't include the first subhalo's particles, since they contain no stars
    x = np.hstack([sub["x"][:,0] for sub in p[1:]])
    y = np.hstack([sub["x"][:,1] for sub in p[1:]])
    m = np.hstack(mp_star[1:])

    # Everything from this point onwards is matplotlib shenanigans.

    n_grid = 400
    cell_area = 4*lim**2 / n_grid**2
    norm = mpl.colors.LogNorm(vmin=0.3, vmax=1e4)
    kwargs = {"extent": [-lim, lim, -lim, lim],
              "norm": norm, "cmap": "inferno", "gridsize": n_grid}
    poly_coll = ax.hexbin(x, y, m/cell_area, **kwargs)

    ax.set_ylabel(r"$Y\ ({\rm kpc})$")
    ax.set_xlabel(r"$X\ ({\rm kpc})$")

    symlib.plot_circle(
        ax, rs["x"][0,-1,0], rs["x"][0,-1,1],
        rs["rvir"][0,-1], c="white"
    )
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_title(r"${\rm Stellar\ surface\ density}\ " +
                 "(M_\odot\,kpc^{-2})$")
    fig.colorbar(poly_coll)
    fig.savefig("%s/stellar_streams_%s_%d.png" %
                (plot_dir, suite, i_host))

def main():
    n_hosts = symlib.n_hosts(suite)
    for i_host in range(n_hosts):
        print("%d/%d" % (i_host+1, n_hosts))
        plot_host(i_host)

if __name__ == "__main__": main()
