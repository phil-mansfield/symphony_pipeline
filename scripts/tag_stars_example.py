import symlib
import numpy as np
import scipy.stats as stats
from colossus.cosmology import cosmology
import matplotlib.pyplot as plt
import palette
from palette import pc

def main():
    palette.configure(False)

    # Get location of the halo data.
    base_dir = "/sdf/home/p/phil1/ZoomIns"
    suite, i_host = "SymphonyMilkyWayHR", 0
    sim_dir = symlib.get_host_directory(base_dir, suite, i_host)

    # You need to specify the galaxy-halo model that you're using. Two fiducial
    # dwarf-scale models are provided, DWARF_GALAXY_HALO_MODEL, which uses
    # UniverseMachine, and DWARF_GALAXY_HALO_MODEL_NO_UM, which uses fits to
    # UniverseMachine. Use the former when you can and the latter for
    # simulations without UniverseMachine data.
    gal_halo = symlib.DWARF_GALAXY_HALO_MODEL_NO_UM

    # Call this function to get star properties. stars is a list of star
    # properites, with each element giving stars for a different subhalo. gals
    # gives information on each galaxy, and ranks stores a bunch of intermediate
    # tagging data (referred to as ranks) whcih is expensive to compute and can
    # be used if re-tagging stars.
    stars, gals, ranks = symlib.tag_stars(sim_dir, gal_halo)
    
    # Each element of the stars list is an array of type symlib.STAR_DTYPE. It
    # has the following fields:
    #   "mp" - the mass of the star particle in Msun
    #   "Fe_H" - the [Fe/H] of the star particle
    #   "a_form" - the scale factor at which this star was formed

    # gals is an array of type symlib.GALAXY_HISTORY_DTYPE. It has the following
    # fields:
    #   "m_star_i" - the initial mass of the galaxy in Msun
    #   "r_half_2d_i" - the initial 2d half-mass radius in kpc
    #   "r_half_3d_i" - the initial 3d half-mass radius in kpc
    #   "Fe_H_i" - the intial mean metallicity, [Fe/H]
    #   "sigma_Fe_H_i" - the initial sqrt(variance) of the MDF
    #   "a50" - the scale factor at which 50% of the stellar mass was formed
    #   "a90" - the scale factor at which 90% of the stellar mass was formed

    # You probably don't want to do anything with ranks directly ¯\_(ツ)_/¯

    # tag_stars only return the unchanging properties of stars, so you'll need
    # separately load positions and velocities.

    # symlib.Particles is a class that manages reading in simulation data. When
    # you intialize it, it loads and analyzes several header files.
    part = symlib.Particles(sim_dir)
    # The read method reads particle properties at a given snapshot. The mode
    # argument determines which particles get read in and how many times they
    # are replicated across subhalos. Setting it to "stars" reads in the
    # particles that get tagged as stars. (For veterans, this is the same as
    # as the "smooth" mode). The formatting of p is the same as the formatting
    # of the stars list, so p[i][j] is the jth particle in the ith halo and
    # stars[i][j] is the same particle. When read in through the "stars" mode,
    # p[i][j] will be the same particle regardless of snapshot.
    p = part.read(235, mode="stars")

    # Each element of p is an array of type symlib.PARTICLE_DTYPE. It has the
    # following fields:
    #   "x" - the 3d position of the particle in (p)kpc relative to the host
    #   "v" - the 3d velocity of the paritcle in km/s relative to the host
    #   "id" - A unique integer ID for the particle
    #   "snap" - the snapshot when the particle first fell into a halo. This is
    #            also the first snapshot where the particle is tracked.
    #   "ok" - true if the particle is being tracked and false otherwise.

    # Beyond this point, I make the density and [Fe/H] plots. You don't
    # neccessarily need to read through this.

    p, stars = np.hstack(p[1:]), np.hstack(stars[1:])

    rs, hist = symlib.read_rockstar(sim_dir)
    fig_density, fig_Fe_H = plot_images(rs, p, stars)

    fig_density.savefig("../plots/stellar_halo/rho_image_%s_%d.png" %
                        (suite, i_host))
    fig_Fe_H.savefig("../plots/stellar_halo/Fe_H_image_%s_%d.png" %
                     (suite, i_host))

def plot_images(rs, p, stars):
    lim = rs[0,-1]["rvir"]
    lim = 80
    figsize=(8, 8)

    ok = (p["x"][:,2] < lim) & (p["x"][:,2] > -lim)

    fig1, ax = plt.subplots(figsize=figsize)
    r_bins = np.linspace(-lim, lim, 200)
    rho, _, _ = np.histogram2d(p["x"][ok,0], p["x"][ok,1], bins=r_bins,
                               weights=stars["mp"][ok])
    area = (r_bins[1] - r_bins[0])**2
    im = ax.imshow(np.log10(rho.T / area), extent=[-lim, +lim, -lim, +lim],
                   cmap="inferno", origin="lower", vmin=0, vmax=8)
    fig1.suptitle(r"$\log_{10}(\rho)\,(M_\odot\,{\rm kpc}^2)$")
    ax.set_xlabel(r"$X\,({\rm kpc})$")
    ax.set_ylabel(r"$Y\,({\rm kpc})$")
    plt.colorbar(im)

    fig2, ax = plt.subplots(figsize=figsize)
    r_bins = np.linspace(-lim, lim, 100)
    rho, _, _ = np.histogram2d(p["x"][ok,0], p["x"][ok,1], bins=r_bins,
                               weights=stars["mp"][ok])
    rho_Fe_H, _, _ = np.histogram2d(p["x"][ok,0], p["x"][ok,1], bins=r_bins,
                                    weights=stars["mp"][ok]*stars["Fe_H"][ok])
    mean, std = np.mean(stars["Fe_H"]), np.std(stars["Fe_H"])
    Fe_H = rho_Fe_H/rho
    mean = np.median(Fe_H)
    im = ax.imshow(Fe_H, extent=[-lim, +lim, -lim, +lim],
                   cmap="RdBu", origin="lower",
                   vmin=-2.5, vmax=-0.5)
    fig2.suptitle(r"${\rm [Fe/H]}$")
    ax.set_xlabel(r"$X\,({\rm kpc})$")
    ax.set_ylabel(r"$Y\,({\rm kpc})$")
    plt.colorbar(im)

    return fig1, fig2

if __name__ == "__main__": main()
