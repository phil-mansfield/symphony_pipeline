import symlib
import numpy as np
import scipy.stats as stats
from colossus.cosmology import cosmology
import matplotlib.pyplot as plt
import palette
from palette import pc
import time

def old_main():
    palette.configure(False)

    base_dir = "/sdf/home/p/phil1/ZoomIns"
    suite = "SymphonyMilkyWay"
    i_host = 0
    i_sub = 1

    sim_dir = symlib.get_host_directory(base_dir, suite, i_host)

    scale = symlib.scale_factors(sim_dir)
    rs, hist = symlib.read_rockstar(sim_dir)
    sf, _ = symlib.read_symfind(sim_dir)

    gal_halo = symlib.DWARF_GALAXY_HALO_MODEL

    stars, gal_hists, ranks = symlib.tag_stars(
        sim_dir, gal_halo, target_subs=np.array([i_sub]))

    print("m_{p,star}")
    print(stars[i_sub]["mp"][:30])
    print("[Fe/H]")
    print(stars[i_sub]["Fe_H"][:30])
    print("a_form")
    print(stars[i_sub]["a_form"][:30])
    print()
    print("history")
    print(gal_hists[i_sub])

    print(np.mean(stars[i_sub]["Fe_H"]))
    print(np.std(stars[i_sub]["Fe_H"]))
    print(stats.spearmanr(stars[i_sub]["Fe_H"], stars[i_sub]["a_form"])[0])

    param = symlib.simulation_parameters(suite)
    mp = param["mp"]/param["h100"]
    cosmo = cosmology.setCosmology('', symlib.colossus_parameters(param))

    z_form = 1/stars[i_sub]["a_form"] - 1
    t_form = cosmo.age(z_form)
    star_age = cosmo.age(1/scale[-1] - 1) - t_form

    plt.figure()
    plt.plot(star_age, stars[i_sub]["Fe_H"], ".", alpha=0.05, c="k")
    plt.xlabel(r"${\rm age\ (Gyr)}$")
    plt.ylabel(r"${\rm [Fe/H]}$")
    plt.xlim(9, 14)
    
    plt.savefig("../plots/stellar_halo/Fe_age_relation_%d_%d.png" %
                (i_host, i_sub))

    um = symlib.read_um(sim_dir)
    print(um[i_sub]["sfr"])

def main():
    palette.configure(False)

    base_dir = "/sdf/home/p/phil1/ZoomIns"
    suite = "SymphonyMilkyWayHR"
    i_host = 3

    sim_dir = symlib.get_host_directory(base_dir, suite, i_host)
    gal_halo = symlib.DWARF_GALAXY_HALO_MODEL_NO_UM
    stars, gal_hists, ranks = symlib.tag_stars(sim_dir, gal_halo)

    part = symlib.Particles(sim_dir)
    p = part.read(235, mode="stars")
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
    im = ax.imshow(Fe_H.T, extent=[-lim, +lim, -lim, +lim],
                   cmap="RdBu", origin="lower",
                   vmin=-2.5, vmax=-0.5)
    fig2.suptitle(r"${\rm [Fe/H]}$")
    ax.set_xlabel(r"$X\,({\rm kpc})$")
    ax.set_ylabel(r"$Y\,({\rm kpc})$")
    plt.colorbar(im)

    return fig1, fig2

if __name__ == "__main__": old_main()
