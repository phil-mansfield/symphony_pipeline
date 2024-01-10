import numpy as np
import symlib
import matplotlib.pyplot as plt
import palette
from palette import pc

palette.configure(False)
base_dir = "/sdf/home/p/phil1/ZoomIns"

def main():
    suite = "SymphonyMilkyWay"
    n_host = symlib.n_hosts(suite)

    fig1, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()
    fig3, ax3 = plt.subplots()
    fig4, ax4 = plt.subplots()

    for i_host in range(n_host):
        print(i_host)
        sim_dir = symlib.get_host_directory(base_dir, suite, i_host)

        gal, _  = symlib.read_galaxies(sim_dir)
        sf, hist = symlib.read_symfind(sim_dir)
        
        ok = sf["ok"][:,-1]
        ax1.plot(gal["v_disp_3d_star"][ok,-1],
                 gal["v_disp_3d_dm"][ok,-1], ".",
                 alpha=0.2, c=pc("r"))
        ax2.plot(gal["v_disp_3d_star"][ok,-1],
                 sf["vmax"][ok,-1], ".",
                 alpha=0.2, c=pc("r"))
        for snap in range(sf.shape[1]):
            if snap != 235: continue
            ok = sf["ok"][:,snap] & (gal["vmax_dm"][:,snap] > 0)
            ax3.plot(sf["m"][ok,snap]/hist["mpeak"][ok],
                     gal["vmax_dm_debias"][ok,snap]/gal["vmax_dm"][ok,snap],
                     ".", alpha=0.1, c=pc("r"))
            ax4.plot(sf["vmax"][ok,snap],
                     gal["vmax_dm_debias"][ok,snap]/gal["vmax_dm"][ok,snap],
                     ".", alpha=0.1, c=pc("r"))

    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlabel(r"$\sigma_{\rm 3d,\star}\ ({\rm km\,s^{-1}})$")
    ax1.set_ylabel(r"$\sigma_{\rm 3d,dm}\ ({\rm km\,s^{-1}})$")
    xlo, xhi = ax1.get_xlim()
    ax1.set_xlim(xlo, xhi)
    ax1.plot([xlo, xhi], [xlo, xhi], "--", c="k", lw=1.5)
    fig1.savefig("../plots/stellar_halo/v_disp_star_dm.png")

    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlabel(r"$\sigma_{\rm 3d,\star}\ ({\rm km\,s^{-1}})$")
    ax2.set_ylabel(r"$v_{\rm max}\ (\rm km\,s^{-1})$")
    xlo, xhi = ax2.get_xlim()
    ax2.set_xlim(xlo, xhi)
    ax2.plot([xlo, xhi], [xlo, xhi], "--", c="k", lw=1.5)
    fig2.savefig("../plots/stellar_halo/v_disp_vmax.png")

    ax3.set_xscale("log")
    ax3.set_xlabel(r"$m/m_{\rm peak}$")
    ax3.set_ylabel(r"$v_{\rm max,debias}/v_{\rm max}$")
    ax3.set_ylim(0.95, 1.25)
    ax3.set_xlim(1e-2, 1)
    fig3.savefig("../plots/stellar_halo/vratio_bias.png")

    ax4.set_xscale("log")
    ax4.set_xlabel(r"$v_{\rm max}$")
    ax4.set_ylabel(r"$v_{\rm max,debias}/v_{\rm max}$")
    ax4.set_ylim(0.95, 1.25)
    ax4.set_xlim(8, 150)
    fig4.savefig("../plots/stellar_halo/vmax_debias.png")


if __name__ == "__main__": main()
