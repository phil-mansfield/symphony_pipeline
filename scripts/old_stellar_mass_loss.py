import numpy as np
import matplotlib.pyplot as plt
import palette
from palette import pc
import symlib

palette.configure(False)

def main():
    fig1, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()
    fig3, ax3 = plt.subplots()
    fig4, ax4 = plt.subplots()

    suites = ["SymphonyLMC", "SymphonyMilkyWay", "SymphonyGroup",
              "SymphonyLCluster"][::-1]
    colors = [pc("r"), pc("o"), pc("b"), pc("p")]

    for i_suite in range(len(suites)):
        suite, color = suites[i_suite], colors[i_suite]
        n_hosts = symlib.n_hosts(suite)
        for i_host in range(n_hosts):
            print(i_host)
            dir_name = symlib.get_host_directory(
                "/sdf/home/p/phil1/ZoomIns", suite, i_host)

            param = symlib.simulation_parameters(dir_name)
            mp = param["mp"]/param["h100"]

            try:
                gal, gal_hist = symlib.read_galaxies(dir_name)
            except:
                print("Error reading galaxy files for:\n%s" % dir_name)
                continue
            sf, hist = symlib.read_symfind(dir_name)

            sub_ok = (hist["merger_ratio"] < 0.15) & (hist["mpeak_pre"]/mp > 3e4)

            for i in range(1, len(sf)):
                snap = hist["first_infall_snap"][i]
                if gal["m_dyn"][i,snap] <= 0 or hist["mpeak_pre"][i]/mp < 3e3:
                    continue

                ok = gal["ok"][i] & (gal["m_star"][i] > 0)
                if sub_ok[i]:
                    ax1.plot(sf["m"][i,ok]/hist["mpeak_pre"][i],
                             gal["m_star"][i,ok]/gal_hist["m_star_i"][i],
                             ".", c="k", alpha=0.1)
                    ax2.plot(sf["m"][i,ok]/hist["mpeak_pre"][i],
                             gal["r_half"][i,ok]/gal_hist["r_half_i"][i],
                             ".", c="k", alpha=0.1)
                    ax3.plot(sf["m"][i,ok]/hist["mpeak_pre"][i],
                             (gal["m_star"][i,ok]/2/gal["m_dyn"][i,ok]) / 
                             (gal["m_star"][i,snap]/2/gal["m_dyn"][i,snap]),
                             ".", c="k", alpha=0.1)
                else:
                    if hist["merger_ratio"][i] < 0.15:
                        ax4.plot(gal_hist["m_star_i"][i],
                                 gal["m_star"][i,snap]/2/gal["m_dyn"][i,snap],
                                 ".", c=color)
                    else:
                        ax4.plot(gal_hist["m_star_i"][i],
                                 gal["m_star"][i,snap]/2/gal["m_dyn"][i,snap],
                                 "x", c=color)

    ax1.set_ylabel(r"$m_\star/m_{\rm \star,infall}$")
    ax1.set_xlabel(r"$m/m_{\rm infall}$")
    ax2.set_ylabel(r"$r_{50}/m_{\rm 50,infall}$")
    ax2.set_xlabel(r"$m/m_{\rm infall}$")
    ax3.set_ylabel(r"$m_\star(<r_{50})/m_{\rm dm}(<r_{50})$")
    ax3.set_xlabel(r"$m/m_{\rm infall}$")
    ax4.set_xlabel(r"$M_\star\,(M_\odot)$")
    ax4.set_ylabel(r"$g_\star(<r_{50})/g_{\rm dm}(< r_{50})$")

    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax3.set_xscale("log")
    ax3.set_yscale("log")
    ax4.set_xscale("log")
    ax4.set_yscale("log")

    ax1.set_ylim(1e-2, 2)
    ax2.set_ylim(0.3, 10)
    ax1.set_xlim(1e-3, 1)
    ax2.set_xlim(1e-3, 1)
    ax3.set_xlim(1e-3, 1)
    ax3.set_ylim(3e-4, 1)
    ax4.set_xlim(1e2, 1e12)
    ax4.set_ylim(1e-6, None)

    for i_suite in range(len(suites)):
        ax4.plot([], [], ".", c=colors[i_suite], label=suites[i_suite])
    ax4.legend(loc="upper left", fontsize=17)
    
    xlo, xhi = ax4.get_xlim()
    ax4.plot([xlo, xhi], [1, 1], "--", c="k")

    fig1.savefig("../plots/stellar_halo/fstar_fdm.png")
    fig2.savefig("../plots/stellar_halo/r50_fdm.png")
    fig3.savefig("../plots/stellar_halo/mratio_fdm.png")
    fig4.savefig("../plots/stellar_halo/mratio_mstar.png")

if __name__ == "__main__": main()
