import numpy as np
import matplotlib.pyplot as plt
import symlib
from palette import pc
import palette
import mock_disruption
import survival
from colossus.cosmology import cosmology
from colossus.halo import mass_so

SUITE_COMPARE = False

def numerical_limit(n_peak):
    ln = np.log10(n_peak)
    return 10**(1.6587 + 0.3861*ln - 0.01853*ln**2)/n_peak

def main():
    palette.configure(True)

    base_dir = "/sdf/home/p/phil1/ZoomIns"
    if not SUITE_COMPARE:
        suites = ["SymphonyMilkyWayHR", "SymphonyMilkyWayLR"]
        colors = [pc("r"), pc("b")]
        labels = [r"${\rm high{-}res}\ (n_{\rm peak}\times 8)$",
                  r"${\rm fiducial}$"]
    
    else:
        suites = ["SymphonyLMC", "SymphonyMilkyWay", "SymphonyGroup",
                  "SymphonyLCluster"]
        colors = [pc("r"), pc("o"), pc("b"), pc("p"), pc("k")]
        labels = [
            r"${\rm LMC}$", r"${\rm MilkyWay}$", r"${\rm Group}$",
            r"${\rm L{-}Cluster}$",
        ]

    bin_edges = [
        (10**2.5, 10**3),
        (10**3, 10**3.5),
        (10**3.5, 10**4),
        (10**4, 10**4.5),
        (10**4.5, 10**5),
    ]
    titles = [
        r"$10^{2.5} < n_{\rm peak} < 10^{3}$",
        r"$10^{3} < n_{\rm peak} < 10^{3.5}$",
        r"$10^{3.5} < n_{\rm peak} < 10^{4}$",
        r"$10^{4} < n_{\rm peak} < 10^{4.5}$",
        r"$10^{4.5} < n_{\rm peak} < 10^{5}$",
    ]

    old_medians = [None]*len(bin_edges)

    for i_suite in range(len(suites)):
        suite = suites[i_suite]
        ts = [[] for _ in range(len(bin_edges))]
        ms = [[] for _ in range(len(bin_edges))]
        npeaks = [[] for _ in range(len(bin_edges))]
        infall_masses = [[] for _  in range(len(bin_edges))]

        n_hosts = symlib.n_hosts(suite)
        param = symlib.simulation_parameters(suite)
        mp = param["mp"]/param["h100"]
        # All this does is make sure that HR and LR are put on the same axes
        # Gets corrected for later on.
        if suite == "SymphonyMilkyWayHR" and len(suites) == 2:
            mp *= 8
        cosmo = cosmology.setCosmology('', symlib.colossus_parameters(param))

        for i_host in range(n_hosts):
            if i_host == 3 and suite in ["SymphonyMilkyWayHR", "SymphonyMilkyWayLR"]:
                continue
            print(suite, i_host)
            sim_dir = symlib.get_host_directory(base_dir, suite, i_host)
            h, hist = symlib.read_subhalos(sim_dir)
            c = symlib.read_cores(sim_dir)
            scale = symlib.scale_factors(sim_dir)
            z = 1/scale - 1
            age = cosmo.age(z)
            t_dyn = mass_so.dynamicalTime(z, "vir", "crossing")

            for j, (bin_low, bin_high) in enumerate(bin_edges):
                npeak = hist["mpeak_pre"]/mp
                idx = np.where((npeak < bin_high) & (npeak > bin_low))[0]
                for i in idx:
                    ok = c["ok"][i]
                    if np.sum(ok) < 5 or hist["merger_ratio"][i] > 0.1:
                        continue
                    snap = hist["first_infall_snap"][i]
                    ms[j].append(c["m_bound"][i, ok]/hist["mpeak_pre"][i])
                    ms[j][-1] = np.minimum(ms[j][-1], 1)
                    ts[j].append((age[ok] - age[snap])/t_dyn[snap])
                    npeaks[j].append(npeak[i])
                    infall_masses[j].append(c["m_bound"][i,snap]/hist["mpeak_pre"][i])

        for j in range(len(ms)):
            plt.figure(j)
            t_eval = np.linspace(0, 30, 100)
            med, err_low, err_high, _, _, ok = mock_disruption.kaplan_meier_mass_loss(ts[j], ms[j], t_eval)
            im = np.median(infall_masses[j])
            med[0], err_low[0], err_high[0] = im, im, im
            print(im)
            t_eval, med = t_eval[ok], med[ok]
            err_low, err_high = err_low[ok], err_high[ok]

            n_med = np.median(npeaks[j])
            if suite == "SymphonyMilkyWayHR": n_med *= 8
            print(n_med)
            n_ratio_lim = numerical_limit(n_med)
            n8_ratio_lim = numerical_limit(8*n_med)
            print(n_ratio_lim, n8_ratio_lim)

            c = colors[i_suite]
            conv = med > n_ratio_lim
            #plt.plot(t_eval[conv], med[conv], c=c,
            #         label=labels[i_suite])
            if np.sum(med > n_ratio_lim) == 0:
                ratio_conv, t_conv = 0, 0
            else:
                ratio_conv, t_conv = np.min(med[conv]), np.max(t_eval[conv])
            conv8 = med > n8_ratio_lim
            if np.sum(med > n8_ratio_lim) == 0:
                ratio8_conv, t8_conv = 0, 0
            else:
                ratio8_conv, t8_conv = np.min(med[conv8]),np.max(t_eval[conv8])
            #plt.plot([t_conv], [ratio_conv], "*", ms=15, color=c)
            #plt.plot([t8_conv], [ratio8_conv], "o", ms=10, color=c)
            plt.plot(t_eval, med, c=c, label=labels[i_suite])
            plt.fill_between(t_eval, err_low, err_high, color=c, alpha=0.2)

            if suite == "SymphonyMilkyWayLR":
                plt.plot([t_conv], [2*ratio_conv], marker=r"$\downarrow$",
                         color=pc("b"), markersize=20)
                plt.plot([t8_conv], [2*ratio8_conv], marker=r"$\downarrow$",
                         color=pc("r"), markersize=20)

    for j in range(len(bin_edges)):
        plt.figure(j)
        if SUITE_COMPARE:
            plt.legend(loc="lower left")
        else:
            #plt.plot([], [], marker=r"$\downarrow$", c=pc("b"), linestyle=None,
            #         label=r"$t_{\rm conv,ideal,fid}$")
            #plt.plot([], [], marker=r"$\downarrow$", c=pc("r"), linestyle=None,
            #         label=r"$t_{\rm conv,ideal,HR}$")
            if j == 0:
                plt.legend(loc="lower left", fontsize=22)
        plt.yscale("log")
        plt.ylim(3e-4, 3)
        plt.xlim(-0.5, 15)
        if not SUITE_COMPARE:
            plt.title(titles[j], fontsize=37)
            plt.xlabel(r"$(t-t_{\rm infall})/t_{\rm cross}$", fontsize=37)
            plt.ylabel(r"$m/m_{\rm peak}$", fontsize=37)
        else:
            plt.title(titles[j])
            plt.xlabel(r"$(t-t_{\rm infall})/t_{\rm cross}$")
            plt.ylabel(r"$m/m_{\rm peak}$")
        if SUITE_COMPARE:
            print("../plots/core_tracking/mass_loss_conv_comp_%d.pdf" % j)
            plt.savefig("../plots/core_tracking/mass_loss_conv_comp_%d.pdf" % j)
        else:
            print("../plots/core_tracking/mass_loss_conv_%d.pdf" % j)
            plt.savefig("../plots/core_tracking/mass_loss_conv_%d.pdf" % j)

if __name__ == "__main__": main()
