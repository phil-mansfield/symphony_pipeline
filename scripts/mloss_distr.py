import numpy as np
import matplotlib.pyplot as plt
import palette
from palette import pc
import symlib

palette.configure(True)

def main():
    base_dir = "/sdf/home/p/phil1/ZoomIns"
    suite = "SymphonyMilkyWay"
    
    for i_host in range(symlib.n_hosts(suite)):
        sim_dir = symlib.get_host_directory(base_dir, suite, i_host)
        param = symlib.simulation_parameters(sim_dir)
        mp = param["mp"]/param["h100"]

        sf, hist = symlib.read_symfind(sim_dir)

        ok = sf["ok"][:,-1]
        
        n.append(sf["m"][ok,-1]/mp)
        npeak.append(hist["mpeak_pre"][ok]/mp)

        #break

    n, npeak = np.hstack(n), np.hstack(npeak)
    
    #bins = np.linspace(-3, 0, 30)
    bins = np.linspace(0, 1, 30)
    plt.figure()

    bin1 = (10**3 < n) & (n < 10**3.5)
    #bin2 = (10**4 < n) & (n < 10**4.5)
    bin3 = (10**3 < npeak) & (npeak < 10**3.5)
    #bin4 = (10**4 < npeak) & (npeak < 10**4.5)
    #ratio = np.log10(n/npeak)
    ratio = n/npeak
    
    plt.hist(ratio[bin1], bins=bins, histtype="step", lw=3, color=pc("r"),
             density=True, label=r"$10^3 < n < 10^{3.5}$")
    #plt.hist(ratio[bin2], bins=bins, histtype="step", lw=3, color=pc("r"),
    #         ls="--", density=True)
    plt.hist(ratio[bin3], bins=bins, histtype="step", lw=3, color=pc("b"),
             density=True, label=r"$10^3 < n_{\rm peak} < 10^{3.5}$")
    #plt.hist(ratio[bin2], bins=bins, histtype="step", lw=3, color=pc("b"),
    #         ls="--", density=True)
    plt.xlabel(r"$\mu \equiv m/m_{\rm peak}$")
    plt.ylabel(r"${\rm Pr}(\mu)$")
    plt.legend(loc="upper right")

    plt.savefig("../plots/core_tracking/mloss_distr.pdf")

    
if __name__ == "__main__": main()
