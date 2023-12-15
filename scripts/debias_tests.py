import numpy as np
import symlib
import palette
from palette import pc
import matplotlib.pyplot as plt

palette.configure(False)

base_dir = "/sdf/home/p/phil1/ZoomIns"
suite = "SymphonyMilkyWay"
out_dir = "../plots/debias_model"

v_bins = 10**np.linspace(0.5, 2.5, 100)

def f_wrong(sf, gal, v_bins=v_bins):
    n_all, _ = np.histogram(sf["vmax"], bins=v_bins)
    n_wrong, _ = np.histogram(sf["vmax"][~gal["m23_m_conv"]], bins=v_bins)
    n_all = np.maximum(n_all, 1)

    n_all = np.cumsum(n_all[::-1])[::-1]
    n_wrong = np.cumsum(n_wrong[::-1])[::-1]

    print(np.sum(gal["m23_m_conv"]), len(gal))

    return n_wrong/n_all

def main():
    sim_dir = symlib.get_host_directory(base_dir, suite, 0)
    print(sim_dir)

    param = symlib.simulation_parameters(sim_dir)
    mp = param["mp"]/param["h100"]

    sf, hist = symlib.read_symfind(sim_dir)
    gal, gal_hist = symlib.read_galaxies(sim_dir)
    sf, gal = sf[:,-1], gal[:,-1]

    r = np.sqrt(np.sum(sf["x"]**2, axis=1))
    hi_r, lo_r = (r < 250) & sf["ok"], (r < 50) & sf["ok"]

    #hi_n, hi_edges = np.histogram(np.log10(sf["vmax"][hi_r]))
    
    npeak = hist["mpeak"]/mp
    m_ratio = sf["m"]/hist["mpeak"]
    print(npeak[sf["ok"]])
    print(m_ratio[sf["ok"]])
    print(gal["m23_m_conv"][sf["ok"]])
    
    f_wrong_hi = f_wrong(sf[hi_r], gal[hi_r])
    f_wrong_lo = f_wrong(sf[lo_r], gal[lo_r])

    print(f_wrong_hi)
    print(f_wrong_lo)

    fig, ax = plt.subplots()
    ax.plot(v_bins[1:], f_wrong_hi, "--", c=pc("r"))
    ax.plot(v_bins[1:], f_wrong_lo, "--", c=pc("b"))

    ax.set_xscale("log")

    #plt.plot(hi_edges[1:], np.cumsum(hi_n[::-1])[::-1], pc("r"))
    fig.savefig("%s/f_wrong.png" % out_dir)

if __name__ == "__main__": main()
