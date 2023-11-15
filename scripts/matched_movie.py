import numpy as np
import matplotlib.pyplot as plt
import symlib
import palette
from palette import pc
import os.path as path

COMOVING_MATCH = True
PLOT_SET = np.arange(0, 41, dtype=int)

def main():
    global PLOT_SET

    palette.configure(False)
    base_dir = "/sdf/group/kipac/g/cosmo/ki21/phil1/simulations/ZoomIns"
    suite_1, suite_2 = "SymphonyMilkyWay", "SymphonyMilkyWayDiskDMO"
    name_1, name_2 = r"${\rm 8K}$", r"${\rm 4K}$"
    host_i = 13
    out_dir = "../plots/disk_tests/match_h%d" % host_i

    sim_dir_1 = symlib.get_host_directory(base_dir, suite_1, host_i)
    sim_dir_2 = symlib.get_host_directory(base_dir, suite_2, host_i)
    
    h_1, hist_1 = symlib.read_rockstar(sim_dir_1)
    h_2, hist_2 = symlib.read_rockstar(sim_dir_2)

    n_min = min(len(h_1), len(h_2))
    PLOT_SET = PLOT_SET[:n_min]

    print(np.sum(h_1["ok"][:,-1]), np.sum(h_2["ok"][:,-1]))

    if COMOVING_MATCH:
        h_cmov_1, hist_cmov_1 = symlib.read_rockstar(sim_dir_1, comoving=True)
        h_cmov_2, hist_cmov_2 = symlib.read_rockstar(sim_dir_2, comoving=True)
        match_1, match_2 = symlib.match_subhalos(
            h_cmov_1, h_cmov_2, hist_cmov_1, hist_cmov_2)
    else:
        match_1, match_2 = symlib.match_subhalos(h_1, h_2, hist_1, hist_2)

    fig, ax = plt.subplots()
    r_max = h_1[0,-1]["rvir"]*1.25
    scale = symlib.scale_factors(sim_dir_1)
    for snap in range(h_1.shape[1]):
        if snap != h_1.shape[1] - 1: continue
        print(snap)
        ax.cla()
        ax.plot([], [], "-", c="tab:gray", label=name_1)
        ax.plot([], [], "--", c="tab:gray", label=name_2)
        ax.set_title("$a= %.4f$" % scale[snap])
        ax.legend(frameon=True)
        symlib.plot_matched_subhalos(ax, h_1, h_2,
                                     match_1, match_2,
                                     snap, r_max,
                                     plot_set_1=PLOT_SET, plot_set_2=PLOT_SET)
        fig.savefig(path.join(out_dir, "match_%03d.png" % snap))


if __name__ == "__main__": main()
