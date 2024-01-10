import symlib
import numpy as np
import matplotlib.pyplot as plt
import palette
from palette import pc
import os.path as path

SUITE = "SymphonyMilkyWay"
OUT_DIR = "../plots/um"
BASE_DIR = "/sdf/home/p/phil1/ZoomIns"

colors = [pc("r"), pc("o"), pc("b"), pc("g"), pc("p")]
colors_ext = [pc("r", 0.35), pc("o", 0.35), pc("b", 0.35), pc("g", 0.35), pc("p", 0.35), pc("r", 0.65), pc("o", 0.65), pc("b", 0.65), pc("g", 0.65), pc("p", 0.65)]

def comp_rs_um():
    n_host = symlib.n_hosts(SUITE)

    fig, ax = plt.subplots()

    for i_host in range(n_host):
        sim_dir = symlib.get_host_directory(BASE_DIR, SUITE, i_host)

        h, hist = symlib.read_subhalos(sim_dir)
        um = symlib.read_um(sim_dir)

        a = symlib.scale_factors(sim_dir)
                
        lim_1 = (1, len(colors_ext)+1)
        lim_2 = (45, len(colors_ext)+45)

        for k, lim in enumerate([lim_1, lim_2]):
            ax.cla()

            low, high = lim
            var_h = "mvir"
            var_um = "mvir"
            
            ok_h = h["ok"][0]
            ok_um = um["ok"][0]
            no_um = um["ok"][0] & (~um["is_orphan"][0])
            ax.plot(a[ok_h], h[var_h][0,ok_h], c="k")
            ax.plot(a[ok_um], um[var_um][0,ok_um], "--", lw=1.5, c="k")
            ax.plot(a[ok_um], um[var_um][0,no_um], "--", lw=1.5, c="k")

            for i in range(high - low):
                j = low + i
                ok_h = h["ok"][j]
                ok_um = um["ok"][j]
                no_um = um["ok"][j] & (~um["is_orphan"][j])
                c = colors_ext[i]
                ax.plot(a[ok_h], h[var_h][j,ok_h], c=c)
                ax.plot(a[ok_um], um[var_um][j,ok_um], "--", lw=1.5, c=c)
                ax.plot(a[no_um], um[var_um][j,no_um], "--", lw=1.5, c=c)

            ax.set_yscale("log")
            ax.set_xscale("log")
            ax.set_ylabel(r"$M_{\rm vir}$")
            ax.set_xlabel(r"$a(t)$")

            fig.savefig(path.join(OUT_DIR, "comp_rs_um_%d_%02d.png" %
                                  (k, i_host)))
        break


def main():
    palette.configure(False)

    comp_rs_um()

if __name__ == "__main__": main()
