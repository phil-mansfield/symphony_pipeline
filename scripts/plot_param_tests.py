import numpy as np
import matplotlib.pyplot as plt
import palette
from palette import pc
import os.path as path

plot_dir = "/home/users/phil1/code/src/github.com/phil-mansfield/symphony_pipeline/plots/core_plots"

def get_k(s):
    return int(s[1:3])
def get_n(s):
    return int(s[5:7])

def main():
    palette.configure(True)

    npeak = [3e2, 1e3, 1e4, 1e5]
    npeak_names = ["300", "10^{3}", "10^{4}", "10^{5}"]

    all_k = [64, 32, 16, 8]
    k_colors = [pc("r"), pc("o"), pc("b"), pc("p")]
    

    n_rows = 19
    rows = np.loadtxt("tables/param_tests.txt",
                      usecols=(1, 2, 3, 4, 5, 6), dtype=int)
    names = np.loadtxt("tables/param_tests.txt", usecols=(0,), dtype=str)

    for i_block in range(4):
        if i_block != 0 and i_block != 2: continue

        names_i = names[i_block*n_rows:(i_block+1)*n_rows]
        cols_i = rows[i_block*n_rows:(i_block+1)*n_rows].T
        n_either, n_both, n_err_c, n_err_rs, n_miss_rs, n_miss_c = cols_i
        
        f_err_c = (n_miss_c + n_err_c) / (n_both + n_miss_c)
        f_err_rs = (n_miss_rs + n_err_rs) / (n_both + n_miss_rs)

        k = np.fromiter(map(get_k, names_i), int)
        n = np.fromiter(map(get_n, names_i), int)

        fig, ax = plt.subplots()
        ax.set_title(r"$n_{\rm peak} > %s$" % npeak_names[i_block])

        for i, k_i in enumerate(all_k):
            idx = np.where(k == k_i)[0]
            ax.plot(n[idx], f_err_rs[idx],
                     c=k_colors[i], label=r"$k=%d$"%k_i)
            ax.plot(n[idx], f_err_c[idx], "--", c=k_colors[i])

            if k_i == 16:
                ax.plot([n[idx[0]]], [f_err_rs[idx[0]]], "*", c=pc("b"), ms=20)
                ax.plot([n[idx[0]]], [f_err_c[idx[0]]], "*", c=pc("b"), ms=20)

        ax.plot([], [], "-", c=pc("a"), label=r"${\rm Rockstar\ errors}$")
        ax.plot([], [], "--", c=pc("a"), label=r"${\rm Particle{-}tracking\ errors}$")

        if i_block == 0:
            ax.legend(fontsize=16, loc="upper left")
        ax.set_xlabel(r"$N_{\rm core}$")
        ax.set_ylabel(r"$f_{\rm err}$")
        ax.set_xscale("log")
        ax.set_ylim(0, 0.23)
        ax.set_xlim(1, 32*1.05)

        fig.savefig(path.join(plot_dir, "param_comp_%d.pdf" % i_block))

if __name__ == "__main__": main()
