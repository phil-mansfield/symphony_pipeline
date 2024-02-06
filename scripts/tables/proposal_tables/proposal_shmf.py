import numpy as np
import matplotlib.pyplot as plt
import palette
from palette import pc

palette.configure(True)

def main():
    plt.figure()
    m_s1, n_s1 = np.loadtxt("Symfind_Rvir.csv").T
    m_r1, n_r1 = np.loadtxt("Rockstar_Rvir.csv").T
    m_s4, n_s4 = np.loadtxt("Symfind_Rvir4.csv").T
    m_r4, n_r4 = np.loadtxt("Rockstar_Rvir4.csv").T

    plt.plot(10**m_s1, 10**n_s1, c=pc("r"),
             label=r"${\rm Symfind\ }(r_{\rm sub}<R_{\rm halo})$")
    plt.plot(10**m_r1, 10**n_r1, "--", c=pc("r"),
             label=r"${\rm Rockstar\ }(r_{\rm sub}<R_{\rm halo})$")
    plt.plot(10**m_s4, 10**n_s4, c=pc("b"),
             label=r"${\rm Symfind\ }(r_{\rm sub}<R_{\rm halo}/4)$")
    plt.plot(10**m_r4, 10**n_r4, "--", c=pc("b"),
             label=r"${\rm Rockstar\ }(r_{\rm sub}<R_{\rm halo}/4)$")

    plt.xscale("log")
    plt.yscale("log")

    plt.xlabel(r"$M_{\rm sub}/M_{\rm halo}$")
    plt.ylabel(r"$\langle{\rm Cumulative\ number\ of\ subhalos}\rangle$")

    plt.legend(loc="upper right", fontsize=17)

    plt.savefig("shmf.png")

if __name__ == "__main__": main()
