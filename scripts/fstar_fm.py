import numpy as np
import matplotlib.pyplot as plt
import symlib
import palette
from palette import pc

base_dir = "/oak/stanford/orgs/kipac"
suite = "SymphonyMilkyWay"


def m_lim_eps(param, scale, h, c, ok, infall_snap):
    # Factor os 1.284 comes from Mansfield & Avestruz
    eps = scale*param["eps"]/param["h100"]
    c_infall = h["cvir"][infall_snap]
    rs_infall = h["rvir"][infall_snap]/h["cvir"][infall_snap]
    def f(x):
        return np.log(1 + x) - x/(1 + x)
    return 1.79/1.284 * ((eps[ok]*c["r50_bound"][ok])/
                         (f(c_infall)*rs_infall**2))

def m_lim_n(param, scale, h, c, ok, infall_snap):
    mp = param["mp"]/param["h100"]
    n_infall = h["mvir"][infall_snap]/mp
    return 0.32*(n_infall/1e3)**-0.8 * np.ones(np.sum(ok))

def plot_vdb_o_limits():
    sim_dir = symlib.get_host_directory(base_dir, suite, 0)
    h, hist = symlib.read_subhalos(sim_dir)
    c = symlib.read_cores(sim_dir, "fid")
    scale = symlib.scale_factors(sim_dir)
    param = symlib.simulation_parameters(sim_dir)

    ok = hist["merger_ratio"] < 1/15
    target_halos = np.where(ok)[0][1:21]

    for ih in target_halos:
        plt.cla()

        ok = c["ok"][ih]
        infall_snap = hist["first_infall_snap"][ih]
        m_infall = h["mvir"][ih,infall_snap]
        lim_1 = m_lim_eps(param, scale, h[ih], c[ih], ok, infall_snap)
        lim_2 = m_lim_n(param, scale, h[ih], c[ih], ok, infall_snap)

        plt.plot(scale[ok], c["m_bound"][ih,ok]/m_infall, c=pc("r"),
                 label=r"$m$")
        plt.plot(scale[ok], lim_1, c=pc("o"), label=r"$m_{{\rm lim},\epsilon}$")
        plt.plot(scale[ok], lim_2, c=pc("b"), label=r"$m_{{\rm lim},n}$")
        plt.legend(loc="lower left", fontsize=16, frameon=True)
        plt.xlabel(r"$a$")
        plt.ylabel(r"$m/m_{\rm infall}$")
        plt.yscale("log")
        plt.savefig("../plots/fstar_fm/limit_examples/vdb_lims_%03d.png" % ih)


def main():

    palette.configure(False)
    
    plot_vdb_o_limits()

if __name__ == "__main__": main()
