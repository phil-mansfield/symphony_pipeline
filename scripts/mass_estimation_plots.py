import numpy as np
import matplotlib.pyplot as plt
import symlib
import gravitree
import palette
from palette import pc

SUITE = "SymphonyMilkyWay"
BASE_DIR = "/oak/stanford/orgs/kipac/users/phil1/simulations/ZoomIns/"

def print_mass_table():
    i_halo = 0
    f = open("tables/mass_estimates.txt", "a+")
    
    for i_halo in range(symlib.n_hosts(SUITE)):
        print(i_halo)

        sim_dir = symlib.get_host_directory(BASE_DIR, SUITE, i_halo)

        c = symlib.read_cores(sim_dir)
        h, hist = symlib.read_subhalos(sim_dir)
        h_cmov, hist_cmov = symlib.read_subhalos(sim_dir, comoving=True)
        param = symlib.simulation_parameters(SUITE)
        mp = param["mp"]/param["h100"]
        eps = param["eps"]/param["h100"]

        host = h[0]
        h, c, hist = h[1:], c[1:], hist[1:]
        m_vir, m_peak = h["mvir"][:,-1], hist["mpeak"]

        dr_hc = np.sqrt(np.sum((h["x"][:,-1] - c["x"][:,-1])**2, axis=1))
        r_h = np.sqrt(np.sum(h["x"][:,-1]**2, axis=1))
        r_c = np.sqrt(np.sum(c["x"][:,-1]**2, axis=1))
        r_half = c["r50_bound"][:,-1]
        m_bound = c["m_bound"][:,-1]
        m_tidal = c["m_tidal"][:,-1]

        bins = [1e8, 1e9, 1e10, 1e11]
        info = symlib.ParticleInfo(sim_dir)
        x = symlib.read_particles(info, sim_dir, 235, "x")
        v = symlib.read_particles(info, sim_dir, 235, "v")
        valid = symlib.read_particles(info, sim_dir, 235, "valid")
        owner = symlib.read_particles(info, sim_dir, 235, "ownership")
        
        is_err = dr_hc > r_half
        ok = (h["ok"][:,-1] & (m_vir > 32*mp) &
              c["ok"][:,-1] & (m_bound > 32*mp) &
              (r_c > r_half) & (~is_err))
        
        subs = np.where(ok)[0]

        for i in range(len(subs)):
            if i < 10: continue
            i_sub = subs[i]
            x_i, v_i, valid_i = x[i_sub+1], v[i_sub+1], valid[i_sub+1]
            x_i = symlib.set_units_x(x_i, h_cmov[i_sub+1,-1], 1, param)
            v_i = symlib.set_units_v(v_i, h_cmov[i_sub+1,-1], 1, param)
            
            E_1 = gravitree.binding_energy(x_i, v_i, mp, eps)
            E_2 = gravitree.binding_energy(x_i, v_i, mp, eps, 2)
            E_3 = gravitree.binding_energy(x_i, v_i, mp, eps, 3)
            E_4 = gravitree.binding_energy(x_i, v_i, mp, eps, 4)
            E_5 = gravitree.binding_energy(x_i, v_i, mp, eps, 5)
            E_6 = gravitree.binding_energy(x_i, v_i, mp, eps, 6)
            E_7 = gravitree.binding_energy(x_i, v_i, mp, eps, 7)
            E_8 = gravitree.binding_energy(x_i, v_i, mp, eps, 8)
            E_9 = gravitree.binding_energy(x_i, v_i, mp, eps, 9)
            E_10 = gravitree.binding_energy(x_i, v_i, mp, eps, 10)

            m_tree_1 = np.sum(E_1 < 0)*mp
            m_tree_2 = np.sum(E_2 < 0)*mp
            m_tree_3 = np.sum(E_3 < 0)*mp
            m_tree_4 = np.sum(E_4 < 0)*mp
            m_tree_5 = np.sum(E_5 < 0)*mp
            m_tree_6 = np.sum(E_6 < 0)*mp
            m_tree_7 = np.sum(E_7 < 0)*mp
            m_tree_8 = np.sum(E_8 < 0)*mp
            m_tree_9 = np.sum(E_9 < 0)*mp
            m_tree_10 = np.sum(E_10 < 0)*mp
            print("%2d %4d %.3g %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f"%(
                i_halo, i_sub+1,
                m_tree_10,
                r_c[i_sub]/host["rvir"][-1],
                m_vir[i_sub]/m_peak[i_sub],
                m_bound[i_sub]/m_peak[i_sub],
                m_tree_1/m_peak[i_sub],
                m_tree_2/m_peak[i_sub],
                m_tree_3/m_peak[i_sub],
                m_tree_4/m_peak[i_sub],
                m_tree_5/m_peak[i_sub],
                m_tree_6/m_peak[i_sub],
                m_tree_7/m_peak[i_sub],
                m_tree_8/m_peak[i_sub],
                m_tree_9/m_peak[i_sub],
                m_tree_10/m_peak[i_sub]
            ), file=f)

def read_mass_table():
    fname = "tables/mass_estimates.txt"
    cols = np.loadtxt(fname).T
    cols = cols[2:]
    dtype = [
        ("mass", "f4"),
        ("dist", "f4"),
        ("mu_rs", "f4"),
        ("mu_sph", "f4"),
        ("mu_iter", ("f4", (9,))),
        ("mu", "f4")
    ]
    x = np.zeros(len(cols[0]), dtype=dtype)
    x["mass"] = cols[0]
    x["dist"] = cols[1]
    x["mu_rs"] = cols[2]
    x["mu_sph"] = cols[3]
    x["mu_iter"][:,0] = cols[4]
    x["mu_iter"][:,1] = cols[5]
    x["mu_iter"][:,2] = cols[6]
    x["mu_iter"][:,3] = cols[7]
    x["mu_iter"][:,4] = cols[8]
    x["mu_iter"][:,5] = cols[9]
    x["mu_iter"][:,6] = cols[10]
    x["mu_iter"][:,7] = cols[11]
    x["mu_iter"][:,8] = cols[12]
    x["mu"] = cols[13]

    return x

def required_iters():
    t = read_mass_table()
    iters = np.arange(1, 11)
    n = np.zeros(10, dtype=int)
    n[0] = np.sum(t["mu_iter"][:,0] == t["mu"])
    for i in range(1, 10):
        if i == 9:
            n[i] = np.sum(t["mu_iter"][:,i-1] != t["mu"])
        else:
            n[i] = np.sum((t["mu_iter"][:,i] == t["mu"]) & 
                          (t["mu_iter"][:,i-1] != t["mu"]))

    n, iters = n[:-1], iters[:-1]
    err = np.sqrt(n) / len(t)
    n = n / len(t)

    p = np.polyfit(iters[1:], np.log(n[1:]), 1)

    plt.plot(iters, n, c=pc("r"))
    plt.plot(iters[1:], np.exp(p[0]*iters[1:] + p[1]), "--", lw=2, c="k",
             label=r"$f(N) = %.3f\cdot e^{%.3f\cdot N}$" % (p[1], p[0]))
    plt.fill_between(iters, n-err, n+err, color=pc("r"), alpha=0.3)
    plt.xlim(0.95, 9.05)
    plt.yscale("log")
    plt.legend(loc="upper right")
    plt.xlabel(r"$N_{\rm iters,convergence}$")
    plt.ylabel(r"$f(N_{\rm iters,convergence})$")

    plt.savefig("../plots/core_plots/iter_covergence.png")

def mass_comparison():
    t = read_mass_table()
    ok = t["mu"] > 0
    t = t[ok]
    
    bins = np.linspace(0, 2, 400)
    plt.hist(1/(t["mu"]/t["mu_rs"]), color=pc("k"), histtype="step", lw=3, bins=bins, cumulative=True, density=True, label=r"$m_X = m_{\rm Rockstar}$")
    plt.hist(1/(t["mu"]/t["mu_iter"][:,0]), color=pc("r"), histtype="step", lw=3, bins=bins, cumulative=True, density=True, label=r"$m_X = m_{\rm bound,iter=1}$")
    plt.hist(1/(t["mu"]/t["mu_iter"][:,1]), color=pc("r"), histtype="step", lw=3, ls="--", bins=bins, cumulative=True, density=True, label=r"$m_X = m_{\rm bound,iter=2}$")
    plt.hist(1/(t["mu"]/t["mu_sph"]), color=pc("b"), histtype="step", lw=3, bins=bins, cumulative=True, density=True, label=r"$m_X = m_{\rm bound,spherical}$")
    plt.xlim(0.90, 1.25)

    plt.xlabel(r"$m_X/m_{\rm bound}$")
    plt.ylabel(r"${\rm CDF}(<m_X/m_{\rm bound})$")
    plt.legend(loc="lower right")
    plt.savefig("../plots/core_plots/mass_comp.pdf")

    print_ratio_stats("rs", t["mu_rs"]/t["mu"])
    print_ratio_stats("iter=1", t["mu_iter"][:,0]/t["mu"])
    print_ratio_stats("iter=2", t["mu_iter"][:,1]/t["mu"])
    print_ratio_stats("iter=3", t["mu_iter"][:,2]/t["mu"])
    print_ratio_stats("sph", t["mu_sph"]/t["mu"])
    
def print_ratio_stats(name, ratio):
    print("%s: med = %.3f std = %.3f log std = %.3f, 68%% = (%.3f - %.3f)" %
          (name, np.median(ratio), np.std(ratio), np.std(np.log10(ratio)),
           np.quantile(ratio, 0.5-0.68/2), np.quantile(ratio, 0.5+0.68/2)))

def radial_trend():
    plt.figure()
    t = read_mass_table()
    ok = t["mu"] > 0

    plt.plot(t["mu"], t["mu_sph"]/t["mu"], ".", c=pc("b"), alpha=0.2)

    plt.xscale("log")
    plt.ylim(0.85, 1.15)

    plt.savefig("../plots/core_plots/mass_comp_radial.pdf")

def main():
    palette.configure(True)
    #print_mass_table()
    #required_iters()
    mass_comparison()
    radial_trend()

if __name__ == "__main__": main()
