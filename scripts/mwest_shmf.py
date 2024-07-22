import symlib
import palette
from palette import pc
import numpy as np
import matplotlib.pyplot as plt

palette.configure(False)

base_dir = "/sdf/home/p/phil1/ZoomIns"
heavy_hosts_dmo = [
    "Halo119", "Halo199", "Halo364", "Halo416", "Halo460",
    "Halo530", "Halo675", "Halo800", "Halo829", "Halo852"
]
heavy_hosts = [name + "_m" for name in heavy_hosts_dmo]
light_hosts_dmo = np.arange(symlib.n_hosts("SymphonyMilkyWay"), dtype=int)
light_hosts = np.copy(light_hosts_dmo)

""" Caching isn't neccessary, it just makes it less slow to reread the halos
in if you keep the mass func funciton simple.
"""
rs_cache, sym_cache = {}, {}
def cache_read(sim_dir, use_rockstar):    
    if use_rockstar:
        if sim_dir in rs_cache:
            h, hist = rs_cache[sim_dir]
        else:
            h, hist = symlib.read_rockstar(sim_dir)
            rs_cache[sim_dir] = (h, hist)
    else:
        if sim_dir in sym_cache:
            h, hist = sym_cache[sim_dir]
        else:
            h, hist = symlib.read_symfind(sim_dir)
            sym_cache[sim_dir] = (h, hist)
        
    return h, hist


def mass_func(ax, suite, host_names, r_cut, use_rockstar=False, mass_def="m",
              c=pc("r"), ls="-"):
    n_host = 0
    m_all = []
    for i in range(len(host_names)):
        #if i > 5: continue
        name = host_names[i]
        sim_dir = symlib.get_host_directory(base_dir, suite, name)
        print("   ", name, sim_dir)
        
        h, hist = cache_read(sim_dir, use_rockstar)
        
        ok = h["ok"][:,-1]
        x = h["x"][:,-1,:]
        r = np.sqrt(np.sum(x**2, axis=1))
        ok = ok & (r < r_cut)
        
        if use_rockstar:
            sf, _ = cache_read(sim_dir, False)
            #ok = ok & sf["ok_rs"][:,-1]
            ok[0] = False

        m = h[mass_def][:,-1] if mass_def in ["m", "vmax"] else hist[mass_def]
        idx = np.where(m[ok] > 2e11)[0]
        m_all.append(m[ok])
        n_host += 1
        
    m = np.hstack(m_all)

    log_low, log_high = np.log10(np.min(m)), np.log10(np.max(m))
    bins = 10**np.linspace(log_low, log_high, 300)

    n, _ = np.histogram(m, bins=bins)
    n = np.cumsum(n[::-1])[::-1] / n_host
    
    ax.plot(bins[1:], n, c=c, ls=ls)
    
def main():
    r_cut = [250, 50]

    for i in range(len(r_cut)):
        fig, ax = plt.subplots()
        print("r_cut =", r_cut)
        mass_func(ax, "SymphonyMilkyWay", light_hosts_dmo, r_cut[i],
                  mass_def="mpeak", c=pc("r"), ls="-", use_rockstar=False)
        mass_func(ax, "EDEN_MilkyWay_8K", light_hosts, r_cut[i],
                  mass_def="mpeak", c=pc("b"), ls="-", use_rockstar=False)
        mass_func(ax, "SymphonyMilkyWay", light_hosts_dmo, r_cut[i],
                  mass_def="mpeak", c=pc("r"), ls="--", use_rockstar=True)
        mass_func(ax, "EDEN_MilkyWay_8K", light_hosts, r_cut[i],
                  mass_def="mpeak", c=pc("b"), ls="--", use_rockstar=True)

        ax.plot([], [], "-", c=pc("r"),
                label=r"${\rm Symfind,\ DMO}$")
        ax.plot([], [], "-", c=pc("b"),
                label=r"${\rm Symfind,\ Disk}$")
        ax.plot([], [], "--", c=pc("r"),
                label=r"${\rm Rockstar,\ DMO}$")
        ax.plot([], [], "--", c=pc("b"),
                label=r"${\rm Rockstar,\ Disk}$")

        ax.legend(loc="upper right")
        ax.set_ylim(0.03, 300)
        ax.set_xlim(1e8, 3e11)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(r"$m_{\rm peak}$")
        ax.set_ylabel(r"$N(>m_{\rm peak})$")
    
        fig.savefig("../plots/m_func_%d.png" % i)
    
    
if __name__ == "__main__": main()
