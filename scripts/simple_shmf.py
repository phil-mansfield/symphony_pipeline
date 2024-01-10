import numpy as np
import matplotlib.pyplot as plt
import palette
import symlib
from palette import pc

base_dir = "/sdf/home/p/phil1/ZoomIns"
suite = "SymphonyMilkyWayHR"

r_cuts = [0.25, 0.5, 1]
colors = [pc("b"), pc("o"), pc("r")]

vmax_range = (2, 200)
mass_range = (1e-5, 0.1)
mpeak_range = (1e-5, 0.1)

def shmf(m, n_host, m_range):
    bins = 10**np.linspace(np.log10(m_range[0]), np.log10(m_range[1]), 200)
    N, _ = np.histogram(m, bins)
    N = np.cumsum(N[::-1])[::-1]
    ok = N > 0
    return bins[1:][ok], N[ok]

def main():
    palette.configure(False)

    vmax_ratio = []
    mass_ratio = []
    mass_rs, mass_sf = [], []
    vmax_rs, vmax_sf = [], []
    mpeak_rs, mpeak_sf = [], []
    r_rs, r_sf = [], []
    
    n_host = symlib.n_hosts(suite)
    for i_host in range(n_host):
        print(i_host)
        if i_host >= 1: break
        sim_dir = symlib.get_host_directory(base_dir, suite, i_host)
        
        sf, hist = symlib.read_symfind(sim_dir)
        rs, _ = symlib.read_rockstar(sim_dir)

        ok_sf, ok_rs = sf["ok"][:,-1], rs["ok"][:,-1]
        ok_sf[0], ok_rs[0] = False, False
        shared = ok_sf & ok_rs
        vmax_ratio.append(sf["vmax"][shared,-1]/rs["vmax"][shared,-1])
        mass_ratio.append(sf["m"][shared,-1]/rs["m"][shared,-1])
        mass_sf.append(sf["m"][ok_sf, -1]/rs["m"][0,-1])
        mass_rs.append(rs["m"][ok_rs, -1]/rs["m"][0,-1])
        vmax_sf.append(sf["vmax"][ok_sf, -1])
        vmax_rs.append(rs["vmax"][ok_rs, -1])
        mpeak_sf.append(hist["mpeak"][ok_sf]/rs["m"][0,-1])
        mpeak_rs.append(hist["mpeak"][ok_rs]/rs["m"][0,-1])
        r_sf.append(np.sqrt(np.sum(sf["x"][:,-1]**2, axis=1))[ok_sf]
                    / rs[0,-1]["rvir"])
        r_rs.append(np.sqrt(np.sum(rs["x"][:,-1]**2, axis=1))[ok_rs]
                    / rs[0,-1]["rvir"])

    vmax_ratio = np.hstack(vmax_ratio)
    mass_ratio = np.hstack(mass_ratio)
    mass_sf, mass_rs = np.hstack(mass_sf), np.hstack(mass_rs)
    vmax_sf, vmax_rs = np.hstack(vmax_sf), np.hstack(vmax_rs)
    mpeak_sf, mpeak_rs = np.hstack(mpeak_sf), np.hstack(mpeak_rs)
    r_sf, r_rs = np.hstack(r_sf), np.hstack(r_rs)

    fig_mass, ax_mass = plt.subplots()
    fig_vmax, ax_vmax = plt.subplots()
    fig_mpeak, ax_mpeak = plt.subplots()
    
    vmax_corr, mass_corr = np.median(vmax_ratio), np.median(mass_ratio)
    print(vmax_corr, mass_corr)

    for i_cut in range(len(r_cuts)):
        ok_sf, ok_rs = r_sf < r_cuts[i_cut], r_rs < r_cuts[i_cut]

        m, N = shmf(mass_sf[ok_sf], n_host, mass_range)
        ax_mass.plot(m, N, c=colors[i_cut])
        m, N = shmf(mass_rs[ok_rs]*mass_corr, n_host, mass_range)
        ax_mass.plot(m, N, "--", c=colors[i_cut])

        m, N = shmf(vmax_sf[ok_sf], n_host, vmax_range)
        ax_vmax.plot(m, N, c=colors[i_cut])
        m, N = shmf(vmax_rs[ok_rs]*vmax_corr, n_host, vmax_range)
        ax_vmax.plot(m, N, "--", c=colors[i_cut])

        m, N = shmf(mpeak_sf[ok_sf], n_host, mpeak_range)
        ax_mpeak.plot(m, N, c=colors[i_cut])
        m, N = shmf(mpeak_rs[ok_rs], n_host, mpeak_range)
        ax_mpeak.plot(m, N, "--", c=colors[i_cut])
        
    ax_mass.set_xscale("log")
    ax_mass.set_yscale("log")
    ax_mass.set_xlabel(r"$m/M_{\rm vir}$")
    ax_mass.set_ylabel(r"$N(>m)$")

    ax_vmax.set_xscale("log")
    ax_vmax.set_yscale("log")
    ax_vmax.set_xlabel(r"$v_{\rm max}\ ({\rm km\,s^{-1}})$")
    ax_vmax.set_ylabel(r"$N(>v_{\rm max})$")

    ax_mpeak.set_xscale("log")
    ax_mpeak.set_yscale("log")
    ax_mpeak.set_xlabel(r"$m_{\rm peak}/M_{\rm vir}$")
    ax_mpeak.set_ylabel(r"$N(>m_{\rm peak})$")

    fig_mass.savefig("../plots/shmf_mass.pdf")
    fig_vmax.savefig("../plots/shmf_vmax.pdf")
    fig_mpeak.savefig("../plots/shmf_mpeak.pdf")
        

if __name__ == "__main__": main()
