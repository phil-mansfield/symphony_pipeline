import symlib
import numpy as np
import matplotlib.pyplot as plt
import palette
from palette import pc

palette.configure(False)

base_dir = "/Users/phil/data/ZoomIns"
eden_suite = "EDEN_MilkyWay_8K"
sym_suite = "SymphonyMilkyWay"

error_list = [9, 10, 28, 49]
m_list = [45, 46, 47, 48, 49, 50]
m_match = {45: 10, 46: 14, 47: 19, 48: 26, 49: 33, 50: 34}

def get_dirs(i):
    if i in error_list: return None, None, None
    eden_dir = symlib.get_host_directory(base_dir, eden_suite, i)
    j = m_match[i] if i in m_list else i
    sym_dir = symlib.get_host_directory(base_dir, sym_suite, j)

    return eden_dir, sym_dir, i in m_list

def is_orbiting_subhalo(sf, hist):
    r_dot_v = np.sum(sf["x"]*sf["v"], axis=2)    
    past_peri = np.any(r_dot_v > 0, axis=1)
    is_ok = sf["ok"][:,-1]
    small_enough = hist["merger_ratio"] < 0.15
    return past_peri & is_ok & small_enough

def pericenters():
    m_flag_e, m_flag_s = [], []
    r_peri_e, r_peri_s = [], []
    mpeak_e, mpeak_s = [], []
    
    n_host, n_host_m = 0, 0
    
    for i in range(symlib.n_hosts(eden_suite)):        
        eden_dir, sym_dir, is_m = get_dirs(i)
        if eden_dir is None: continue
        print(i)
        
        sf_e, hist_e = symlib.read_symfind(eden_dir)
        sf_s, hist_s = symlib.read_symfind(sym_dir)
        ok_e = is_orbiting_subhalo(sf_e, hist_e)
        ok_s = is_orbiting_subhalo(sf_s, hist_s)
        
        if is_m:
            m_flag_e.append(np.ones(np.sum(ok_e), dtype=bool))
            m_flag_s.append(np.ones(np.sum(ok_s), dtype=bool))
            n_host_m += 1
        else:
            m_flag_e.append(np.zeros(np.sum(ok_e), dtype=bool))
            m_flag_s.append(np.zeros(np.sum(ok_s), dtype=bool))
            n_host += 1

        r_e = np.sqrt(np.sum(sf_e["x"]**2, axis=2))
        r_s = np.sqrt(np.sum(sf_s["x"]**2, axis=2))
        for j in range(len(sf_e)):
            r_e[j,~sf_e["ok"][j]] = np.inf
        for j in range(len(sf_s)):
            r_s[j,~sf_s["ok"][j]] = np.inf

        # We need to do this more carefully.
        r_peri_e.append(np.min(r_e, axis=1)[ok_e])
        r_peri_s.append(np.min(r_s, axis=1)[ok_s])
        mpeak_e.append(hist_e["mpeak"][ok_e])
        mpeak_s.append(hist_s["mpeak"][ok_s])
        
    param = symlib.simulation_parameters(eden_dir)
    mp = param["mp"]/param["h100"]
        
    m_flag_e, m_flag_s = np.hstack(m_flag_e), np.hstack(m_flag_s)
    mpeak_e, mpeak_s = np.hstack(mpeak_e), np.hstack(mpeak_s)
    r_peri_e, r_peri_s = np.hstack(r_peri_e), np.hstack(r_peri_s)

    print("%g" % (mp*4e3))
    
    fig, ax = plt.subplots(1, 2, sharey=True, figsize=(14, 6.5))
    bins = np.linspace(0, 120, 25)
    dx = bins[1] - bins[0]
    ok_e = (mpeak_e/mp > 4e3) & (m_flag_e == 0)
    ok_s = (mpeak_s/mp > 4e3) & (m_flag_s == 0)
    ax[0].hist(r_peri_e[ok_e], histtype="step", lw=3, color=pc("r"),
               bins=bins, weights=np.ones(np.sum(ok_e))/n_host/dx, density=True)
    ax[0].hist(r_peri_s[ok_s], histtype="step", lw=3, color=pc("b"),
               bins=bins, weights=np.ones(np.sum(ok_s))/n_host/dx, density=True)
    
    bins = np.linspace(0, 120, 12)
    dx = bins[1] - bins[0]
    ok_e = (mpeak_e/mp > 4e3) & (m_flag_e == 1)
    ok_s = (mpeak_s/mp > 4e3) & (m_flag_s == 1)
    ax[1].hist(r_peri_e[ok_e], histtype="step", lw=3, color=pc("r"),
               bins=bins, weights=np.ones(np.sum(ok_e))/n_host_m/dx, density=True)
    ax[1].hist(r_peri_s[ok_s], histtype="step", lw=3, color=pc("b"),
               bins=bins, weights=np.ones(np.sum(ok_s))/n_host_m/dx, density=True)

    ax[0].set_ylabel(r"${\rm PDF}$")
    ax[0].set_xlabel(r"$r_{\rm peri}$")
    ax[1].set_xlabel(r"$r_{\rm peri}$")
    ax[0].set_title(r"${\rm fiducial-mass}$")
    ax[1].set_title(r"${\rm high-mass}$")
    ax[0].set_xlim(0, 120)
    ax[1].set_xlim(0, 120)

    ax[0].plot([], [], c=pc("r"), label=r"${\rm EDEN\ (disk)}$")
    ax[0].plot([], [], c=pc("b"), label=r"${\rm Symphony\ (DMO)}$")
    ax[0].legend(loc="upper right", fontsize=16)

def infall_time():
    m_flag_e, m_flag_s = [], []
    a_infall_e, a_infall_s = [], []
    mpeak_e, mpeak_s = [], []
    
    n_host, n_host_m = 0, 0
    
    for i in range(symlib.n_hosts(eden_suite)):        
        eden_dir, sym_dir, is_m = get_dirs(i)
        if eden_dir is None: continue
        print(i)
        
        sf_e, hist_e = symlib.read_symfind(eden_dir)
        sf_s, hist_s = symlib.read_symfind(sym_dir)
        ok_e = sf_e["ok"][:,-1]
        ok_s = sf_s["ok"][:,-1]
        
        if is_m:
            m_flag_e.append(np.ones(np.sum(ok_e), dtype=bool))
            m_flag_s.append(np.ones(np.sum(ok_s), dtype=bool))
            n_host_m += 1
        else:
            m_flag_e.append(np.zeros(np.sum(ok_e), dtype=bool))
            m_flag_s.append(np.zeros(np.sum(ok_s), dtype=bool))
            n_host += 1
        scale = symlib.scale_factors(sym_dir)
            
        a_infall_e.append(scale[hist_e["merger_snap"][ok_e]])
        a_infall_s.append(scale[hist_s["merger_snap"][ok_s]])
        mpeak_e.append(hist_e["mpeak"][ok_e])
        mpeak_s.append(hist_s["mpeak"][ok_s])
        
    param = symlib.simulation_parameters(eden_dir)
    mp = param["mp"]/param["h100"]
        
    m_flag_e, m_flag_s = np.hstack(m_flag_e), np.hstack(m_flag_s)
    mpeak_e, mpeak_s = np.hstack(mpeak_e), np.hstack(mpeak_s)
    a_infall_e, a_infall_s = np.hstack(a_infall_e), np.hstack(a_infall_s)

    print("%g" % (mp*4e3))
    
    fig, ax = plt.subplots(1, 2, sharey=True, figsize=(14, 6.5))
    bins = np.linspace(0, 1, 25)
    dx = bins[1] - bins[0]
    ok_e = (mpeak_e/mp > 4e3) & (m_flag_e == 0)
    ok_s = (mpeak_s/mp > 4e3) & (m_flag_s == 0)
    ax[0].hist(a_infall_e[ok_e], histtype="step", lw=3, color=pc("r"),
               bins=bins, weights=np.ones(np.sum(ok_e))/n_host/dx, density=True)
    ax[0].hist(a_infall_s[ok_s], histtype="step", lw=3, color=pc("b"),
               bins=bins, weights=np.ones(np.sum(ok_s))/n_host/dx, density=True)
    
    bins = np.linspace(0, 1, 12)
    dx = bins[1] - bins[0]
    ok_e = (mpeak_e/mp > 4e3) & (m_flag_e == 1)
    ok_s = (mpeak_s/mp > 4e3) & (m_flag_s == 1)
    ax[1].hist(a_infall_e[ok_e], histtype="step", lw=3, color=pc("r"),
               bins=bins, weights=np.ones(np.sum(ok_e))/n_host_m/dx, density=True)
    ax[1].hist(a_infall_s[ok_s], histtype="step", lw=3, color=pc("b"),
               bins=bins, weights=np.ones(np.sum(ok_s))/n_host_m/dx, density=True)

    ax[0].set_ylabel(r"${\rm PDF}$")
    ax[0].set_xlabel(r"$a_{\rm infall}$")
    ax[1].set_xlabel(r"$a_{\rm infall}$")
    ax[0].set_title(r"${\rm fiducial-mass}$")
    ax[1].set_title(r"${\rm high-mass}$")
    ax[0].set_xlim(0, 1)
    ax[1].set_xlim(0, 1)
    ax[0].set_ylim
    
    ax[0].plot([], [], c=pc("r"), label=r"${\rm EDEN\ (disk)}$")
    ax[0].plot([], [], c=pc("b"), label=r"${\rm Symphony\ (DMO)}$")
    ax[0].legend(loc="upper right", fontsize=16)
    
def shmf():
    vmax_e, vmax_s = [], []
    mpeak_e, mpeak_s = [], []
    mu_e, mu_s = [], []
    r_e, r_s = [], []
    m_flag_e, m_flag_s = [], []
    
    n_host, n_host_m = 0, 0
    
    for i in range(symlib.n_hosts(eden_suite)):        
        eden_dir, sym_dir, is_m = get_dirs(i)
        if eden_dir is None: continue
        print(i)
        
        sf_e, hist_e = symlib.read_symfind(eden_dir)
        sf_s, hist_s = symlib.read_symfind(sym_dir)
        rs_e, _ = symlib.read_rockstar(eden_dir)
        rs_s, _ = symlib.read_rockstar(sym_dir)
        ok_e, ok_s = sf_e["ok"][:,-1], sf_s["ok"][:,-1]
        rvir_e, rvir_s = rs_e["rvir"][0,-1], rs_s["rvir"][0,-1]
        
        if is_m:
            m_flag_e.append(np.ones(np.sum(sf_e["ok"][:,-1]), dtype=bool))
            m_flag_s.append(np.ones(np.sum(sf_s["ok"][:,-1]), dtype=bool))
            n_host_m += 1
        else:
            m_flag_e.append(np.zeros(np.sum(sf_e["ok"][:,-1]), dtype=bool))
            m_flag_s.append(np.zeros(np.sum(sf_s["ok"][:,-1]), dtype=bool))
            n_host += 1
        
        vmax_e.append(sf_e["vmax"][ok_e,-1])
        vmax_s.append(sf_s["vmax"][ok_s,-1])
        mpeak_e.append(hist_e["mpeak"][ok_e])
        mpeak_s.append(hist_s["mpeak"][ok_s])
        mu_e.append(sf_e["m"][ok_e,-1]/hist_e["mpeak"][ok_e])
        mu_s.append(sf_s["m"][ok_s,-1]/hist_s["mpeak"][ok_s])
        r_e.append(np.sqrt(np.sum(sf_e["x"][ok_e,-1]**2, axis=1))/rvir_e)
        r_s.append(np.sqrt(np.sum(sf_s["x"][ok_s,-1]**2, axis=1))/rvir_s)
        
    m_flag_e, m_flag_s = np.hstack(m_flag_e), np.hstack(m_flag_s)
    vmax_e, vmax_s = np.hstack(vmax_e), np.hstack(vmax_s)
    mpeak_e, mpeak_s = np.hstack(mpeak_e), np.hstack(mpeak_s)
    mu_e, mu_s = np.hstack(mu_e), np.hstack(mu_s)
    r_e, r_s = np.hstack(r_e), np.hstack(r_s)
    
    for mode in range(2):
        fig, ax = plt.subplots(1, 2, sharey=True, figsize=(14, 6.5))
        if mode == 0:
            ok_e, ok_s = ~m_flag_e, ~m_flag_s
            weight_e = np.ones(np.sum(ok_e))/n_host_m
            weight_s = np.ones(np.sum(ok_s))/n_host_m
        else:
            ok_e, ok_s = m_flag_e, m_flag_s
            weight_e = np.ones(np.sum(ok_e))/n_host
            weight_s = np.ones(np.sum(ok_s))/n_host
    
        for sim in range(2):
            if sim == 0:
                ls = "-"
                weight, ok = weight_e, ok_e
                vmax, mpeak, mu, r = vmax_e[ok], mpeak_e[ok], mu_e[ok], r_e[ok]
            else:
                ls = "--"
                weight, ok = weight_s, ok_s
                vmax, mpeak, mu, r = vmax_s[ok], mpeak_s[ok], mu_s[ok], r_s[ok]

            vmax_order = np.argsort(vmax)
            mpeak_order = np.argsort(mpeak)

            r_cuts = [0.25, 0.5, 100.0]
            colors = [pc("b"), pc("o"), pc("r")]

            for r_i in range(3):
                ok = r[vmax_order] < r_cuts[r_i]
                ax[0].plot(vmax[vmax_order][ok],
                           np.arange(np.sum(ok))[::-1],
                           c=colors[r_i], ls=ls)
                ok = r[mpeak_order] < r_cuts[r_i]
                ax[1].plot(mpeak[mpeak_order][ok],
                           np.arange(np.sum(ok))[::-1],
                           c=colors[r_i], ls=ls)

        ax[0].set_xlabel(r"$v_{\rm max}\ ({\rm km\,s^{-1}})$")
        ax[1].set_xlabel(r"$m_{\rm peak}\ (M_\odot)$")
        ax[0].set_ylabel(r"$N(>v_{\rm max};\,>m_{\rm peak})$")
        ax[0].set_xscale("log")
        ax[1].set_xscale("log")
        ax[0].set_yscale("log")
        ax[1].set_yscale("log")
        ax[0].set_xlim(10, 200)
        ax[1].set_xlim(1e8, 1e12)
        ax[0].plot([], [], pc("r"), label=r"${\rm all\ subhalos}$")
        ax[0].plot([], [], pc("o"), label=r"$r < R_{\rm vir}/2$")
        ax[0].plot([], [], pc("b"), label=r"$r < R_{\rm vir}/4$")
        ax[0].plot([], [], pc("a"), label=r"${\rm EDEN\ (disk)}$")
        ax[0].plot([], [], pc("a"), ls="--", label=r"${\rm Symphony\ (CDM)}$")

        ax[0].legend(loc="upper right", fontsize=16)
        if mode == 0:
            plt.suptitle(r"${\rm Fiducial-mass\ disks}$")
        else:
            plt.suptitle(r"${\rm High-mass\ disks}$")
        
def main():
    shmf()
    pericenters()
    infall_time()
    
    plt.show()
    
if __name__ == "__main__": main()
