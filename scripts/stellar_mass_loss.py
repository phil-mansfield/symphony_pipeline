import numpy as np
import matplotlib.pyplot as plt
import palette
from palette import pc
import cache_stars
import symlib
from colossus.cosmology import cosmology
import survival

base_dir = "/sdf/home/p/phil1/ZoomIns"
HYDRO = False

def t_relax_t_orbit(Nr, r, eps):
    return Nr/4 * (np.log(r**2/eps**2 + 1) +
                   eps**2-2*r**2/(3*(eps**2+r**2)) -
                   np.log(3/2))**-1

def t_orbit(Mr, r):
    return 7.5 * (r/40)**1.5 * (1e10/Mr)**0.5


def mixing_times(p, mp, eps, hydro=False):
    r = np.sqrt(np.sum(p["x"]**2, axis=1))
        
    order = np.argsort(r)
    Nr = np.zeros(len(r))
    Nr_sort = np.arange(len(Nr)) + 1
    Nr[order] = Nr_sort

    t_orb = t_orbit(Nr*mp, r)
    t_relax = t_relax_t_orbit(Nr, r, eps)*t_orb

    if HYDRO:
        f_baryon = 0.02212 / (0.1206 + 0.02212)
        t_relax *= f_baryon

    return t_relax * 0.3

def example_stellar_mass_loss(suite, i_host, i_sub):
    plt.figure()

    sim_dir = symlib.get_host_directory(base_dir, suite, i_host)
    param = symlib.simulation_parameters(sim_dir)
    scale = symlib.scale_factors(sim_dir)
    mp, eps = param["mp"]/param["h100"], param["eps"]/param["h100"]*scale

    cosmo = cosmology.setCosmology('', symlib.colossus_parameters(param))
    age = cosmo.age(1/scale - 1)

    sf, hist = symlib.read_symfind(sim_dir)
    
    r50_mults = [0.005, 0.008, 0.015, 0.025, 0.05]
    method_names = ["r=0.005", "r=0.008", "r=0.015", "r=0.025", "r=0.05"]
    labels = [r"$r_{\star,50}/r_{\rm vir} = 0.005$",
              r"$r_{\star,50}/r_{\rm vir} = 0.008$",
              r"$r_{\star,50}/r_{\rm vir} = 0.015$",
              r"$r_{\star,50}/r_{\rm vir} = 0.025$",
              r"$r_{\star,50}/r_{\rm vir} = 0.050$"]
    colors = [pc("r"), pc("o"), pc("g"), pc("b"), pc("p")]

    snap_i = hist["first_infall_snap"][i_sub]
    part = symlib.Particles(sim_dir)
    p = part.read(snap_i, halo=i_sub, mode="smooth")
    p["x"] -= sf["x"][i_sub, snap_i]
    p["v"] -= sf["v"][i_sub, snap_i]

    t_mix = mixing_times(p, mp, eps[snap_i], hydro=HYDRO)
    dt = age - age[snap_i]

    for i in range(len(r50_mults)):
        print(method_names[i])
        gal, gal_hist = symlib.read_galaxies(sim_dir, method_names[i])
        ok = sf["ok"][i_sub]
        conv = gal["m23_m_conv"][i_sub] & ok
        m_frac = sf["m"][i_sub]/sf["m"][i_sub,snap_i]
        m_star_frac = gal["m_star"][i_sub]/gal_hist["m_star_i"][i_sub]

        plt.plot(m_frac[ok], m_star_frac[ok], c=colors[i],
                 label=labels[i])

        mp_star = cache_stars.read_stars(suite, i_host, method_names[i])[0][i_sub]
        m_mixed = np.zeros(len(age))
        p_ok = p["ok"]
        m_star_i = np.sum(mp_star[p_ok])
        for snap in range(snap_i, len(age)):
            if not sf["ok"][i_sub,snap]: break
            m_mixed[snap] = np.sum(mp_star[p_ok & (t_mix < dt[snap])])
        m_mixed_frac = m_mixed/m_star_i
        plt.plot(m_frac[ok], m_mixed_frac[ok], "--", c=colors[i], lw=1.5)

    plt.xscale("log")
    #plt.yscale("log")
    plt.ylim(0, 1.4)

    plt.xlabel(r"$m/m_{\rm infall}$")
    plt.ylabel(r"$m_\star/m_{\star, \rm infall}$")
    plt.legend(loc="upper left", frameon=True, fontsize=16)

    if HYDRO:
        plt.savefig("../plots/stellar_halo/examples/m_star_loss_%s_%d_%d_hydro.png" %
                    (suite, i_host, i_sub))
    else:
        plt.savefig("../plots/stellar_halo/examples/m_star_loss_%s_%d_%d.png" %
                    (suite, i_host, i_sub))

def print_cutoffs(suite, i_host, n_peak_min, n_peak_max):
    sim_dir = symlib.get_host_directory(base_dir, suite, i_host)
    param = symlib.simulation_parameters(sim_dir)
    scale = symlib.scale_factors(sim_dir)
    mp, eps = param["mp"]/param["h100"], param["eps"]/param["h100"]*scale

    cosmo = cosmology.setCosmology('', symlib.colossus_parameters(param))
    age = cosmo.age(1/scale - 1)

    sf, hist = symlib.read_symfind(sim_dir)
    
    r50_mults = [0.005, 0.008, 0.015, 0.025, 0.05]
    method_names = ["r=0.005", "r=0.008", "r=0.015", "r=0.025", "r=0.05"]

    gals, gal_hists = [None]*5, [None]*5
    for i in range(len(gals)):
        try:
            gals[i], gal_hists[i] = symlib.read_galaxies(sim_dir, method_names[i])
        except:
            pass

    part = symlib.Particles(sim_dir)
    mp_stars = cache_stars.read_stars(suite, i_host, method_names[i])[0]
    for i_sub in range(1, len(sf)):
        if not n_peak_min < hist["mpeak"][i_sub]/mp < n_peak_max or hist["merger_ratio"][i_sub] > 0.1:
            continue

        snap_i = hist["first_infall_snap"][i_sub]
        p = part.read(snap_i, halo=i_sub, mode="smooth")
        p["x"] -= sf["x"][i_sub, snap_i]
        p["v"] -= sf["v"][i_sub, snap_i]
        mp_star = mp_stars[i_sub]

        t_mix = mixing_times(p, mp, eps[snap_i], hydro=HYDRO)
        dt = age - age[snap_i]

        for i in range(len(r50_mults)):
            gal, gal_hist = gals[i], gal_hists[i]
            if gal is None: continue
            ok = sf["ok"][i_sub]
            conv = gal["m23_m_conv"][i_sub] & ok
            m_frac = sf["m"][i_sub]/sf["m"][i_sub,snap_i]
            m_star_frac = gal["m_star"][i_sub]/gal_hist["m_star_i"][i_sub]

            m_mixed = np.zeros(len(age))
            p_ok = p["ok"]
            m_star_i = np.sum(mp_star[p_ok])
            for snap in range(snap_i, len(age)):
                if not sf["ok"][i_sub,snap]: break
                m_mixed[snap] = np.sum(mp_star[p_ok & (t_mix < dt[snap])])
            m_mixed_frac = m_mixed/m_star_i

            if np.sum(m_mixed_frac < m_star_frac) == 0:
                last_frac = 1.0
                censored = True
            else:
                last_snap = np.where(m_mixed_frac < m_star_frac)[0][-1]
                last_frac = m_star_frac[last_snap]
                censored = (last_snap == len(m_mixed_frac) - 1 or
                            not sf["ok"][i_sub][last_snap+1])

            print(i_host, i_sub, i, last_frac, 1 if censored else 0)

def m_star_frac_lim_dist():
    in_names = ["tables/m_star_frac_limits_1e5.txt",
                "tables/m_star_frac_limits_1e4.txt",
                "tables/m_star_frac_limits_1e5_hydro.txt",
                "tables/m_star_frac_limits_1e4_hydro.txt"]
    titles = [r"$10^5 < n_{\rm peak} < 10^6\ ({\rm DMO})$",
              r"$10^4 < n_{\rm peak} < 10^5\ ({\rm DMO})$",
              r"$10^5 < n_{\rm peak} < 10^6\ ({\rm Hydro\ Model})$",
              r"$10^4 < n_{\rm peak} < 10^5\ ({\rm Hydro\ Model})$"]
    out_names = ["../plots/stellar_halo/m_star_frac_lim_1e5.png",
                 "../plots/stellar_halo/m_star_frac_lim_1e4.png",
                 "../plots/stellar_halo/m_star_frac_lim_1e5_hydro.png",
                 "../plots/stellar_halo/m_star_frac_lim_1e4_hydro.png"]

    colors = [pc("r"), pc("o"), pc("g"), pc("b"), pc("p")]

    for i_part in range(len(in_names)):
        plt.cla()

        _, _, i_mult, m_frac, cen = np.loadtxt(in_names[i_part]).T
        cen = cen > 0

        bad_data = np.isnan(m_frac) | (m_frac > 1) | (m_frac < 0)
        m_frac[bad_data] = 1
        cen[bad_data] = True

        m_eval = 10**np.linspace(-2, 1, 200)
        m_eval = np.linspace(0, 1, 200)
        for i in range(5):
            ok = i == i_mult
            S, S_err = survival.kaplan_meier(m_frac[ok], cen[ok], m_eval,
                                             decreasing=True)
            
            plt.plot(m_eval, S, c=colors[i])
            plt.fill_between(m_eval, S+S_err, S-S_err,
                             color=colors[i], alpha=0.2)

        if m_eval[0] != 0:
            plt.xscale("log")
        plt.xlim(m_eval[0], m_eval[-1])

        plt.title(titles[i_part])
        plt.xlabel(r"$\mu_{\star,\rm lim} \equiv m_{\star, \rm lim}/m_{\star,\rm infall}$")
        plt.ylabel(r"${\rm Pr}(< \mu_{\star,\rm lim})$")
        plt.savefig(out_names[i_part])

def main():
    palette.configure(False)
    #"""
    suite = "SymphonyMilkyWayHR"
    print("# 1e5")
    for i_host in range(symlib.n_hosts(suite)):
        print_cutoffs(suite, i_host, 1e5, 1e6)
    print()
    print()
    print("# 1e4")
    for i_host in range(symlib.n_hosts(suite)):
        print_cutoffs(suite, i_host, 1e4, 1e5)
    #"""
    #m_star_frac_lim_dist()    

if __name__ == "__main__": main()
