import numpy as np
import matplotlib.pyplot as plt
import symlib
from colossus.cosmology import cosmology
from colossus.halo import mass_defs, mass_so
import gravitree
import palette
from palette import pc

base_dir = "/sdf/home/p/phil1/ZoomIns"

def crossing_count(age, scale, hist, h):
    z = 1/scale - 1
    snap_i = hist["first_infall_snap"]
    snap_f = np.zeros(snap_i.shape, dtype=int)

    for i in range(1,len(snap_f)):
        snap_f[i] = np.max(np.where(h["ok"][i])[0])
    snap_f[0] = snap_i[0]

    T_cross = mass_so.dynamicalTime(z[snap_i], "vir", "crossing")
    i0 = 8

    return (age[snap_f]-age[snap_i])/T_cross, T_cross, snap_i, snap_f

def star_centric_energy(x, v, mp_star, mp_dm, eps, ok, iter=5):
    n = len(x)
    idx = np.arange(n, dtype=int)
    E_out = np.ones(n)*np.inf

    x, v, idx, mp_star = x[ok], v[ok], idx[ok], mp_star[ok]
    dx0, dv0 = 0.0, 0.0

    for i in range(5):
        E = gravitree.binding_energy(x, v, mp_dm, eps, n_iter=1)
        ok = E < 0
        if np.sum(ok) < 0: return None, None, None, False

        v0 = np.average(v[ok], axis=0, weights=mp_star[ok])
        x0 = np.average(x[ok], axis=0, weights=mp_star[ok])
        v = v - v0
        x = x - x0
        dx0, dv0 = dx0 + x0, dv0 + v0

        # Don't actually remove any particles based on that initial
        # velocity guess.
        if i > 0:
            x, v, idx = x[ok], v[ok], idx[ok]
            mp_star = mp_star[ok]

    E_out[idx] = E[ok]

    return E_out, dx0, dv0, True

def energy_distribution():
    out_dir = "../plots/stellar_halo/examples"
    suite = "SymphonyMilkyWayHR"
    #suite = "SymphonyMilkyWay"

    i_host = 0
    sim_dir = symlib.get_host_directory(base_dir, suite, i_host)

    param = symlib.simulation_parameters(suite)
    cosmo = cosmology.setCosmology('', symlib.colossus_parameters(param))
    scale = symlib.scale_factors(sim_dir)
    mp = param["mp"]/param["h100"]
    eps = param["eps"]/param["h100"] * scale
    age = cosmo.age(1/scale - 1)

    rs, hist = symlib.read_rockstar(sim_dir)
    sf, _ = symlib.read_symfind(sim_dir)

    dt_T, T_cross, snap_i, snap_f = crossing_count(age, scale, hist, sf)
    n_peak = hist["mpeak_pre"]/mp
    #ok = (dt_T > 10) & (n_peak > 1e5) & (hist["merger_ratio"] < 0.1)
    ok = (dt_T > 10) & (n_peak > 1e5) & (hist["merger_ratio"] < 0.1)
    #ok = (dt_T > 6) & (n_peak > 1e4) & (hist["merger_ratio"] < 0.1)
    print(suite, i_host, "->", np.sum(ok))
    if np.sum(ok) == 0: return

    r_ratio = 0.01
    ratio_name = "10"
    prof_name = "plummer"

    gal_halo = symlib.GalaxyHaloModel(
        symlib.UniverseMachineMStarFit(),
        symlib.FixedRHalf(r_ratio),
        symlib.PlummerProfile() if prof_name == "plummer" else None,
        symlib.Kirby2013Metallicity(),
        no_scatter=True
    )

    plt.figure()

    part = symlib.Particles(sim_dir)
    mp_star, ranks, mstar, r_half, Fe_H = symlib.tag_stars(sim_dir, gal_halo)
    for i in np.where(ok)[0]:
        if i != 16: continue
        plt.cla()

        dt_T_i = (age - age[snap_i[i]])/T_cross[i]
        #eval_snap = np.searchsorted(dt_T_i, [0, 3, 6])
        eval_snap = np.searchsorted(dt_T_i, [0, 2.5, 5, 7.5])

        m_tot_initial = 0.0
        m_star_initial = 0.0
        E_norm_initial = 0.0

        print("sub_i", i)
        for i_snap, snap  in enumerate(eval_snap):
            p = part.read(snap, halo=i, mode="all")
            mp_dm_i = np.ones(len(p))*mp
            mp_star_i = np.zeros(len(p))

            print("s %3d M %.2g" % (snap, sf["m"][i,snap]))

            mp_star_i[p["smooth"]] = mp_star[i]

            E, dx0, dv0, is_bound = star_centric_energy(
                p["x"]-sf["x"][i,snap], p["v"]-sf["v"][i,snap],
                mp_star_i, mp, eps[snap], p["ok"])
            ok = E < 0
            if m_tot_initial == 0.0:
                m_tot_initial = np.sum(ok)*mp
                m_star_initial = np.sum(mp_star_i[ok])

            E_bins = np.linspace(-7, 0, 75)
            #E_bins = np.linspace(-1.5, 0, 75)
            dE = E_bins[1] - E_bins[0]
            E = E / sf["vmax"][i,eval_snap[0]]**2

            f_norm = 0.001
            E_norm = np.quantile(E, f_norm)
            if E_norm_initial == 0.0:
                E_norm_initial = E_norm
            else:
                #pass
                E = E + (E_norm_initial - E_norm)
            print(E_norm, E_norm_initial)

            #E = E/np.abs(E_norm)
            """
            n_ok = np.sum(ok)
            q = n_norm/n_ok
            if q >= 1: continue
            E_norm = np.quantile(E[ok], q)
            E /= np.abs(E_norm)/2
            """
            #E_norm = np.sort(E)[5000]
            #E /= np.abs(E_norm)
            #E = E / sf["vmax"][i,snap]**2

            print("%.2g %.3g %.3g" % (np.sum(mp_star_i[ok])/m_tot_initial,
                                      np.sum(mp_star_i[ok]), m_tot_initial))

            red_min, red_max = 0.25, 0.75
            black_min, black_max = 0.5, 0.8
            f_norm = len(eval_snap)-1
            f_red = red_min + (red_max - red_min)*i_snap/f_norm
            f_black = black_min + (black_max - black_min)*i_snap/f_norm

            is_cumu = False
            """
            plt.hist(E, bins=E_bins, weights=mp_star_i/m_star_initial/dE,
                     color=pc("r", f_red), lw=3, histtype="step",
                     cumulative=is_cumu)
            plt.hist(E, bins=E_bins, weights=mp_dm_i/m_tot_initial/dE,
                     color=pc("k", f_black), lw=3, histtype="step",
                     cumulative=is_cumu)
            """
            n, _ = np.histogram(E[ok], bins=E_bins)
            n_star, _ = np.histogram(E[ok], bins=E_bins,
                                   weights=mp_star_i[ok]/m_tot_initial/dE)
            n_dm, _ = np.histogram(E[ok], bins=E_bins,
                                     weights=mp_dm_i[ok]/m_tot_initial/dE)

            E_mid = (E_bins[1:] + E_bins[:-1])/2
            ok = n > 4
            #plt.plot(E_mid[ok], np.cumsum(n_dm[ok])*dE, c=pc("k", f_black))
            #plt.plot(E_mid[ok], np.cumsum(n_star[ok])*dE, c=pc("r", f_red))
            plt.plot(E_mid[ok], n_dm[ok], c=pc("k", f_black))
            plt.plot(E_mid[ok], n_star[ok], c=pc("r", f_red))
            #plt.xlim((E_bins[0], E_bins[-1]))

            """
            print("          ", snap)
            print("    dv0   ", dv0)
            print("    masses %.3f %.3g %.3g" % (np.sum(mp_dm_i[ok])/sf["m_raw"][i,snap], sf["m_raw"][i,snap], np.sum(mp_star_i[ok])))

            print(np.quantile(E, 0.001))
            print(np.quantile(E, 0.01))
            print(np.quantile(E, 0.1))
            """
        
        plt.yscale("log")
        plt.ylim(1e-8, 1)
        plt.xlabel(r"$\varepsilon = (E - E_0)/V_{\rm max,infall}^2$")
        plt.ylabel(r"$d(m/m_{\rm dm,infall})/d\varepsilon$")
        plt.savefig("%s/E_distr_%s_h%d_s%d_r%s_%s.png" %
                    (out_dir, suite, i_host, i, ratio_name, prof_name))
        #break

def convergence_mass_comp():
    out_dir = "../plots/stellar_halo/examples"
    suite = "SymphonyMilkyWayHR"
    #suite = "SymphonyMilkyWay"

    i_host = 0
    sim_dir = symlib.get_host_directory(base_dir, suite, i_host)

    param = symlib.simulation_parameters(suite)
    cosmo = cosmology.setCosmology('', symlib.colossus_parameters(param))
    scale = symlib.scale_factors(sim_dir)
    mp = param["mp"]/param["h100"]
    eps = param["eps"]/param["h100"] * scale
    age = cosmo.age(1/scale - 1)

    rs, hist = symlib.read_rockstar(sim_dir)
    sf, _ = symlib.read_symfind(sim_dir)

    dt_T, T_cross, snap_i, snap_f = crossing_count(age, scale, hist, sf)
    n_peak = hist["mpeak_pre"]/mp
    ok = (dt_T > 6) & (n_peak > 1e4) & (hist["merger_ratio"] < 0.1)
    print(suite, i_host, "->", np.sum(ok))
    if np.sum(ok) == 0: return

    r_ratio = 0.01
    ratio_name = "10"
    prof_name = "plummer"

    gal_halo = symlib.GalaxyHaloModel(
        symlib.UniverseMachineMStarFit(),
        symlib.FixedRHalf(r_ratio),
        symlib.PlummerProfile() if prof_name == "plummer" else None,
        symlib.Kirby2013Metallicity(),
        no_scatter=True
    )

    plt.figure()

    part = symlib.Particles(sim_dir)
    mp_star, ranks, mstar, r_half, Fe_H = symlib.tag_stars(sim_dir, gal_halo)
    for i in np.where(ok)[0]:
        if i != 16: continue
        plt.cla()

        eval_snap = np.arange(snap_i[i], snap_f[i]+1, 5)

        m_tot_initial = 0.0
        m_star_initial = 0.0
        E_norm_initial = 0.0

        print("sub_i", i)
        for i_snap, snap  in enumerate(eval_snap):
        
        plt.yscale("log")
        plt.ylim(1e-8, 1)
        plt.xlabel(r"$\varepsilon = (E - E_0)/V_{\rm max,infall}^2$")
        plt.ylabel(r"$d(m/m_{\rm dm,infall})/d\varepsilon$")
        plt.savefig("%s/E_distr_%s_h%d_s%d_r%s_%s.png" %
                    (out_dir, suite, i_host, i, ratio_name, prof_name))
        #break

def main():
    palette.configure(False)
    #energy_distribution()
    convergence_mass_comp()

if __name__ == "__main__": main()
