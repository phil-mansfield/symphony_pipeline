import numpy as np
import matplotlib.pyplot as plt
import palette
from palette import pc
import cache_stars
import symlib
from colossus.cosmology import cosmology

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

def main():
    palette.configure(False)
    example_stellar_mass_loss("SymphonyMilkyWayHR", 0, 16)

if __name__ == "__main__": main()
