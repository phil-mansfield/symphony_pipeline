import numpy as np
import symlib
import convert_core_catalogue
import numpy.random as ranomd
import gravitree
import os
import sys
import scipy.signal as signal
import scipy.special as special

def make_model(use_um, r_name):
    if r_name == "jiang":
        radius_model = Jiang2019RHalf()
    elif len(r_name) > 2 and r_name[:2] == "r=":
        radius_model = symlib.FixedRHalf(float(r_name[2:]))
    else:
        raise ValueError("Unrecognized r_name, %s" % r_name)

    if use_um:
        m_star_model = symlib.StellarMassModel(
            symlib.UniverseMachineMStar(),
            symlib.UniverseMachineSFH()
        )
    else:
        m_star_model = symlib.StellarMassModel(
            symlib.UniverseMachineMStarFit(),
            symlib.DarkMatterSFH()
        )
    return symlib.GalaxyHaloModel(
        m_star_model,
        symlib.ProfileModel(
            radius_model,
            symlib.PlummerProfile()
        ),
        symlib.MetalModel(
            symlib.Kirby2013Metallicity(),
            symlib.Kirby2013MDF(model_type="gaussian"),
            symlib.FlatFeHProfile(),
	    symlib.GaussianCoupalaCorrelation()
        )   
    )


models = {
    "um": make_model(True, "jiang"),
    "um_fit": make_model(False, "jiang"),
    "r=0.005": make_model(False, "r=0.005")
    "r=0.008": make_model(False, "r=0.008"),
    "r=0.015": make_model(False, "r=0.015"),
    "r=0.025": make_model(False, "r=0.025"),
    "r=0.05": make_model(False, "r=0.05"),
    "r=0.1": make_model(False, "r=0.1"),
    "r=0.2": make_model(False, "r=0.2"),
    "r=0.4": make_model(False, "r=0.4")
}

model_names = sorted(models.keys())

def v_disp(v, mp):
    mp = mp * np.ones(len(v))
    sigma_sq_3d = 0.0
    for dim in range(3):
        sigma_sq_3d += np.sum(v[:,dim]**2*mp)
    sigma_sq_3d /= (3*np.sum(mp))
    
    return np.sqrt(sigma_sq_3d)

def capped_rel_max(x, x_min, debug=False):
    i = signal.argrelextrema(x, np.greater)[0]
    if len(i) <= 1:
        return x_min
    else:
        return np.max(x[i[1:]])


def vmax(x, mp, eps):
    if len(x) == 0: return 0, 0
    r = np.sqrt(np.sum(x**2, axis=1))
    order = np.argsort(r)
    r = r[order]

    m = np.cumsum(mp*np.ones(len(r)))
    v = 2.074e-3 * m**0.5 * (r**2 + eps**2)**-0.25
    vmax_min = 2.074e-3 * mp**0.5 * eps**-0.5
    vmax = capped_rel_max(v, vmax_min)

    h = eps / 0.357
    A, beta = 0.172, -0.522
    v_debias = v / (1 - np.exp(-(A*h/np.sqrt(r**2 + eps**2))**beta))

    vmax_debias = capped_rel_max(v_debias, vmax_min)
    return vmax, vmax_debias

def m23_S_moments(n_peak):
    z90 = 1.2816
    log_n = np.log10(n_peak)
    p9 = -0.3473 -0.3756*log_n
    p5 = -0.5054 -0.5034*log_n
    p1 = 0.0526 - 0.8121*log_n

    return p1, p5, p9

def m23_S(n_peak, mu):
    p1, p5, p9 = m23_S_moments(n_peak)

    log_mu = np.log10(mu)
    S = np.zeros(len(mu))
    low = log_mu < p5

    d_high = (log_mu - p5)/(p9 - p5)
    d_low = (log_mu - p5)/(p5 - p1)

    S[low] = (1+special.erf(d_low[low]*1.2816/np.sqrt(2)))/2
    S[~low] = (1+special.erf(d_high[~low]*1.2816/np.sqrt(2)))/2

    return S

def get_m23_weight(n_peak, mu):
    return 1/m23_S(n_peak, mu)

def get_m23_v_conv_lim(npeak):
    x = np.log10(npeak)
    b2, b1, b0 = -0.01853, 0.3861, 1.6597
    return 10**(x*x*b2 + x*b1 + b0)
    
def galaxy_catalog(sim_dir, i_host, model_names):

    n_model

    gal_dir = os.path.join(sim_dir, "galaxies")
    if not os.path.exists(gal_dir):
        os.makedirs(gal_dir)
        
    scale = symlib.scale_factors(sim_dir)
    param = symlib.simulation_parameters(sim_dir)
    eps = param["eps"]/param["h100"] * scale
    mp = param["mp"]/param["h100"]

    print("Starting I/O")
    part = symlib.Particles(sim_dir)
    sf, hist = symlib.read_symfind(sim_dir)

    print("Starting to tag")
    mp_star, _, m_star_i, r_half_i, _ = symlib.tag_stars(
        sim_dir, models[model_name])
    print("Tagging done")

    r_half_i *= (1/0.5**(2/3) - 1)**-0.5
    r_half_i[0] = 0

    r_half = np.zeros(sf.shape, dtype=np.float32)
    m_star = np.zeros(sf.shape, dtype=np.float32)
    x0_star = np.zeros(sf.shape + (3,), dtype=np.float32)
    v0_star = np.zeros(sf.shape + (3,), dtype=np.float32)
    m_dyn = np.zeros(sf.shape, dtype=np.float32)
    v_disp_3d_star = np.zeros(sf.shape, dtype=np.float32)
    v_disp_3d_dm = np.zeros(sf.shape, dtype=np.float32)
    vmax_dm = np.zeros(sf.shape, dtype=np.float32)
    vmax_dm_debias = np.zeros(sf.shape, dtype=np.float32)
    m23_weight = np.zeros(sf.shape, dtype=np.float32)
    m23_m_conv = np.zeros(sf.shape, dtype=bool)
    m23_v_conv = np.zeros(sf.shape, dtype=bool)

    m_star_i = np.asarray(m_star_i, dtype=np.float32)
    r_half_i = np.asarray(r_half_i, dtype=np.float32)    

    for snap in range(len(scale)):
        if np.sum(sf["ok"][1:,snap]) == 0: continue
        if snap % 10 == 0: print("   ", snap)

        p = part.read(snap, mode="all")

        ok = sf["ok"][:,snap]
        m23_weight[ok,snap] = get_m23_weight(
            hist["mpeak"][ok]/mp,
            sf["m"][ok,snap]/hist["mpeak"][ok]
        )

        m, npeak = sf["m"][:,snap], hist["mpeak"]/mp
        m23_m_conv[:,snap] = m > get_m23_v_conv_lim(8*npeak)*mp
        m23_v_conv[:,snap] = m > get_m23_v_conv_lim(npeak)*mp

        for i in range(1, len(sf)):
            if not sf["ok"][i,snap]: continue

            x, v, ok = p[i]["x"], p[i]["v"], p[i]["ok"]
            x_all = x - sf["x"][i,snap]
            x, v = x - sf["x"][i,snap], v - sf["v"][i, snap]
            smooth = p[i]["smooth"]
            is_bound = np.zeros(len(x), dtype=bool)
            idx = np.arange(len(x))[ok]

            E = gravitree.binding_energy(x[ok], v[ok], mp, eps[snap],
                                         n_iter=4)
            is_bound[idx] = E < 0

            mp_star_i = mp_star[i][is_bound[smooth]]
            x_star = x[is_bound & smooth]
            v_star = v[is_bound & smooth]
            v_star_w, x_star_w = np.copy(v_star), np.copy(x_star)
            for dim in range(3):
                x_star_w[:,dim] *= mp_star_i
                v_star_w[:,dim] *= mp_star_i

            if np.sum(mp_star_i) > 0:
                x0_star[i,snap] = np.sum(x_star_w, axis=0)/np.sum(mp_star_i)
                v0_star[i,snap] = np.sum(v_star_w, axis=0)/np.sum(mp_star_i)

            r = np.sqrt(np.sum((x - x0_star[i,snap])**2, axis=1))
            r = r[is_bound & smooth]
            order = np.argsort(r)
            r = r[order]
            mp_star_i = mp_star_i[order]

            x_all = x_all - x0_star[i,snap]
            r_all = np.sum(x_all**2, axis=1)

            m_star[i,snap] = np.sum(mp_star_i)
            c_mass = np.cumsum(mp_star_i)
            if len(c_mass) != 0:
                r_half[i,snap] = r[np.searchsorted(c_mass, c_mass[-1]/2)]
                m_dyn[i,snap] = np.sum(r_all < r_half[i,snap])*mp

            x0_star[i,snap] += sf["x"][i,snap]
            v0_star[i,snap] += sf["v"][i,snap]

            v_disp_3d_dm[i,snap] = v_disp(v[is_bound], mp)
            v_disp_3d_star[i,snap] = v_disp(v_star, mp_star_i)
            vmax_dm[i,snap], vmax_dm_debias[i,snap] = vmax(
                x[is_bound], mp, eps[snap])        


    file_name = os.path.join(gal_dir, "%s.dat" % model_name)
    n_snap, n_halo = sf.shape
    with open(file_name, "w+") as fp:
        m_star_i.tofile(fp)
        r_half_i.tofile(fp)
        m_star.tofile(fp)
        r_half.tofile(fp)
        m_dyn.tofile(fp)
        x0_star.tofile(fp)
        v0_star.tofile(fp)
        v_disp_3d_dm.tofile(fp)
        v_disp_3d_star.tofile(fp)
        vmax_dm.tofile(fp)
        vmax_dm_debias.tofile(fp)
        m23_weight.tofile(fp)
        m23_m_conv.tofile(fp)
        m23_v_conv.tofile(fp)
        sf["ok"].tofile(fp)

def main():
    config_name, idx_str, model_names = sys.argv[1], sys.argv[2], sys.argv[3:]
    print("Running models %s on config %s, galaxy %d" %
          (model_names, config_name, idx_str))
    target_idx = int(idx_str)

    sim_dirs = convert_core_catalogue.get_sim_dirs(config_name)
    n_host = len(sim_dirs)
    
    for i_host in range(n_host):
        if target_idx == -1 or i_host == target_idx:
            galaxy_catalog(sim_dirs[i_host], i_host, model_names)

if __name__ == "__main__": main()
