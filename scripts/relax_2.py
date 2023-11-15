import symlib
import numpy as np
import matplotlib.pyplot as plt
import palette
from palette import pc
import gravitree
from colossus.cosmology import cosmology
from colossus.halo import mass_so
import scipy.stats as stats
import scipy.interpolate as interp

SUITE = "SymphonyMilkyWay"
OUT_DIR = "../plots/fm_fstar"
BASE_DIR = "/sdf/home/p/phil1/ZoomIns"

params = { 'flat': True, 'H0': 70.0, 'Om0': 0.286,
           'Ob0': 0.049, 'sigma8': 0.82, 'ns': 0.95 }
cosmo = cosmology.setCosmology('myCosmo', params)
suffix = "fid3"
        
class MassProfile(object):
    def __init__(self, dx, mp, eps):
        dr = np.sqrt(np.sum(dx**2, axis=1))
        order = np.argsort(dr)
        self.dr = dr[order]
        self.mp = mp
        self.eps = eps

    def Nr(self, r):
        i = np.searchsorted(self.dr, r)
        return i + 1

    def Mr(self, r):
        return self.Nr(r)*self.mp

    def t_orbit(self, r):
        return 7.5 * (r/40)**1.5 * (1e10/self.Mr(r))**0.5

    def t_relax_t_orbit(self, r):
        """ Ludlow et al. (2018)
        """
        N = self.Nr(r)
        eps = self.eps
        #N_term = np.sqrt(200)*N/4
        eps_term = (np.log(r**2/eps**2 + 1) +
                    eps**2-2*r**2/(3*(eps**2+r**2)) -
                    np.log(3/2))**-1
        #rho_term = (N*mp/(4*np.pi/3 * r**3) / self.rho_crit)**(-0.5)
        #return N_term*eps_term*rho_term
        return N/4.0 * eps_term

    def t_relax(self, r):
        return self.t_orbit(r)*self.t_relax_t_orbit(r)

def energy_radii(dx, E, E_eval):
    r = np.sqrt(np.sum(dx**2, axis=1))
    order = np.argsort(E)
    r_E, E_E = r[order], E[order]
    r_r = np.sort(r)

    eval_ok = (E_eval > np.min(E)) & (E_eval < np.max(E))
    r50, rq = np.zeros(len(E_eval)), np.zeros(len(E_eval))

    i_E = np.searchsorted(E_E, E_eval)
    for i in range(len(E_eval)):
        if not eval_ok[i]: continue

        rq[i] = r_r[i_E[i]]
        r50 = np.median(r_E[:i_E[i]])

    return rq, r50, eval_ok

def m_lim_eps(param, scale, h, c, ok, infall_snap):
    # returns a mass ratio, not a mass
    # Factor of 1.284 comes from Mansfield & Avestruz (2021)
    eps = scale*param["eps"]/param["h100"]
    c_infall = h["cvir"][infall_snap]
    rs_infall = h["rvir"][infall_snap]/h["cvir"][infall_snap]
    def f(x):
        return np.log(1+x) - x/(1+x)
    return 1.79/1.284 * ((eps[ok]*c["r50_bound"][ok])/
                         (f(c_infall)*rs_infall**2))

def m_lim_n(param, scale, h, c, ok, infall_snap):
    # Returns a mass ratio, not a mass
    mp = param["mp"]/param["h100"]
    n_infall = h["mvir"][infall_snap]/mp
    return 0.32*(n_infall/1e3)**-0.8 * np.ones(np.sum(ok))

def timescale_comparison(i0):
    i_host = 0
    sim_dir = symlib.get_host_directory(BASE_DIR, SUITE, i_host)

    param = symlib.simulation_parameters(sim_dir)
    scale = symlib.scale_factors(sim_dir)
    z = 1/scale - 1
    age = cosmo.age(z)

    t_dyn = mass_so.dynamicalTime(z, "vir", "crossing")
    mp, eps = param["mp"]/param["h100"], param["eps"]/param["h100"]*scale

    h, hist = symlib.read_subhalos(sim_dir)
    c = symlib.read_cores(sim_dir, suffix=suffix)

    snap0 = hist["first_infall_snap"][i0]
    dt = age - age[snap0]

    part = symlib.Particles(sim_dir)
    p = part.read(snap0, mode="all")[i0]
    dx = p["x"] - h["x"][i0,snap0]
    dr = np.sqrt(np.sum(dx**2, axis=1))
    dv = p["v"] - h["v"][i0,snap0]

    E = gravitree.binding_energy(
        p["x"], dv, mp, eps[snap0], n_iter=3, ok=p["ok"])
    E /= h["vmax"][i0,snap0]**2
    
    bound_infall = (E < 0) & p["ok"]
    prof = MassProfile(dx[bound_infall], mp, eps[snap0])

    t_vir = prof.t_orbit(h["rvir"][i0,snap0])
    t_relax = prof.t_relax(dr[bound_infall])
    
    t_relax_t_orbit = prof.t_relax_t_orbit(dr[bound_infall])

    relaxed_mass = np.zeros(len(scale))
    relaxed_mass_hydro = np.zeros(len(scale))
    f_particle = params["Ob0"] / (params["Om0"] - params["Ob0"])
    for snap in range(snap0, len(scale)):
        if not c["ok"][i0,snap0]: continue
        p_age = age[snap] - age[p["snap"][bound_infall]]
        relaxed_mass[snap] = np.sum(t_relax < p_age)*mp
        relaxed_mass_hydro[snap] = np.sum(t_relax*f_particle < dt[snap])*mp

    fig, ax = plt.subplots()

    ok = c["ok"][i0]
    m_infall = h["mvir"][i0,hist["first_infall_snap"][i0]]
    m_lim_1 = m_lim_eps(param, scale, h[i0], c[i0], ok, snap0)*m_infall
    m_lim_2 = m_lim_n(param, scale, h[i0], c[i0], ok, snap0)*m_infall
    dt /= t_dyn[hist["first_infall_snap"][i0]]

    ax.plot(dt[ok], c["m_bound"][i0,ok], c="k",
            label=r"${\rm subhalo\ mass}$")
    #ax.plot(dt[ok], relaxed_mass[ok], c=pc("r"),
    #        label=r"$m_{\rm mix}\ ({\rm L18})$")
    #ax.plot(dt[ok], relaxed_mass_hydro[ok], "--", c=pc("r"),
    #        label=r"$m_{\rm mix,hydro}\ ({\rm L18})$")
    #ax.plot(dt[ok], m_lim_1, c=pc("o"),
    #        label=r"$m_{\rm lim,\epsilon}\ ({\rm vdBO18})$")
    #ax.plot(dt[ok], m_lim_2, c=pc("b"),
    #        label=r"$m_{{\rm lim},n}\ ({\rm vdBO18})$")
    ax.plot(dt[ok], m_lim_1, c=pc("o"),
            label=r"${\rm force{-}softening\ limit}$")
    ax.plot(dt[ok], m_lim_2, c=pc("b"),
            label=r"${\rm discreteness\ limit}$")
    ax.plot(dt[ok], relaxed_mass[ok], c=pc("r"),
            label=r"${\rm two{-}body\ relaxation\ limit}$")
    
    ax.set_ylabel(r"$m\ (M_\odot)$")
    ax.set_xlabel(r"$(t - t_{\rm infall})/t_{\rm cross}$")
    ax.legend(loc="upper right", fontsize=16)
    ax.set_yscale("log")
    m_min = min(np.min(c["m_bound"][i0,ok])/10,
                np.min(m_lim_1), np.min(m_lim_2))
    if i0 == 15: 
        m_min = 3e6
    ax.set_ylim(m_min, None)
    
    print("%s/mass_limits_%02d_%03d.pdf" % (OUT_DIR, i_host, i0))
    fig.savefig("%s/mass_limits_%02d_%03d.pdf" % (OUT_DIR, i_host, i0))

def last_snap(ok):
    if np.sum(ok) < 2: return -1
    return np.where(ok)[0][-1]

def stack_mass_loss():
    suites = ["SymphonyMilkyWay"]

    all_m, all_t = [], []
    m_min = []
    m_max = []
    t_max = []
    m_disrupt = []
    mt = []
    t_eval = np.linspace(0, 2, 100)
    
    for suite in suites:
        n_host = symlib.n_hosts(suite)

        fig, ax = plt.subplots()

        for i_host in range(n_host):
            sim_dir = symlib.get_host_directory(BASE_DIR, suite, i_host)
            
            param = symlib.simulation_parameters(sim_dir)
            scale = symlib.scale_factors(sim_dir)
            age = cosmo.age(1/scale - 1)
            
            h, hist = symlib.read_subhalos(sim_dir)
            c = symlib.read_cores(sim_dir, suffix=suffix)

            conv_snap = np.zeros(len(h), dtype=int)
            final_snap = np.zeros(len(h), dtype=int)
            conv_snap[0], final_snap[0] = -1, -1

            for i in range(1, len(h)):
                ok = c["ok"][i]
                
                infall_snap = hist["first_infall_snap"][i]
                
                m_infall = c["m_bound"][i,infall_snap]
                m = c["m_bound"][i]/m_infall
                m_lim_1 = m_lim_eps(param, scale, h[i], c[i], ok, infall_snap)
                m_lim_2 = m_lim_n(param, scale, h[i], c[i], ok, infall_snap)

                is_conv = np.zeros(len(m), dtype=bool)
                is_conv[ok] = (m[ok] > m_lim_1) & (m[ok] > m_lim_2)
                is_conv[ok] = m[ok] > m_lim_2
                conv_snap[i] = last_snap(is_conv & ok)
                final_snap[i] = last_snap(ok)

            print(i_host)
            max_snap = h.shape[1] - 1
            
            n_conv = conv_snap - hist["first_infall_snap"]
            ok = ((conv_snap != -1) & (final_snap != max_snap) &
                  (n_conv > 40) & (conv_snap < final_snap) &
                  (hist["merger_ratio"] < 0.1))

            m_conv = c["m_bound"][conv_snap]
            t_conv = age[conv_snap]
            t_start = age[hist["first_infall_snap"]]
            delta_t = t_conv - t_start
            
            for i in range(len(h)):
                if not ok[i] and not (i in [5, 7, 8, 15, 25] and i_host == 0): continue

                m_infall = c["m_bound"][i][hist["first_infall_snap"][i]]
                m_conv = c["m_bound"][i][conv_snap[i]]
                if m_infall/m_conv < 10: continue

                print("%g %g" % (m_infall, m_conv))

                log_m_infall = np.log10(m_infall)
                log_m_conv = np.log10(m_conv)
                delta_log_m = log_m_infall - log_m_conv

                t = (age[c["ok"][i]] - t_start[i]) / delta_t[i]
                m = (np.log10(c["m_bound"][i,c["ok"][i]]) - log_m_infall) / delta_log_m
                color = pc("k") if final_snap[i] < max_snap else pc("r")


                snap = np.arange(c.shape[1], dtype=int)[c["ok"][i]]
                m_max.append(np.max(m[snap < conv_snap[i]]))
                m_min.append(np.min(m[snap < conv_snap[i]]))
                t_max.append(np.max(t))
                m_disrupt.append(m[-1])

                #ax.plot(t, m, c=color, lw=1, alpha=0.1)
                if i == 15 and i_host == 0:
                    ax.plot(t, m, pc("k"))
                
                all_m.append(m)
                all_t.append(t)
                mt.append(eval_mt(m, t, t_eval))
                
            #if i_host > 3: break

    m = np.hstack(all_m)
    t = np.hstack(all_t)
    mt = np.array(mt)


    print(m_disrupt)
    print(np.median(m_disrupt))
    #print(m_max)
    #print(m_min)
    #print(np.median(m_min))
    #print(np.median(t_max))

    m_edges = np.linspace(-2, 0, 30)
    t_edges = np.linspace(0, 2, 30)

    m_mid = (m_edges[1:] + m_edges[:-1])/2 
    t_mid = (t_edges[1:] + t_edges[:-1])/2 

    m_med, _, _ = stats.binned_statistic(t, m, "median", bins=t_edges)
    t_med, _, _ = stats.binned_statistic(m, t, "median", bins=m_edges)

    ax.set_xlim(0, 2)
    ax.set_ylim(-1.9, 0)
    #ax.plot(t_med, m_mid, "o", color=pc("r"))
    #ax.plot(t_mid, m_med, "o", color=pc("b"))


    mt_low = np.quantile(mt, 0.5-0.68/2, axis=0)
    mt_high = np.quantile(mt, 0.5+0.68/2, axis=0)
    mt_med = np.median(mt, axis=0)
    print(mt_med.shape)


    ax.plot(t_eval, mt_med, pc("r"))
    ax.plot([0, 2], [0, -2], "--", lw=1.5, c=pc("r"))
    ax.plot([0, 2], [-1, -1], "--", lw=1, c=pc("a"))
    ax.plot([1, 1], [0, -2], "--", lw=1, c=pc("a"))

    ax.set_xlabel(r"$\Delta t/t_{\rm conv}$")
    ax.set_ylabel(r"$\Delta \log_{10}\,m / (\log_{10}\,m_{\rm infall} - \log_{10}\,m_{\rm conv})$")

    fig.savefig("%s/stacked_mass_loss.pdf" % OUT_DIR)


def eval_mt(m, t, t_eval):
    f = interp.interp1d(t, m, bounds_error=False, fill_value=(m[0], -100))
    return f(t_eval)

def main():
    palette.configure(True)

    #for i in range(1, 40):
    #    timescale_comparison(i)

    timescale_comparison(15)
    #stack_mass_loss()

if __name__ == "__main__": main()
