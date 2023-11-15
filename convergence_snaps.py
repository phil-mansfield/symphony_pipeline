import numpy as np
import find_infall_cores
from colossus.cosmology import cosmology
import sys
import symlib
import gravitree
import os.path as path
import struct

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
        eps_term = (np.log(r**2/eps**2 + 1) +
                    eps**2-2*r**2/(3*(eps**2+r**2)) -
                    np.log(3/2))**-1
        return N/4.0 * eps_term

    def t_relax(self, r):
        return self.t_orbit(r)*self.t_relax_t_orbit(r)

def m_lim_eps(param, scale, h, c, infall_snap):
    # returns a mass ratio, not a mass
    # Factor of 1.284 comes from Mansfield & Avestruz (2021)
    eps = scale*param["eps"]/param["h100"]
    c_infall = h["cvir"][infall_snap]
    rs_infall = h["rvir"][infall_snap]/h["cvir"][infall_snap]
    def f(x):
        return np.log(1+x) - x/(1+x)
    return 1.79/1.284 * ((eps*c["r50_bound"])/
                         (f(c_infall)*rs_infall**2))

def m_lim_n(param, scale, h, c, infall_snap):
    # Returns a mass ratio, not a mass
    mp = param["mp"]/param["h100"]
    n_infall = h["mvir"][infall_snap]/mp
    return 0.32*(n_infall/1e3)**-0.8 * np.ones(len(c))

def mpeak_loss(m, ok):
    out = np.zeros(len(m))
    idx = np.where(ok)[0][::-1]
    mpeak = m[idx[0]]
    for i in idx:
        if m[i] > mpeak: mpeak = m[i]
        out[i] = mpeak
    return out

def convergence_snap(m, m_lim, ok):
    snaps = np.where(ok & (m < m_lim))[0]
    if len(snaps) == 0:
        return -1
    else:
        return snaps[0]

def convergence_snaps(sim_dir):
    param = symlib.simulation_parameters(sim_dir)
    c_param = symlib.colossus_parameters(param)
    cosmo = cosmology.setCosmology("", c_param)

    scale = symlib.scale_factors(sim_dir)
    mp, eps = param["mp"]/param["h100"], param["eps"]/param["h100"]*scale
    age = cosmo.age(1/scale - 1)

    h, hist = symlib.read_subhalos(sim_dir)
    c = symlib.read_cores(sim_dir)

    conv_disc = np.ones(len(c), dtype=int)*-1
    conv_eps = np.ones(len(c), dtype=int)*-1
    conv_relax = np.ones(len(c), dtype=int)*-1
    conv_relax_hydro = np.ones(len(c), dtype=int)*-1
    snap_disrupt = np.ones(len(c), dtype=int)*-1
    snap_disrupt_rs = np.ones(len(c), dtype=int)*-1

    part = symlib.Particles(sim_dir)
    f_particle = param["Ob0"] / (param["Om0"] - param["Ob0"])

    for snap in range(c.shape[1]):
        idx = np.where(hist["first_infall_snap"] == snap)[0]
        if len(idx) == 0: continue
        print(snap, idx)

        p = part.read(snap)

        for i in idx:
            dx = p[i]["x"] - h["x"][i,snap]
            dr = np.sqrt(np.sum(dx**2, axis=1))
            dv = p[i]["v"] - h["v"][i,snap]
            E = gravitree.binding_energy(
                dx, dv, mp, eps[snap], n_iter=3, ok=p[i]["ok"])
            bound = E < 0

            prof = MassProfile(dx[bound], mp, eps[snap])
            t_relax = prof.t_relax(dr[bound])

            m_relax, m_relax_hydro = np.zeros(len(scale)), np.zeros(len(scale))
            for s in range(snap, len(scale)):
                p_age = age[s] - age[p[i]["snap"][bound]]
                m_relax[s] = np.sum(t_relax < p_age)*mp
                m_relax_hydro[s] = np.sum(t_relax*f_particle<p_age)*mp
                
            m_infall = c["m_bound"][i,snap]
            m = c["m_bound"][i]
            ok = c["ok"][i]
            m_disc = m_infall * m_lim_n(param, scale, h[i], c[i], snap)
            m_eps = m_infall * m_lim_eps(param, scale, h[i], c[i], snap)
            
            conv_disc[i] = convergence_snap(m, m_disc, ok)
            conv_eps[i] = convergence_snap(m, m_eps, ok)
            conv_relax[i] = convergence_snap(m, m_relax, ok)
            conv_relax_hydro[i] = convergence_snap(m, m_relax_hydro, ok)
            
            ok_snaps = np.where(c["ok"][i])[0]
            if len(ok_snaps) == 0:
                pass
            else:
                snap_disrupt[i] = ok_snaps[-1]+1

            ok_snaps_rs = np.where(h["ok"][i] & c["ok_rs"][i])[0]
            if len(ok_snaps_rs) == 0:
                pass
            else:
                snap_disrupt_rs[i] = ok_snaps_rs[-1]+1

    snap_disrupt[snap_disrupt == len(scale)] = -1
    snap_disrupt_rs[snap_disrupt_rs == len(scale)] = -1

    write_convergence_limits(sim_dir, conv_disc, conv_eps,
                             conv_relax, conv_relax_hydro,
                             snap_disrupt, snap_disrupt_rs)
    h, hist = symlib.read_subhalos(sim_dir)
    print(hist["disrupt_snap"][:20])
    print(hist["disrupt_snap_rs"][:20])
    print(hist["conv_snap_discrete"][:20])
    print(hist["conv_snap_eps"][:20])
    print(hist["conv_snap_relax"][:20])
    print(hist["conv_snap_relax_hydro"][:20])

def write_convergence_limits(sim_dir, conv_disc, conv_eps, conv_relax, conv_relax_hydro, snap_disrupt, snap_disrupt_rs):
    file_name = path.join(sim_dir, "convergence_limits.dat")
    with open(file_name, "wb+") as fp:
        fp.write(struct.pack("q", len(conv_disc)))
        np.asarray(conv_disc, dtype=np.int32).tofile(fp)
        np.asarray(conv_eps, dtype=np.int32).tofile(fp)
        np.asarray(conv_relax, dtype=np.int32).tofile(fp)
        np.asarray(conv_relax_hydro, dtype=np.int32).tofile(fp)
        np.asarray(snap_disrupt, dtype=np.int32).tofile(fp)
        np.asarray(snap_disrupt_rs, dtype=np.int32).tofile(fp)

def main():
    config_name, idx_str = sys.argv[1], sys.argv[2]
    target_idx = int(idx_str)

    sim_dirs = find_infall_cores.get_sim_dirs(config_name)

    for host_i in range(len(sim_dirs)):
        if host_i != target_idx and target_idx != -1: continue
        sim_dir = sim_dirs[host_i]
        if sim_dir[-1] == "/": sim_dir = sim_dir[:-1]

        convergence_snaps(sim_dir)

if __name__ == "__main__": main()
