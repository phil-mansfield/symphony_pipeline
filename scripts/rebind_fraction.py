import numpy as np
import symlib
import matplotlib.pyplot as plt
import palette
import gravitree
from palette import pc

def hbt_plus(x, v, ok_hbt, mp, eps, c):
    idx = np.arange(len(x), dtype=int)
    idx, x, v = idx[ok_hbt], x[ok_hbt], v[ok_hbt]
    x0, v0 = x, v

    dv = np.zeros(3)
    prev_len = len(x)
    for i in range(10):
        mean_v = np.mean(v, axis=0)
        v = v - mean_v
        dv += mean_v
        ok = gravitree.binding_energy(x, v, mp, eps, n_iter=1) < 0
        x, v = x[ok], v[ok]
        if len(x) == prev_len: break
        prev_len = len(x)

    n_bound = len(x)
    if n_bound == 0:
        return 0.0, np.array([0, 0, 0]), ok_hbt, False

    mass = mp*n_bound
    mean_r = np.sqrt(np.sum(np.mean(x, axis=0)**2))
    dx = x - np.mean(x, axis=0)
    r = np.sqrt(np.sum(dx**2, axis=1))
    r_half = np.median(r)
    
    r_c = np.sqrt(np.sum((np.mean(x, axis=0) - c["x"])**2))

    is_intact = r_half < mean_r and (not c["ok"] or r_c < r_half)

    if n_bound*3 >= len(x0):
        return mass, dv, ok_hbt, is_intact
    else:
        cutoff_frac = 3*n_bound/len(x0)
        phi = gravitree.potential_energy(x0, mp, eps)
        cutoff_phi = np.quantile(phi, cutoff_frac)
        ok_idx = idx[phi <= cutoff_phi]
        
        ok_hbt = np.zeros(len(ok_hbt), dtype=bool)
        ok_hbt[ok_idx] = True

        return mass, dv, ok_hbt, is_intact
    

def print_rebind_fractions(i_host):
    suite = "SymphonyMilkyWay"
    base_dir = "/sdf/home/p/phil1/ZoomIns"

    sim_dir = symlib.get_host_directory(base_dir, suite, i_host)

    part_info = symlib.ParticleInfo(sim_dir)
    param = symlib.simulation_parameters(sim_dir)
    h, hist = symlib.read_subhalos(sim_dir)
    h_cmov, hist_cmov = symlib.read_subhalos(sim_dir, comoving=True)
    c = symlib.read_cores(sim_dir)
    scale = symlib.scale_factors(sim_dir)

    #i_subs = [14, 16]
    #i_subs = [11, 15, 18, 19, 30]
    i_subs = [11, 14, 16]
    
    for i in i_subs:
        print("sub", i)
        fig, ax = plt.subplots()
        ok_c, ok_rs = c["ok"][i], c["ok_rs"][i] & c["ok"][i]
        s_ok = np.where(ok_c)[0]
        
        snaps = np.arange(len(scale), dtype=int)
        s_ok = snaps[snaps >= np.min(s_ok)]

        msub_1 = np.zeros(len(scale))
        msub_prev = np.zeros(len(scale))
        prev_ok = None
        mh, ok_h = np.zeros(len(scale)), np.zeros(len(scale), dtype=bool)

        for s in s_ok:
            print("    snap", s)
            ok = symlib.read_particles(part_info, sim_dir, s, "valid")[i]
            if s == s_ok[0]:
                ok_hbt = ok

            x = symlib.read_particles(part_info, sim_dir, s, "x")[i]
            v = symlib.read_particles(part_info, sim_dir, s, "v")[i]
            x = symlib.set_units_x(x, h_cmov[0,s], scale[s], param)
            v = symlib.set_units_v(v, h_cmov[0,s], scale[s], param)

            x -= c["x"][i,s]
            v -= c["v"][i,s]

            a = scale[s]
            mp, eps = param["mp"]/param["h100"], param["eps"]/param["h100"]*a

            is_bound_1 = gravitree.binding_energy(
                x[ok], v[ok], mp, eps, n_iter=1) < 0
            m_hbt, dv_hbt, ok_hbt, is_intact = hbt_plus(
                x+c["x"][i,s], v+c["v"][i,s], ok_hbt, mp, eps, c[i,s])
            mh[s], ok_h[s] = m_hbt, is_intact

            msub_1[s] = mp*np.sum(is_bound_1)

        ok_rs = c["ok_rs"][i]
        ax.plot(scale[ok_rs], h["mvir"][i,ok_rs], pc("r"), label=r"${\rm Rockstar\ +\ consistent-trees}$")
        ax.plot(scale[ok_h], mh[ok_h], pc("o"), label=r"${\rm HBT+}$")
        ax.plot(scale[ok_c], c["m_bound"][i,ok_c], pc("b"), label=r"${\rm This\ paper}$")

        ax.legend(loc="lower right", fontsize=17)
        ax.set_yscale("log")
        ax.set_xlabel(r"$a$")
        ax.set_ylabel(r"$m_{\rm sub}$")
        fig.savefig("../plots/core_tracking/rebind_%d.png" % i)


def main():
    palette.configure(True)
    print_rebind_fractions(0)

if __name__ == "__main__": main()
