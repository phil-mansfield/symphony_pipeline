import numpy as np
import matplotlib.pyplot as plt
import symlib
import palette
import gravitree
import time
from palette import pc
import struct
import numpy.random as random
import sys

base_dir = "/sdf/home/p/phil1/ZoomIns"
#base_dir = "/oak/stanford/orgs/kipac"
suite = "SymphonyMilkyWay"

def plot_basic_mass_loss():
    sim_dir = symlib.get_host_directory(base_dir, suite, 0)
    param = symlib.simulation_parameters(sim_dir)
    scale = symlib.scale_factors(sim_dir)
    mp, eps = param["mp"]/param["h100"], param["eps"]*scale/param["h100"]

    f_rvir = [0.01, 0.02, 0.03, 0.04, 0.05]
    gal_halo_models = [None]*len(f_rvir)
    mp_stars = [None]*len(f_rvir)
    m_stars = [None]*len(f_rvir)
    r_halfs = [None]*len(f_rvir)

    target_subs = np.array([i0])

    for i in range(len(gal_halo_models)):
        gal_halo_models[i] = symlib.GalaxyHaloModel(
            symlib.UniverseMachineMStarFit(),
            symlib.FixedRHalf(f_rvir[i]),
            symlib.PlummerProfile(),
            symlib.Kirby2013Metallicity(),
            no_scatter=False
        )

        mp_stars[i], _, m_stars[i], r_halfs[i], _ = symlib.tag_stars(
            sim_dir, gal_halo_models[i], target_subs=target_subs)

        mp_stars[i], m_stars[i], r_halfs[i] = mp_stars[i][i0], m_stars[i][i0], r_halfs[i][i0]

    h, hist = symlib.read_subhalos(sim_dir)
    c = symlib.read_cores(sim_dir)

    part = symlib.Particles(sim_dir)
    mstar = np.zeros((len(h[i0]), len(f_rvir)))
    for snap in range(len(h[i0])):
        if snap % 10 == 0: print(snap)
        if not c["ok"][i0,snap]: continue

        p = part.read(snap, mode="all")[i0]
        dx = p["x"] - c["x"][i0,snap]
        dv = p["v"] - c["v"][i0,snap]

        E = gravitree.binding_energy(dx, dv, mp, eps[i0],
                                     n_iter=10, ok=p["ok"])
        is_bound = E <= 0
        star_is_bound = is_bound[p["smooth"]]

        for i in range(len(f_rvir)):
            mstar[snap,i] = np.sum(mp_stars[i][star_is_bound])

    ok = c["ok"][i0]
    colors = [pc("r"), pc("o"), pc("g"), pc("b"), pc("p")]
    i_start = np.where(ok)[0][0]
    for i in range(len(f_rvir)):
        for fig in range(3):
            plt.figure(fig)
            plt.plot(c["m_bound"][i0,ok]/c["m_bound"][i0,i_start],
                     mstar[ok,i]/mstar[i_start,i],
                     colors[i],
                     label=r"$f_{\rm Rvir}=%.02f$" % f_rvir[i])
    plt.figure(1)
    plt.xscale("log")
    plt.figure(2)
    plt.xscale("log")
    plt.yscale("log")

    for fig in range(3):
        plt.figure(fig)
        plt.xlabel(r"$m/m_{\rm infall}$")
        plt.ylabel(r"$m_\star/m_{\rm \star,infall}$")
        plt.legend(loc="lower right")
        plt.savefig("../plots/fm_fstar/example_fm_fstar_%d.png" % fig)

def rhalf(dx, mp_star, is_bound):
    dx, mp_star = dx[is_bound], mp_star[is_bound]
    dr = np.sqrt(np.sum(dx**2, axis=1))
    order = np.argsort(dr)
    m_sum = np.cumsum(mp_star[order])
    i_mid = np.searchsorted(m_sum, m_sum[-1]/2)
    return dr[order][i_mid]

def mass_loss_table(sim_dir, gal_halo_model):
    param = symlib.simulation_parameters(sim_dir)
    scale = symlib.scale_factors(sim_dir)
    mp = param["mp"]/param["h100"]
    eps = param["eps"]*scale/param["h100"]*scale

    mp_star, _, m_stars, r_half, _ = symlib.tag_stars(sim_dir, gal_halo_model)

    h, hist = symlib.read_subhalos(sim_dir)
    c = symlib.read_cores(sim_dir)

    t0 = time.time()
    part = symlib.Particles(sim_dir)
    mstar = np.zeros(c.shape)
    rhalf_star = np.zeros(c.shape)
    for snap in range(len(scale)):
        if np.sum(c["ok"][:,snap]) == 0: continue
        print("%3d %5d %.2f" % (snap, np.sum(c["ok"][:,snap]), time.time()-t0))

        p = part.read(snap, mode="all")

        for i in range(1, len(c)):
            if not c[i,snap]["ok"]: continue

            dx = p[i]["x"] - c["x"][i,snap]
            dv = p[i]["v"] - c["v"][i,snap]

            E = gravitree.binding_energy(dx, dv, mp, eps[snap],
                                         n_iter=3, ok=p[i]["ok"])
            is_bound = E <= 0
            star_is_bound = is_bound[p[i]["smooth"]]
            mstar[i,snap] = np.sum(mp_star[i][star_is_bound])
            rhalf_star[i,snap] = rhalf(
                dx[p[i]["smooth"]], mp_star[i], star_is_bound)
            
    return  c["ok"], mstar, rhalf_star, c["m_bound"], c["r50_bound"]
    
            
def write_mass_loss_table(file_name, ok, m_star, rhalf_star, m_dm, rhalf_dm):
    with open(file_name, "wb+") as fp:
        fp.write(struct.pack("ii", ok.shape[0], ok.shape[1]))
        ok.tofile(fp)
        np.asarray(m_star, dtype=np.float32).tofile(fp)
        np.asarray(rhalf_star, dtype=np.float32).tofile(fp)
        np.asarray(m_dm, dtype=np.float32).tofile(fp)
        np.asarray(rhalf_dm, dtype=np.float32).tofile(fp)

def read_mass_loss_table(file_name):
    with open(file_name, "rb") as fp:
        n_halo = struct.unpack("i", fp.read(4))[0]
        n_snap = struct.unpack("i", fp.read(4))[0]
        n_tot = n_halo*n_snap
        shape = (n_halo, n_snap)
        
        ok = np.fromfile(fp, np.bool, n_tot)
        m_star = np.fromfile(fp, np.float32, n_tot)
        rhalf_star = np.fromfile(fp, np.float32, n_tot)
        m_dm = np.fromfile(fp, np.float32, n_tot)
        rhalf_dm = np.fromfile(fp, np.float32, n_tot)

        ok = ok.reshape(shape)
        m_star = m_star.reshape(shape)
        rhalf_star = rhalf_star.reshape(shape)
        m_dm = m_dm.reshape(shape)
        rhalf_dm = rhalf_dm.reshape(shape)

    return ok, m_star, rhalf_star, m_dm, rhalf_dm

def main():
    palette.configure(False)
    if len(sys.argv) > 1:
        target = int(sys.argv[1])
    else:
        target = 0

    gal_halo_model = symlib.GalaxyHaloModel(
        symlib.UniverseMachineMStarFit(),
        symlib.Jiang2019RHalf(),
        symlib.PlummerProfile(),
        symlib.Kirby2013Metallicity(),
        no_scatter=False
    )
    
    for host_i in range(symlib.n_hosts(suite)):
        if host_i != target and target != -1: continue
        random.seed(host_i)

        print(host_i)
        continue

        sim_dir = symlib.get_host_directory(base_dir, suite, host_i)
        ok, m_star, rhalf_star, m_dm, rhalf_dm = mass_loss_table(
            sim_dir, gal_halo_model)

        if suite == "SymphonyMilkyWay":
            file_name = "/sdf/group/kipac/g/cosmo/ki21/phil1/data/fstar_fm_tables/mass_loss_%02d.dat" % host_i
        elif suite == "SymphonyMilkyWayHR":
            file_name = "/sdf/group/kipac/g/cosmo/ki21/phil1/data/fstar_fm_tables/mass_loss_hr_%02d.dat" % host_i
        write_mass_loss_table(file_name, ok, m_star, rhalf_star, m_dm, rhalf_dm)
        ok, m_star, rhalf_star, m_dm, rhalf_dm = read_mass_loss_table(file_name)



if __name__ == "__main__": main()
