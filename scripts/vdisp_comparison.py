import numpy as np
import gravitree
import symlib
import numpy.random as random

try:
    import palette
    palette.configure(False)
    from palette import pc
except:
    pc = lambda x: x

random.seed(1337)
gal_model = symlib.GalaxyHaloModel(
    symlib.UniverseMachineMStarFit(),
    symlib.Jiang2019RHalf(),
    symlib.PlummerProfile(),
    symlib.Kirby2013Metallicity(),
    no_scatter=False
)

sim_dir = symlib.get_host_directory(
    "/sdf/home/p/phil1/ZoomIns", "SymphonyMilkyWayHR", 0
)

def expand_mp_star(mp_star, p):
    out = [None]*len(p)
    for i in range(len(out)):
        out[i] = np.zeros(len(p[i]))
        out[i][p[i]["smooth"]] = mp_star[i]
    return out

def center_of_mass(x, mp):
    w = np.zeros(x.shape)
    for dim in range(3): w[:,dim] = x[:,dim]*mp
    return np.sum(w, axis=0)/np.sum(mp)

def velocity_dispersion(v, mp):
    v0 = center_of_mass(v, mp)
    dv = np.zeros(v.shape)
    for dim in range(3): dv[:, dim] = v[:,dim] - v0[dim]
    
    sigma_sq = 0.0
    for dim in range(3): sigma_sq += np.sum(mp*dv[:,dim]**2)
    sigma_sq /= (3*np.sum(mp))

    return np.sqrt(sigma_sq)

def main():
    # Simulation parameters
    scale = symlib.scale_factors(sim_dir)
    param = symlib.simulation_parameters(sim_dir)
    eps = param["eps"]/param["h100"] * scale
    mp = param["mp"]/param["h100"]

    # I/O
    part = symlib.Particles(sim_dir)
    sf, hist = symlib.read_symfind(sim_dir)
    snap = sf.shape[1]-1
    p = part.read(snap, mode="all")

    # Tag stars
    mp_star, _, _, _, _ = symlib.tag_stars(
        sim_dir, gal_model)

    # Particle masses
    mp_star = expand_mp_star(mp_star, p)
    mp_dm = [mp*np.ones(len(p[i])) for i in range(len(p))]

    for i in range(10, 20+1):
        if not sf["ok"][i,-1]: continue
        print("Subhalo %d" % i)
        print("M*,tot = %.3g  Mdm,tot = %.3g" % (np.sum(mp_star[i]),
                                                 np.sum(mp_dm[i])))

        # Unbound quantities
        x_com_unbound_dm = center_of_mass(p[i]["x"], mp_dm[i])
        v_com_unbound_dm = center_of_mass(p[i]["v"], mp_dm[i])
        v_disp_unbound_dm = velocity_dispersion(p[i]["v"], mp_dm[i])
        x_com_unbound_star = center_of_mass(p[i]["x"], mp_star[i])
        v_com_unbound_star = center_of_mass(p[i]["v"], mp_star[i])
        v_disp_unbound_star = velocity_dispersion(p[i]["v"], mp_star[i])

        # Unbinding
        ok, is_bound = p[i]["ok"], np.zeros(len(p[i]), dtype=bool)
        idx = np.arange(len(p[i]), dtype=int)
        E = gravitree.binding_energy(
            p[i]["x"][ok], p[i]["v"][ok]-sf["v"][i,snap], mp, eps[snap],
            n_iter=10
        )
        is_bound[idx] = E < 0

        # Bound quantities
        x_com_bound_dm = center_of_mass(p[i]["x"][is_bound],
                                        mp_dm[i][is_bound])
        v_com_bound_dm = center_of_mass(p[i]["v"][is_bound],
                                        mp_dm[i][is_bound])
        v_disp_bound_dm = velocity_dispersion(p[i]["v"][is_bound],
                                                mp_dm[i][is_bound])
        x_com_bound_star = center_of_mass(p[i]["x"][is_bound],
                                          mp_star[i][is_bound])
        v_com_bound_star = center_of_mass(p[i]["v"][is_bound],
                                          mp_star[i][is_bound])
        v_disp_bound_star = velocity_dispersion(p[i]["v"][is_bound],
                                                mp_star[i][is_bound])

        print("M*,bound = %.3f, Mdm,bound = %.3g" % 
              (np.sum(mp_star[i][is_bound]), np.sum(mp_dm[i][is_bound])))

        # Half-mass radii
        dx = p[i]["x"] - sf["x"][i,snap]
        r = np.sqrt(np.sum(dx**2, axis=1))
        print("r_50,dm,all = %.2f r_50,dm,bound = %.2f" % 
              (np.median(r), np.median(r[is_bound])))

        print("Position comparison")
        print("  Symfind     ", sf["x"][i,snap])
        print("  bound dm    ", x_com_bound_dm)
        print("  bound star  ", x_com_bound_star)
        print("  unbound DM  ", x_com_unbound_dm)
        print("  unbound star", x_com_unbound_star)
        print("Velocity comparison")
        print("  Symfind     ", sf["v"][i,snap])
        print("  bound dm    ", v_com_bound_dm)
        print("  bound star  ", v_com_bound_star)
        print("  unbound DM  ", v_com_unbound_dm)
        print("  unbound star", v_com_unbound_star)
        print("Velocity dispersion comparison")
        print("  Symfind vmax", sf["vmax"][i,snap])
        print("  bound dm    ", v_disp_bound_dm)
        print("  bound star  ", v_disp_bound_star)
        print("  unbound DM  ", v_disp_unbound_dm)
        print("  unbound star", v_disp_unbound_star)

        print()
        print()


if __name__ == "__main__": main()
