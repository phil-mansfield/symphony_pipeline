import numpy as np
import symlib

R_BINS = 10**np.linspace(-4, 0, 121)
out_dir = "/fs/ddn/sdf/group/kipac/g/cosmo/ki21/phil1/data/subhalo_profiles"

def save_table(suite, host):
    print(suite, host)
    
    base_dir = "/sdf/home/p/phil1/ZoomIns"

    sim_dir = symlib.get_host_directory(base_dir, suite, host)
    param = symlib.simulation_parameters(sim_dir)
    
    sf, hist = symlib.read_symfind(sim_dir)
    rs, hist = symlib.read_rockstar(sim_dir)
    profs_u = np.zeros((len(hist), len(R_BINS)-1), dtype=int)
    profs_b = np.zeros((len(hist), len(R_BINS)-1), dtype=int)
    scale = symlib.scale_factors(sim_dir)
    
    zs = np.array([0, 1, 2])
    snaps = [np.argmin(np.abs(scale - (1/(z+1)))) for z in zs]

    part = symlib.Particles(sim_dir)
    for i_snap, snap in enumerate(snaps):
        p = part.read(snap, mode="current")

        for i in range(len(hist)):
            if i > 0 and not sf["ok"][i,snap]: continue

            dv = p[i]["v"] - sf["v"][i,snap]
            dx = p[i]["x"] - sf["x"][i,snap]
            r = np.sqrt(np.sum(dx**2, axis=1))/rs["rvir"][i,snap]
            order = np.argsort(r)
            
            ke = np.sum(dv**2, axis=1)/2
            ok = np.ones(len(ke), dtype=bool)
            for _ in range(3):
                _, vmax, pe, _ = symlib.profile_info(
                    param, dx, ok=ok)
                E = ke + pe*vmax**2
                ok = E < 0
            
            profs_u[i,:] = np.histogram(r, bins=R_BINS)[0]
            profs_b[i,:] = np.histogram(r[ok], bins=R_BINS)[0]
            
        np.save("%s/r_bins.npy" % out_dir, R_BINS)
        np.save("%s/prof_%s_h%02d_z%d_u.npy" %
                (out_dir, suite, host, zs[i_snap]), profs_u)
        np.save("%s/prof_%s_h%02d_z%d_b.npy" %
                (out_dir, suite, host, zs[i_snap]), profs_b)
        

def main():
    suites = [
        "MWest"
        #"SymphonyLMC",
        #"SymphonyMilkyWay",
        #"SymphonyGroup",
        #"SymphonyLCluster"
    ]
    for suite in suites:
        for i in range(symlib.n_hosts(suite)):
            save_table(suite, i)
    
if __name__ == "__main__": main()
