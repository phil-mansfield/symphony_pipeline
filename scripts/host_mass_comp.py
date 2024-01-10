import numpy as np
import symlib
import gravitree

def main():
    for i_host in range(45):
        sim_dir = symlib.get_host_directory("/sdf/home/p/phil1/ZoomIns", "SymphonyMilkyWay", i_host)

        r, _ = symlib.read_rockstar(sim_dir)
        part = symlib.Particles(sim_dir)
        p = part.read(235)[0]

        p = p[p["ok"]]
        dist = np.sqrt(np.sum(p["x"]**2, axis=1))

        param = symlib.simulation_parameters(sim_dir)
        mp = param["mp"] / param["h100"]
        eps = param["eps"] / param["h100"]

        #is_bound = gravitree.binding_energy(p["x"], p["v"], mp, eps,n_iter=3)<0
        in_rvir = dist < r[0,-1]["rvir"]

        #print("%g %g %g" % (r[0,-1]["m"], mp*np.sum(in_rvir),
        #                    mp*np.sum(in_rvir & is_bound)))
        print("%2d %8.3g %8.3g %.3f" % (i_host, r[0,-1]["m"], mp*np.sum(in_rvir),
                                  mp*np.sum(in_rvir)/r[0,-1]["m"]))

if __name__ == "__main__": main()
