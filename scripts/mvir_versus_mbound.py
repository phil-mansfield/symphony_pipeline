import numpy as np
import matplotlib.pyplot as plt
import gravitree
import symlib

OUT_DIR = "../plots/core_tracking"
BASE_DIR = "/sdf/home/p/phil1/ZoomIns"

def main():
    ratio = []
    for i_host in range(symlib.n_hosts("SymphonyMilkyWay")):
        sim_dir = symlib.get_host_directory(
            BASE_DIR, "SymphonyMilkyWay", i_host)

        param = symlib.simulation_parameters(sim_dir)
        mp, eps = param["mp"]/param["h100"], param["eps"]/param["h100"]
        sy, hist = symlib.read_symfind(sim_dir)
        rock, _ = symlib.read_rockstar(sim_dir)
        part = symlib.Particles(sim_dir)

        host = rock[0,-1]
        sub = sy[1:,-1]
        p_sub = part.read(235)[1:]

        rho_vir = host["m"] / (4*np.pi/3 * host["rvir"]**3)

        for i in range(len(sub)):
            if not sub["ok"][i]: continue

            p = p_sub[i]
            p = p[p["ok"]]

            p["x"] -= sub["x"][i]
            r = np.sqrt(np.sum(p["x"]**2, axis=1))
            vol = 4*np.pi/3 * r**3
            p["v"] -= sub["v"][i]
    
            E = gravitree.binding_energy(
                p["x"], p["v"], mp, eps, n_iter=3,
            )
            ok = E < 0

            vol = np.sort(vol[ok])
            rho = np.arange(len(vol))*mp / vol
            mvir = np.sum(rho > rho_vir)*mp

            ratio.append(mvir/(mp*np.sum(ok)))
        
    print("%.3f" % np.quantile(ratio, 0.5-0.68/2))
    print("%.3f" % np.quantile(ratio, 0.5))
    print("%.3f" % np.quantile(ratio, 0.5+0.68/2))
        
    
if __name__ == "__main__": main()
