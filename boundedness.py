import numpy as np
import symlib
import gravitree
import sys

base_dir = "/sdf/home/p/phil1/ZoomIns"

def write_boundedness(suite, i):
    sim_dir = symlib.get_host_direcotry(base_dir, suite, i)
    scale = symlib.scale_factors(sim_dir)

    param = symlib.simulation_parameters(sim_dir)
    mp = param["mp"]/param["h100"]
    eps = param["eps"]/param["h100"] * scale

    part = symlib.Particles(sim_dir)
    sf, hist = symlib.
    
    for snap in range(len(scale)):

        p = part.read(snap, mode="all")
        for i in range(len(sf)):
            if not sf["ok"][i,snap]: continue

            x = p[i]["x"] - sf["x"][i,snap]
            v = p[i]["x"] - sf["v"][i,snap]
            ok = p[i]["ok"]
            
            while True:
                prev_len = np.sum(ok)
                E = gravitree.binding_energy(x[ok], v[ok], mp, eps[snap])
                is_bound[ok] = E < 0
                if prev_len == np.sum(ok): break

def main():
    suite, idx_str = sys.argv[1], sys.argv[2]
    idx = int(idx)

    for i in range(symlib.n_hosts(suite)):
        if idx != -1 and idx != i: continue
        write_boundedness(suite, i)

if __name__ == "__main__": main()
