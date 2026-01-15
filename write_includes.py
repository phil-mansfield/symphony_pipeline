import numpy as np
import symlib
import find_infall_cores
import sys
import gravitree
import time

def write_includes(sim_dir):
    print(sim_dir)

    t0 = time.time()
    
    sf, hist = symlib.read_symfind(sim_dir)
    rs, hist = symlib.read_rockstar(sim_dir)
    scale = symlib.scale_factors(sim_dir)
    n_snap = len(scale)
    param = symlib.simulation_parameters(sim_dir)
    mp = param["mp"]/param["h100"]
    eps = param["eps"]/param["h100"]*scale
    
    symlib.write_include_header(sim_dir)
    part = symlib.Particles(sim_dir)

    t1 = time.time()
    print("Setup time: %5.2f s" % (t1 - t0))
    
    for snap in range(n_snap):
        t0 = time.time()
        
        print("    snap %4d" % snap)
        p = part.read(snap, mode="current")
        t2 = time.time()
        
        E = [None]*len(p)
        E_sph = [None]*len(p)
        m_enc = [None]*len(p)
        
        for i in range(len(sf)):
            if len(p[i]) == 0 or (not sf["ok"][i,snap] and not rs["ok"][i,snap]):
                E[i] = np.zeros(len(p[i]))
                E_sph[i] = np.zeros(len(p[i]))
                m_enc[i] = np.zeros(len(p[i]))
                continue

            if sf["ok"][i,snap]:
                dx = p[i]["x"] - sf["x"][i,snap]
                dv = p[i]["v"] - sf["v"][i,snap]
            else:
                dx = p[i]["x"] - rs["x"][i,snap]
                dv = p[i]["v"] - rs["v"][i,snap]

            # This can be sped up with binsort if needed
            rmax, vmax, pe_sph_vmax2, order = symlib.profile_info(param, dx)
                
            pe_sph = pe_sph_vmax2*vmax*vmax
            ke = np.sum(dv*dv, axis=1)/2
            E_sph[i] = pe_sph + ke
            m_enc[i] = mp*(order+1)
            
            if np.sum(E_sph[i] < 0) > 1 and i != 0:
                t = gravitree.Tree(dx, eps[snap], mp, gravitree.G_COSMO)    
                E[i], d = gravitree.unbind(t, dv, return_diagnostics=True)
                dt_tot = sum([x[2] for x in d])
                if i < 10:
                    print("        sub %2d: np = %8d n_iter = %2d dt = %8.2f s"
                          % (i, len(E[i]), len(d), dt_tot))
            else:
                E[i] = np.zeros(len(p[i]))
        
        t3 = time.time()
                
        symlib.write_include_variable(sim_dir, snap, "E", E)
        symlib.write_include_variable(sim_dir, snap, "E_sph", E_sph)
        symlib.write_include_variable(sim_dir, snap, "m_enc", m_enc)

        t1 = time.time()
        
        print("    snap %4d time: %8.2f s (I/O: %8.2f s)" %
              (snap, t1 - t0, t2 - t0 + t1 - t3))
            
def main():
    base_dir, suite, idx_str = sys.argv[1],  sys.argv[2], sys.argv[3]
    target_idx = int(idx_str)

    for host_i in range(symlib.n_hosts(suite)):
        if host_i != target_idx and target_idx != -1: continue
        sim_dir = symlib.get_host_directory(base_dir, suite, host_i)
        write_includes(sim_dir)
        
if __name__ == "__main__": main()
