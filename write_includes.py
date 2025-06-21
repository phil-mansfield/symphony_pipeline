import numpy as np
import symlib
import find_infall_cores
import sys
import gravitree
import time

def write_includes(sim_dir):
    sf, hist = symlib.read_symfind(sim_dir)
    scale = symlib.scale_factors(sim_dir)
    n_snap = len(scale)
    param = symlib.simulation_parameters(sim_dir)
    mp = param["mp"]/param["h100"]
    eps = param["eps"]/param["h100"]*scale
    
    symlib.write_include_header(sim_dir)

    part = symlib.Particles(sim_dir)
    
    for snap in range(n_snap):
        if snap % 10 == 0: print("    %d" % snap)
        p = part.read(snap, mode="current")
        
        E = [None]*len(p)
        E_sph = [None]*len(p)
        m_enc = [None]*len(p)

        r_host_0 = None
        
        for i in range(len(sf)):
            if len(p[i]) == 0 or not sf["ok"][i,snap]:
                E[i] = np.zeros(len(p[i]))
                E_sph[i] = np.zeros(len(p[i]))
                m_enc[i] = np.zeros(len(p[i]))
                continue
            
            dx = p[i]["x"] - sf["x"][i,snap]
            dv = p[i]["v"] - sf["v"][i,snap]

            # This can be sped up with binsort if needed
            rmax, vmax, pe_sph_vmax2, order = symlib.profile_info(param, dx)

            if i == 0:
                r_host_0 = np.sqrt(np.sum(dx*dx, axis=1))
                r_host_0 = r_host_0[order]
                
            pe_sph = pe_sph_vmax2*vmax*vmax
            ke = np.sum(dv*dv, axis=1)/2
            E_sph[i] = pe_sph + ke
            m_enc[i] = mp*(order+1)
            
            if i != 0 and np.sum(E_sph[i] < 0) > 1:
                t = gravitree.Tree(dx, eps[snap], mp, gravitree.G_COSMO)
                E[i] = gravitree.unbind(t, dv) + ke
            else:
                E[i] = np.zeros(len(p[i]))

        symlib.write_include_variable(sim_dir, snap, "E", E)
        symlib.write_include_variable(sim_dir, snap, "E_sph", E_sph)
        symlib.write_include_variable(sim_dir, snap, "m_enc", m_enc)
            
def main():
    base_dir, suite, idx_str = sys.argv[1],  sys.argv[2], sys.argv[3]
    target_idx = int(idx_str)

    for host_i in range(symlib.n_hosts(suite)):
        if host_i != target_idx and target_idx != -1: continue
        sim_dir = symlib.get_host_directory(base_dir, suite, host_i)
        write_includes(sim_dir)
        
if __name__ == "__main__": main()
