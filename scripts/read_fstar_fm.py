import numpy as np
import matplotlib.pyplot as plt
import symlib
import palette
import gravitree
import time
from palette import pc
import struct
import numpy.random as random
            
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

def m_lim_eps(param, scale, h, c, ok, infall_snap):
    # returns a mass ratio, not a mass
    # Factor of 1.284 comes from Mansfield & Avestruz (2021) 
    eps = scale*param["eps"]/param["h100"]
    c_infall = h["cvir"][infall_snap]
    rs_infall = h["rvir"][infall_snap]/h["cvir"][infall_snap]
    def f(x):
        return np.log(1+x) - x/(1+x)
    return 1.79/1.284 * ((eps[ok]*c["r50_bound"][ok])/
                         (f(c_infall)*rs_infall**2))

def m_lim_n(param, scale, h, c, ok, infall_snap):
    # Returns a mass ratio, not a mass
    mp = param["mp"]/param["h100"]
    n_infall = h["mvir"][infall_snap]/mp
    return 0.32*(n_infall/1e3)**-0.8 * np.ones(np.sum(ok))

def is_converged_vdb18(sim_dir):
    param = symlib.simulation_parameters(sim_dir)
    scale = symlib.scale_factors(sim_dir)
    h, hist = symlib.read_subhalos(sim_dir)
    c = symlib.read_cores(sim_dir)

    is_conv = np.zeros(h.shape, dtype=bool)
    snap = np.arange(h.shape[1], dtype=int)

    for i in range(1, len(hist)):
        ok = c["ok"][i]
        m_lim_1 = m_lim_eps(
            param, scale, h[i], c[i], ok, hist["first_infall_snap"][i])
        m_lim_2 = m_lim_n(
            param, scale, h[i], c[i], ok, hist["first_infall_snap"][i])
        m_lim = np.maximum(m_lim_1, m_lim_2)
        mass_lim = m_lim * h["mvir"][i][hist["first_infall_snap"][i]]
        
        conv_snaps = snap[ok][c["m_bound"][i,ok] > mass_lim]
        if len(conv_snaps) == 0: continue

        is_conv[i] = ((snap >= hist["first_infall_snap"][i]) & 
                      (snap <= np.max(conv_snaps)))
                      
    return is_conv & c["ok"]

def main():
    base_dir = "/sdf/home/p/phil1/ZoomIns"
    suite = "SymphonyMilkyWay"
    sim_dir = symlib.get_host_directory(base_dir, suite, 0)

    file_name = "tables/example_mass_loss_0.dat"
    ok, m_star, rhalf_star, m_dm, rhalf_dm = read_mass_loss_table(file_name)
    print("Properties for subhalo 14:")
    print("m_star\n", m_star[14][ok[14]])
    print("m_dm\n", m_dm[14][ok[14]])
    print("r_half,star\n", rhalf_star[14][ok[14]])
    print("r_half,dm\n", rhalf_dm[14][ok[14]])
    
    is_conv = is_converged_vdb18(sim_dir)
    print("Within convergence limits:")
    print("m_star\n", m_star[14][is_conv[14]])
    print("m_dm\n", m_dm[14][is_conv[14]])
    print("r_half,star\n", rhalf_star[14][is_conv[14]])
    print("r_half,dm\n", rhalf_dm[14][is_conv[14]])

if __name__ == "__main__": main()
