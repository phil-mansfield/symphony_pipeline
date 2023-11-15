import numpy as np
import symlib

SUITE = "SymphonyMilkyWay"
OUT_DIR = "../plots/core_tracking"
BASE_DIR = "/sdf/home/p/phil1/ZoomIns"

suffixes = [
    "commit_2", "commit_3", "k16_n32"
]

ref = "fid"

def main():
    sim_dir = symlib.get_host_directory(BASE_DIR, SUITE, 0)
    param = symlib.simulation_parameters(sim_dir)

    h, hist = symlib.read_subhalos(sim_dir)
    c_ref = symlib.read_cores(sim_dir, suffix=ref)

    n_ok_rs, n_ok_pt = np.sum(h["ok"][:,-1]), np.sum(c_ref["ok"][:,-1])
    
    print("n_ok_rs = %d" % n_ok_rs)
    print("n_ok_pt = %d" % n_ok_pt)
    print()

    for suffix in suffixes:
        c = symlib.read_cores(sim_dir, suffix=suffix)
        n_diff = np.sum(c["ok"][:,-1]) - n_ok_pt

        both_ok = c["ok"][:,-1] & c_ref["ok"][:,-1]
        c_ok, c_ref_ok = c[both_ok,-1], c_ref[both_ok,-1]

        m_frac_diff = c_ok["m_bound"]/c_ref_ok["m_bound"] - 1
        mt_frac_diff = c_ok["m_tidal"]/c_ref_ok["m_tidal"] - 1
        rt_frac_diff = c_ok["r_tidal"]/c_ref_ok["r_tidal"] - 1

        print(suffix)
        print("f_diff_ok = %.3f" % (n_diff / n_ok_pt))
        print("m = %.3f +/- %.3f" % (np.mean(m_frac_diff),
                                      np.std(m_frac_diff)))
        print("mt = %.3f +/- %.3f" % (np.mean(mt_frac_diff),
                                      np.std(mt_frac_diff)))
        print("rt = %.3f +/- %.3f" % (np.mean(rt_frac_diff),
                                      np.std(rt_frac_diff)))
        
        print()
        

if __name__ == "__main__": main()
