import numpy as np
import symlib

base_dir = "/sdf/group/kipac/g/cosmo/ki21/phil1/simulations/ZoomIns/"
suite = "SymphonyMilkyWay"

def main():
    n_hosts = symlib.n_hosts(suite)
    
    n_err, n_tot, n_last = 0, 0, 0

    for i_host in range(n_hosts):
        sim_dir = symlib.get_host_directory(base_dir, suite, i_host)
        param = symlib.simulation_parameters(sim_dir)
        mp = param["mp"]/param["h100"]

        h, hist = symlib.read_subhalos(sim_dir, include_false_selections=True)
        c1 = symlib.read_cores(sim_dir, include_false_selections=True,
                               suffix="fid")
        c2 = symlib.read_cores(sim_dir, include_false_selections=True,
                               suffix="fid2")

        idx = np.arange(len(h), dtype=int)
        ok = hist["mpeak_pre"]/mp > 300
        ok[0] = False
        
        err1 = (~c1["ok"][idx, hist["first_infall_snap"]]) & ok
        err2 = (~c2["ok"][idx, hist["first_infall_snap"]]) & ok

        err_last = (~c2["ok"][:,-1]) & (h["ok"][:,-1] & c2["ok_rs"][:,-1]) & ok & err2

        n_err, n_tot = n_err+np.sum(err2), np.sum(ok)+n_tot
        n_last += np.sum(err_last)

        print(i_host, np.sum(err2), np.sum(err_last), np.sum(ok),
              "%.4f" % (n_err/n_tot), "%.4f" % (n_last/n_tot))


if __name__ == "__main__": main()
