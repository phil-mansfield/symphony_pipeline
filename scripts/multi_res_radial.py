import numpy as np
import palette
from palette import pc
import symlib
import os.path as path

def main():
    base_dir = "/sdf/home/p/phil1/ZoomIns/"
    suites = ["SymphonyMilkyWayLR", "SymphonyMilkyWayHR"]

    mass_lim = 4e5 * 8 * 300
    print("mpeak cutoff: %g" % mass_lim)

    n_rs, n_rs2, n_sy = 0, 0, 0

    for suite in suites:
        for i_host in range(symlib.n_hosts(suite)):
            sim_dir = symlib.get_host_directory(base_dir, suite, i_host)
            
            rs, hist = symlib.read_rockstar(sim_dir)
            sy, _ = symlib.read_symfind(sim_dir)
            host = rs[0]

            ok = hist["mpeak_pre"] > mass_lim
            rs, hist, sy = rs[ok], hist[ok], sy[ok]

            r_rs = np.sqrt(np.sum(rs["x"][:,-1,:]**2, axis=1))
            r_sy = np.sqrt(np.sum(sy["x"][:,-1,:]**2, axis=1))

            ok_rs  = (r_rs < host["rvir"][-1]) & rs["ok"][:,-1]
            ok_rs2 = (r_rs < host["rvir"][-1]) & sy["ok_rs"][:,-1]
            ok_sy  = (r_sy < host["rvir"][-1]) & sy["ok"][:,-1]

            suite_name = suite
            if suite_name == "SymphonyMilkyWayLR":
                suite_name = "SymphonyMilkyWay"

            print(suite_name, path.basename(sim_dir),
                  "%2d %2d %2d %.2f %.2f %.2f" %(
                      np.sum(ok_rs),
                      np.sum(ok_rs2),
                      np.sum(ok_sy),
                      np.min(r_rs[ok_rs]),
                      np.min(r_rs[ok_rs2]),
                      np.min(r_sy[ok_sy])))


if __name__ == "__main__": main()
