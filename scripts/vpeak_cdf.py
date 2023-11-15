import numpy as np
import matplotlib.pyplot as plt
import palette
import symlib
from palette import pc

base_dir = "/sdf/home/p/phil1/ZoomIns"

def main():
    suites = [
        #"SymphonyLMC", "SymphonyMilkyWay", "SymphonyGroup",
        #"SymphonyLCluster", 
        "SymphonyCluster"
    ]

    for suite in suites:
        print(suite)
        n_host = symlib.n_hosts(suite)

        vmax_all, vpeak_all = [], []

        for i_host in range(n_host):
            if i_host != 13: continue
            print(i_host)
            sim_dir = symlib.get_host_directory(base_dir, suite, i_host)
            h, _ = symlib.read_subhalos(sim_dir)
            I = np.arange(len(h), dtype=int)
            i_peak = np.argmax(h["mvir"], axis=1)
            vmax_peak = h["vmax"][I,i_peak]
            vmax_all.append(vmax_peak)
            
        vmax = np.hstack(vmax_all)
        for p in [0.5, 0.1, 0.05, 0.01]:
            print("    %.2f: %6.2f" % (p, np.quantile(vmax, p)))
            

if __name__ == "__main__": main()
