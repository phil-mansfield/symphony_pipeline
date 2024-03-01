import numpy as np
import symlib

def main():
    suites = ["SymphonyCluster"]

    base_dir = "/oak/stanford/orgs/kipac/users/phil1/simulations/ZoomIns/"

    for i_suite in range(len(suites)):
        suite = suites[i_suite]
        print(suite)

        mc = 0
        n_host = symlib.n_hosts(suite)
        for i_host in range(n_host):
            sim_dir = symlib.get_host_directory(base_dir, suite, i_host)
            h, hist = symlib.read_subhalos(sim_dir)
            ratio = hist["mpeak"][1:]/hist["mpeak"][0]
            mc += np.sum(ratio > 0.01)
        print(mc/n_host)


if __name__ == "__main__": main()
