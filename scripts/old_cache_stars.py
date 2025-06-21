import numpy as np
import symlib
import pickle
import sys

cache_dir = "/sdf/home/p/phil1/data/data/star_cache"
base_dir = "/sdf/home/p/phil1/ZoomIns"

def read_stars(suite, i_host, method_name=""):
    if method_name == "":
        fname = "%s/%s_%s.dat" % (cache_dir, suite, i_host)
    else:
        fname = "%s/%s_%s_%s.dat" % (cache_dir, suite, method_name, i_host)

    with open(fname, "rb") as fp:
        return pickle.load(fp)

def main():
    #suites = ["SymphonyLMC", "SymphonyMilkyWay", "MWest",
    #          "SymphonyGroup", "SymphonyCluster"]
    suites = ["SymphonyMilkyWayHR"]

    if len(sys.argv) == 1:
        method_name = ""
    else:
        method_name = sys.argv[1]

    if method_name == "":
        gal_halo = symlib.GalaxyHaloModel(
            symlib.UniverseMachineMStarFit(),
            symlib.Jiang2019RHalf(),
            symlib.PlummerProfile(),
            symlib.Kirby2013Metallicity(),
            no_scatter=False
        )
    elif len(method_name) > 2 and method_name[:2] == "r=":
        size = {
            "r=0.005": symlib.FixedRHalf(0.005),
            "r=0.008": symlib.FixedRHalf(0.008),
            "r=0.015": symlib.FixedRHalf(0.015),
            "r=0.025": symlib.FixedRHalf(0.025),
            "r=0.05": symlib.FixedRHalf(0.05),
        }[method_name]
        gal_halo = symlib.GalaxyHaloModel(
            symlib.UniverseMachineMStarFit(),
            size,
            symlib.PlummerProfile(),
            symlib.Kirby2013Metallicity(),
            no_scatter=False
        )

    for suite in suites:
        n_host = symlib.n_hosts(suite)

        for i_host in range(n_host):
            print(suite, i_host+1, "/", n_host)
            sim_dir = symlib.get_host_directory(base_dir, suite, i_host)

            r, _ = symlib.read_rockstar(sim_dir)
            target_subs = np.arange(1, len(r), dtype=int)

            mp_star, _, m_star, r_half, Fe_H = symlib.tag_stars(
                sim_dir, gal_halo, target_subs=target_subs)
            
            stars = (mp_star, Fe_H, m_star, r_half)

            if method_name == "":
                fname = "%s/%s_%s.dat" % (cache_dir, suite, i_host)
            else:
                fname = "%s/%s_%s_%s.dat" % (cache_dir, suite, method_name, i_host)


            with open(fname, "wb+") as fp:
                pickle.dump(stars, fp)

if __name__ == "__main__": main()
