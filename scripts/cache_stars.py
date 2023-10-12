import numpy as np
import symlib
import pickle

cache_dir = "/sdf/home/p/phil1/data/data/star_cache"
base_dir = "/sdf/home/p/phil1/ZoomIns"

def read_stars(suite, i_host):
    with open("%s/%s_%s.dat" % (cache_dir, suite, i_host), "rb") as fp:
        return pickle.load(fp)

def main():
    #suites = ["SymphonyLMC", "SymphonyMilkyWay", "MWest",
    #          "SymphonyGroup", "SymphonyCluster"]
    suites = ["SymphonyMilkyWayHR"]

    gal_halo = symlib.GalaxyHaloModel(
        symlib.UniverseMachineMStarFit(),
        symlib.Jiang2019RHalf(),
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

            with open("%s/%s_%s.dat" % (cache_dir, suite, i_host), "wb+") as fp:
                pickle.dump(stars, fp)

if __name__ == "__main__": main()
