import cProfile
import symlib
import numpy as np

base_dir = "/sdf/home/p/phil1/ZoomIns"
suite = "SymphonyMilkyWay"
i_host = 0

sim_dir = symlib.get_host_directory(base_dir, suite, i_host)

gal_halo = symlib.DWARF_GALAXY_HALO_MODEL

stars, gals, ranks = symlib.tag_stars(sim_dir, gal_halo)
for i in range(len(stars)):
    print(i, "%g" % np.sum(stars[i]["mp"]))

print("profiling")

def main():
    symlib.retag_stars(sim_dir, gal_halo, ranks)

if __name__ == "__main__": cProfile.run("main()", sort="cumulative")
