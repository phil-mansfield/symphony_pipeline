import symlib
import numpy as np
import palette
from palette import pc
import time

def main():
    palette.configure(False)

    base_dir = "/home/users/phil1/oak/simulations/ZoomIns/"
    suite = "MWest"
    i_host = 4
    i_sub = 1

    sim_dir = symlib.get_host_directory(base_dir, suite, i_host)

    scale = symlib.scale_factors(sim_dir)

    gal_halo = symlib.DWARF_GALAXY_HALO_MODEL_NO_UM

    stars, gal_hists, ranks = symlib.tag_stars(sim_dir, gal_halo)

if __name__ == "__main__": main()
