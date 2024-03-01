import cprofile

def main():
    palette.configure(False)

    base_dir = "/sdf/home/p/phil1/ZoomIns"
    suite = "SymphonyMilkyWay"
    i_host = 0

    sim_dir = symlib.get_host_directory(base_dir, suite, i_host)

    scale = symlib.scale_factors(sim_dir)

    gal_halo = symlib.DWARF_GALAXY_HALO_MODEL

    cProfile.run("symlib.tag_stars(sim_dir, gal_halo)")
    

if __name__ == "__main__": main()
