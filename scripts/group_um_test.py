import symlib

base_dir = "/sdf/home/p/phil1/ZoomIns"

def main():
    suite = "SymphonyGroup"
    sim_dir = symlib.get_host_directory(base_dir, suite, 9)
    um = symlib.read_um(sim_dir)
    i = 235
    print(um["m_star"][0,-i:])
    print()
    print(um["x"][0,-i:])
    print()
    print(um["v"][0,-i:])
    print()
    print(um["is_orphan"][0,-i:])
    print()
    print("----")
    print()
    print(um["m_star"][1,-i:])
    print()
    print(um["x"][1,-i:])
    print()
    print(um["v"][1,-i:])
    print()
    print(um["is_orphan"][1,-i:])
    print()
    
if __name__ == "__main__": main()
