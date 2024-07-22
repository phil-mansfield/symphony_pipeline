import numpy as np
import symlib

missing_ids = [
    70818252,
    37842198,
    62207651
]

halo_names = [
    "Halo015",
    "Halo175",
    "Halo909"
]

def main():
    base_dir = "/sdf/home/p/phil1/ZoomIns"
    suite = "SymphonyGroup"
    for i_host in range(3):
        id0 = missing_ids[i_host]
        sim_dir = symlib.get_host_directory(base_dir, suite, halo_names[i_host])
        b = symlib.read_branches(sim_dir)
        upids, ids, snaps = symlib.read_tree(sim_dir, ["upid", "id", "snap"])
        h, hist = symlib.read_rockstar(sim_dir)
        
        print(id0 in ids)
        print(id0 in h["id"])
        for i_sub in range(h.shape[0]):
            if missing_ids[i_host] in h["id"][i_sub]:
                print(i_sub, np.where(h["id"][i_sub] == id0)[0])

        print(h["id"][3,181], h["id"][0,181], h["id"][17,181])

if __name__ == "__main__": main()
