import numpy as np
import symlib
import sys
import os.path

def main():
    config_name, idx_str = sys.argv[1], sys.argv[2]
    target_idx = int(idx_str)

    n_halo = len(np.loadtxt(config_name, dtype=str))

    for i in range(n_halo):
        if target_idx == -1 or i == target_idx:
            write_um_file(config_name, i)
            

def write_um_file(config_file, target_idx):
    config_cols = np.loadtxt(config_file, dtype=str).T
    sim_dir, um_fmt = config_cols[7][target_idx], config_cols[9][target_idx]

    scale = symlib.scale_factors(sim_dir)
    n_snap = len(scale)
    um_files = [None]*n_snap
    n = 0
    for i in range(n_snap):
        um_file = um_fmt % scale[i]
        if os.path.exists(um_file):
            um_files[i] = um_file
            n += 1
    print("%2d" % target_idx, n)

    
#def get_sim_dirs(config_name):
#    with open(config_name, "r") as fp: text = fp.read()
#    lines = [line for line in text.split("\n") if len(line) > 0]
#    for i in range(len(lines)):
#        lines[i] = [tok for tok in lines[i].split(" ") if len(tok) > 0]
#    return [line[7] for line in lines]

if __name__ == "__main__": main()
