####################
# file_nonsense.py #
####################

""" This file handles various tasks related to setting files for data transfer.

Run as python3 file_nonsense.py <mode> <suite> <halo_num>

It can be run in several modes:
move_trees - move tree files from their old location HaloXYZ/halos to
             HaloXYZ/trees
tar_halos - tar the relevant portions of the halos subdirectory
tar_trees - tar the trees subdirectory
tar_particles - tar the particle subdirectory
tar_full_snapshots - tar the full_snapshot subdirectory
tar - Run all four of the tar modes

<suite> can be any of the suites in ZoomIns, or it can be all to iterate over 
all suites.  halo_num can be -1 to loop thorugh all halos or to a specific
number to perform a task on just one
"""

import sys
import os.path as path
import os
import symlib
import numpy as np
import glob
import tarfile

cd_dir = "/fs/ddn/sdf/group/kipac/g/cosmo/ki21/phil1/simulations/ZoomIns"
base_dir = "."
tar_dir = "/fs/ddn/sdf/group/kipac/g/cosmo/ki21/phil1/simulations/ZoomIns/tar_files"

if cd_dir != os.getcwd():
    print("Need to run file_nonsense.py from %s" % cd_dir)
    exit(1)

def main():
    try:
        mode, suite, halo = sys.argv[1:]
        halo = int(halo)
    except:
        print("""Usage:
python3 file_nonsense.py <mode> <suite> <halo>

See file docstring for additional details.""")
        exit(1)

    if mode == "move_trees":
        move_trees(suite, halo)
    if mode in ["tar_halos", "tar"]:
        tar_halos(suite, halo)
    if mode in ["tar_trees", "tar"]:
        tar_trees(suite, halo)
    if mode in ["tar_particles", "tar"]:
        tar_particles(suite, halo)
    if mode in ["tar_full_snapshots", "tar"]:
        tar_full_snapshots(suite, halo)

def suites_and_halos(suite, halo):
    if suite == "all":
        suites = ["SymphonyLMC", "SymphonyMilkyWay", "SymphonyGroup",
                  "SymphonyLCluster", "SymphonyCluster", "MWest"]
    else:
        suites = [suite]

    halos = [None]*len(suites)
    for i in range(len(suites)):
        if halo >= 0:
            halos[i] = np.array([halo], dtype=int)
        elif halo == -1:
            halos[i] = np.arange(symlib.n_hosts(suites[i]), dtype=int)

    return suites, halos
        
def move_trees(suite, halo):
    suites, halos = suites_and_halos(suite, halo)
    for j, suite in enumerate(suites):
        print(suite)
        for i in halos[j]:
            print("   ", i)
            sim_dir = symlib.get_host_directory(base_dir, suite, i)
            tree_dir = path.join(sim_dir, "trees")
            if not path.exists(tree_dir): os.mkdir(tree_dir)
            tree_files = glob.glob(path.join(sim_dir, "halos", "tree_*"))
            for tree_file in tree_files:
                base = path.basename(tree_file)
                target = path.join(tree_dir, base)
                os.rename(tree_file, target)
            os.rename(
                path.join(tree_dir, "tree_header.dat"),
                path.join(sim_dir, "halos", "tree_header.dat")
            )
                

def tar_directory(suite, halo, dir_name):
    suites, halos = suites_and_halos(suite, halo)
    for j, suite in enumerate(suites):
        print(suite)
        for i in halos[j]:
            print("   ", i)
            sim_dir = symlib.get_host_directory(base_dir, suite, i)
            halo_name = path.basename(sim_dir)
            
            data_dir = path.join(sim_dir, dir_name)
            tar_target = path.join(tar_dir, "%s_%s_%s.tar" % (suite, halo_name, dir_name))
            
            if path.exists(tar_target):
                os.remove(tar_target)
                
            f = tarfile.open(tar_target, "w")
            if dir_name == "halos":
                if suite == "SymphonyCluster":
                    file_names = ["subhalos.dat",
                                  "snap_scale.dat",
                                  "tree_header.dat"]
                else:
                    if path.exists(path.join(data_dir, "cores_fid4.dat")):
                        file_names = ["infall_cores.dat", "subhalos.dat",
                                      "cores_fid4.dat", "snap_scale.dat",
                                      "tree_header.dat"]
                    else:
                        file_names = ["infall_cores.dat", "subhalos.dat",
                                      "cores_fid3.dat", "snap_scale.dat",
                                      "tree_header.dat"]
                        
                files = [path.join(data_dir, file_name) for 
                         file_name in file_names]
                for file_name in files: f.add(file_name)
            else:
                f.add(data_dir)

            f.close()

def tar_halos(suite, halo):
    tar_directory(suite, halo, "halos")                

def tar_trees(suite, halo):
    tar_directory(suite, halo, "trees")                

def tar_particles(suite, halo):
    tar_directory(suite, halo, "particles")

def tar_full_snapshots(suite, halo):
    tar_directory(suite, halo, "full_snapshots")

if __name__ == "__main__": main()
