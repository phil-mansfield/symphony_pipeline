import numpy as np

ids = np.loadtxt("config_richie.txt", usecols=(0,), dtype=int)
n_snap = 235
eps = 0.00036
mp = 0.0002255849403297816 * 1e10
n_file = 4
particle_fmt = "/sdf/home/p/phil1/data/data/SymphonyMilkyWayDisk/%s/snapshot_%%03d.%%d"
tree_fmt = "/sdf/group/kipac/u/ycwang/MW4K_disk/%s/output_disk/rockstar_all/trees/"
out_fmt = "/sdf/home/p/phil1/ZoomIns/SymphonyMilkyWayDisk/%s/"
tree_type = "ct_rvmax"
um_fmt = "/sdf/group/kipac/u/ycwang/MWmass_4K/%s/output/rockstar/groupcat/sfr_catalog_%%s.txt"

raw_tree_dirs = np.loadtxt("config_richie.txt", usecols=(6,), dtype=str)

for i in range(len(ids)):
    halo = raw_tree_dirs[i].split("/")[7]
    particle = particle_fmt % halo
    tree = tree_fmt % halo
    out = out_fmt % halo
    um = um_fmt % halo

    print("%d %d %f %g %d %s %s %s %s %s" %
          (ids[i], n_snap, eps, mp, n_file, particle,
           tree, out, tree_type, um))
