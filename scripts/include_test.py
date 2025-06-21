import sys
sys.path.append("/sdf/home/p/phil1/code/src/github.com/phil-mansfield/nimbus_plots")

import symlib
import lib
import cache_stars
import numpy as np

i_host = 0
suite = "MWest"
snap = 235

# This should point to wherever your simulation directory is
sim_dir = symlib.get_host_directory("/sdf/home/p/phil1/ZoomIns", suite, i_host)

# Read in subhalos and particles. The include flag includes additional particle
# properties.
sf, hist = symlib.read_symfind(sim_dir)
part = symlib.Particles(sim_dir, include=["E"])
p = part.read(235, mode="stars")

print(p[1]["E"])
print(np.sum(p[1]["v"]**2, axis=1))

# I'm doing this just to make the test script fast, but in practice, you'd actually
# get your stars from symlib.tag_stars() or symlib.retag_stars().
stars, _ = cache_stars.read_stars("fid_dwarf", "SymphonyMilkyWay", i_host)

x, v, mp_star = [], [], []
for i in range(1, len(sf)):
    # For now there's a bug where I forgot to add in the kinetic energy to "E".
    # Will fix that soon.
    dv = p[i]["v"] - sf["v"][i,snap]
    is_bound = np.sum(dv*dv, axis=1)/2 + p[i]["E"] < 0
    
    x.append(p[i]["x"][~is_bound])
    v.append(p[i]["v"][~is_bound])
    mp_star.append(stars[i]["mp"][~is_bound])
    
    if i < 20:
        print("%3d: f_dm = %.3f, f_star = %.3f" %
              (i, sf["m"][i,snap]/hist["mpeak_pre"][i],
               np.sum(stars[i]["mp"][is_bound])/np.sum(stars[i]["mp"])))
        
x, v, mp_star = np.concatenate(x), np.concatenate(v), np.concatenate(mp_star)

print("Stellar halo particles:", len(x))
print("Satellite particles:", sum(map(len, p[1:])) - len(x))
print("Unbound stellar halo mass: %.3g" % (np.sum(mp_star),))
print("Bound satelite mass: %.3g" % (sum(map(lambda x: np.sum(x["mp"]), stars)) - np.sum(mp_star),))
