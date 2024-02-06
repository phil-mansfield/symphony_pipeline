import numpy as np
import symlib

base_dir = "/sdf/home/p/phil1/ZoomIns"
suite = "SymphonyMilkyWay"

vmax_vpeak = []
vmax_vmax = []
mvir_vpeak = []
mvir_vmax = []
a_vpeak = []
a_vmax = []

for i_host in range(symlib.n_hosts(suite)):
    sim_dir = symlib.get_host_directory(base_dir, suite, i_host)
    scale = symlib.scale_factors(sim_dir)    

    rs, hist = symlib.read_rockstar(sim_dir)

    ok = rs["ok"][:,-1]
    rs, hist = rs[ok][1:], hist[ok][1:]
    i1, i2 = np.argmax(hist["vpeak"]), np.argmax(rs["vmax"][:,-1])
    
    vmax_vpeak.append(rs["vmax"][i1,-1])
    vmax_vmax.append(rs["vmax"][i2,-1])
    mvir_vpeak.append(rs["m"][i1,-1])
    mvir_vmax.append(rs["m"][i2,-1])
    a_vpeak.append(scale[hist["merger_snap"][i1]])
    a_vmax.append(scale[hist["merger_snap"][i2]])

print(np.median(vmax_vpeak))
print(np.median(vmax_vmax))
print("%g" % np.median(mvir_vpeak))
print("%g" % np.median(mvir_vmax))
print(np.median(a_vpeak))
print(np.median(a_vmax))

print()

print(np.mean(vmax_vpeak))
print(np.mean(vmax_vmax))
print("%g" % np.mean(mvir_vpeak))
print("%g" % np.mean(mvir_vmax))
print(np.mean(a_vpeak))
print(np.mean(a_vmax))
