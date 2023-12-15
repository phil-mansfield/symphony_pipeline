import symlib
import numpy as np
import time

base_dir = "/sdf/home/p/phil1/ZoomIns"
suite = "SymphonyMilkyWayHR"
i_host = 0

sim_dir = symlib.get_host_directory(base_dir, suite, i_host)

snap =  200
part = symlib.Particles(sim_dir)
ps = part.read(snap, mode="smooth")
pa = part.read(snap, mode="all")
pc = part.read(snap, mode="current")

ids = [p["id"] for p in ps]
ida = [p["id"] for p in pa]
idc = [p["id"] for p in pc]
ok = [p["ok"] for p in pa]

_ids = symlib.transform_smooth_particles(ids, part.part_info, "smooth")
t0 = time.time()
_ida = symlib.transform_smooth_particles(ids, part.part_info, "all")
t1 = time.time()
_idc = symlib.transform_smooth_particles(ids, part.part_info, "current", ok)

print(t1 - t0)

#for i in range(len(ids)):
#    print(np.all(ids[i] == _ids[i]),
#          np.all(ida[i] == _ida[i]),
#          np.all(idc[i] == _idc[i]))
