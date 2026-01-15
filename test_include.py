import symlib

print("Testing include=... readers")

base_dir = "/sdf/home/p/phil1/ZoomIns/"
suite = "SymphonyMilkyWay"
snap = 235
i_sub = 3

print("Reading headers")
sim_dir = symlib.get_host_directory(base_dir, suite, 0)
part = symlib.Particles(sim_dir, include=["E_sph", "E"])

print("Reading particles in 'current' mode")
p = part.read(235, mode="current")
print("E")
print(p[i_sub]["E"])
print("E_sph")
print(p[i_sub]["E_sph"])
print("Reading particles in 'stars' mode")
p = part.read(235, mode="stars")
print("E")
print(p[i_sub]["E"])
print("E_sph")
print(p[i_sub]["E_sph"])
print("Completed reading without errors")
