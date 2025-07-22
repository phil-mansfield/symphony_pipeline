import symlib

print("Testing include=... readers")

print("Reading headers")
part = symlib.Particles("/oak/stanford/orgs/kipac/users/phil1/simulations/ZoomIns/MWest/Halo004/", include=["E_sph", "E"])
print("Reading particles in 'all' mode")
p = part.read(235, mode="current")
print("Reading particles in 'stars' mode")
p = part.read(235, mode="stars")
print("Completed reading without errors")
