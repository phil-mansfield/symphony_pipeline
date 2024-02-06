import symlib

# Specify the halo you want to read. This bit would probably be done inside a loop
base_dir = "/sdf/group/kipac/g/cosmo/ki21/phil1/simulations/ZoomIns"
suite = "SymphonyGroup"
# you can choose the halo name either through its index within the halo list
# or its full name, e.g., "Halo015". Whichever is easiest for you.
halo = 0 

sim_dir = symlib.get_host_directory(base_dir, suite, halo)
scale = symlib.scale_factors(sim_dir)
um = symlib.read_um(sim_dir)

print(scale)
print(um["m_star"][0])
