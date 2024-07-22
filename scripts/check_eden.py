import symlib

base_dir = "/sdf/home/p/phil1/ZoomIns"
suite = "EDEN_MilkyWay_8K"

for i in range(symlib.n_hosts(suite)):
    sim_dir = symlib.get_host_directory(base_dir, suite, i)
    try:
        symlib.read_symfind(sim_dir)
    except:
        print("Failure with %d" % i)

print("All halos checked")
