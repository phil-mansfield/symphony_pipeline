import numpy as np
import symlib
import gravitree

def main():
    base_dir = "/oak/stanford/orgs/kipac/users/phil1/simulations/ZoomIns"
    suite = "SymphonyMilkyWay"
    i_halo = 0
    sim_dir = symlib.get_host_directory(base_dir, suite, i_halo)

    param = symlib.simulation_parameters(suite)
    h, hist = symlib.read_subhalos(sim_dir)

    info = symlib.ParticleInfo(sim_dir)

    x = symlib.read_particles(info, sim_dir, 235, "x")
    v = symlib.read_particles(info, sim_dir, 235, "v")
    ok = symlib.read_particles(info, sim_dir, 235, "valid")

    for i in range(10, 20):
        if not h[i,-1]["ok"]: continue
        print(i)

        x_i = symlib.set_units_x(x[i], h[i,-1] 1.0, param)[ok]
        v_i = symlib.set_units_x(x[i], h[i,-1] 1.0, param)[ok]

        ok_i = gravitree.binding_energy() < 0

        rmax_1, vmax_2, _, _ = symlib.profile_info(param, )
        

if __name__ == "__main__": main()
