import numpy as np
import symlib
import palette
from palette import pc
import matplotlib.pyplot as plt

suite = "SymphonyCluster"

def get_particles():
    base_dir = "/sdf/home/p/phil1/ZoomIns"
    sim_dir = symlib.get_host_directory(base_dir, suite, 0)
    
    a = symlib.scale_factors(sim_dir)
    n_snap = len(a)

    idxs = np.array([
        96929, 1317931, 266911, 201679, 189202,
        492554, 417363, 483668, 467444, 371843,
    ])
    #snap = [141, 130, 146, 144, 144,
    #        153, 151, 153, 152, 149]

    info = symlib.ParticleInfo(sim_dir)

    x, v = np.zeros((len(idxs), n_snap, 3)), np.zeros((len(idxs), n_snap, 3))
    a = symlib.scale_factors(sim_dir)

    id = symlib.read_particles(info, sim_dir, 0, "id", owner=0)
    print(id[idxs])

    h, hist = symlib.read_subhalos(sim_dir, comoving=True)
    for snap in [140, 141, 148, 149]:
        print(h["id"][0,snap], h["x"][0,snap], h["rvir"][0,snap])

    for snap in range(n_snap):
        if snap % 20 == 0: print(snap)
        if snap == 0:
            s = symlib.read_particles(info, sim_dir, snap, "snap", owner=0)
            print(s[idxs])

        ok  = symlib.read_particles(info, sim_dir, snap, "valid", owner=0)
        xx = symlib.read_particles(info, sim_dir, snap, "x", owner=0)
        vv = symlib.read_particles(info, sim_dir, snap, "v", owner=0)

        for j, i in enumerate(idxs):
            x[j,snap], v[j,snap] = xx[i]*1e3, vv[i]

        if snap in [140, 141, 148, 149]:
            print(snap)
            print(h["x"][0,snap])
            print(1, xx[96929])
            print(2, xx[371843])
        
    return x, v

def main():
    palette.configure(False)

    x, v = get_particles()

    base_dir = "/sdf/home/p/phil1/ZoomIns"
    sim_dir = symlib.get_host_directory(base_dir, "SymphonyCluster", 0)

    a = symlib.scale_factors(sim_dir)
    h, hist = symlib.read_subhalos(sim_dir, comoving=True)

    colors = [pc("r"), pc("o"), pc("g"), pc("b"), pc("p")]

    for i in range(len(x)):
        if i % 5 == 0:
            plt.figure(i//5)
            plt.xlabel(r"$a$")
            plt.ylabel(r"$r\ (ckpc)$")

            rvir = h["rvir"][0][h["ok"][0]]*1e3
            plt.plot(a, rvir, "--", c="k")

        ok = ~np.isnan(x[i,:,0])
        dx = x[i,ok] - h[0,ok]["x"]*1e3
        r = np.sqrt(np.sum(dx**2, axis=1))

        plt.plot(a[ok], r, colors[i % 5])

    plt.figure(0)
    plt.savefig("../plots/core_tracking/p_tracks_1.png")
    plt.figure(1)
    plt.savefig("../plots/core_tracking/p_tracks_2.png")


if __name__ == "__main__": main()
