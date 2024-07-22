import symlib
import numpy as np
import os.path as path

base_dir = "/sdf/home/p/phil1/ZoomIns"
suite = "SymphonyMilkyWay"
out_dir = "../tmp_data"

def main():
    sim_dir = symlib.get_host_directory(base_dir, suite, 0)
    a = symlib.scale_factors(sim_dir)
    rs, hist = symlib.read_rockstar(sim_dir)
    
    with open(path.join(out_dir, "scale.txt"), "w+") as fp:
        for i in range(len(a)): print("%f" % a[i], file=fp)
    
    part = symlib.Particles(sim_dir)
    bins_per_decade = 64
    r_bins = 10**np.linspace(np.log10(0.02), np.log10(200), 4*bins_per_decade)
    n_bins = len(r_bins)
    for snap in range(0, len(a), 10):
        with open(path.join(out_dir, "prof_%03d.txt" % snap), "w+") as fp:
            for i in range(n_bins):
                print("%f" % r_bins[i], end=" ", file=fp)
            print("", file=fp)

            p = part.read(snap)

            for i_sub in range(len(rs)):
                if not rs["ok"][i_sub,snap]:
                    n = np.ones(n_bins, dtype=int)*-1
                else:
                    bins = np.hstack([[0], r_bins])
                    dx = p[i_sub]["x"] - rs["x"][i_sub,snap]
                    r = np.sqrt(np.sum(dx**2, axis=1))
                    n, _ = np.histogram(r, bins)

                for i in range(n_bins):
                    print("%d" % n[i], end=" ", file=fp)
                print("", file=fp)
                    
if __name__ == "__main__": main()
