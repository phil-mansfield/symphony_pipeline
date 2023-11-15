import numpy as np
import symlib

base_dir = "/sdf/home/p/phil1/ZoomIns"
suite = "SymphonyCluster"
i_host = 0

def main():
    sim_dir = symlib.get_host_directory(base_dir, suite, i_host)

    scale = symlib.scale_factors(sim_dir)
    h, hist = symlib.read_subhalos(sim_dir)
    h_cmov, hist_cmov = symlib.read_subhalos(sim_dir, comoving=True)
    part = symlib.Particles(sim_dir)

    for snap in range(199, 200):
        rvir = h[0,snap]["rvir"]
        p = part.read(snap, mode="smooth")[0]
        r = np.sqrt(np.sum(p["x"]**2, axis=1))
        vr = np.sum(p["x"]*p["v"], axis=1)/r

        ok = (p["snap"] == snap) & p["smooth"]
        too_soon = (r > rvir) & ok
        too_late = (vr > 0) & ok
        h100 = symlib.simulation_parameters(sim_dir)["h100"]
        print("%3d %5d %5d %5d %5d" % (snap, np.sum(ok), np.sum(too_soon), np.sum(too_late), np.sum(too_soon & too_late)))
        print(np.where(too_soon))
        print("id", p["id"][4:6])
        print("rvir", rvir/scale[snap]/1e3*h100, h[0,snap]["id"], scale[snap])
        print("dr", r[4:6]/scale[snap]/1e3*h100)
        print(r[too_soon]/scale[snap]/1e3*h100)
        print("dx", p["x"][4:6]/scale[snap]/1e3*h100)
        print("hx", h_cmov["x"][snap,snap])
        print("x", (p["x"]/scale[snap]/1e3*h100 + h_cmov["x"][0,snap])[4:6])
        break

if __name__ == "__main__": main()
