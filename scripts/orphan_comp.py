import numpy as np
import palette
from palette import pc
import symlib
import matplotlib.pyplot as plt
import matplotlib as mpl

palette.configure(False)

def main():
    base_dir = "/sdf/home/p/phil1/ZoomIns"
    suite = "SymphonyGroup"

    info = [
        (181, "Halo015", 70817991, 2),
        (155, "Halo444", 56285063, 13),
        (129, "Halo962", 21551142, 17),
        (142, "Halo399", 21519041, 42),
        (175, "Halo581", 62832734, 19)
    ]

    colors = [pc("r"), pc("o"), pc("g"), pc("b"), pc("p")]

    fig1, ax1 = plt.subplots(figsize=(7,7))
    fig2, ax2 = plt.subplots(figsize=(7,7))
    
    for i in range(5):
        snap, name, _, i_sub = info[i]
        sim_dir = symlib.get_host_directory(base_dir, suite, name)
        rs, hist = symlib.read_rockstar(sim_dir)
        um = symlib.read_um(sim_dir)
        sf, _ = symlib.read_symfind(sim_dir)
        a = symlib.scale_factors(sim_dir)

        rs_ok = rs["ok"][i_sub]
        um_ok = um["ok"][i_sub]
        sf_ok = sf["ok"][i_sub]
        
        ax1.plot(a[rs_ok], rs["m"][i_sub,rs_ok], c=colors[i])
        ax1.plot(a[um_ok], um["mvir"][i_sub,um_ok], lw=1.5, c=colors[i])
        ax1.plot(a[sf_ok], sf["m"][i_sub,sf_ok], "--", c=colors[i])

        ax2.plot(a[um_ok], um["m_star"][i_sub,um_ok], c=colors[i])

        fig0, ax0 = plt.subplots()
        snap0 = hist["first_infall_snap"][i_sub]+1
        host, sub, sub_sf = rs[0,snap0], rs[i_sub,snap0], sf[i_sub,snap0]
        lim = 1.5*host["rvir"]
        ax0.set_xlim(-lim, lim)
        ax0.set_ylim(-lim, lim)

        print(name, "i_sub =", i_sub, "snap =", snap)
        print("x =", sub["x"])
        print("v =", sub["v"])
        print("m = %.4g" % sub["m"])
        print()

        ix, iy = 0, 1
        c = pc("r")
        symlib.plot_circle(ax0, host["x"][ix], host["x"][iy], host["rvir"],
                           ls=":", c=c, lw=3)
        if sub["ok"]:
            symlib.plot_circle(ax0, sub["x"][ix], sub["x"][iy], sub["rvir"],
                               ls="-", c=c, lw=3)
        if sub_sf["ok"]:
            symlib.plot_circle(ax0, sub_sf["x"][ix], sub_sf["x"][iy],
                               sub_sf["r_half"], ls="--", c=c, lw=3)

        part = symlib.Particles(sim_dir)
        p = part.read(snap0, mode="smooth", halo=i_sub)
        ci = part.core_indices(mode="smooth", halo=i_sub)
        cores = p[ci]

        cr = np.sqrt(np.sum((cores["x"]-sub["x"])**2, axis=1))
        #print(cr)
        if np.sum(cr > 2) > 0:
            print(ci)
            print("problem children")
            print(ci[cr > 2])
            print(cores["x"][cr > 2])
            print(cr[cr > 2])
            print()
        
        norm = mpl.colors.LogNorm(vmin=1, vmax=10000)
        kwargs = {"extent": [-lim, lim, -lim, lim],
                  "norm": norm, "cmap": "inferno", "gridsize": 200}
        ax0.hexbin(p["x"][:,ix], p["x"][:,iy], **kwargs)
            
        ax0.set_xlabel(r"$%s\ ({\rm kpc})$" % ["x", "y", "z"][ix])
        ax0.set_ylabel(r"$%s\ ({\rm kpc})$" % ["x", "y", "z"][ix])
        ax0.plot(cores["x"][:,ix], cores["x"][:,iy], ".", c="white")
        fig0.savefig("../plots/um_comp/part_%s_%d_%d.png" % (name, i_sub, snap))
            
    for ax in [ax1, ax2]:
        ax.set_yscale("log")
        ax.set_xlabel(r"$a(r)$")
        ax.set_ylabel(r"$m\ (M_\odot)$")

        ylo, yhi = ax.get_ylim()
        ax.set_ylim(ylo, yhi)
        for i in range(5):
            snap = info[i][0]
            ax.plot([a[snap], a[snap]], [ylo, yhi],
                    "--", c=colors[i], lw=1)

    ax2.set_ylabel(r"$m_\star\ (M_\odot)$")
        
    fig1.savefig("../plots/um_comp/um_comp1_new.png")
    fig2.savefig("../plots/um_comp/um_comp2_new.png")
    
    
if __name__ == "__main__": main()
