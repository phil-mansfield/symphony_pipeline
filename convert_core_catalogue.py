import numpy as np
import struct
import symlib
import sys
import os
import os.path as path
import scipy.interpolate as interpolate
import scipy.optimize as optimize

def get_sim_dirs(config_name):
    with open(config_name, "r") as fp: text = fp.read()
    lines = [line for line in text.split("\n") if len(line) > 0]
    for i in range(len(lines)):
        lines[i] = [tok for tok in lines[i].split(" ") if len(tok) > 0]
    return [line[7] for line in lines]

def parse_sim_dir(sim_dir):
    suite_dir, halo = path.split(sim_dir)
    base_dir, suite = path.split(suite_dir)
    return base_dir, suite, halo

def get_matching_file_names(file_fmt):
    file_names = []
    while True:
        file_name = file_fmt % len(file_names)
        if not path.exists(file_name): return file_names
        file_names.append(file_name)

def parse_flags(flags):
    out = { }
    for flag in flags:
        if len(flag) < 3 or flag[:2] != "--": continue
        tok = flag[2:].split("=")
        if len(tok) == 1:
            out[tok[0]] = ""
        else:
            out[tok[0]] = tok[1]
    return out

def set_ok_flags(c, h, hist):
    # You get a gold star if you can follow all the cuts here. They're
    # described mroe cleanly in the paper.
    c["ok"] = c["m_bound"] > 0

    r_core = np.sqrt(np.sum(c["x"]**2, axis=2))
    r_rs = np.sqrt(np.sum(h["x"]**2, axis=2))
    r_core_rs = np.sqrt(np.sum((h["x"] - c["x"])**2, axis=2))

    c["is_err"] = (c["f_core"] == 0) | (c["r50_bound"] > r_core)
    c["is_err_rs"] = (c["f_core_rs"] == 0)|(c["r50_bound_rs"] > r_rs)

    both_ok = ((c["ok"] & (~c["is_err"])) &
               (h["ok"] & (~c["is_err_rs"])))
    disagree = (both_ok & (r_core_rs > c["r50_bound"]) &
                (r_core_rs > c["r50_bound_rs"]))

    disagree_rs_wins = disagree & (c["f_core_rs"] > c["f_core"])
    disagree_core_wins = disagree & (~disagree_rs_wins)
    c["is_err"] = (c["is_err"] | disagree_rs_wins) & c["ok"]
    c["is_err_rs"] = (c["is_err_rs"] | disagree_core_wins) & h["ok"]

    c["ok"] = c["ok"] & (~c["is_err"])
    c["ok_rs"] = h["ok"] & (~c["is_err_rs"])

    has_neighbor = np.zeros(c.shape, dtype=bool)
    has_neighbor[:,:-1] = c["ok"][:,1:]
    has_neighbor[:,1:] = has_neighbor[:,1:] | c["ok"][:,:-1]
    recent_infall = hist["first_infall_snap"] == c.shape[1]-1
    has_neighbor[:,-1] |= recent_infall
    c["ok"] = c["ok"] & has_neighbor

def interpolate_cores(scales, c, h, hist, mp):
    log_a = np.log10(scales)
    # This funciton is kind of cursed, and I'm sorry
    for i in range(1, len(c)):
        first_snap = hist["first_infall_snap"][i]
        if np.sum(c["ok"][i]) == 0:
            last_snap = hist["first_infall_snap"][i]
        else:
            last_snap = np.max(np.where(c["ok"][i])[0])

        h_infall = h[i,first_snap]

        # This almost never happens (0.7% of subhaloes). Only needed to keep
        # interpolation from blowing up. In general, particle-tracking and
        # rockstar get pretty similar halo properties immediately after infall,
        # and the particle tracking is essentially initialized to start out
        # with the same phase space position as the halo because the core
        # particles are calculated based on Rockstar's location.
        # Error and interpolation flags get set to make sure that this still
        # counts as a "loss" for particle-tracking.
        if not c["ok"][i,first_snap]:
            r = np.sqrt(np.sum(h[i,first_snap]["x"]**2))
            rvir = h[i,first_snap]["rvir"]
            cvir = h[i,first_snap]["cvir"]
            m_ratio = h[i,first_snap]["mvir"]/h[0,first_snap]["mvir"]
            r50_rs = _half_mass_nfw(cvir, 0.5)/cvir * rvir
            r_tidal_rs = r*(m_ratio/3)**(1/3.0)
            x_tidal_rs = r_tidal_rs/(rvir/cvir)
            m_tidal_rs = _m_enc_nfw(x_tidal_rs)/_m_enc_nfw(cvir)
            m_tidal_rs = h[i,first_snap]["mvir"]*min(1, m_tidal_rs)
            c[i,first_snap]["x"] = h[i,first_snap]["x"]
            c[i,first_snap]["v"] = h[i,first_snap]["v"]
            c[i,first_snap]["r_tidal"] = r_tidal_rs
            c[i,first_snap]["r50_bound"] = r50_rs
            c[i,first_snap]["r50_bound_rs"] = r50_rs
            c[i,first_snap]["m_tidal"] = m_tidal_rs
            c[i,first_snap]["m_tidal_bound"] = m_tidal_rs
            c[i,first_snap]["m_bound"] = h[i,first_snap]["mvir"]
            c[i,first_snap]["vmax"] = h[i,first_snap]["vmax"]
            c[i,first_snap]["f_core"] = 0
            c[i,first_snap]["f_core_rs"] = 1
            c[i,first_snap]["d_core_mbp"] = 0
            c[i,first_snap]["ok"] = True
            c[i,first_snap]["ok_rs"] = True
            c[i,first_snap]["is_err"] = True
            c[i,first_snap]["is_err_rs"] = False
            c[i,first_snap]["interp"] = True

        n_interp = last_snap-first_snap + 1

        bad_t = c["ok"][i] & (c["m_tidal"][i] <= 0)
        c["m_tidal"][i,bad_t] = mp
        c["r_tidal"][i,bad_t] = 0
        bad_bt = c["ok"][i] & (c["m_tidal_bound"][i] <= 0)
        c["m_tidal_bound"][i,bad_bt] = mp
        
        # If there are no errors, don't bother interpolating
        if np.sum(c["ok"][i]) == n_interp: continue

        snap = np.arange(len(log_a), dtype=int)
        skipped = snap[(~c["ok"][i]) & (snap >= first_snap) &
                       (snap <= last_snap)]

        # Interpolate x with the highest order method you can use
        if n_interp-len(skipped) > 3:
            x_kind = "cubic"
        elif n_interp-len(skipped) > 2:
            x_kind = "quadratic"
        elif n_interp-len(skipped) > 1:
            x_kind = "linear"
        else:
            continue

        # Actually interpolate
        ok = c["ok"][i]
        def intr(x, kind="linear"):
            return interpolate.interp1d(
                log_a[ok], x, kind=kind)(log_a[skipped])

        for dim in range(3):
            c["x"][i,skipped,dim] = intr(c["x"][i,ok,dim], kind=x_kind)
            c["v"][i,skipped,dim] = intr(c["v"][i,ok,dim])
        names = [
            "r_tidal", "r50_bound", "r50_bound_rs", "m_tidal", "m_tidal_bound",
            "m_bound", "vmax", "f_core", "d_core_mbp"
        ]
        for name in names:
            if name in ["m_bound", "m_tidal", "m_tidal_bound"]:
                m = c[name][i,ok]
                m[m <= 0] = mp
                c[name][i,skipped] = 10**intr(np.log10(m))
            else:
                c[name][i,skipped] = intr(c[name][i,ok])
            
        c["ok"][i,skipped] = True
        c["is_err"][i,skipped] = True
        c["interp"][i,skipped] = True

def _f(x): return np.log(1+x) - x/(1+x)
def _x_max_nfw(): return 2.1626
def _v_vmax_nfw(x): return 2.1506 * np.sqrt(_f(x) / x)
def _m_enc_nfw(x): return _f(x)
def _alpha_nfw(x): return -1 - 2*x/(1 + x)

def _half_mass_nfw(x, mass_fraction):
    def f_mass_ratio(xx):
        return _m_enc_nfw(xx) / _m_enc_nfw(x) - mass_fraction
    sol = optimize.root_scalar(f_mass_ratio, bracket=[1e-4*x, x])
    return sol.root

def main():
    config_name, idx_str = sys.argv[1], sys.argv[2]
    flags = parse_flags(sys.argv[3:])

    target_idx = int(idx_str)
    sim_dirs = get_sim_dirs(config_name)

    for i_host in range(len(sim_dirs)):
        if target_idx != -1 and target_idx != i_host: continue

        sim_dir = sim_dirs[i_host]
        if sim_dir[-1] == "/": sim_dir = sim_dir[:-1]

        print("Converting host %d" % i_host, sim_dir)
        convert(sim_dir, flags)

def convert(sim_dir, flags): 
    if "suffix" not in flags:
        file_fmt = path.join(sim_dir, "halos", "core.%d.txt")
    else:
        file_fmt = path.join(sim_dir, "halos", "core_%s.%%d.txt" %
                             flags["suffix"])
    file_names = get_matching_file_names(file_fmt)

    base_dir, suite_name, halo_name = parse_sim_dir(sim_dir)
    param = symlib.parameter_table[suite_name]
    h, hist = symlib.read_subhalos(sim_dir, include_false_selections=True)
    mp = param["mp"]/param["h100"]
    scales = symlib.scale_factors(sim_dir)

    out = np.ones(h.shape, dtype=symlib.CORE_DTYPE)
    out["x"], out["v"] = -1, -1
    out["r_tidal"], out["r50_bound"], out["r50_bound_rs"] = -1, -1, -1
    out["m_tidal"], out["m_tidal_bound"], out["m_bound"] = -1, -1, -1

    for file_name in file_names:
        cols = np.loadtxt(file_name, dtype=int, usecols=(0,1)).T
        if cols.shape == (0,): continue

        snaps, subs = cols
        cols = np.loadtxt(file_name).T
        x, v = cols[2:5].T, cols[5:8].T
        r_tidal, r50_bound, r50_bound_rs = cols[8:11]
        m_tidal, m_tidal_bound, m_bound = cols[11:14]
        vmax, f_core, f_core_rs, d_core_mbp = cols[14:18]

        for i in range(len(subs)):
            sub, snap = subs[i], snaps[i]
            out[sub,snap]["x"], out[sub,snap]["v"] = x[i], v[i]
            out[sub,snap]["r_tidal"] = r_tidal[i]
            out[sub,snap]["m_tidal"] = m_tidal[i]
            out[sub,snap]["r50_bound"] = r50_bound[i]
            out[sub,snap]["r50_bound_rs"] = r50_bound_rs[i]
            out[sub,snap]["m_tidal"] = m_tidal[i]
            out[sub,snap]["m_tidal_bound"] = m_tidal_bound[i]
            out[sub,snap]["m_bound"] = m_bound[i]
            out[sub,snap]["vmax"] = vmax[i]
            out[sub,snap]["f_core"] = f_core[i]
            out[sub,snap]["f_core_rs"] = f_core_rs[i]
            out[sub,snap]["d_core_mbp"] = d_core_mbp[i]

    set_ok_flags(out, h, hist)
    interpolate_cores(scales, out, h, hist, mp)

    if "suffix" not in flags:
        file_root = "cores.dat"
    else:
        file_root = "cores_%s.dat"
        
    file_root = ("cores.dat" if "suffix" not in flags else
                 "cores_%s.dat" % flags["suffix"])
    out_file_name = path.join(sim_dir, "halos", file_root)
    out.reshape(out.shape[0]*out.shape[1])
    with open(out_file_name, "wb") as fp:
        fp.write(struct.pack("qq", out.shape[0], out.shape[1]))
        out["x"].tofile(fp)
        out["v"].tofile(fp)
        out["r_tidal"].tofile(fp)
        out["r50_bound"].tofile(fp)
        out["r50_bound_rs"].tofile(fp)
        out["m_tidal"].tofile(fp)
        out["m_tidal_bound"].tofile(fp)
        out["m_bound"].tofile(fp)
        out["vmax"].tofile(fp)
        out["f_core"].tofile(fp)
        out["f_core_rs"].tofile(fp)
        out["d_core_mbp"].tofile(fp)
        out["ok"].tofile(fp)
        out["ok_rs"].tofile(fp)
        out["is_err"].tofile(fp)
        out["is_err_rs"].tofile(fp)
        out["interp"].tofile(fp)
        
    if "suffix" in flags:
        c = symlib.read_cores(sim_dir, suffix=flags["suffix"])
    else:
        c = symlib.read_cores(sim_dir)

if __name__ == "__main__": main()
