import symlib
import numpy as np

def get_all_halos(sim_dir, snap):
    param = symlib.simulation_parameters(sim_dir)
    h100 = param["h100"]
    
    rs, hist = symlib.read_rockstar(sim_dir)
    sf, _ = symlib.read_symfind(sim_dir)
    rs_cmov, _ = symlib.read_rockstar(sim_dir, comoving=True)
    b = symlib.read_branches(sim_dir)

    ids, snaps, x, mvir = symlib.read_tree(
        sim_dir, ["id", "snap", "x", "mvir"])
    mvir /= h100
    
    mpeak = np.zeros(len(b))
    start = b["start"]
    start_snap = snaps[start]
    infall = start_snap - b["first_infall_snap"] + start
    has_infall = b["first_infall_snap"] != -1
    start[has_infall] = infall[has_infall]

    mpeak_pre = np.zeros(b.shape)
    for i in range(len(b)):
        mpeak_pre[i] = np.max(mvir[start[i]: b["end"][i]])

    dx = (x - rs_cmov["x"][0,snap])/h100 * symlib.scale_factors(sim_dir)[snap]
    ok = snaps == snap
    ok_tree = ok
    idx = np.where(ok)[0]
    #branch_idx = np.searchsorted(b["start"], idx, side="right") - 1
    branch_idx = np.searchsorted(b["end"], idx, side="right")
    branch_ok = np.zeros(b.shape, dtype=bool)
    branch_ok[branch_idx] = True

    problem_idx = np.where(branch_idx[1:] == branch_idx[:-1])[0]
    ids, dx, mvir, mpeak_pre = ids[ok], dx[ok], mvir[ok], mpeak_pre[branch_ok]
    
    dtype = [("x", "f4", (3,)), ("m", "f4"), 
             ("mpeak", "f4"), ("is_tracked", "?")]
    out_1 = np.zeros(len(sf), dtype=dtype)
    out_1["x"] = sf["x"][:,snap]
    out_1["m"] = sf["m"][:,snap]
    out_1["mpeak"] = hist["mpeak_pre"]
    out_1["is_tracked"] = True
    out_1["x"][0] = rs["x"][0,snap]
    out_1["m"][0] = rs["m"][0,snap]
    out_1["mpeak"][0] = hist["mpeak_pre"][0]

    ok = sf["ok"][:,snap]
    ok[0] = True
    out_1 = out_1[ok]

    tracked_ids = np.sort(rs["id"][rs["ok"][:,snap],snap])

    out_2 = np.zeros(len(ids), dtype=dtype)
    out_2["x"], out_2["m"], out_2["mpeak"] = dx, mvir, mpeak_pre
    idx = np.searchsorted(tracked_ids, ids)
    idx[idx >= len(tracked_ids)] = len(tracked_ids) - 1
    ok = ~(tracked_ids[idx] == ids)
    out_2 = out_2[ok]
    
    return np.hstack([out_1, out_2])
    
def main():
    base_dir = "/sdf/home/p/phil1/ZoomIns"
    sim_dir = symlib.get_host_directory(base_dir, "SymphonyMilkyWay", 0)

    h = get_all_halos(sim_dir, 235)
    r = np.sqrt(np.sum(h["x"]**2, axis=1))
    print(r[:10])
    print(h["mpeak"][:10])
    print(h["m"][:10])
    print(np.sum(h["is_tracked"]), np.sum(~h["is_tracked"]))

if __name__ == "__main__": main()
