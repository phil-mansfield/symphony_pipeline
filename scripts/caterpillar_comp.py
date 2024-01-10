import numpy as np
import symlib

base_dir = "/sdf/home/p/phil1/ZoomIns"
suite = "SymphonyMilkyWayHR"

n_hosts = symlib.n_hosts(suite)

n_rct, n_sym = 0, 0
np_rct, np_sym = 0, 0
na_rct, na_sym = 0, 0

for i_host in range(n_hosts):
    print(suite, i_host)
    sim_dir = symlib.get_host_directory(base_dir, suite, i_host)
    param = symlib.simulation_parameters(sim_dir)
    mp = param["mp"]/param["h100"]
    
    h, hist = symlib.read_subhalos(sim_dir)
    c = symlib.read_cores(sim_dir)
    
    n_rct += np.sum((h["mvir"][1:,-1] > 1e8) & c["ok_rs"][1:,-1] & h["ok"][1:,-1])
    n_sym += np.sum((c["m_bound"][1:,-1] > 1e8) & c["ok"][1:,-1])
    np_rct += np.sum((hist["mpeak_pre"][1:] > 1e8) & c["ok_rs"][1:,-1])
    np_sym += np.sum((hist["mpeak_pre"][1:] > 1e8) & c["ok"][1:,-1])
    na_rct += np.sum(c["ok_rs"][1:,-1] & h["ok"][1:,-1])
    na_sym += np.sum(c["ok"][1:,-1])



print(n_rct/n_hosts)
print(n_sym/n_hosts)
print(np_rct/n_hosts)
print(np_sym/n_hosts)
print(na_rct/n_hosts)
print(na_sym/n_hosts)
