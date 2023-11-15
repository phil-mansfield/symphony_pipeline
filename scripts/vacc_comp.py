import numpy as np
import matplotlib.pyplot as plt
import symlib

def main():
    suite = "SymphonyMilkyWay"
    base_dir = "/sdf/home/p/phil1/ZoomIns"

    vacc, vpeak, vmpeak = [], [], []
    macc, mpeak = [], []

    n_hosts = symlib.n_hosts(suite)
    for i_host in range(n_hosts):
        sim_dir = symlib.get_host_directory(base_dir, suite, i_host)
        h, hist = symlib.read_subhalos(sim_dir)

        for i_sub in range(1, len(h)):
            s0 = hist["first_infall_snap"][i_sub]
            macc.append(h["mvir"][i_sub,s0])
            vacc.append(h["vmax"][i_sub,s0])
            mpeak.append(np.max(h["mvir"][i_sub,:s0+1]))
            vpeak.append(np.max(h["vmax"][i_sub,:s0+1]))
            i_mpeak = np.argmax(h["mvir"][i_sub,:s0+1])
            vmpeak.append(h["vmax"][i_sub,i_mpeak])
            
    macc = np.array(macc)
    vacc = np.array(vacc)
    vpeak = np.array(vpeak)
    mpeak = np.array(mpeak)
    vmpeak = np.array(vmpeak)

    print(np.median(vpeak/vmpeak))
    print(np.median(vpeak/vacc))
    print(np.median(mpeak/macc))

if __name__ == "__main__": main()
