import numpy as np
import matplotlib.pyplot as plt
import symlib
import palette
from palette import pc
from colossus.cosmology import cosmology
from colossus.halo import mass_so

class Matches(object):
    def __init__(self, file_name):
        host, suite, i_lr, i_hr = np.loadtxt(file_name, dtype=int).T

        self.lr_to_hr = []
        self.hr_to_lr = []

        host_max = np.max(host)
        for i_host in range(host_max+1):
            is_lr = (host == i_host) & (suite == 0)
            is_hr = (host == i_host) & (suite == 1)

            self.lr_to_hr.append(i_hr[is_lr])
            self.hr_to_lr.append(i_lr[is_hr])

def pre_infall_mpeak(h, hist):
    mpeak = np.zeros(len(hist))
    for i in range(len(hist)):
        #print(hist["first_infall_snap"][i][)
        mpeak[i] = np.max(h["mvir"][i,:hist["first_infall_snap"][i]])
    return mpeak

def main():
    match = Matches("tables/lr_hr_match.txt")

    base_dir = "/oak/stanford/orgs/kipac/users/phil1/simulations/ZoomIns/"
    suite_lr, suite_hr = "SymphonyMilkyWayLR", "SymphonyMilkyWayHR"

    n_host = symlib.n_hosts(suite_lr)

    min_ratio = 0

    dt = []
    mass_frac = []
    for i in range(n_host):
        dir_lr = symlib.get_host_directory(base_dir, suite_lr, i)
        dir_hr = symlib.get_host_directory(base_dir, suite_hr, i)
        h_lr, hist_lr = symlib.read_subhalos(dir_lr)
        h_hr, hist_hr = symlib.read_subhalos(dir_hr)

        ok = match.lr_to_hr[i] != -1
        i_hr = match.lr_to_hr[i][ok]
        h_lr, hist_lr = h_lr[ok], hist_lr[ok]
        h_hr, hist_hr = h_hr[i_hr], hist_hr[i_hr]
        mpeak_lr = pre_infall_mpeak(h_lr, hist_lr)
        mpeak_hr = pre_infall_mpeak(h_hr, hist_hr)

        ok_lr = mpeak_lr/h_lr["mvir"][0,-1] > min_ratio
        h_lr, hist_lr, mpeak_lr = h_lr[ok_lr], hist_lr[ok_lr], mpeak_lr[ok_lr]
        h_hr, hist_hr, mpeak_hr = h_hr[ok_lr], hist_hr[ok_lr], mpeak_hr[ok_lr]

        scale_lr = symlib.scale_factors(dir_lr)
        scale_hr = symlib.scale_factors(dir_hr)
        param = symlib.simulation_parameters(dir_lr)
        cosmo = cosmology.setCosmology("", symlib.colossus_parameters(param))
        
        t_lr = cosmo.age(1/scale_lr - 1)
        t_hr = cosmo.age(1/scale_hr - 1)
        T_orbit_lr = mass_so.dynamicalTime(1/scale_lr - 1, "vir", "orbit")

        snap_lr = hist_lr["first_infall_snap"]
        snap_hr = hist_hr["first_infall_snap"]
        dt.append((t_hr[snap_hr] - t_lr[snap_lr]) / T_orbit_lr[snap_lr])

        surv = h_lr["ok"][:,-1] & h_hr["ok"][:,-1]
        mass_frac_lr = h_lr["mvir"][:,-1]/mpeak_lr 
        mass_frac_hr = h_hr["mvir"][:,-1]/mpeak_hr 
        mass_frac.append(mass_frac_hr[surv]/mass_frac_lr[surv])

    dt = np.hstack(dt)
    mass_frac = np.hstack(mass_frac)

    print("dt/T_orbit:")
    print("%.4f" % np.percentile(dt, 50))
    print("(68%%: %.4f to %.4f)" % (np.percentile(dt, 50-68/2), np.percentile(dt, 50+68/2)))
    print("(95%%: %.4f to %.4f)" % (np.percentile(dt, 2.5), np.percentile(dt, 97.5)))

    print("mass_frac_hr/mass_frac_lr:")
    print("%.4f" % np.percentile(mass_frac, 50))
    print("(68%%: %.4f to %.4f)" % (np.percentile(mass_frac, 50-68/2), np.percentile(mass_frac, 50+68/2)))
    print("(95%%: %.4f to %.4f)" % (np.percentile(mass_frac, 2.5), np.percentile(mass_frac, 97.5)))
        

if __name__ == "__main__": main()
