import numpy as np
import matplotlib.pyplot as plt
import palette
from palette import pc
import symlib
import os.path as path
from colossus.cosmology import cosmology
from colossus.halo import mass_so
import scipy.signal as signal
import scipy.stats as stats
import os.path as path

ratio_limit = 0.1
suite = "SymphonyMilkyWay"
n_hr_max = 5
suite_lr = "SymphonyMilkyWayLR"
suite_hr = "SymphonyMilkyWayHR"
base_dir = "/sdf/home/p/phil1/ZoomIns"
plot_dir_fmt = "/scratch/phil1/plots/tracking/mah/h%d"

def plot_percentiles(ax, ls, color, x, y, label=None):
    bins, n_min = 20, 50
    n, edges, _ = stats.binned_statistic(np.log10(x), y, "count", bins=bins)
    hi, _, _ = stats.binned_statistic(
        np.log10(x), y, lambda xx: np.percentile(xx, 50+68/2), bins=bins)
    mi, _, _ = stats.binned_statistic(
        np.log10(x), y, lambda xx: np.percentile(xx, 50), bins=bins)
    lo, _, _ = stats.binned_statistic(
        np.log10(x), y, lambda xx: np.percentile(xx, 50-68/2), bins=bins)
    ok = n > n_min

    mid = 10**((edges[1:] + edges[:-1])/2)

    ax.plot(mid[ok], mi[ok], ls, c=color, label=label)
    #ax.plot(mid[ok], hi[ok], "--", c=color, lw=1.5)
    #ax.plot(mid[ok], lo[ok], "--", c=color, lw=1.5)

def mass_loss_vars():
    n_method = 2
    method_colors = [pc("r"), pc("b")]

    both_exist = [[] for _ in range(n_method)]
    dm_dt = [[] for _ in range(n_method)]
    z_sub = [[] for _ in range(n_method)]
    msub_mhost = [[] for _ in range(n_method)]
    msub_mpeak = [[] for _ in range(n_method)]
    csub_chost = [[] for _ in range(n_method)]
    rsub_rhost = [[] for _ in range(n_method)]
    lnm_sub = [[] for _ in range(n_method)]

    for i_host in range(symlib.n_hosts(suite)):
        if "HR" in suite:
            cut_mult = 8
        else:
            cut_mult = 1
        plot_dir = plot_dir_fmt % i_host
        host_dir = symlib.get_host_directory(base_dir, suite, i_host)
        param = symlib.simulation_parameters(host_dir)
        cosmo = cosmology.setCosmology('', symlib.colossus_parameters(param))
        h, hist = symlib.read_subhalos(host_dir)

        c = symlib.read_cores(host_dir)
        scale = symlib.scale_factors(host_dir)
        z = 1/scale - 1
        t = cosmo.age(z)
        T = mass_so.dynamicalTime(z, "vir", "crossing")
        
        mp = param["mp"]/param["h100"]
        ok = hist["mpeak"]/mp > 300*cut_mult
        c, h, hist = c[ok], h[ok], hist[ok]
        
        smooth_order = 4
        smooth_window = 15

        for i in range(1, len(h)):
            R_host, m_host, ok_host = h["rvir"][0], h["mvir"][0], h["ok"][0]
            m_rs, m_c = h["mvir"][i], c["m_bound"][i]
            ok1_rs, ok2_rs, ok_c = h["ok"][i], c["ok_rs"][i], c["ok"][i]
            ok1_rs, ok2_rs = ok1_rs & h["ok"][0,:], ok2_rs & h["ok"][0,:]
            r_rs = np.sqrt(np.sum(h["x"][i]**2, axis=1))
            r_c = np.sqrt(np.sum(c["x"][i]**2, axis=1))

            snap = np.arange(h.shape[1], dtype=int)
            i_start = hist["merger_snap"][i]
            is_sub_rs = ok2_rs & (snap >= i_start)

            is_sub_method = [
                ok2_rs & (snap >= i_start),
                ok_c,
            ]
            m_method = [m_rs, m_c]
            r_method = [r_rs, r_c]

            for i_method in range(n_method):
                r, m = r_method[i_method], m_method[i_method]
                is_sub = is_sub_method[i_method]

                if (np.sum(is_sub) < smooth_window or
                    hist["merger_ratio"][i] > ratio_limit): continue
            
                lnm_smooth = signal.savgol_filter(
                    np.log10(m[is_sub]), smooth_window, smooth_order,
                    mode="nearest")[1:-1]
                dlnm_smooth = signal.savgol_filter(
                    np.log(m[is_sub]), smooth_window, smooth_order,
                    deriv=1, delta=1.0, mode="nearest")[1:-1]
                dt = (t[is_sub][2:] - t[is_sub][:-2])/2
                dm_dt_i = dlnm_smooth/dt * T[is_sub][1:-1]

                both_exist_i = ok2_rs[is_sub][1:-1]
                z_sub_i = z[is_sub][1:-1]
                msub_mhost_i = (m[is_sub]/m_host[is_sub])[1:-1]
                msub_mpeak_i = (m[is_sub]/np.max(m_rs[:i_start]))[1:-1]
                csub_chost_i = h[i,i_start]["cvir"]*np.ones(len(msub_mhost_i))
                rsub_rhost_i = (r[is_sub]/h[0,is_sub]["rvir"])[1:-1]
                lnm_sub_i = lnm_smooth

                both_exist[i_method].append(both_exist_i)
                dm_dt[i_method].append(dm_dt_i)
                z_sub[i_method].append(z_sub_i)
                msub_mhost[i_method].append(msub_mhost_i)
                msub_mpeak[i_method].append(msub_mpeak_i)
                csub_chost[i_method].append(csub_chost_i)
                rsub_rhost[i_method].append(rsub_rhost_i)
                lnm_sub[i_method].append(lnm_sub_i)

                #if np.sum(z_sub_i > 3) > 10:
                #    print("   ", i_host, i)

    for i in range(n_method):
        both_exist[i] = np.hstack(both_exist[i])
        dm_dt[i] = np.hstack(dm_dt[i])
        z_sub[i] = np.hstack(z_sub[i])
        msub_mhost[i] = np.hstack(msub_mhost[i])
        msub_mpeak[i] = np.hstack(msub_mpeak[i])
        csub_chost[i] = np.hstack(csub_chost[i])
        rsub_rhost[i] = np.hstack(rsub_rhost[i])
        lnm_sub[i] = np.hstack(lnm_sub[i])

        ok = (msub_mhost[i] < 0.1) & (msub_mpeak[i] < 1.0) #& (z_sub[i] > 3) & both_exist[i]
        dm_dt[i] = dm_dt[i][ok]
        z_sub[i] = z_sub[i][ok]
        msub_mhost[i] = msub_mhost[i][ok]
        msub_mpeak[i] = msub_mpeak[i][ok]
        csub_chost[i] = csub_chost[i][ok]
        rsub_rhost[i] = rsub_rhost[i][ok]
        lnm_sub[i] = lnm_sub[i][ok]

    fig, axs = plt.subplots(2, 3, figsize=(24, 16))
    ax_z, ax_msub_mhost, ax_msub_mpeak = axs[0,0], axs[0,1], axs[0,2]
    ax_csub_chost, ax_rsub_rhost = axs[1,0], axs[1,1]
    fig.delaxes(axs[1,2])

    for i in range(n_method):
        c = method_colors[i]
        plot_percentiles(ax_z, "-", c, z_sub[i]+1, dm_dt[i], None)
        plot_percentiles(ax_msub_mhost, "-", c, msub_mhost[i], dm_dt[i], None)
        plot_percentiles(ax_msub_mpeak, "-", c, msub_mpeak[i], dm_dt[i], None)
        
        plot_percentiles(ax_csub_chost, "-", c, csub_chost[i], dm_dt[i], None)
        plot_percentiles(ax_rsub_rhost, "-", c, rsub_rhost[i], dm_dt[i], None)

    axs[0,0].set_ylabel(r"$\dot{m}_{\rm sub}/(m_{\rm sub}/\tau_{\rm crossing})$")
    axs[1,0].set_ylabel(r"$\dot{m}_{\rm sub}/(m_{\rm sub}/\tau_{\rm crossing})$")

    ax_z.set_xlabel(r"$1+z$")
    ax_msub_mhost.set_xlabel(r"$m_{\rm sub}/M_{\rm vir}$")
    ax_msub_mpeak.set_xlabel(r"$m_{\rm sub}/m_{\rm sub,peak}$")
    ax_csub_chost.set_xlabel(r"$c_{\rm sub,infall}$")
    ax_rsub_rhost.set_xlabel(r"$r_{\rm sub}/R_{\rm vir}$")

    ax_z.set_xscale("log")
    ax_msub_mhost.set_xscale("log")
    ax_msub_mpeak.set_xscale("log")
    ax_csub_chost.set_xscale("log")
    ax_rsub_rhost.set_xscale("log")


    ax_z.set_xlim()
    ax_msub_mhost.set_ylim((-3, 0))
    ax_msub_mpeak.set_ylim((-3, 0))
    ax_csub_chost.set_ylim((-3, 0))
    ax_rsub_rhost.set_ylim((-3, 0))

    fig.savefig("../plots/core_tracking/mass_loss_vars_%s.pdf" % suite)

def mass_loss_msub_mpeak():
    n_method = 4
    colors = [pc("r"), pc("b"), pc("o"), pc("b")]
    panel = [0, 0, 1, 1]
    labels = [r"${\rm Rockstar}$", r"${\rm Particle{-}tracking}$",
              r"{\rm Fiducial}", r"${\rm High{-}res}$"]
    dm_dt = [[] for _ in range(n_method)]
    msub_mpeak = [[] for _ in range(n_method)]
    msub_mhost = [[] for _ in range(n_method)]

    for i_host in range(symlib.n_hosts(suite)):
        print(i_host, "out of", symlib.n_hosts(suite))
        plot_dir = plot_dir_fmt % i_host
        host_dir = symlib.get_host_directory(base_dir, suite, i_host)
        param = symlib.simulation_parameters(host_dir)
        cosmo = cosmology.setCosmology('', symlib.colossus_parameters(param))
        h, hist = symlib.read_subhalos(host_dir)

        c = symlib.read_cores(host_dir)
        scale = symlib.scale_factors(host_dir)
        z = 1/scale - 1
        t = cosmo.age(z)
        T = mass_so.dynamicalTime(z, "vir", "crossing")
        
        mp = param["mp"]/param["h100"]
        ok = hist["mpeak"]/mp > 300*8
        c, h, hist = c[ok], h[ok], hist[ok]

        if i_host < n_hr_max:
            host_dir = symlib.get_host_directory(base_dir, suite_hr, i_host)
            param = symlib.simulation_parameters(host_dir)
            cosmo = cosmology.setCosmology(
                '', symlib.colossus_parameters(param))
            h_hr, hist_hr = symlib.read_subhalos(host_dir)
            
            c_hr = symlib.read_cores(host_dir)
            scale = symlib.scale_factors(host_dir)

            mp = param["mp"]/param["h100"]
            ok = hist_hr["mpeak"]/mp > 300*8
            c_hr, h_hr, hist_hr = c_hr[ok], h_hr[ok], hist_hr[ok]

        smooth_order = 4
        smooth_window = 15

        for i in range(1, len(h)):
            R_host, m_host, ok_host = h["rvir"][0], h["mvir"][0], h["ok"][0]
            m_rs, m_c = h["mvir"][i], c["m_bound"][i]
            ok1_rs, ok2_rs, ok_c = h["ok"][i], c["ok_rs"][i], c["ok"][i]
            ok1_rs, ok2_rs = ok1_rs & h["ok"][0,:], ok2_rs & h["ok"][0,:]

            if i_host < n_hr_max:
                ok_c_hr = c_hr["ok"][i]
                m_hr = c_hr["m_bound"][i]
            else:
                ok_c_hr = None
                m_hr = None

            snap = np.arange(h.shape[1], dtype=int)
            i_start = hist["merger_snap"][i]

            is_sub_method = [
                ok2_rs & (snap >= i_start),
                ok_c, ok_c, ok_c_hr
            ]
            m_method = [m_rs, m_c, m_c, m_hr]

            for i_method in range(n_method):
                if ((i_method == 3 or i_method == 2) and
                    i_host >= n_hr_max): continue

                m = m_method[i_method]
                is_sub = is_sub_method[i_method]

                if (np.sum(is_sub) < smooth_window or
                    hist["merger_ratio"][i] > ratio_limit): continue
            
                dlnm_smooth = signal.savgol_filter(
                    np.log(m[is_sub]), smooth_window, smooth_order,
                    deriv=1, delta=1.0, mode="interp")[:-1]
                #dlnm_smooth = signal.savgol_filter(
                #    np.log(m[is_sub]), smooth_window, smooth_order,
                #    deriv=1, delta=1.0)[:-1]
                dt = t[is_sub][1:] - t[is_sub][:-1]
                dm_dt_i = dlnm_smooth/dt * T[is_sub][:-1]
                msub_mpeak_i = (m[is_sub]/np.max(m_rs[:i_start]))[:-1]
                msub_mhost_i = (m[is_sub]/m_host[is_sub])[:-1]

                dm_dt[i_method].append(dm_dt_i)
                msub_mhost[i_method].append(msub_mhost_i)
                msub_mpeak[i_method].append(msub_mpeak_i)

        #break

    for i in range(n_method):
        dm_dt[i] = np.hstack(dm_dt[i])
        msub_mpeak[i] = np.hstack(msub_mpeak[i])
        msub_mhost[i] = np.hstack(msub_mhost[i])

        ok = (msub_mhost[i] < 0.1) & (msub_mpeak[i] < 1.0)
        dm_dt[i] = dm_dt[i][ok]
        msub_mpeak[i] = msub_mpeak[i][ok]
        msub_mhost[i] = msub_mhost[i][ok]

    fig, ax = plt.subplots(1, 2, figsize=(16, 8), sharey=True)

    plot_percentiles(ax[0], pc("r"), msub_mpeak[0], dm_dt[0], labels[0])
    plot_percentiles(ax[0], pc("b"), msub_mpeak[1], dm_dt[1], labels[1])
    plot_percentiles(ax[1], pc("r"), msub_mpeak[2], dm_dt[2], labels[2])
    plot_percentiles(ax[1], pc("b"), msub_mpeak[3], dm_dt[3], labels[3])
    
    ax[0].set_ylabel(r"$\dot{m}_{\rm sub}/(m_{\rm sub}/\tau_{\rm crossing})$")
    ax[0].set_xscale("log")
    ax[1].set_xscale("log")
    ax[0].set_xlabel(r"$m_{\rm sub}/m_{\rm peak}$")
    ax[1].set_xlabel(r"$m_{\rm sub}/m_{\rm peak}$")

    ax[0].set_xlim(1e-3, 1)
    ax[1].set_xlim(1e-3, 1)
    ax[0].set_ylim(-3, 0)
    ax[1].set_ylim(-3, 0)

    ax[0].set_ylim((-3, 0))
    ax[1].set_ylim((-3, 0))

    ax[0].legend(loc="lower right", fontsize=17)
    ax[1].legend(loc="lower right", fontsize=17)

    fig.savefig("../plots/core_tracking/mass_loss_msub_mpeak.pdf")

def mass_loss_msub_mpeak_2():
    n_method = 4
    colors = [pc("r"), pc("b"), pc("o"), pc("b")]
    panel = [0, 0, 1, 1]
    labels = [r"${\rm Rockstar}$", r"${\rm Particle{-}tracking}$",
              r"{\rm Fiducial}", r"${\rm High{-}res}$"]
    dm_dt = [[] for _ in range(n_method)]
    msub_mpeak = [[] for _ in range(n_method)]
    mpeak = [[] for _ in range(n_method)]

    for i_host in range(symlib.n_hosts(suite)):
        print(i_host, "out of", symlib.n_hosts(suite))
        plot_dir = plot_dir_fmt % i_host
        host_dir = symlib.get_host_directory(base_dir, suite, i_host)
        param = symlib.simulation_parameters(host_dir)
        cosmo = cosmology.setCosmology('', symlib.colossus_parameters(param))
        h, hist = symlib.read_subhalos(host_dir)

        c = symlib.read_cores(host_dir)
        scale = symlib.scale_factors(host_dir)
        z = 1/scale - 1
        t = cosmo.age(z)
        T = mass_so.dynamicalTime(z, "vir", "crossing")

        if i_host < n_hr_max:
            host_dir = symlib.get_host_directory(base_dir, suite_hr, i_host)
            param = symlib.simulation_parameters(host_dir)
            cosmo = cosmology.setCosmology(
                '', symlib.colossus_parameters(param))
            h_hr, hist_hr = symlib.read_subhalos(host_dir)
            c_hr = symlib.read_cores(host_dir)

        smooth_order = 4
        smooth_window = 15

        for i in range(1, len(h)):
            R_host, m_host, ok_host = h["rvir"][0], h["mvir"][0], h["ok"][0]
            m_rs, m_c = h["mvir"][i], c["m_bound"][i]
            ok1_rs, ok2_rs, ok_c = h["ok"][i], c["ok_rs"][i], c["ok"][i]
            ok1_rs, ok2_rs = ok1_rs & h["ok"][0,:], ok2_rs & h["ok"][0,:]

            if i_host < n_hr_max:
                ok_c_hr = c_hr["ok"][i]
                m_hr = c_hr["m_bound"][i]
            else:
                ok_c_hr = None
                m_hr = None

            snap = np.arange(h.shape[1], dtype=int)
            i_start = hist["first_infall_snap"][i]

            is_sub_method = [
                ok2_rs & (snap >= i_start),
                ok_c, ok_c, ok_c_hr
            ]
            m_method = [m_rs, m_c, m_c, m_hr]

            for i_method in range(n_method):
                if ((i_method == 3 or i_method == 2) and
                    i_host >= n_hr_max): continue

                m = m_method[i_method]
                is_sub = is_sub_method[i_method]

                if (np.sum(is_sub) < smooth_window or
                    hist["merger_ratio"][i] > ratio_limit): continue
            
                dlnm_smooth = signal.savgol_filter(
                    np.log(m[is_sub]), smooth_window, smooth_order,
                    deriv=1, delta=1.0, mode="interp")[1:-1]

                dt = (t[is_sub][2:] - t[is_sub][:-2]) / 2
                # Note that this is also scaled by msub now
                dm_dt_i = dlnm_smooth/dt * T[is_sub][1:-1]
                msub_mpeak_i = (m[is_sub]/np.max(m_rs[:i_start]))[1:-1]
                mpeak_i = np.max(m_rs[:i_start])*np.ones(len(msub_mpeak_i))

                ok = msub_mpeak_i <= 1
                
                dm_dt[i_method].append(dm_dt_i[ok])
                msub_mpeak[i_method].append(msub_mpeak_i[ok])
                mpeak[i_method].append(mpeak_i[ok])

    for i in range(n_method):
        dm_dt[i] = np.hstack(dm_dt[i])
        msub_mpeak[i] = np.hstack(msub_mpeak[i])
        mpeak[i] = np.hstack(mpeak[i])
        print(mpeak[i])

    fig, ax = plt.subplots(1, 2, figsize=(16, 8), sharey=True)

    colors = [pc("r"), pc("o"), pc("b"), pc("p")]
    #lims = [(10**2.5, 10**3), (10**3, 10**3.5), (10**3.5, 10**4), (10**4, 10**4.5)]
    lims = [(10**2.5, 10**6), (10**3, 10**6), (10**3.5, 10**6), (10**4, 10**6)]

    axs = [0, 0, 1, 1]
    ls = ["-", "--", "--", "-"]

    mp = param["mp"]/param["h100"]
    for i in range(len(lims)):
        c = colors[i]
        for j in range(len(msub_mpeak)):
            ok = (mpeak[j]/mp >= lims[i][0]) & (mpeak[j]/mp < lims[i][1])
            if j == 3:
                ok = (mpeak[j]/mp >= 8*lims[i][0]) & (mpeak[j]/mp < 8*lims[i][1])
            plot_percentiles(ax[axs[j]], ls[j], c, msub_mpeak[j][ok],
                             dm_dt[j][ok], None)

    
    ax[0].set_ylabel(r"$\dot{m}_{\rm sub}/(m_{\rm sub}/\tau_{\rm crossing})$")
    ax[0].set_xscale("log")
    ax[1].set_xscale("log")
    ax[0].set_xlabel(r"$m_{\rm sub}/m_{\rm peak}$")
    ax[1].set_xlabel(r"$m_{\rm sub}/m_{\rm peak}$")

    ax[0].set_xlim(1e-3, 1)
    ax[1].set_xlim(1e-3, 1)
    ax[0].set_ylim(-3, 0)
    ax[1].set_ylim(-3, 0)

    ax[0].set_ylim((-3, 0))
    ax[1].set_ylim((-3, 0))

    #ax[0].legend(loc="lower right", fontsize=17)
    #ax[1].legend(loc="lower right", fontsize=17)

    fig.savefig("../plots/core_tracking/mass_loss_msub_mpeak_2.pdf")

def polyfit(x, y, deg):
    p = np.polyfit(x, y, deg)
    pd = np.copy(p)
    for i in range(len(pd)):
        pd[i] *= len(pd) - i - 1
    return p, pd[:-1]

def f_poly(x, p):
    y = np.zeros(x.shape)
    deg = len(p) - 1
    for i in range(deg+1):
        y += p[i]*x**(deg-i)
    return y

def dlnm_dlna_to_m_dot(dlnm_dlna, a, cosmo):
    z = 1/a - 1
    t = cosmo.age(z)
    z_dot = cosmo.age(t, derivative=1, inverse=True)
    return -a * dlnm_dlna * z_dot

def m_min(m):
    m = np.copy(m)
    m_lim = m[-1]
    for i in range(len(m)-2,-1,-1):
        m_lim = max(m_lim, m[i])
        m[i] = m_lim

    return m

def indiv_mass_loss():
    targets = [
        (0, [38, 43, 45, 54, 60, 61, 64, 68, 75, 97])
    ]
    out_dir = "../plots/core_tracking/mass_loss_tests"
    
    for i_host, idxs in targets:
        host_dir = symlib.get_host_directory(base_dir, suite, i_host)
        param = symlib.simulation_parameters(host_dir)
        cosmo = cosmology.setCosmology('', symlib.colossus_parameters(param))
        h, hist = symlib.read_subhalos(host_dir)

        c = symlib.read_cores(host_dir)
        scale = symlib.scale_factors(host_dir)
        z = 1/scale - 1
        t = cosmo.age(z)
        T = mass_so.dynamicalTime(z, "vir", "crossing")
        
        fig, ax = plt.subplots(2, 2, sharex=True, figsize=(16, 16))
        m_ax, dm_raw_ax = ax[0,0], ax[0,1]
        dm_m_ax, dm_m_t_ax = ax[1,0], ax[1,1]
        log_a = np.log10(scale)

        smooth_window=15
        smooth_order=5

        for i in idxs:
            m_ax.cla()
            dm_raw_ax.cla()
            dm_m_ax.cla()
            dm_m_t_ax.cla()

            m_rs, m_c = h["mvir"][i], c["m_bound"][i]
            ok_rs, ok_c = c["ok_rs"][i], c["ok"][i]
            ok_rs, ok_c = ok_rs & h["ok"][0,:], ok_c & h["ok"][0,:]

            #m_rs, m_c = m_min(m_rs), m_min(m_c)

            snap = np.arange(h.shape[1], dtype=int)
            i_start = hist["first_infall_snap"][i]

            is_sub_rs = ok_rs & (snap >= i_start)
            is_sub_c = ok_c

            zp1_rs = (z+1)[is_sub_rs][1:-1]
            zp1_c = (z+1)[is_sub_c][1:-1]
            lnm_smooth_rs = signal.savgol_filter(
                np.log10(m_rs[is_sub_rs]), smooth_window, smooth_order,
                mode="nearest")[1:-1]
            lnm_smooth_c = signal.savgol_filter(
                np.log10(m_c[is_sub_c]), smooth_window, smooth_order,
                mode="nearest")[1:-1]
            dlnm_smooth_rs = signal.savgol_filter(
                np.log(m_rs[is_sub_rs]), smooth_window, smooth_order,
                deriv=1, delta=1.0, mode="nearest")[1:-1]
            dlnm_smooth_c = signal.savgol_filter(
                np.log(m_c[is_sub_c]), smooth_window, smooth_order,
                deriv=1, delta=1.0, mode="nearest")[1:-1]

            dt_rs = (t[is_sub_rs][2:] - t[is_sub_rs][:-2]) / 2
            dt_c = (t[is_sub_c][2:] - t[is_sub_c][:-2]) / 2            

            da_rs = (log_a[is_sub_rs][2:] - log_a[is_sub_rs][:-2]) / 2
            da_c = (log_a[is_sub_c][2:] - log_a[is_sub_c][:-2]) / 2
            dm_da_rs = dlnm_smooth_rs/da_rs
            dm_da_c = dlnm_smooth_c/da_c
            dm_dt_rs = dlnm_smooth_rs/dt_rs
            dm_dt_c = dlnm_smooth_c/dt_c

            log_a_rs, log_a_c = log_a[is_sub_rs][1:-1], log_a[is_sub_c][1:-1]

            deg = 3
            p_rs, pd_rs = polyfit(log_a_rs,np.log10(m_rs[is_sub_rs][1:-1]),deg)
            p_c, pd_c = polyfit(log_a_c,np.log10(m_c[is_sub_c][1:-1]),deg)

            m_poly_rs = 10**f_poly(log_a_rs, p_rs)
            dm_poly_rs = f_poly(log_a_rs, pd_rs)
            m_poly_c = 10**f_poly(log_a_c, p_c)
            dm_poly_c = f_poly(log_a_c, pd_c)

            m_ax.plot(log_a[is_sub_rs], m_rs[is_sub_rs], pc("r", 0.3))
            m_ax.plot(log_a[is_sub_c], m_c[is_sub_c], pc("b", 0.3))

            m_ax.plot(log_a_rs, 10**lnm_smooth_rs, pc("r"))
            m_ax.plot(log_a_c, 10**lnm_smooth_c, pc("b"))

            m_ax.plot(log_a_rs, m_poly_rs, pc("r", 0.7))
            m_ax.plot(log_a_c, m_poly_c, pc("b", 0.7))

            dm_m_ax.plot(log_a_rs, dm_da_rs/2, pc("r"))
            dm_m_ax.plot(log_a_c, dm_da_c/2, pc("b"))
            dm_m_ax.plot(log_a_rs, dm_poly_rs, pc("r", 0.7))
            dm_m_ax.plot(log_a_c, dm_poly_c, pc("b", 0.7))

            m_dot_poly_rs = dlnm_dlna_to_m_dot(
                dm_poly_rs, 10**log_a_rs, cosmo)
            m_dot_poly_c = dlnm_dlna_to_m_dot(
                dm_poly_c, 10**log_a_c, cosmo)

            dm_raw_ax.plot(log_a_rs, dm_dt_rs, pc("r"))
            dm_raw_ax.plot(log_a_c, dm_dt_c, pc("b"))
            dm_raw_ax.plot(log_a_rs, m_dot_poly_rs, pc("r", 0.7))
            dm_raw_ax.plot(log_a_c, m_dot_poly_c, pc("b", 0.7))

            T_rs = T[is_sub_rs][1:-1]
            T_c = T[is_sub_c][1:-1]
            dm_m_t_ax.plot(log_a_rs, dm_dt_rs*T_rs, pc("r"))
            dm_m_t_ax.plot(log_a_c, dm_dt_c*T_c, pc("b"))
            dm_m_t_ax.plot(log_a_rs, m_dot_poly_rs*T_rs, pc("r", 0.7))
            dm_m_t_ax.plot(log_a_c, m_dot_poly_c*T_c, pc("b", 0.7))


            m_ax.set_ylabel(r"$m_{\rm sub}$")
            dm_m_ax.set_ylabel(r"$d\,\log{(m_{\rm sub})}/d\,\log{(a)}$")
            dm_m_ax.set_xlabel(r"$\log_{10}(a)$")
            dm_raw_ax.set_ylabel(r"$\dot{m}_{\rm sub}/m_{\rm sub}\ ({\rm Gyr}^{-1})$")
            dm_m_t_ax.set_ylabel(r"$\dot{m}_{\rm sub}/(m_{\rm sub}/\tau_{\rm crossing})\ ({\rm Gyr}^{-1})$")
            dm_m_t_ax.set_xlabel(r"$\log_{10}(a)$")

            m_ax.set_yscale("log")
            fig.suptitle(r"$N_{\rm rs}=%d,\ N_c=%d$" %
                         (np.sum(is_sub_rs), np.sum(is_sub_c)))
            fig.savefig(path.join(out_dir, "sub_%03d.png" % i))

if __name__ == "__main__":
    palette.configure(True)
    #mass_loss_vars()
    #mass_loss_msub_mpeak()
    mass_loss_msub_mpeak_2()
    #indiv_mass_loss()
