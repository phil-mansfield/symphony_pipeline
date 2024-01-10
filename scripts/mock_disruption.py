import numpy as np
import matplotlib.pyplot as plt
import palette
from palette import pc
import numpy.random as random
import scipy.interpolate as interpolate
import survival

def kaplan_meier_mass_loss(ts, ms, t_eval,
                           m_eval=10**np.linspace(-5, 0, 250)):
    """ ms needs to be m/m0 along with whatever smoothing
    you want to do.
    """

    for i in range(len(ts)):
        if len(ms[i]) != len(ts[i]) or len(ms) < 2:
            raise ValueError("len(ts[%d]) = %d, len(ms[%d]) = %d" % 
                             (i, len(ts), i, len(ms)))

    fs = [interpolate.interp1d(t, m, fill_value=(m[0], m[-1]),
                               bounds_error=False) for t, m in zip(ts, ms)]
    
    cen = np.zeros((len(ts), len(t_eval)), dtype=bool)
    m_final = np.zeros((len(ts), len(t_eval)))

    for i in range(len(ts)):
        cen[i,:] = np.max(ts[i]) < t_eval
        m_final[i,:] = fs[i](t_eval)

    med = np.zeros(len(t_eval))
    err_low = np.zeros(len(t_eval))
    err_high = np.zeros(len(t_eval))
    low = np.zeros(len(t_eval))
    high = np.zeros(len(t_eval))            
    ok = np.ones(len(t_eval), dtype=bool)

    for i in range(len(t_eval)):
        if np.sum(~cen[:,i]) < 10:
            ok[i] = False
            continue
        S, err = survival.kaplan_meier(
            m_final[:,i], cen[:,i], m_eval, decreasing=True)
        Si = S[~np.isnan(S)]
        if np.min(Si) > 0.5:
            ok[i] = False
            continue

        survival.remove_lower_NaNs(S)
        survival.remove_upper_NaNs(S)
        survival.remove_lower_NaNs(err)
        survival.remove_upper_NaNs(err, val=0.0)

        err_low[i], med[i], err_high[i] = survival.survival_quantile(
            m_eval, S, err, 0.5) 
        _, low[i], _ = survival.survival_quantile(
            m_eval, S, err, 0.5 - 0.68/2) 
        _, high[i], _ = survival.survival_quantile(
            m_eval, S, err, 0.5 + 0.68/2) 


    return med, err_low, err_high, low, high, ok

def f(x, y, alpha):
    return -alpha*x, -y/alpha

def main():
    palette.configure(True)
    
    x = np.linspace(0, 10, 200)
    y = np.linspace(0, -3, 200)
    n = 1000
    
    alpha = random.uniform(0.25, 1, n)
    cutoff_y = random.uniform(-2, -1, n)
    cutoff_y = np.minimum(random.randn(n) - 1.5, -0.5)

    lower_lim = -4

    all_y = np.zeros((n, len(x)))
    cut_y = np.ones((n, len(x)))*lower_lim
    cut_x = np.ones((n, len(y)))*-1

    x_cuts = np.zeros(n)
    
    for i in range(n):
        y_i, x_i = f(x, y, alpha[i])
        ok = y_i > cutoff_y[i]
        if i < 30:
            plt.plot(x[ok], 10**y_i[ok], lw=1, c="k", alpha=0.5)

        all_y[i,:] = y_i
        cut_y[i,ok] = y_i[ok]
        x_cuts[i] = np.max(x[ok])

        ok = y > cutoff_y[i]
        cut_x[i,ok] = x_i[ok]
        cut_x[i,~ok] = np.max(x_i[ok])

    med_y_all = np.median(all_y, axis=0) 
    med_y_cut = np.median(cut_y, axis=0)
    med_x_cut = np.median(cut_x, axis=0)
    
    med_y_skip = np.zeros(len(x))
    for i in range(len(x)):
        ok = cut_y[:,i] > lower_lim
        if np.sum(ok) > 0:
            med_y_skip[i] = np.median(cut_y[ok,i])
        else:
            med_y_skip[i] = lower_lim
    
    y_med_cut = -1.5
    x_med_cut = np.median(x_cuts)

    #plt.plot([0, 4], 2*[10**y_med_cut], "--", c="k")
    #plt.plot(2*[x_med_cut], [1, 0.001], "--", c="k")
    
    plt.plot([], [], pc("r"), label=r"${\rm True\ median}$")

    plt.plot(x, 10**med_y_skip, c=pc("o"), label=r"${\rm Stacking\ (Method\ A)}$")
    plt.plot(x, 10**med_y_cut, c=pc("b"), label=r"${\rm Stacking\ (Method\ B)}$")
    plt.plot(x, 10**med_y_all, c=pc("r"))

    t_eval = np.linspace(0, 6, 40)[1:]
    ts = [None for _ in range(len(all_y))]
    ms = [None for _ in range(len(all_y))]
    for i in range(len(all_y)):
        ok = x < x_cuts[i]
        ts[i], ms[i] = x[ok], 10**all_y[i,ok]

    med, err_low, err_high, low, high, ok = kaplan_meier_mass_loss(ts, ms, t_eval)

    plt.title(r"${\rm Toy\ Model}$")
    plt.plot(t_eval, med, c=pc("p"), label=r"${\rm Kaplan{-}Meier}$")
    plt.fill_between(t_eval, err_low, err_high, color=pc("p"), alpha=0.2)
    
    plt.xlabel(r"$\tilde{t}$")
    plt.ylabel(r"$\tilde{m}(\tilde{t})$")
    
    plt.legend(loc="upper right", fontsize=17)

    plt.yscale("log")
    plt.ylim(0.001, 1)
    plt.xlim(0, 7)

    print("../plots/core_tracking/mock_disruption.pdf")
    plt.savefig("../plots/core_tracking/mock_disruption.pdf")
    
if __name__ == "__main__":
    main()
