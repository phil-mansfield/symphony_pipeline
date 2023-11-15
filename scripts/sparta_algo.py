import numpy as np
import matplotlib.pyplot as plt
import symlib

try:
    import palette
    palette.configure(False)
except:
    pass

class ParticleBuffer(object):
    def __init__(self, sim_dir, start_snap, mode="all", halo=0):
        """ Initialize buffer. start_snap is the first loaded snapshot, not the
        first snapshot that you can find pericentres at.
        
        You will need to manually load the first four snapshots into the buffer
        yourself.
        """
        # Information about the central
        rs, hist = symlib.read_rockstar(sim_dir)
        self.rvir = rs["rvir"][0,:]
        self.halo, self.mode = halo, mode

        # Set up internal rotating particle buffers
        self.part = symlib.Particles(sim_dir)
        n_p = self.part.count(mode=mode, halo=halo)
        self.n_p = n_p

        self.r, self.vr = np.zeros((4, n_p)), np.zeros((4, n_p))
        self.ok = np.zeros((4, n_p), dtype=bool)

        self.infall_snap = np.ones(n_p, dtype=int)*-1
        self.start_snap = start_snap

        # Set up two snapshot counters and initialize rotating buffer
        self.leading_snap = start_snap-1
        self.snap = self.leading_snap-1
        
    def load_next_snap(self):
        """ load_next_snap increments the snapshot by 1 and updates buffer
        internal state accordingly
        """
        self.snap += 1
        self.leading_snap += 1
        print(self.leading_snap)
        
        p = self.part.read(self.leading_snap, mode=self.mode, halo=self.halo)
        # Get rid of values that are going to cause math function warnings
        p["x"][~p["ok"]], p["v"][~p["ok"]] = 1, 1

        r = np.sqrt(np.sum(p["x"]**2, axis=1))
        vr = np.sum(p["x"]*p["v"], axis=1)/r
        
        # Rotate buffers
        self.r[:-1,:] = self.r[1:,:]
        self.vr[:-1,:], self.ok[:-1,:] = self.vr[1:,:], self.ok[1:,:]
        self.r[-1,:], self.vr[-1,:], self.ok[-1,:] = r, vr, p["ok"]
        
        # Flag the snapshot of first infall for particles.
        pre_infall = self.infall_snap == -1
        in_rvir = (r < self.rvir[self.leading_snap]) & p["ok"]
        self.infall_snap[pre_infall & in_rvir] = self.leading_snap

class Trajectories(object):
    def __init__(self, buf, targets):
        """ Initialize trajectories for a set of indices, targets. buf
        is an initialized ParticleBuffer
        """
        self.buf = buf
        self.targets = targets

        n_snap = len(buf.rvir)
        self.r = np.zeros((len(targets), n_snap))
        self.vr = np.zeros((len(targets), n_snap))
        self.ok = np.zeros((len(targets), n_snap), dtype=bool) 

    def update(self):
        """ update trajectories unsing the current leading_snap of a
        ParticleBuffer
        """
        idx = np.arange(len(self.targets), dtype=int)
        self.r[idx,self.buf.leading_snap] = self.buf.r[-1,self.targets]
        self.vr[idx,self.buf.leading_snap] = self.buf.vr[-1,self.targets]
        self.ok[idx,self.buf.leading_snap] = self.buf.ok[-1,self.targets]

class Pericentres(object):
    def __init__(self, buf, n_peri):
        self.buf = buf
        self.n_peri = n_peri

        self.idx = np.arange(buf.n_p, dtype=int)
        self.peri_count = np.zeros(buf.n_p, dtype=int)
        self.peri_snap = np.ones((buf.n_p, n_peri), dtype=int)*-1

    def update(self):
        """ update updates the current pericentre snapshots using the current
        state of the ParticleBuffer
        """
        is_peri = self.is_peri()
        self.peri_count[is_peri] += 1
        do_update = is_peri & (self.peri_count <= self.n_peri)
        i0, i1 = self.idx[do_update], self.peri_count[do_update]-1
        self.peri_snap[i0, i1] = self.buf.snap

    def is_peri(self):
        """ is_peri returns true if a particle has a pericentre between
        buf.snap and buf.snap -1 and false otherwise.
        """
        always_ok = self.buf.ok[0]
        neg_before = (self.buf.vr[0] < 0) & (self.buf.vr[1] < 0)
        pos_after = (self.buf.vr[2] > 0) & (self.buf.vr[3] > 0)
        past_infall = self.buf.infall_snap <= self.buf.snap
        infall_seen = self.buf.start_snap < self.buf.infall_snap
        return always_ok & neg_before & pos_after & past_infall & infall_seen

def initialize(buf, traj, peri):
    for _ in range(4):
        buf.load_next_snap()
        traj.update()
    peri.update()

def plot_trajectories(sim_dir, a, targets, traj, peri):
    rs, hist = symlib.read_rockstar(sim_dir)
    host_rvir, host_ok = rs["rvir"][0,:], rs["ok"][0,:]

    fig, ax = plt.subplots(2, sharex=True, figsize=(8, 16))
    fig.subplots_adjust(hspace=0.05)
    ax[1].set_xlabel(r"$a(t)$")
    ax[1].set_ylabel(r"$r\ ({\rm physical\ kpc})$")
    ax[0].set_ylabel(r"$v_r\ ({\rm km\,s^{-1}})$")

    colors = ["tab:red", "tab:orange", "tab:green", "tab:blue", "tab:purple"]
    for i, j in enumerate(targets):
        ok, r, vr = traj.ok[i], traj.r[i], traj.vr[i]
        ax[1].plot(a[ok], r[ok], c=colors[i])
        ax[1].plot(a[ok], r[ok], "o", c=colors[i])
        ax[0].plot(a[ok], vr[ok], c=colors[i])
        ax[0].plot(a[ok], vr[ok], "o", c=colors[i])
        
        snap1, snap2 = peri.peri_snap[j,0], peri.peri_snap[j,1]
        snap3 = peri.peri_snap[j,2]
        ms = 15
        if snap1 != -1:
            ax[1].plot(a[snap1], r[snap1], "o", c="k", ms=ms)
            ax[0].plot(a[snap1], vr[snap1], "o", c="k", ms=ms)
        if snap2 != -1:
            ax[1].plot(a[snap2], r[snap2], "*", c="k", ms=ms)
            ax[0].plot(a[snap2], vr[snap2], "*", c="k", ms=ms)
        if snap3 != -1:
            ax[1].plot(a[snap3], r[snap3], "^", c="k", ms=ms)
            ax[0].plot(a[snap3], vr[snap3], "^", c="k", ms=ms)

    xlo, xhi = ax[0].get_xlim()
    ax[0].set_xlim(xlo, xhi)
    ax[0].plot([xlo, xhi], [0, 0], "--", c="k", lw=2)
    ax[1].plot(a[host_ok], host_rvir[host_ok], "--", c="k", lw=2)
    fig.savefig("../plots/sparta_test.png")

def main():
    halo, mode = 0, "all"

    base_dir = "/sdf/home/p/phil1/ZoomIns"
    sim_dir = symlib.get_host_directory(base_dir, "SymphonyCluster", 0)
    a = symlib.scale_factors(sim_dir)
    n_snap = len(a)

    # Select target particles for trajectories.
    p = symlib.Particles(sim_dir).read(n_snap - 1, mode=mode, halo=halo)
    targets = np.where(p["snap"] == 125)[0][:5]
    print(targets)

    # Initialize
    buf = ParticleBuffer(sim_dir, 124, halo=halo, mode=mode)
    peri = Pericentres(buf, 3)
    traj = Trajectories(buf, targets)
    initialize(buf, traj, peri)

    # Move through all snapshots
    while buf.leading_snap < n_snap-1:
        buf.load_next_snap()
        traj.update()
        peri.update()

    # Plot
    plot_trajectories(sim_dir, a, targets, traj, peri)

if __name__ == "__main__": main()
