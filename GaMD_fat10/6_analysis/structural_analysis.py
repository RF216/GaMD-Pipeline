import os
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    "font.size": 10,
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "axes.linewidth": 0.8,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "xtick.major.size": 4,
    "ytick.major.size": 4,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": False,
    "ytick.right": False,
    "legend.frameon": True,
    "legend.framealpha": 0.9,
    "legend.edgecolor": "0.8",
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.05,
    "axes.facecolor": "white",
    "figure.facecolor": "white",
    "axes.grid": False,
})

# One distinct color per replica, plus a red for mean/reference lines
REP_COLORS = {
    1: "#0072B2",  # blue
    2: "#E69F00",  # orange
    3: "#009E73",  # green
}
MEAN_COLOR = "#D55E00"   # red/orange
BAND_COLOR = "#CC79A7"   # purple

REPLICAS = [1, 2, 3]
TOTAL_TIME_NS = 1000.0  # each replica is an independent 1 us production run
DATA_DIR = "structural"
RESULTS_DIR = "results"
os.makedirs(RESULTS_DIR, exist_ok=True)


# ── Helpers ─────────────────────────────────────────────
def load(fname):
    x, y = [], []
    with open(fname) as f:
        for line in f:
            if not line.strip() or line.startswith("#") or line.startswith("@"):
                continue
            parts = line.split()
            if len(parts) >= 2:
                x.append(float(parts[0]))
                y.append(float(parts[1]))
    return np.array(x), np.array(y)


def to_ns(frames):
    if len(frames) < 2:
        return frames
    return (frames - frames[0]) / (frames[-1] - frames[0]) * TOTAL_TIME_NS


# ── 1. RMSD (protein only, all 3 replicas overlaid) ────
def plot_rmsd_overlaid():
    fig, ax = plt.subplots(figsize=(6.5, 3.5))

    for rep in REPLICAS:
        frames, rmsd = load(f"{DATA_DIR}/rmsd_rep{rep}.dat")
        time = to_ns(frames)
        ax.plot(
            time, rmsd,
            color=REP_COLORS[rep],
            linewidth=0.7,
            alpha=0.85,
            label=f"Replica {rep} (mean {np.mean(rmsd):.2f} Å)"
        )

    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("RMSD (Å)")
    ax.set_title("Protein RMSD - 3 Replicas", pad=8)
    ax.set_xlim(0, TOTAL_TIME_NS)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend()

    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, "rmsd_plot_overlaid.png"))
    plt.close()


def plot_rmsd_averaged():
    fig, ax = plt.subplots(figsize=(6.5, 3.5))

    time = None
    rmsd_all = []

    for rep in REPLICAS:
        frames, rmsd = load(f"{DATA_DIR}/rmsd_rep{rep}.dat")
        if time is None:
            time = to_ns(frames)
        rmsd_all.append(rmsd)

    rmsd_all = np.array(rmsd_all)
    mean_rmsd = rmsd_all.mean(axis=0)
    std_rmsd = rmsd_all.std(axis=0, ddof=1)

    ax.fill_between(
        time,
        mean_rmsd - std_rmsd,
        mean_rmsd + std_rmsd,
        color=BAND_COLOR,
        alpha=0.2,
        label="±1 SD across replicas"
    )

    ax.plot(
        time,
        mean_rmsd,
        color=MEAN_COLOR,
        linewidth=1.1,
        label=f"Mean RMSD (3 replicas, {np.mean(mean_rmsd):.2f} Å)"
    )

    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("RMSD (Å)")
    ax.set_title("Protein RMSD - Mean of 3 Replicas", pad=8)
    ax.set_xlim(0, TOTAL_TIME_NS)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend()

    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, "rmsd_plot_averaged.png"))
    plt.close()


# ── 2. RMSF (mean ± stdev across 3 replicas, per residue) ─
def plot_rmsf():
    fig, ax = plt.subplots(figsize=(6.5, 3.5))

    residues = None
    rmsf_all = []

    for rep in REPLICAS:
        res, rmsf = load(f"{DATA_DIR}/rmsf_rep{rep}.dat")
        if residues is None:
            residues = res
        rmsf_all.append(rmsf)

    rmsf_all = np.array(rmsf_all)
    mean_rmsf = rmsf_all.mean(axis=0)
    std_rmsf = rmsf_all.std(axis=0, ddof=1)

    # Individual replicas, light lines for context
    for rep, rmsf in zip(REPLICAS, rmsf_all):
        ax.plot(
            residues,
            rmsf,
            color=REP_COLORS[rep],
            linewidth=0.5,
            alpha=0.35
        )

    # Mean ± stdev band
    ax.fill_between(
        residues,
        mean_rmsf - std_rmsf,
        mean_rmsf + std_rmsf,
        color=BAND_COLOR,
        alpha=0.2,
        label="±1 SD across replicas"
    )

    ax.plot(
        residues,
        mean_rmsf,
        color=MEAN_COLOR,
        linewidth=1.1,
        label="Mean RMSF (3 replicas)"
    )

    ax.set_xlabel("Residue Number")
    ax.set_ylabel("RMSF (Å)")
    ax.set_title("Per-Residue RMSF - Mean of 3 Replicas", pad=8)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend()

    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, "rmsf_plot.png"))
    plt.close()


# ── 3. Radius of gyration - OVERLAID ────────────────────
def plot_rog_overlaid():
    fig, ax = plt.subplots(figsize=(6.5, 3.5))

    for rep in REPLICAS:
        frames, rog = load(f"{DATA_DIR}/rog_rep{rep}.dat")
        time = to_ns(frames)

        ax.plot(
            time,
            rog,
            color=REP_COLORS[rep],
            linewidth=0.7,
            alpha=0.85,
            label=f"Replica {rep} (mean {np.mean(rog):.2f} Å)"
        )

    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Rg (Å)")
    ax.set_title("Radius of Gyration of GaMD FAT10 - 3 Replicas", pad=8)
    ax.set_xlim(0, TOTAL_TIME_NS)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend()

    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, "rog_plot_overlaid.png"))
    plt.close()


# ── 3b. Radius of gyration - AVERAGED ───────────────────
def plot_rog_averaged():
    fig, ax = plt.subplots(figsize=(6.5, 3.5))

    time = None
    rog_all = []

    for rep in REPLICAS:
        frames, rog = load(f"{DATA_DIR}/rog_rep{rep}.dat")

        if time is None:
            time = to_ns(frames)

        rog_all.append(rog)

    rog_all = np.array(rog_all)

    mean_rog = rog_all.mean(axis=0)
    std_rog = rog_all.std(axis=0, ddof=1)

    # Mean ± SD
    ax.fill_between(
        time,
        mean_rog - std_rog,
        mean_rog + std_rog,
        color=BAND_COLOR,
        alpha=0.2,
        label="±1 SD across replicas"
    )

    ax.plot(
        time,
        mean_rog,
        color=MEAN_COLOR,
        linewidth=1.1,
        label=f"Mean Rg (3 replicas, {np.mean(mean_rog):.2f} Å)"
    )

    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Rg (Å)")
    ax.set_title("Radius of Gyration of GaMD FAT10 - Mean of 3 Replicas", pad=8)
    ax.set_xlim(0, TOTAL_TIME_NS)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend()

    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, "rog_plot_averaged.png"))
    plt.close()


# ── 4. SASA - OVERLAID ─────────────────────────────────
def plot_sasa_overlaid():
    fig, ax = plt.subplots(figsize=(6.5, 3.5))

    for rep in REPLICAS:
        frames, sasa = load(f"{DATA_DIR}/surf_rep{rep}.dat")
        time = to_ns(frames)

        ax.plot(
            time,
            sasa,
            color=REP_COLORS[rep],
            linewidth=0.7,
            alpha=0.85,
            label=f"Replica {rep} (mean {np.mean(sasa):.0f} Å²)"
        )

    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("SASA (Å²)")
    ax.set_title("Solvent Accessible Surface Area - 3 Replicas", pad=8)
    ax.set_xlim(0, TOTAL_TIME_NS)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend()

    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, "sasa_plot_overlaid.png"))
    plt.close()


# ── 4b. SASA - AVERAGED ────────────────────────────────
def plot_sasa_averaged():
    fig, ax = plt.subplots(figsize=(6.5, 3.5))

    time = None
    sasa_all = []

    for rep in REPLICAS:
        frames, sasa = load(f"{DATA_DIR}/surf_rep{rep}.dat")

        if time is None:
            time = to_ns(frames)

        sasa_all.append(sasa)

    sasa_all = np.array(sasa_all)

    mean_sasa = sasa_all.mean(axis=0)
    std_sasa = sasa_all.std(axis=0, ddof=1)

    # Mean ± SD
    ax.fill_between(
        time,
        mean_sasa - std_sasa,
        mean_sasa + std_sasa,
        color=BAND_COLOR,
        alpha=0.2,
        label="±1 SD across replicas"
    )

    ax.plot(
        time,
        mean_sasa,
        color=MEAN_COLOR,
        linewidth=1.1,
        label=f"Mean SASA (3 replicas, {np.mean(mean_sasa):.0f} Å²)"
    )

    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("SASA (Å²)")
    ax.set_title("Solvent Accessible Surface Area - Mean of 3 Replicas", pad=8)
    ax.set_xlim(0, TOTAL_TIME_NS)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend()

    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, "sasa_plot_averaged.png"))
    plt.close()


# ── RUN ────────────────────────────────────────────────
if __name__ == "__main__":

    # RMSD
    plot_rmsd_overlaid()
    plot_rmsd_averaged()

    # RMSF
    plot_rmsf()

    # Radius of gyration
    plot_rog_overlaid()
    plot_rog_averaged()

    # SASA
    plot_sasa_overlaid()
    plot_sasa_averaged()

    print("All 3-replica plots saved!")
