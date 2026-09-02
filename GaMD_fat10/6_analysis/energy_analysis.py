# plot_md_analysis.py
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

# One distinct color per replica, plus red for reference lines
REP_COLORS = {
    1: "#0072B2",  # blue
    2: "#E69F00",  # orange
    3: "#009E73",  # green
}
REF_COLOR = "#D55E00"  # red

REPLICAS = [1, 2, 3]
DATA_DIR = "energy"  # expects energy/rep1, energy/rep2, energy/rep3
RESULTS_DIR = "results"
os.makedirs(RESULTS_DIR, exist_ok=True)


# ── Helper ─────────────────────────────────────────────────────────────────
def load(fname):
    times, vals = [], []
    with open(fname) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("@"):
                continue
            parts = line.split()
            if len(parts) >= 2:
                times.append(float(parts[0]))
                vals.append(float(parts[1]))
        return np.array(times) / 1000, np.array(vals)  # ps → ns, values


def rep_path(rep, fname):
    return f"{DATA_DIR}/rep{rep}/{fname}"


# ── 1. Energy (ETOT, EPTOT, EKTOT), 3 replicas each ─────────────────────────
def plot_energy():
    components = {
        "E$_{tot}$": "summary.ETOT",
        "E$_{pot}$": "summary.EPTOT",
        "E$_{kin}$": "summary.EKTOT",
    }

    fig, axes = plt.subplots(3, 1, figsize=(6.5, 8.0), sharex=True)

    for ax, (label, fname) in zip(axes, components.items()):
        for rep in REPLICAS:
            t, e = load(rep_path(rep, fname))
            ax.plot(t, e,
                    color=REP_COLORS[rep],
                    linewidth=0.7,
                    alpha=0.85,
                    label=f"Replica {rep}")
        ax.set_ylabel(f"{label}\n(kcal mol$^{{-1}}$)")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.minorticks_off()

    axes[0].set_title("MD Simulation Energy Components - 3 Replicas", pad=8)
    axes[0].legend(loc="upper right")
    axes[-1].set_xlabel("Time (ns)")

    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, "energy_plot.png"))
    plt.close()
    print("Saved energy_plot.png")


# ── 2. Temperature, 3 replicas overlaid ─────────────────────────────────────
def plot_temperature(target_temp=300):
    fig, ax = plt.subplots(figsize=(6.5, 3.5))

    for rep in REPLICAS:
        t, temp = load(rep_path(rep, "summary.TEMP"))
        ax.plot(t, temp,
                color=REP_COLORS[rep],
                linewidth=0.7,
                alpha=0.85,
                label=f"Replica {rep}")

    ax.axhline(target_temp, color=REF_COLOR, linewidth=0.8, linestyle="--",
               alpha=0.7, label=f"Target ({target_temp} K)")

    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Temperature (K)")
    ax.set_title("MD Simulation Temperature - 3 Replicas", pad=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.minorticks_off()
    ax.legend(loc="upper right")
    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, "temperature_plot.png"))
    plt.close()
    print("Saved temperature_plot.png")


# ── 3. Density, 3 replicas overlaid ─────────────────────────────────────────
def plot_density():
    fig, ax = plt.subplots(figsize=(6.5, 3.5))

    for rep in REPLICAS:
        t, dens = load(rep_path(rep, "summary.DENSITY"))
        ax.plot(t, dens,
                color=REP_COLORS[rep],
                linewidth=0.7,
                alpha=0.85,
                label=f"Replica {rep}")

    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Density (g cm$^{-3}$)")
    ax.set_title("MD Simulation Density - 3 Replicas", pad=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.minorticks_off()
    ax.legend(loc="upper right")
    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, "density_plot.png"))
    plt.close()
    print("Saved density_plot.png")


# ── 4. Volume, 3 replicas overlaid ──────────────────────────────────────────
def plot_volume():
    fig, ax = plt.subplots(figsize=(6.5, 3.5))

    for rep in REPLICAS:
        t, vol = load(rep_path(rep, "summary.VOLUME"))
        ax.plot(t, vol,
                color=REP_COLORS[rep],
                linewidth=0.7,
                alpha=0.85,
                label=f"Replica {rep}")

    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Volume (Å$^{3}$)")
    ax.set_title("MD Simulation Volume - 3 Replicas", pad=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.minorticks_off()
    ax.legend(loc="upper right")
    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, "volume_plot.png"))
    plt.close()
    print("Saved volume_plot.png")


# ── 5. Pressure, 3 replicas (running average only - raw pressure per replica
#      would put 6+ noisy lines on one axis and be unreadable) ──────────────
def plot_pressure():
    fig, ax = plt.subplots(figsize=(6.5, 3.5))

    for rep in REPLICAS:
        t, pres = load(rep_path(rep, "summary.PRES"))
        window = max(1, len(pres) // 50)
        pres_smooth = np.convolve(pres, np.ones(window) / window, mode="same")
        ax.plot(t, pres_smooth,
                color=REP_COLORS[rep],
                linewidth=1.1,
                alpha=0.9,
                label=f"Replica {rep} (running avg, n={window})")

    ax.axhline(1.0, color=REF_COLOR, linewidth=0.8, linestyle="--",
               alpha=0.7, label="1 atm")

    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Pressure (bar)")
    ax.set_title("MD Simulation Pressure - 3 Replicas", pad=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.minorticks_off()
    ax.legend(loc="upper right")
    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, "pressure_plot.png"))
    plt.close()
    print("Saved pressure_plot.png")


# ── Run all ────────────────────────────────
if __name__ == "__main__":
    plot_energy()
    plot_temperature(target_temp=310)  # ← change to match your temp0
    plot_density()
    plot_volume()
    plot_pressure()