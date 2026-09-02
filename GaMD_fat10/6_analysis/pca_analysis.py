# plot_pca.py
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter
import os

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

COLORS = {
    "blue": "#0072B2",
    "orange": "#E69F00",
    "green": "#009E73",
    "purple": "#CC79A7",
    "red": "#D55E00",
}

# Per-replica colors, used for the "colored by replica" overlay plots
REPLICA_COLORS = {
    1: COLORS["blue"],
    2: COLORS["orange"],
    3: COLORS["green"],
}

# ══════════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ══════════════════════════════════════════════════════════════════════════
TOTAL_TIME_NS = 1000.0       # length of EACH replica, in ns
DATA_DIR = "pca"             # folder containing the input .dat files
OUTPUT_DIR = "results"       # folder where output PNGs are saved
REPLICAS = [1, 2, 3]         # replica indices, matches *_rep{N}.dat naming
SMOOTH_SIGMA = 2.5            # Gaussian smoothing (increase → smoother basins)
ENERGY_CAP = 5.0               # kcal/mol - cap ΔG to focus color range on basins
NBINS = 100                    # histogram bins per axis
TEMPERATURE = 300.0            # K - match your simulation temperature


def projection_file(rep):
    return os.path.join(DATA_DIR, f"pca_projections_rep{rep}.dat")


def eigenvector_file(rep):
    return os.path.join(DATA_DIR, f"pca_eigenvectors_rep{rep}.dat")


def output_path(fname):
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    return os.path.join(OUTPUT_DIR, fname)


# ── Helpers ────────────────────────────────────────────────────────────────
def load_projections(fname):
    """Load cpptraj projection file → (frames, PC1, PC2, PC3)."""
    frames, pc1, pc2, pc3 = [], [], [], []
    with open(fname) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("@"):
                continue
            parts = line.split()
            if len(parts) >= 4:
                frames.append(float(parts[0]))
                pc1.append(float(parts[1]))
                pc2.append(float(parts[2]))
                pc3.append(float(parts[3]))
    return np.array(frames), np.array(pc1), np.array(pc2), np.array(pc3)


def to_ns(frames, total_time_ns=TOTAL_TIME_NS):
    """Scale frame indices → 0 to total_time_ns."""
    if len(frames) < 2:
        return frames
    return (frames - frames[0]) / (frames[-1] - frames[0]) * total_time_ns


def load_all_replicas():
    """Load projections for every replica in REPLICAS.

    Returns a dict: {rep: {"frames": ..., "pc1": ..., "pc2": ..., "pc3": ..., "time_ns": ...}}
    """
    data = {}
    for rep in REPLICAS:
        fname = projection_file(rep)
        frames, pc1, pc2, pc3 = load_projections(fname)
        data[rep] = {
            "frames": frames,
            "pc1": pc1,
            "pc2": pc2,
            "pc3": pc3,
            "time_ns": to_ns(frames),
        }
        print(f"Loaded {fname}: {len(frames)} frames")
    return data


def pooled(data, xkey, ykey):
    """Concatenate a given PC pair across all replicas (for combined free energy)."""
    x = np.concatenate([data[rep][xkey] for rep in REPLICAS])
    y = np.concatenate([data[rep][ykey] for rep in REPLICAS])
    return x, y


# ── 1. PC1 vs PC2 Scatter - one panel per replica, coloured by time ────────
def plot_pc1_vs_pc2():
    data = load_all_replicas()

    fig, axes = plt.subplots(1, len(REPLICAS), figsize=(5.0 * len(REPLICAS), 5.0),
                              sharex=True, sharey=True)
    if len(REPLICAS) == 1:
        axes = [axes]

    sc = None
    for ax, rep in zip(axes, REPLICAS):
        d = data[rep]
        sc = ax.scatter(d["pc1"], d["pc2"], c=d["time_ns"], cmap="viridis",
                         s=4, alpha=0.6, edgecolors="none", rasterized=True,
                         vmin=0, vmax=TOTAL_TIME_NS)
        ax.set_xlabel("PC1 (Å)")
        ax.set_title(f"Replica {rep}", pad=8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    axes[0].set_ylabel("PC2 (Å)")

    fig.suptitle("PCA - PC1 vs PC2 (coloured by time)", y=1.02)
    cbar = fig.colorbar(sc, ax=axes, pad=0.02, shrink=0.9)
    cbar.set_label("Time (ns)")

    plt.savefig(output_path("pca_pc1_pc2.png"))
    plt.close()
    print(f"Saved {output_path('pca_pc1_pc2.png')}")


# ── 2. PC1 vs PC3 Scatter - one panel per replica, coloured by time ────────
def plot_pc1_vs_pc3():
    data = load_all_replicas()

    fig, axes = plt.subplots(1, len(REPLICAS), figsize=(5.0 * len(REPLICAS), 5.0),
                              sharex=True, sharey=True)
    if len(REPLICAS) == 1:
        axes = [axes]

    sc = None
    for ax, rep in zip(axes, REPLICAS):
        d = data[rep]
        sc = ax.scatter(d["pc1"], d["pc3"], c=d["time_ns"], cmap="viridis",
                         s=4, alpha=0.6, edgecolors="none", rasterized=True,
                         vmin=0, vmax=TOTAL_TIME_NS)
        ax.set_xlabel("PC1 (Å)")
        ax.set_title(f"Replica {rep}", pad=8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    axes[0].set_ylabel("PC3 (Å)")

    fig.suptitle("PCA - PC1 vs PC3 (coloured by time)", y=1.02)
    cbar = fig.colorbar(sc, ax=axes, pad=0.02, shrink=0.9)
    cbar.set_label("Time (ns)")

    plt.savefig(output_path("pca_pc1_pc3.png"))
    plt.close()
    print(f"Saved {output_path('pca_pc1_pc3.png')}")


# ── 3. Replica overlay - PC1 vs PC2, coloured by replica ID ────────────────
def plot_replica_overlay():
    data = load_all_replicas()

    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    for rep in REPLICAS:
        d = data[rep]
        ax.scatter(d["pc1"], d["pc2"], s=4, alpha=0.5, edgecolors="none",
                   rasterized=True, color=REPLICA_COLORS.get(rep, "gray"),
                   label=f"Replica {rep}")

    ax.set_xlabel("PC1 (Å)")
    ax.set_ylabel("PC2 (Å)")
    ax.set_title("PCA - PC1 vs PC2 (coloured by replica)", pad=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    leg = ax.legend(markerscale=3, loc="best")

    plt.tight_layout()
    plt.savefig(output_path("pca_pc1_pc2_by_replica.png"))
    plt.close()
    print(f"Saved {output_path('pca_pc1_pc2_by_replica.png')}")


# ── 4. Free Energy Landscape (PC1 vs PC2), pooled across all replicas ──────
def plot_free_energy():
    data = load_all_replicas()
    pc1, pc2 = pooled(data, "pc1", "pc2")

    fig, ax = plt.subplots(figsize=(6.5, 5.5))

    # 2D histogram → Gaussian smooth → free energy
    H, xedges, yedges = np.histogram2d(pc1, pc2, bins=NBINS, density=True)
    H_smooth = gaussian_filter(H, sigma=SMOOTH_SIGMA)

    kB = 0.001987204  # kcal/(mol·K)
    kT = kB * TEMPERATURE

    # Compute ΔG; unvisited bins → NaN
    H_safe = np.where(H_smooth > 0, H_smooth, np.nan)
    G = -kT * np.log(H_safe / np.nanmax(H_safe))

    # Cap energy and fill unvisited bins at the cap
    G_capped = np.clip(G, 0, ENERGY_CAP)
    G_filled = np.where(np.isnan(G_capped), ENERGY_CAP, G_capped)

    # Heatmap - blue = low ΔG (basins), red = high ΔG / unsampled
    im = ax.pcolormesh(xedges, yedges, G_filled.T,
                       cmap="jet", vmin=0, vmax=ENERGY_CAP,
                       shading="flat", rasterized=True)

    # Contour lines to delineate basins
    xcenters = 0.5 * (xedges[:-1] + xedges[1:])
    ycenters = 0.5 * (yedges[:-1] + yedges[1:])
    X, Y = np.meshgrid(xcenters, ycenters)
    contour_levels = np.arange(0.5, ENERGY_CAP, 0.5)
    ax.contour(X, Y, G_filled.T, levels=contour_levels,
               colors="k", linewidths=0.4, alpha=0.5)

    cbar = plt.colorbar(im, ax=ax, pad=0.02)
    cbar.set_label("ΔG (kcal mol$^{-1}$)")

    ax.set_xlabel("PC1 (Å)")
    ax.set_ylabel("PC2 (Å)")
    n_reps = len(REPLICAS)
    ax.set_title(f"Free Energy Landscape - PC1 vs PC2 (pooled, {n_reps} replicas)", pad=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_path("pca_free_energy.png"))
    plt.close()
    print(f"Saved {output_path('pca_free_energy.png')}")


# ── 5. Per-replica free energy landscapes, side by side (same color scale) ─
def plot_free_energy_per_replica():
    data = load_all_replicas()
    kB = 0.001987204
    kT = kB * TEMPERATURE

    # Shared axis range across all replicas, so every panel covers the same
    # PC1/PC2 extent and uses identical bin edges (no white gaps, no
    # mismatched axes between panels).
    all_pc1, all_pc2 = pooled(data, "pc1", "pc2")
    xrange = (all_pc1.min(), all_pc1.max())
    yrange = (all_pc2.min(), all_pc2.max())

    fig, axes = plt.subplots(1, len(REPLICAS), figsize=(5.5 * len(REPLICAS), 5.0),
                              sharex=True, sharey=True)
    if len(REPLICAS) == 1:
        axes = [axes]

    im = None
    for ax, rep in zip(axes, REPLICAS):
        pc1, pc2 = data[rep]["pc1"], data[rep]["pc2"]
        H, xedges, yedges = np.histogram2d(pc1, pc2, bins=NBINS,
                                           range=[xrange, yrange], density=True)
        H_smooth = gaussian_filter(H, sigma=SMOOTH_SIGMA)
        H_safe = np.where(H_smooth > 0, H_smooth, np.nan)
        G = -kT * np.log(H_safe / np.nanmax(H_safe))
        G_capped = np.clip(G, 0, ENERGY_CAP)
        G_filled = np.where(np.isnan(G_capped), ENERGY_CAP, G_capped)

        im = ax.pcolormesh(xedges, yedges, G_filled.T,
                           cmap="jet", vmin=0, vmax=ENERGY_CAP,
                           shading="flat", rasterized=True)

        xcenters = 0.5 * (xedges[:-1] + xedges[1:])
        ycenters = 0.5 * (yedges[:-1] + yedges[1:])
        X, Y = np.meshgrid(xcenters, ycenters)
        contour_levels = np.arange(0.5, ENERGY_CAP, 0.5)
        ax.contour(X, Y, G_filled.T, levels=contour_levels,
                   colors="k", linewidths=0.4, alpha=0.5)

        ax.set_xlabel("PC1 (Å)")
        ax.set_title(f"Replica {rep}", pad=8)
        ax.set_xlim(xrange)
        ax.set_ylim(yrange)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    axes[0].set_ylabel("PC2 (Å)")
    fig.suptitle("Free Energy Landscape - PC1 vs PC2, per replica", y=1.02)
    cbar = fig.colorbar(im, ax=axes, pad=0.02, shrink=0.9)
    cbar.set_label("ΔG (kcal mol$^{-1}$)")

    plt.savefig(output_path("pca_free_energy_per_replica.png"))
    plt.close()
    print(f"Saved {output_path('pca_free_energy_per_replica.png')}")


# ── Run all ────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    plot_pc1_vs_pc2()
    plot_pc1_vs_pc3()
    plot_replica_overlay()
    plot_free_energy()
    plot_free_energy_per_replica()
    