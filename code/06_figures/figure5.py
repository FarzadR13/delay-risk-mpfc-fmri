"""Figure 5: meta-analytic mPFC mask. A, B, heatmap-model cell means; C, D, GLM3 main effect vs gamma (N = 22)."""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import C, correlation, lattimore_params, save_source  # noqa: E402
from plot_utils import clean_axes, ols_band, save_figure  # noqa: E402

PROB = ["[.05, .10]", "[.25, .50, .75]", "[.90, .95]"]
DIST = ["[-30, -20]", "[-10, CE, +10]", "[+20, +30]"]


def main():
    cells = pd.read_csv(C.ROI / "heatmap_cells.csv")
    means = cells.groupby(["condition", "probability_bin", "distance_bin"]).value.mean()
    lo, hi = means.min(), means.max()
    fig = plt.figure(figsize=(10, 10.5), facecolor="white")
    for i, cond in enumerate(C.CONDITIONS):
        ax = fig.add_axes([0.12 + 0.47 * i, 0.60, 0.33, 0.33])
        m = np.array([[means[(cond, p, d)] for d in (1, 2, 3)] for p in (1, 2, 3)])
        im = ax.imshow(m, cmap="viridis", vmin=lo, vmax=hi)
        for p in range(3):
            for d in range(3):
                ax.text(d, p, f"{m[p, d]:.3f}", ha="center", va="center", color="black" if m[p, d] > lo + 0.6 * (hi - lo) else "white")
        ax.set_xticks(range(3), DIST, fontsize=8); ax.set_yticks(range(3), PROB, fontsize=8)
        ax.set(xlabel="Distance to certainty equivalent (%)", ylabel="Objective probability")
        ax.set_title(f"{cond.capitalize()} choices", fontweight="bold", fontsize=13)
        save_source(f"Fig5{'AB'[i]}", cells[cells.condition == cond])
    cb = fig.colorbar(im, cax=fig.add_axes([0.94, 0.60, 0.015, 0.33]))
    cb.set_label("Contrast estimate (a.u.)", fontsize=8)

    b = pd.read_csv(C.ROI / "bartra_mpfc.csv").set_index("participant_id").join(lattimore_params())
    for i, cond in enumerate(C.CONDITIONS):
        ax = fig.add_axes([0.12 + 0.47 * i, 0.07, 0.36, 0.38])
        x, y = b[f"glm3_{cond}_main_effect"], b[f"gamma_{cond}"]
        ols_band(ax, x, y, color="#222222", band="#D9DCE0")
        ax.scatter(x, y, s=70, color="#5DA5DA", edgecolor="black", lw=0.8, zorder=3)
        c = correlation(x, y)
        ax.text(0.03, 0.97, f"Pearson r = {c['r']:.2f}, p = {c['p']:.3f}\nSpearman rho = {c['rho']:.2f}, p = {c['rho_p']:.3f}",
                transform=ax.transAxes, va="top", fontsize=9)
        ax.set(xlabel="mPFC response during choice (a.u.)", ylabel="Curvature parameter gamma")
        ax.set_title(f"{cond.capitalize()} choices", fontweight="bold", fontsize=13); clean_axes(ax)
        save_source(f"Fig5{'CD'[i]}", pd.DataFrame(dict(participant_id=b.index, mpfc_main_effect=x.to_numpy(), gamma=y.to_numpy())))
    for letter, (x, y) in zip("ABCD", ((0.03, 0.96), (0.50, 0.96), (0.03, 0.48), (0.50, 0.48))):
        fig.text(x, y, letter, fontsize=18, fontweight="bold")
    save_figure(fig, "Figure5")


if __name__ == "__main__":
    main()
