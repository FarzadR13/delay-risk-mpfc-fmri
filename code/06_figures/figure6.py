"""Figure 6: GLM2 delayed negative w(p) map (A) and probability-bin profiles in the three clusters (B), N = 22."""
import sys
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import brain_style as bs  # noqa: E402
from common import C, save_source  # noqa: E402
from plot_utils import BLUE, GREEN  # noqa: E402

PEAK = (-9, 59, 26)
ROIS = [("cluster_mpfc", "mPFC"), ("cluster_pcc", "PCC"), ("cluster_left_ag", "Left AG")]


def main():
    fig = plt.figure(figsize=(15 / 2.54, 12 / 2.54), facecolor="black")
    ext, core = bs.extent_img("glm2_delayed_wp_negative"), bs.core_img("glm2_delayed_wp_negative", "glm2_delayed_wp_negative")
    for i, (mode, c) in enumerate(zip("yxz", (PEAK[1], PEAK[0], PEAK[2]))):
        ax = fig.add_axes([0.02 + 0.32 * i, 0.62, 0.285, 0.33], facecolor="black")
        bs.slice_panel(ax, fig, ext, core, mode, c, bs.COOL_EXT, bs.COOL_CORE)
    bs.split_colorbar(fig.add_axes([0.958, 0.62, 0.016, 0.33], facecolor="black"), bs.COOL_EXT, "#4466FF", "#00CCFF")
    fig.text(0.5, 0.585, "Delayed reward", ha="center", color="white", fontsize=11, fontweight="bold")
    fig.legend(handles=[mpatches.Patch(color="#4400CC", label="p<0.005 uncorrected"),
                        mpatches.Patch(color="#00CCFF", label="p<0.001 FWE cluster-corrected")],
               loc="center", ncol=2, fontsize=8, facecolor="#111111", edgecolor="#444444", labelcolor="white",
               bbox_to_anchor=(0.5, 0.535))

    prof = pd.read_csv(C.RESULTS / "tables" / "wp_bin_profiles.csv")
    w = pd.read_csv(C.BEHAVIOR / "prescan_ce.csv")
    w = w[w.participant_id.isin(pd.read_csv(C.ROI / "wp_bins.csv").participant_id.unique())]
    w["probability_bin"] = w.probability.map(C.PROBABILITY_BINS)
    wsem = w.groupby(["participant_id", "condition", "probability_bin"]).w_empirical.mean().groupby(["condition", "probability_bin"]).sem()
    for i, (roi, title) in enumerate(ROIS):
        ax = fig.add_axes([0.07 + 0.31 * i, 0.10, 0.25, 0.32], facecolor="black")
        for cond, col in (("immediate", BLUE), ("delayed", GREEN)):
            p = prof[(prof.roi == roi) & (prof.condition == cond)].sort_values("probability_bin")
            ax.errorbar(p.mean_w_empirical, p["mean"], yerr=p["sem"], xerr=wsem[cond].to_numpy(), color=col, marker="o",
                        ms=4, capsize=2, lw=1.4, label=cond.capitalize())
        ax.set_title(title, color="white", fontsize=9)
        ax.set_xlabel("Subjective probability weight w(p)", color="white", fontsize=7)
        if i == 0:
            ax.set_ylabel("Parameter estimate (a.u.)", color="white", fontsize=7)
        ax.tick_params(colors="white", labelsize=6)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("bottom", "left"):
            ax.spines[s].set_color("white")
    ax.legend(loc="upper center", bbox_to_anchor=(-0.75, -0.25), ncol=2, frameon=False, labelcolor="white", fontsize=8)
    fig.text(0.0, 0.95, "A", color="white", fontsize=13, fontweight="bold")
    fig.text(0.0, 0.46, "B", color="white", fontsize=13, fontweight="bold")
    bs.save(fig, "Figure6")
    save_source("Fig6B", prof[prof.roi.isin([r for r, _ in ROIS])])


if __name__ == "__main__":
    main()
