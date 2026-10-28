"""Figure 4: GLM1 signed distance to the CE. A, immediate choices; B, immediate > delayed (N = 22)."""
import sys
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
import brain_style as bs  # noqa: E402

ROWS = [("A", "glm1_immediate_distance", "glm1_immediate_distance_positive", (6, 50, 14), "Immediate reward",
         "Parametric effect of signed distance to the certainty equivalent (CE)"),
        ("B", "glm1_immediate_minus_delayed_distance", "glm1_immediate_minus_delayed_distance_positive", (6, 50, 17),
         "Immediate reward > Delayed reward", "Difference in the parametric effect of signed distance to CE")]


def main():
    fig = plt.figure(figsize=(15 / 2.54, 14 / 2.54), facecolor="black")
    for r, (letter, name, table, peak, title, subtitle) in enumerate(ROWS):
        y0 = 0.53 - 0.5 * r
        ext, core = bs.extent_img(name), bs.core_img(name, table)
        for i, (mode, c) in enumerate(zip("yxz", (peak[1], peak[0], peak[2]))):
            ax = fig.add_axes([0.02 + 0.32 * i, y0 + 0.17, 0.285, 0.29], facecolor="black")
            bs.slice_panel(ax, fig, ext, core, mode, c, bs.WARM_EXT, bs.WARM_CORE)
        bs.split_colorbar(fig.add_axes([0.958, y0 + 0.17, 0.016, 0.29], facecolor="black"), bs.WARM_EXT, "#FF6A00", "#FFE600")
        fig.text(0.0, y0 + 0.46, letter, color="white", fontsize=13, fontweight="bold")
        fig.text(0.5, y0 + 0.14, title, ha="center", color="white", fontsize=11, fontweight="bold")
        fig.text(0.5, y0 + 0.105, subtitle, ha="center", color="white", fontsize=8)
        handles = [mpatches.Patch(color="#DD1100", label="p<0.005 uncorrected")]
        if core is not None:
            handles.append(mpatches.Patch(color="#FFE600", label="p<0.001 FWE cluster-corrected"))
        fig.legend(handles=handles, loc="center", ncol=2, fontsize=8, facecolor="#111111", edgecolor="#444444",
                   labelcolor="white", bbox_to_anchor=(0.5, y0 + 0.055))
    bs.save(fig, "Figure4")


if __name__ == "__main__":
    main()
