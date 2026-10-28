"""Figure 7: low-p delta w vs the delayed w(p) slope. A, left AG cluster; B, left IFG pars triangularis (N = 22)."""
import sys
from pathlib import Path

import matplotlib.colorbar as mcolorbar
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import nibabel as nib
import pandas as pd
from nilearn import plotting
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import brain_style as bs  # noqa: E402
from common import C, correlation  # noqa: E402
from plot_utils import ols_band, p_text  # noqa: E402

T_COV = stats.t.isf(C.DISPLAY_P, 20)   # display threshold of the covariate map (df = 20)


def scatter(ax, x, y, title, lines):
    ols_band(ax, x, y, color="white", band="#3A3A3A")
    ax.scatter(x, y, s=50, color="#5DA5DA", edgecolor="white", lw=0.6, zorder=3)
    ax.text(0.97, 0.97, "\n".join(lines), transform=ax.transAxes, ha="right", va="top", color="white", fontsize=8)
    ax.set_title(title, color="white", fontweight="bold")
    ax.set_xlabel("Low-p delta w (Delayed - Immediate)", color="white", fontsize=8)
    ax.set_ylabel("Delayed w(p) slope (a.u.)", color="white", fontsize=8)
    ax.set_facecolor("black"); ax.tick_params(colors="white", labelsize=7)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("bottom", "left"):
        ax.spines[s].set_color("white")


def colorbar(fig, rect, vmin, vmax, ticks, label, core_line=False):
    ax = fig.add_axes(rect)
    cb = mcolorbar.ColorbarBase(ax, cmap=bs.COOL_EXT, norm=mcolors.Normalize(vmin, vmax), orientation="horizontal")
    cb.set_ticks(ticks); cb.ax.tick_params(colors="white", labelsize=7)
    cb.set_label(label, color="white", fontsize=7)
    if core_line:
        ax.axvline(bs.T_CORE, color="white", ls="--", lw=1)


def main():
    a = pd.read_csv(C.SOURCE_DIR / "Fig7A.csv"); b = pd.read_csv(C.SOURCE_DIR / "Fig7B.csv")
    fig = plt.figure(figsize=(9.2, 10), facecolor="black")
    ax = fig.add_axes([0.03, 0.45, 0.44, 0.50], facecolor="black")
    ext, core = bs.extent_img("glm2_delayed_wp_negative"), bs.core_img("glm2_delayed_wp_negative", "glm2_delayed_wp_negative")
    d = bs.slice_panel(ax, fig, ext, core, "z", 26, bs.COOL_EXT, bs.COOL_CORE)
    d.add_contours(nib.load(str(C.OWN_MASKS / "fig6a_left_ag_cluster.nii.gz")), levels=[0.5], colors="white", linewidths=1.5)
    colorbar(fig, [0.07, 0.43, 0.34, 0.012], bs.T_EXTENT, bs.T_MAX, [bs.T_EXTENT, bs.T_CORE, 5.0, bs.T_MAX],
             "T, Delayed w(p) (negative); dashed line: FWE core threshold", core_line=True)
    ax = fig.add_axes([0.53, 0.45, 0.44, 0.50], facecolor="black")
    cov = nib.load(str(C.GROUP_MAPS / "glm2_low_p_delta_w_covariate_spmT.nii.gz"))
    tmax = float(cov.get_fdata().max())
    d = plotting.plot_stat_map(cov, bg_img=bs.ANAT, threshold=T_COV, display_mode="x", cut_coords=[-51], cmap=bs.COOL_EXT,
                               vmin=T_COV, vmax=tmax, black_bg=True, draw_cross=False, colorbar=False, axes=ax, figure=fig)
    d.add_contours(nib.load(str(C.OWN_MASKS / "left_ifg_pars_triangularis.nii.gz")), levels=[0.5], colors="white", linewidths=1.5)
    colorbar(fig, [0.57, 0.43, 0.34, 0.012], T_COV, tmax, [round(T_COV, 2), round(tmax, 2)], "T, low-p delta w covariate (df = 20)")
    fig.text(0.25, 0.96, "Left AG", ha="center", color="white", fontsize=12, fontweight="bold")
    fig.text(0.75, 0.96, "Left IFG, pars triangularis", ha="center", color="white", fontsize=12, fontweight="bold")

    c = correlation(a.low_p_delta_w, a.delayed_wp_slope_left_ag_cluster)
    scatter(fig.add_axes([0.09, 0.06, 0.37, 0.28]), a.low_p_delta_w, a.delayed_wp_slope_left_ag_cluster, "Left AG cluster",
            [f"Pearson r = {c['r']:.2f}, {p_text(c['p'])}", f"Bonferroni p = {min(1, 3 * c['p']):.3f}",
             f"Spearman rho = {c['rho']:.2f}, {p_text(c['rho_p'])}"])
    c = correlation(b.low_p_delta_w, b.delayed_wp_slope_left_ifg_pars_triangularis)
    scatter(fig.add_axes([0.59, 0.06, 0.37, 0.28]), b.low_p_delta_w, b.delayed_wp_slope_left_ifg_pars_triangularis,
            "Left IFG, pars triangularis", [f"Pearson r = {c['r']:.2f}, {p_text(c['p'])}", f"Spearman rho = {c['rho']:.2f}, {p_text(c['rho_p'])}"])
    fig.text(0.01, 0.97, "A", color="white", fontsize=16, fontweight="bold")
    fig.text(0.51, 0.97, "B", color="white", fontsize=16, fontweight="bold")
    bs.save(fig, "Figure7")


if __name__ == "__main__":
    main()
