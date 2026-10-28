"""Brain-panel style shared by Figures 4, 6 and 7 (black background, MNI152 template from nilearn).

Two layers per map: an extent layer (T >= 2.88, i.e. about p < 0.005, clusters of at least 20 voxels,
6-connectivity) and a core layer showing the clusters that survive whole-brain cluster-level FWE
correction (T > 3.53, i.e. p < 0.001, 18-connectivity, as in SPM).
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.colorbar as mcolorbar  # noqa: E402
import matplotlib.colors as mcolors  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import nibabel as nib  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from nilearn import plotting  # noqa: E402
from nilearn.datasets import load_mni152_template  # noqa: E402
from scipy import ndimage as ndi  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import C  # noqa: E402

T_EXTENT, T_CORE, T_MAX, K_EXTENT = 2.88, 3.53, 7.0, 20
FRAC = (T_CORE - T_EXTENT) / (T_MAX - T_EXTENT)
ANAT = load_mni152_template()
DPI = 300

WARM_EXT = mcolors.LinearSegmentedColormap.from_list("warm_ext", [
    (0, "#BB0000"), (FRAC * .5, "#DD1100"), (FRAC * .85, "#EE3300"), (FRAC, "#FF7700"),
    (FRAC + .18, "#FFB300"), (FRAC + .45, "#FFE600"), (FRAC + .70, "#FFFF99"), (1, "#FFFFFF")], N=512)
WARM_CORE = mcolors.LinearSegmentedColormap.from_list("warm_core", [
    (0, "#FF7700"), (.25, "#FFB300"), (.55, "#FFE600"), (.8, "#FFFF99"), (1, "#FFFFFF")], N=256)
COOL_EXT = mcolors.LinearSegmentedColormap.from_list("cool_ext", [
    (0, "#3300AA"), (FRAC * .5, "#4400CC"), (FRAC * .85, "#6600DD"), (FRAC, "#0055FF"),
    (FRAC + .18, "#0099FF"), (FRAC + .45, "#00CCFF"), (FRAC + .70, "#AAEEFF"), (1, "#FFFFFF")], N=512)
COOL_CORE = mcolors.LinearSegmentedColormap.from_list("cool_core", [
    (0, "#0055FF"), (.25, "#0099FF"), (.55, "#00CCFF"), (.8, "#AAEEFF"), (1, "#FFFFFF")], N=256)


def tmap(name):
    img = nib.load(str(C.GROUP_MAPS / f"{name}_spmT.nii.gz"))
    return img, np.nan_to_num(np.asarray(img.dataobj, np.float32))


def extent_img(name, t_thresh=T_EXTENT, k_min=K_EXTENT):
    img, t = tmap(name)
    lab, n = ndi.label(t >= t_thresh)
    sizes = np.bincount(lab.ravel())
    keep = np.isin(lab, np.flatnonzero(sizes >= k_min)) & (lab > 0)
    return nib.Nifti1Image(np.where(keep, t, 0).astype(np.float32), img.affine)


def core_img(name, table):
    """Clusters (T > 3.53, 18-connectivity) that are significant at cluster-level pFWE < 0.05 in the SPM table."""
    img, t = tmap(name)
    tab = pd.read_csv(C.GROUP_MAPS / f"{table}_clusters.csv")
    sig = tab[tab.cluster_pFWE < 0.05]
    if sig.empty:
        return None
    k_min = sig.cluster_k.min()
    lab, _ = ndi.label(t > T_CORE, ndi.generate_binary_structure(3, 2))
    sizes = np.bincount(lab.ravel())
    keep = np.isin(lab, np.flatnonzero(sizes >= k_min)) & (lab > 0)
    return nib.Nifti1Image(np.where(keep, t, 0).astype(np.float32), img.affine)


def slice_panel(ax, fig, ext, core, mode, coord, cmap_ext, cmap_core, annotate=True):
    d = plotting.plot_stat_map(ext, bg_img=ANAT, threshold=T_EXTENT, display_mode=mode, cut_coords=[coord],
                               cmap=cmap_ext, black_bg=True, vmax=T_MAX, draw_cross=False, annotate=annotate,
                               colorbar=False, axes=ax, figure=fig)
    if core is not None:
        d.add_overlay(core, threshold=T_CORE, cmap=cmap_core, vmax=T_MAX, colorbar=False)
    return d


def split_colorbar(ax, cmap_ext, lo_color, hi_color):
    cb = mcolorbar.ColorbarBase(ax, cmap=cmap_ext, norm=mcolors.Normalize(T_EXTENT, T_MAX), orientation="vertical")
    cb.set_ticks([T_EXTENT, T_CORE, 5.0, T_MAX])
    cb.set_ticklabels([f"{T_EXTENT:.2f}", f"{T_CORE:.2f}", "5.0", f"{T_MAX:.0f}"])
    cb.set_label("T-value", color="white", fontsize=7, labelpad=3)
    cb.ax.yaxis.set_tick_params(color="white", labelsize=6)
    plt.setp(cb.ax.yaxis.get_ticklabels(), color="white")
    ax.axhline(FRAC, color="white", lw=1.0, ls="--", alpha=0.85)
    ax.text(1.4, (FRAC + 1) / 2, "FWE\ncore", transform=ax.transAxes, color=hi_color, fontsize=5.5, va="center")
    ax.text(1.4, FRAC / 2, "Extent\np<.005", transform=ax.transAxes, color=lo_color, fontsize=5.5, va="center")


def save(fig, name):
    C.FIGURES.mkdir(parents=True, exist_ok=True)
    fig.savefig(C.FIGURES / f"{name}.png", dpi=DPI, facecolor=fig.get_facecolor(), bbox_inches="tight", pad_inches=0.1)
    fig.savefig(C.FIGURES / f"{name}.pdf", facecolor=fig.get_facecolor(), bbox_inches="tight", pad_inches=0.1)
    plt.close(fig)
