"""Extract every participant-level ROI value from the first-level contrast images.

The participant-level images are available on reasonable request. Expected
layout: <contrasts>/<participant_id>/<participant_id>_<contrast>.nii.gz, with the contrast names used by
code/02_first_level (for example glm1_immediate_distance, glm2_delayed_wp, heatmap_immediate_cell3).
Writes, in --out (default data/roi):
  loso_mpfc.csv        leave-one-subject-out mPFC spheres, GLM1 immediate and delayed distance slopes
  bartra_mpfc.csv      meta-analytic mPFC mask means: GLM1 distance slopes and GLM3 main effects
  heatmap_cells.csv    meta-analytic mPFC mask means of the 18 heatmap-model cells
  wp_slopes.csv        GLM2 w(p) slopes in the Figure 6A clusters and the left AG and left IFG masks
  wp_bins.csv          heatmap-model probability bins (mean over distance bins) in the Figure 6A clusters
Requires the masks built by code/00_masks/build_masks.py.
"""
import argparse
import sys
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import C, participants  # noqa: E402

GRID = {(pb, db): (db - 1) * 3 + pb for pb in (1, 2, 3) for db in (1, 2, 3)}  # heatmap cell number


def load_mask(path):
    return np.asarray(nib.load(str(path)).dataobj) > 0


def image(root, pid, name):
    return nib.load(str(root / pid / f"{pid}_{name}.nii.gz"))


def stack(root, ids, name, ref):
    out = []
    for pid in ids:
        img = image(root, pid, name)
        if img.shape != ref.shape or not np.allclose(img.affine, ref.affine):
            raise AssertionError(f"grid mismatch: {pid} {name}")
        out.append(np.asarray(img.dataobj, dtype=np.float64))
    return np.stack(out)


def roi_mean(arr, mask):
    v = arr[:, mask].mean(axis=1)
    if not np.isfinite(v).all():
        raise AssertionError("non-finite ROI value")
    return v


def world(img):
    ijk = np.indices(img.shape).reshape(3, -1).T
    return nib.affines.apply_affine(img.affine, ijk).reshape((*img.shape, 3))


def loso(root, ids, ref, analysis):
    """For each participant, peak of the immediate distance t map of the other 21 within HO mPFC."""
    search = load_mask(C.DERIVED_MASKS / "ho_mpfc.nii.gz")
    imm = stack(root, ids, "glm1_immediate_distance", ref)
    dly = stack(root, ids, "glm1_delayed_distance", ref)
    xyz = world(ref)
    rows = []
    for i, pid in enumerate(ids):
        others = np.delete(imm, i, axis=0)
        t = others.mean(0) / (others.std(0, ddof=1) / np.sqrt(len(others)))
        peak = np.unravel_index(np.argmax(np.where(search & np.isfinite(t), t, -np.inf)), t.shape)
        sphere = (np.linalg.norm(xyz - xyz[peak], axis=-1) <= C.LOSO_RADIUS_MM) & analysis
        rows.append(dict(participant_id=pid, peak_x=xyz[peak][0], peak_y=xyz[peak][1], peak_z=xyz[peak][2],
                         n_voxels=int(sphere.sum()), immediate_slope=imm[i][sphere].mean(),
                         delayed_slope=dly[i][sphere].mean()))
    return pd.DataFrame(rows)


def main(contrast_root, out):
    out.mkdir(parents=True, exist_ok=True)
    ids = participants(fmri_only=True)
    ref = nib.load(str(C.OWN_MASKS / "analysis_mask.nii.gz"))
    analysis = np.asarray(ref.dataobj) > 0
    bartra = load_mask(C.DERIVED_MASKS / "bartra_mpfc.nii.gz")

    loso(contrast_root, ids, ref, analysis).to_csv(out / "loso_mpfc.csv", index=False)

    b = pd.DataFrame(dict(participant_id=ids))
    for name in ("glm1_immediate_distance", "glm1_delayed_distance",
                 "glm3_immediate_main_effect", "glm3_delayed_main_effect"):
        b[name] = roi_mean(stack(contrast_root, ids, name, ref), bartra)
    b.to_csv(out / "bartra_mpfc.csv", index=False)

    rois = {"cluster_mpfc": load_mask(C.OWN_MASKS / "fig6a_mpfc_cluster.nii.gz"),
            "cluster_pcc": load_mask(C.OWN_MASKS / "fig6a_pcc_cluster.nii.gz"),
            "cluster_left_ag": load_mask(C.OWN_MASKS / "fig6a_left_ag_cluster.nii.gz"),
            "ho_left_ag": load_mask(C.DERIVED_MASKS / "ho_left_ag.nii.gz"),
            "ho_left_ifg_tri": load_mask(C.DERIVED_MASKS / "ho_left_ifg_tri.nii.gz")}

    cells, bins = [], []
    for cond in C.CONDITIONS:
        arrays = {n: stack(contrast_root, ids, f"heatmap_{cond}_cell{n}", ref) for n in range(1, 10)}
        for (pb, db), n in GRID.items():
            for pid, v in zip(ids, roi_mean(arrays[n], bartra)):
                cells.append(dict(participant_id=pid, condition=cond, probability_bin=pb, distance_bin=db, cell=n, value=v))
        for roi, m in rois.items():
            if not roi.startswith("cluster_"):
                continue
            for pb in (1, 2, 3):
                v = np.mean([roi_mean(arrays[GRID[(pb, db)]], m) for db in (1, 2, 3)], axis=0)
                for pid, x in zip(ids, v):
                    bins.append(dict(participant_id=pid, roi=roi, condition=cond, probability_bin=pb, value=x))
    pd.DataFrame(cells).to_csv(out / "heatmap_cells.csv", index=False)
    pd.DataFrame(bins).to_csv(out / "wp_bins.csv", index=False)

    slopes = {n: stack(contrast_root, ids, n, ref) for n in ("glm2_delayed_wp", "glm2_immediate_wp")}
    rows = []
    for roi, m in rois.items():
        vals = {n: roi_mean(a, m) for n, a in slopes.items()}
        for i, pid in enumerate(ids):
            rows.append(dict(participant_id=pid, roi=roi, n_voxels=int(m.sum()), **{n: v[i] for n, v in vals.items()}))
    pd.DataFrame(rows).to_csv(out / "wp_slopes.csv", index=False)

    print("ROI values written to", out)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--contrasts", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=C.ROI)
    a = ap.parse_args()
    main(a.contrasts, a.out)
