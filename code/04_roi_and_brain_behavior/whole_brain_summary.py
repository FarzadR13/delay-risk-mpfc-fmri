"""Whole-brain and small-volume results (Figures 4, 6 and 7B) from the SPM cluster tables and group maps.

The tables in data/group_maps/*_clusters.csv are written by code/03_second_level/group_tables.m
(p < 0.001 uncorrected, cluster-level FWE). This script reads them, and computes the map-level
quantities that need no SPM: the cluster-forming threshold, the voxelwise similarity of the GLM1 and
GLM2 distance maps, the distance between the two mPFC peaks and the overlap of the two mPFC clusters.
"""
import sys
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import C, Recorder, save_source  # noqa: E402


def table(name):
    t = pd.read_csv(C.GROUP_MAPS / f"{name}_clusters.csv")
    t[["cluster_pFWE", "cluster_k"]] = t[["cluster_pFWE", "cluster_k"]].ffill()   # secondary peaks of a cluster
    return t


def nearest(t, xyz):
    d = np.linalg.norm(t[["x", "y", "z"]].to_numpy(float) - np.asarray(xyz, float), axis=1)
    return t.iloc[int(np.argmin(d))]


def xyz(row):
    return "[" + " ".join(f"{int(round(v))}" for v in (row.x, row.y, row.z)) + "]"


def load(name):
    return np.asarray(nib.load(str(C.GROUP_MAPS / f"{name}_spmT.nii.gz")).dataobj, float)


def mask(name):
    return np.asarray(nib.load(str(C.OWN_MASKS / f"{name}.nii.gz")).dataobj) > 0


def main():
    rec = Recorder("whole_brain_summary")
    u = stats.t.isf(C.CLUSTER_P, C.DF_GROUP)
    rec.add("meth_t_threshold", u)

    t = table("glm1_immediate_distance_positive")
    m = nearest(t, C.PEAK_FIG4A_MPFC)
    for k, v in dict(xyz=xyz(m), T=m.peak_T, Z=m.peak_Z, k=m.cluster_k, pfwe=m.cluster_pFWE).items():
        rec.add(f"fig4a_mpfc_{k}", v)
    m = nearest(t, (60, -40, -1))
    for k, v in dict(xyz=xyz(m), k=m.cluster_k, pfwe=m.cluster_pFWE).items():
        rec.add(f"fig4a_mtg_{k}", v)
    svc = pd.read_csv(C.GROUP_MAPS / "glm1_immediate_distance_positive_svc_bartra.csv")
    top = svc.loc[svc.peak_T.idxmax()]
    for k, v in dict(xyz=xyz(top), T=top.peak_T, Z=top.peak_Z, pfwe=top.peak_pFWE).items():
        rec.add(f"fig4a_svc_{k}", v)
    rec.add("fig4_delayed_min_pfwe", min(table(f"glm1_delayed_distance_{s}").cluster_pFWE.min() for s in ("positive", "negative")))
    rec.add("fig4b_xyz", xyz(nearest(table("glm1_immediate_minus_delayed_distance_positive"), C.PEAK_FIG4A_MPFC)))

    t = table("glm2_delayed_wp_negative")
    for key, peak in C.PEAK_FIG6A.items():
        m = nearest(t, peak)
        k2 = "ag" if key == "left_ag" else key
        for k, v in dict(xyz=xyz(m), Z=m.peak_Z, T=m.peak_T, k=m.cluster_k, pfwe=m.cluster_pFWE).items():
            rec.add(f"fig6a_{k2}_{k}", v)
    rec.add("fig6_imm_min_pfwe", min(table(f"glm2_immediate_wp_{s}").cluster_pFWE.min() for s in ("positive", "negative")))
    rec.add("fig6_direct_min_pfwe", table("glm2_delayed_more_negative_wp").cluster_pFWE.min())
    direct = load("glm2_delayed_more_negative_wp")
    rec.add("fig6_direct_opposite", int(np.sum(np.nan_to_num(-direct) > u)))

    p4 = nearest(table("glm1_immediate_distance_positive"), C.PEAK_FIG4A_MPFC)
    p6 = nearest(table("glm2_delayed_wp_negative"), C.PEAK_FIG6A["mpfc"])
    rec.add("fig6_peak_distance", float(np.linalg.norm(p4[["x", "y", "z"]].to_numpy(float) - p6[["x", "y", "z"]].to_numpy(float))))
    rec.add("fig6_overlap", int((mask("fig4a_mpfc_cluster") & mask("fig6a_mpfc_cluster")).sum()))

    t = table("glm2_low_p_delta_w_covariate")
    rec.add("fig7_cov_none", "none" if t.cluster_pFWE.min() >= 0.05 else "some", note=f"smallest cluster pFWE = {t.cluster_pFWE.min():.3f}")
    top = t.loc[t.peak_T.idxmax()]
    rec.add("fig7_cov_peak_xyz", xyz(top)); rec.add("fig7_cov_peak_T", top.peak_T)
    m = nearest(t, (-51, 20, 2))
    rec.add("fig7_ifg_xyz", xyz(m)); rec.add("fig7_ifg_T", m.peak_T); rec.add("fig7_ifg_k", m.cluster_k)

    a = np.asarray(nib.load(str(C.OWN_MASKS / "analysis_mask.nii.gz")).dataobj) > 0
    g1, g2 = load("glm1_immediate_distance"), load("glm2_immediate_distance")
    ok = a & np.isfinite(g1) & np.isfinite(g2)
    rec.add("meth_glm1_glm2_r", stats.pearsonr(g1[ok], g2[ok]).statistic)
    rec.add("meth_n_fmri", 22)
    rec.save()
    for panel, name in (("Fig4A_map", "glm1_immediate_distance_positive"), ("Fig4B_map", "glm1_immediate_minus_delayed_distance_positive"),
                        ("Fig6A_map", "glm2_delayed_wp_negative"), ("Fig7B_map", "glm2_low_p_delta_w_covariate")):
        save_source(f"{panel}_cluster_table", pd.read_csv(C.GROUP_MAPS / f"{name}_clusters.csv"))


if __name__ == "__main__":
    main()
