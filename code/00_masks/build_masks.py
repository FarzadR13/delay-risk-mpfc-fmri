"""Rebuild the masks derived from third-party atlases (not redistributed in this repository).

1. Meta-analytic mPFC mask: Bartra, McGuire and Kable (2013) Figure 3A (positive subjective-value
   effects), restricted to |x| <= 20 mm and y >= 20 mm, minus Harvard-Oxford striatal voxels, keeping
   the component that contains the ventromedial peak [2 46 -8], resampled (nearest neighbor) to the
   3-mm analysis grid and intersected with the group analysis mask. Expected size: 460 voxels.
2. Harvard-Oxford anatomical masks (cort-maxprob-thr25-2mm), resampled to the analysis grid:
   mPFC (Frontal Medial Cortex + Paracingulate Gyrus; search region of the leave-one-subject-out ROI),
   left angular gyrus (x < 0) and left IFG pars triangularis (x < 0).

The only network access of the pipeline happens here (Kable lab download and nilearn atlas fetch).
"""
import hashlib
import sys
import urllib.request
import zipfile

import nibabel as nib
import numpy as np
from nilearn.datasets import fetch_atlas_harvard_oxford
from nilearn.image import resample_to_img
from scipy import ndimage as ndi

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[1]))
from common import C, Recorder  # noqa: E402

CONN = ndi.generate_binary_structure(3, 2)   # 18-connectivity, as in SPM
EXPECTED = {"bartra_mpfc": 460, "ho_mpfc": 736, "ho_left_ag": 245, "ho_left_ifg_tri": 118}


def world(img):
    ijk = np.indices(img.shape).reshape(3, -1).T
    return nib.affines.apply_affine(img.affine, ijk).reshape((*img.shape, 3))


def to_grid(img, ref):
    return np.asarray(resample_to_img(img, ref, interpolation="nearest", force_resample=True,
                                      copy_header=True).dataobj)


def atlas(kind):
    item = fetch_atlas_harvard_oxford(kind, data_dir=str(C.CACHE / "nilearn"))
    img = item.maps if isinstance(item.maps, nib.spatialimages.SpatialImage) else nib.load(str(item.maps))
    return img, list(item.labels)


def bartra_source():
    C.CACHE.mkdir(parents=True, exist_ok=True)
    z = C.CACHE / "nifti_roi-masks_fig03.zip"
    if not z.exists():
        urllib.request.urlretrieve(C.BARTRA_URL, z)
    digest = hashlib.sha256(z.read_bytes()).hexdigest()
    if digest != C.BARTRA_SHA256:
        raise RuntimeError(f"SHA-256 mismatch for {z.name}: {digest}")
    out = C.CACHE / "bartra_fig03"
    with zipfile.ZipFile(z) as f:
        f.extractall(out)
    return nib.load(str(next(out.rglob(C.BARTRA_FILE.split("/")[-1]))))


def save(mask, ref, name):
    C.DERIVED_MASKS.mkdir(parents=True, exist_ok=True)
    nib.save(nib.Nifti1Image(mask.astype(np.uint8), ref.affine), str(C.DERIVED_MASKS / f"{name}.nii.gz"))


def main():
    ref = nib.load(str(C.OWN_MASKS / "analysis_mask.nii.gz"))
    analysis = np.asarray(ref.dataobj) > 0
    xyz = world(ref)
    masks = {}

    # 1. Bartra et al. (2013) Figure 3A, medial prefrontal component
    src = bartra_source()
    fig3a = np.asarray(src.dataobj) > 0
    sxyz = world(src)
    box = fig3a & (np.abs(sxyz[..., 0]) <= 20) & (sxyz[..., 1] >= 20)
    sub, sub_labels = atlas("sub-maxprob-thr25-2mm")
    striatal = [i for i, l in enumerate(sub_labels) if any(k in l for k in ("Caudate", "Putamen", "Pallidum", "Accumbens"))]
    kept = box & ~np.isin(to_grid(sub, src).astype(int), striatal)
    lab, _ = ndi.label(kept, CONN)
    dist = np.where(kept, np.linalg.norm(sxyz - np.array([2, 46, -8]), axis=-1), np.inf)
    component = lab == lab[np.unravel_index(np.argmin(dist), dist.shape)]
    masks["bartra_mpfc"] = (to_grid(nib.Nifti1Image(component.astype(np.uint8), src.affine), ref) > 0) & analysis

    # 2. Harvard-Oxford cortical masks
    cort, labels = atlas("cort-maxprob-thr25-2mm")
    a = to_grid(cort, ref).astype(int)

    def region(*names):
        return np.isin(a, [labels.index(n) for n in names])

    left = xyz[..., 0] < 0
    masks["ho_mpfc"] = region("Frontal Medial Cortex", "Paracingulate Gyrus") & analysis
    masks["ho_left_ag"] = region("Angular Gyrus") & left & analysis
    ifg_full = region("Inferior Frontal Gyrus, pars triangularis") & left
    masks["ho_left_ifg_tri"] = ifg_full & analysis

    rec = Recorder("build_masks")
    for name, m in masks.items():
        save(m, ref, name)
        n = int(m.sum())
        if n != EXPECTED[name]:
            raise AssertionError(f"{name}: {n} voxels, expected {EXPECTED[name]}")
        rec.add(f"mask_{name}_voxels", n)
    rec.add("meth_bartra_voxels", int(masks["bartra_mpfc"].sum()))
    shared = np.asarray(nib.load(str(C.OWN_MASKS / "left_ifg_pars_triangularis.nii.gz")).dataobj) > 0
    if not np.array_equal(shared, ifg_full):
        raise AssertionError("rebuilt left IFG mask differs from data/masks/left_ifg_pars_triangularis.nii.gz")
    rec.save()
    print("masks:", {k: int(v.sum()) for k, v in masks.items()})


if __name__ == "__main__":
    main()
