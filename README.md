# delay-risk-mpfc-fmri

Code and data for a model-based fMRI study of how a 12-month delay changes probability weighting and medial prefrontal value-comparison signals during risky choice.

## Summary

Healthy volunteers chose between a sure amount and a gamble (probability p of winning €20) whose outcome was paid either immediately or after 12 months. Certainty equivalents measured before scanning were used to fit a two-parameter probability-weighting function and to place the scanner offers at fixed distances (−30% to +30%) from each participant's certainty equivalent. This repository reproduces the fMRI results shown in Figures 4 to 7 of the accompanying manuscript, the paired tests of the weighting parameters, and the related statements of the Methods and Discussion from the shared data, and documents how the SPM12 first- and second-level models were built.

## Quick start

```bash
python -m venv .venv && source .venv/bin/activate      # Windows: .venv\Scripts\activate
pip install -r requirements.txt
python run_all.py
```

`run_all.py` reads `data/` (plus a one-time download of third-party atlases, see Masks), writes every statistic to `results/stats_reproduced.csv`, compares them with the printed values in `results/check_against_manuscript.csv`, writes the figures to `results/figures/` and the plotted data to `source_data/Source_Data.xlsx`. It takes under a minute after the download.

## Figure panels, statistics, scripts and data

| Result | Script | Data |
|---|---|---|
| Paired tests of the weighting parameters (Results) | `code/01_behavior/parameter_tests.py` | `behavior/lattimore_params.csv` |
| Weighting changes per probability, elevation/curvature decomposition (Discussion, Methods) | `code/01_behavior/weighting_decomposition.py` | `behavior/prescan_ce.csv`, `behavior/lattimore_params.csv` |
| Fig 4A,B and whole-brain statistics | `code/06_figures/figure4.py`; `code/04_roi_and_brain_behavior/whole_brain_summary.py` | `group_maps/glm1_*` |
| Leave-one-subject-out mPFC ROI (Results, Fig 4) | `code/04_roi_and_brain_behavior/loso_roi.py` | `roi/loso_mpfc.csv` |
| Fig 5A,B | `code/06_figures/figure5.py`; `heatmap_gradients.py` | `roi/heatmap_cells.csv` |
| Fig 5C,D | `code/06_figures/figure5.py`; `main_effect_gamma.py` | `roi/bartra_mpfc.csv`, `behavior/lattimore_params.csv` |
| Fig 6A and whole-brain statistics | `code/06_figures/figure6.py`; `whole_brain_summary.py` | `group_maps/glm2_*` |
| Fig 6B | `code/06_figures/figure6.py`; `wp_bin_profiles.py` | `roi/wp_bins.csv`, `behavior/prescan_ce.csv` |
| Fig 7A,B | `code/06_figures/figure7.py`; `lowp_brain_behavior.py`, `whole_brain_summary.py` | `roi/wp_slopes.csv`, `behavior/prescan_ce.csv`, `group_maps/glm2_low_p_delta_w_covariate_*` |

## Pipeline

| Step | Code | Language | Needs |
|---|---|---|---|
| Masks | `code/00_masks/build_masks.py` | Python | internet once (Kable lab, nilearn) |
| Behavior | `code/01_behavior/` | Python | `data/behavior` |
| First level | `code/02_first_level/` (`make_events.m`, `glm1.m`, `glm2.m`, `glm3.m`, `glm_heatmap.m`, `run_first_level.m`, `run_all_first_level.m`) | MATLAB, SPM12 | preprocessed images |
| Second level | `code/03_second_level/` (`one_sample_ttests.m`, `low_p_covariate_model.m`, `group_tables.m`, `cluster_table.m`) | MATLAB, SPM12 | participant-level contrast images |
| ROI extraction | `code/04_roi_and_brain_behavior/extract_roi_values.py` | Python | participant-level contrast images |
| ROI and brain-behavior statistics | `code/04_roi_and_brain_behavior/` | Python | `data/roi`, `data/group_maps` |
| Figures | `code/06_figures/` | Python | `data/` |
| Check | `code/07_check/check_against_manuscript.py` | Python | `results/stats_reproduced.csv` |

The preprocessed and participant-level images are available on reasonable request (see Contact); with them, the MATLAB scripts and `extract_roi_values.py` regenerate `data/group_maps` and `data/roi`.

### First-level models (SPM12)

All models: one session per run (four runs; delayed runs first, then immediate runs), stimulus-onset stick functions, canonical HRF, units in seconds, TR = 2 s, microtime resolution 42 and onset 21, six realignment parameters per run, 128-s high-pass filter, AR(1), implicit mask threshold 0.8, no global normalisation. Parametric modulators are not orthogonalized and are entered as they are (SPM12 does not mean-centre them when orthogonalization is off; the parametric-modulator estimates are the same as with centred modulators because the condition regressor is in the model).

| Model | Paper name | Regressors per run | Contrasts used |
|---|---|---|---|
| `glm1.m` | GLM1 | condition; × signed distance (% of CE); × objective probability | immediate distance, delayed distance, immediate minus delayed distance |
| `glm2.m` | GLM2 | condition; × signed distance; × w(p) | immediate and delayed w(p), delayed minus immediate w(p), distance slopes |
| `glm3.m` | GLM3 | condition only | immediate and delayed main effects (sum of the two runs) |
| `glm_heatmap.m` | heatmap model | nine cells: 3 probability bins × 3 distance bins | the 18 cell contrasts |

In GLM2, w(p) is the participant's weight w(p) = CE/20 at the probability of the trial, from the pre-scan certainty equivalents (`data/behavior/prescan_ce.csv`). `make_events.m` builds the events from `data/behavior/scan_choices.csv`: GLM2, GLM3 and the heatmap model use the valid responses to scheduled trials (`included_in_glm`); GLM1 uses every valid response (`included_in_glm1`), with re-presented trials entered with distance = 0 and probability = 0. Each model script takes a preprocessed-data root laid out as `<root>/<participant_id>/run-<r>/swr*.nii` and `rp_*.txt`. Rebuilt with these scripts, the design matrices of all participants and models are identical to those of the analysis reported in the paper (maximum absolute difference 2e-14).

### Brain-behavior predictor

low-p Δw (Figure 7) is the mean, over p = 0.05 and 0.10, of the delayed-minus-immediate weight w(p) = CE/20 (`code/common.py`, `delta_w`). The sensitivity predictors add p = 0.25, then p = 0.50; the high-probability predictor uses p = 0.90 and 0.95. The whole-brain covariate model (`low_p_covariate_model.m`) uses the same predictor, mean-centred.

### Masks

`code/00_masks/build_masks.py` downloads the Bartra, McGuire and Kable (2013) Figure 3 maps from the Kable lab (SHA-256 checked) and fetches the Harvard-Oxford atlases through nilearn, then rebuilds: the meta-analytic mPFC mask (Figure 3A positive subjective-value map, |x| ≤ 20 mm, y ≥ 20 mm, minus striatal voxels, component containing [2 46 −8], 3-mm grid, analysis mask; 460 voxels) and the Harvard-Oxford mPFC (736 voxels; search region of the leave-one-subject-out ROI), left angular gyrus (245) and left IFG pars triangularis (118) masks. Masks derived from our own data are in `data/masks/`.

## Data dictionary

`data/participants.tsv`: `participant_id` (sub-01 to sub-24), `in_fmri_sample` (yes/no; two participants were excluded from the fMRI analyses for head motion).

`data/behavior/prescan_ce.csv`: `participant_id`, `condition` (immediate/delayed), `probability`, `ce_eur` (pre-scan certainty equivalent, €), `w_empirical` (= `ce_eur` / 20).

`data/behavior/lattimore_params.csv`: `participant_id`, `condition`, `delta` (elevation), `gamma` (curvature) of w(p) = δp^γ / (δp^γ + (1 − p)^γ), as used in the paper.

`data/behavior/scan_choices.csv` (one row per presented trial; events of the first-level models): `participant_id`, `run` (1-4, acquisition order), `condition`, `trial` (presentation order within the run), `onset_s` (stimulus onset, s from run start), `probability`, `sure_amount_eur`, `distance_level` (sure amount relative to the CE, %: −30 to +30; GLM distance modulator), `signed_distance` (sure amount minus CE, €), `choice` (risky/sure/none), `rt_s`, `is_replay` (1 = re-presentation of a missed trial), `included_in_glm` (1 = modeled in GLM2, GLM3 and the heatmap model), `included_in_glm1` (1 = modeled in GLM1).

`data/roi/` (N = 22; values are means of first-level contrast estimates within the ROI, arbitrary units):
- `loso_mpfc.csv`: `participant_id`, `peak_x/y/z` (leave-one-subject-out sphere center, MNI mm), `n_voxels`, `immediate_slope`, `delayed_slope` (GLM1 distance).
- `bartra_mpfc.csv`: meta-analytic mPFC mask; `glm1_immediate_distance`, `glm1_delayed_distance`, `glm3_immediate_main_effect`, `glm3_delayed_main_effect`.
- `heatmap_cells.csv`: `condition`, `probability_bin` (1: .05/.10, 2: .25/.50/.75, 3: .90/.95), `distance_bin` (1: −30/−20, 2: −10/0/+10, 3: +20/+30), `cell`, `value`.
- `wp_slopes.csv`: `roi` (Figure 6A clusters `cluster_*`; Harvard-Oxford `ho_left_ag`, `ho_left_ifg_tri`), `n_voxels`, `glm2_delayed_wp`, `glm2_immediate_wp`.
- `wp_bins.csv`: heatmap-model responses in the Figure 6A clusters averaged over distance bins: `roi`, `condition`, `probability_bin`, `value`.

`data/group_maps/`: unthresholded group maps (`*_spmT.nii.gz`, `*_con.nii.gz`; MNI, 3 mm, 53 × 63 × 52; N = 22) and SPM tables (`*_clusters.csv`; p < 0.001 uncorrected; columns `cluster_pFWE`, `cluster_k`, `peak_pFWE`, `peak_T`, `peak_Z`, `peak_p_unc`, `x`, `y`, `z`; secondary peaks of a cluster have empty cluster columns). `glm1_immediate_distance_positive_svc_bartra.csv` is the small-volume correction within the meta-analytic mPFC mask. The `*_negative`, `*_more_negative` and `*_less_negative` tables test the negative direction or the opposite contrast.

`data/masks/`: `analysis_mask` (group analysis mask, 43,935 voxels), `fig4a_mpfc_cluster` (48 voxels), `fig6a_mpfc_cluster`, `fig6a_pcc_cluster` (71 each), `fig6a_left_ag_cluster` (111), `left_ifg_pars_triangularis` (Harvard-Oxford, resampled to the analysis grid).

`results/check_against_manuscript.csv`: `stat_id`, `figure_or_section`, `quantity`, `manuscript_value`, `reproduced_value`, `status` (MATCH when the reproduced value, rounded to the printed precision, equals the printed value).

## Data availability

Raw MRI data, preprocessed images and participant-level first-level contrast images are available on reasonable request (see Contact). They are not posted publicly to protect participant privacy.

## Software

- MATLAB R2024a (24.1.0.2537033); SPM12 revision 7771 (`spm('Ver')` = SPM12, `spm('Version')` = SPM12 (7771)). The first-level models of the paper were estimated with an earlier SPM12 release (`spm_spm` v6015); re-estimation with revision 7771 gives contrast images with voxelwise r > 0.99999, and the second-level models re-estimated from the participant-level images reproduce the group maps and tables exactly.
- Python 3.10.11 with the versions pinned in `requirements.txt` (numpy 2.2.6, scipy 1.15.3, pandas 2.3.3, nibabel 5.4.2, nilearn 0.13.1, scikit-learn 1.7.2, matplotlib 3.10.9, openpyxl 3.1.5).

## Licenses

Code: MIT (`LICENSE`). Data in `data/` and `source_data/`: CC BY 4.0 (`LICENSE-DATA`).

## Contact

Farzad Rostami, Institut des Sciences Cognitives Marc Jeannerod, CNRS UMR5229, Lyon, France. For questions or data requests, please open an issue on this repository.
