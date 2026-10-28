"""Figure 7: delay-related change in low-probability weighting vs the delayed w(p) slope (GLM2), N = 22.

low-p delta w = mean over p = 0.05 and 0.10 of the delayed-minus-immediate empirical weight w(p) = CE / 20
(data/behavior/prescan_ce.csv). The neural measure is the mean of the GLM2 delayed w(p) contrast image
within each ROI. Bonferroni correction over the three Figure 6A clusters; the left IFG is reported
uncorrected. Sensitivity predictors add p = 0.25, then p = 0.50; the high-probability predictor uses
p = 0.90 and 0.95.
"""
import sys
from pathlib import Path

import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import C, Recorder, correlation, delta_w, leave_one_out_r, participants, save_source  # noqa: E402

CLUSTERS = {"cluster_mpfc": "mpfc", "cluster_pcc": "pcc", "cluster_left_ag": "ag"}


def main():
    ids = participants(fmri_only=True)
    s = pd.read_csv(C.ROI / "wp_slopes.csv").pivot(index="participant_id", columns="roi", values="glm2_delayed_wp").loc[ids]
    x = delta_w([0.05, 0.10], ids)
    rec = Recorder("lowp_brain_behavior")
    for roi, key in CLUSTERS.items():
        c = correlation(x, s[roi])
        rec.add_many(f"fig7_{key}", c, dict(r="r", ci_low="r_ci_low", ci_high="r_ci_high", p="p", rho="rho", rho_p="rho_p"))
        rec.add(f"fig7_{key}_pbonf", min(1.0, 3 * c["p"]))
    loo = leave_one_out_r(x, s.cluster_left_ag)
    rec.add("fig7_ag_robust_min", loo.min()); rec.add("fig7_ag_robust_max", loo.max())
    ts = stats.theilslopes(s.cluster_left_ag, x, 0.95)
    rec.add("fig7_ag_theilsen_slope", ts.slope, note=f"95% CI [{ts.low_slope:.2f}, {ts.high_slope:.2f}]")
    rec.add_many("fig7_hoag", correlation(x, s.ho_left_ag), dict(r="r", p="p"))
    rec.add_many("fig7_ifg", correlation(x, s.ho_left_ifg_tri), dict(r="r", p="p", rho="rho", rho_p="rho_p"))
    rec.add_many("fig7_sens3", correlation(delta_w([0.05, 0.10, 0.25], ids), s.ho_left_ag), dict(r="r", p="p"))
    rec.add_many("fig7_sens4", correlation(delta_w([0.05, 0.10, 0.25, 0.50], ids), s.ho_left_ag), dict(r="r", p="p"))
    high = delta_w([0.90, 0.95], ids)
    rec.add("fig7_highp_min_pbonf", min(min(1.0, 3 * correlation(high, s[roi])["p"]) for roi in CLUSTERS))
    rec.save()
    src = pd.DataFrame(dict(participant_id=ids, low_p_delta_w=x.to_numpy(), delayed_wp_slope_left_ag_cluster=s.cluster_left_ag.to_numpy(),
                            delayed_wp_slope_left_ifg_pars_triangularis=s.ho_left_ifg_tri.to_numpy()))
    save_source("Fig7A", src[["participant_id", "low_p_delta_w", "delayed_wp_slope_left_ag_cluster"]])
    save_source("Fig7B", src[["participant_id", "low_p_delta_w", "delayed_wp_slope_left_ifg_pars_triangularis"]])

if __name__ == "__main__":
    main()
