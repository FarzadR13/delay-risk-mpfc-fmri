"""Discussion: elevation/curvature decomposition of the fitted w(p) change, and empirical w(p) changes per probability.

For each participant, the delayed-minus-immediate change in fitted w(p) is split into the change obtained by
replacing only delta (elevation) or only gamma (curvature) with its delayed value. Behavioral sample, N = 24.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import C, Recorder, empirical_weights, lattimore, lattimore_params, one_sample, paired, save_table  # noqa: E402


def main():
    w = lattimore_params()
    p = np.array(C.PROBABILITIES)
    dn, dd, gn, gd = (w[c].to_numpy()[:, None] for c in ("delta_immediate", "delta_delayed", "gamma_immediate", "gamma_delayed"))
    base = lattimore(p, dn, gn)
    parts = dict(total=lattimore(p, dd, gd) - base, elevation_only=lattimore(p, dd, gn) - base,
                 curvature_only=lattimore(p, dn, gd) - base)
    rows = []
    for j, pr in enumerate(p):
        row = dict(probability=pr)
        for k, v in parts.items():
            s = one_sample(v[:, j]); row[f"{k}_mean"] = s["mean"]; row[f"{k}_p"] = s["p"]
        rows.append(row)
    dec = pd.DataFrame(rows)
    save_table("weighting_decomposition", dec)

    e = empirical_weights()
    emp = []
    for pr in p:
        s = paired(e[("delayed", pr)], e[("immediate", pr)])
        emp.append(dict(probability=pr, immediate_mean=e[("immediate", pr)].mean(), delayed_mean=e[("delayed", pr)].mean(),
                        diff=s["mean"], ci_low=s["ci_low"], ci_high=s["ci_high"], t=s["t"], p=s["p"]))
    emp = pd.DataFrame(emp)
    save_table("empirical_weight_changes", emp)

    rec = Recorder("weighting_decomposition")
    rec.add("beh_wemp_p090_p", emp.loc[emp.probability == 0.90, "p"].item())
    elev_up = bool((dec.elevation_only_mean > 0).all())
    curv_shape = bool((dec.curvature_only_mean[p <= 0.25] > 0).all() and (dec.curvature_only_mean[p >= 0.75] < 0).all())
    rec.add("beh_decomp", "consistent" if elev_up and curv_shape else "not consistent",
            note=f"elevation-only change > 0 at all p: {elev_up}; curvature-only > 0 for p <= .25 and < 0 for p >= .75: {curv_shape}")
    rec.save()


if __name__ == "__main__":
    main()
