"""Fit the Lattimore weighting function to the empirical weights w(p) = CE / 20 and compare with the parameters used.

Least squares on the seven probabilities, 100 random starting values from a fixed seed, best residual
sum of squares retained (w(0) = 0 and w(1) = 1 hold by construction of the function). Reports, per
participant and condition, the difference from data/behavior/lattimore_params.csv, which all downstream
analyses use.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import C, Recorder, empirical_weights, lattimore, lattimore_params, save_table  # noqa: E402

SEED, N_STARTS, TOLERANCE = 2024, 100, 0.01


def fit(p, w, rng):
    best = None
    for start in rng.uniform([0.05, 0.05], [3.0, 2.0], size=(N_STARTS, 2)):
        r = least_squares(lambda x: lattimore(p, x[0], x[1]) - w, start, bounds=([1e-6, 1e-6], [50, 20]))
        if best is None or r.cost < best.cost:
            best = r
    return best.x, 2 * best.cost


def main():
    rng = np.random.default_rng(SEED)
    w = empirical_weights()
    orig = lattimore_params()
    p = np.array(C.PROBABILITIES)
    rows = []
    for pid in w.index:
        for cond in C.CONDITIONS:
            (delta, gamma), rss = fit(p, w.loc[pid, cond].to_numpy(float), rng)
            rows.append(dict(participant_id=pid, condition=cond, delta_refit=delta, gamma_refit=gamma, rss=rss,
                             delta_original=orig.loc[pid, f"delta_{cond}"], gamma_original=orig.loc[pid, f"gamma_{cond}"]))
    t = pd.DataFrame(rows)
    t["max_abs_diff"] = np.maximum(abs(t.delta_refit - t.delta_original), abs(t.gamma_refit - t.gamma_original))
    save_table("lattimore_refit", t)
    rec = Recorder("fit_lattimore")
    rec.add("refit_max_abs_diff", t.max_abs_diff.max())
    rec.add("refit_n_above_tolerance", int((t.max_abs_diff > TOLERANCE).sum()),
            note="participant x condition fits differing from the original parameters by more than 0.01")
    rec.save()
    off = t[t.max_abs_diff > TOLERANCE]
    print(off.round(4).to_string(index=False) if len(off) else f"all refits within {TOLERANCE} of the original parameters")


if __name__ == "__main__":
    main()
