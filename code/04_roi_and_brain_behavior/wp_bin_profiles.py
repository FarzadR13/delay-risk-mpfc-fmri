"""Figure 6B: heatmap-model responses per probability bin in the Figure 6A clusters, N = 22.

Responses are averaged over the three distance bins and plotted against the mean empirical weight
w(p) = CE / 20 of each probability bin and condition (descriptive).
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import C, empirical_weights, participants, save_table  # noqa: E402


def bin_weights(ids):
    w = empirical_weights().loc[ids]
    rows = []
    for cond in C.CONDITIONS:
        for b in (1, 2, 3):
            ps = [p for p, k in C.PROBABILITY_BINS.items() if k == b]
            rows.append(dict(condition=cond, probability_bin=b, mean_w_empirical=np.mean([w[(cond, p)] for p in ps], axis=0).mean()))
    return pd.DataFrame(rows)


def main():
    b = pd.read_csv(C.ROI / "wp_bins.csv")
    prof = b.groupby(["roi", "condition", "probability_bin"]).value.agg(["mean", "sem", "count"]).reset_index()
    prof = prof.merge(bin_weights(participants(fmri_only=True)), on=["condition", "probability_bin"])
    save_table("wp_bin_profiles", prof)


if __name__ == "__main__":
    main()
