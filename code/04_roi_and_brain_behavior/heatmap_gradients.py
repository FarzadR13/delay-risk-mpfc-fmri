"""Figure 5A,B: 3 x 3 probability-by-distance cell means in the meta-analytic mPFC mask and distance gradients.

Gradient = highest minus lowest distance bin, averaged over the three probability bins, tested against
zero per condition (two-tailed one-sample t-tests), N = 22.
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import C, Recorder, one_sample, paired, save_table  # noqa: E402


def gradients(cells):
    out = {}
    for cond in C.CONDITIONS:
        w = cells[cells.condition == cond].pivot_table(index="participant_id", columns=["probability_bin", "distance_bin"], values="value")
        out[cond] = sum(w[(pb, 3)] - w[(pb, 1)] for pb in (1, 2, 3)) / 3
    return pd.DataFrame(out)


def main():
    cells = pd.read_csv(C.ROI / "heatmap_cells.csv")
    g = gradients(cells)
    rec = Recorder("heatmap_gradients")
    rec.add_many("heat_imm", one_sample(g.immediate), dict(grad="mean", ci_low="ci_low", ci_high="ci_high", t="t", p="p"))
    rec.add_many("heat_del", one_sample(g.delayed), dict(grad="mean", t="t", p="p"))
    rec.add_many("heat_diff", paired(g.immediate, g.delayed), dict(mean="mean", t="t", p="p"))
    rec.save()
    means = cells.groupby(["condition", "probability_bin", "distance_bin"]).value.agg(["mean", "sem", "count"]).reset_index()
    save_table("heatmap_cell_means", means)


if __name__ == "__main__":
    main()
