"""Figure 4B: leave-one-subject-out mPFC ROI, immediate vs delayed signed-distance slopes (GLM1), N = 22."""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import C, Recorder, one_sample, paired, save_source  # noqa: E402


def main():
    t = pd.read_csv(C.ROI / "loso_mpfc.csv")
    rec = Recorder("loso_roi")
    keys = dict(mean="mean", ci_low="ci_low", ci_high="ci_high", t="t", p="p", dz="dz")
    rec.add_many("loso_imm", one_sample(t.immediate_slope), keys)
    rec.add_many("loso_del", one_sample(t.delayed_slope), keys)
    rec.add_many("loso_diff", paired(t.immediate_slope, t.delayed_slope), keys)
    rec.save()
    save_source("Fig4B_LOSO_values", t[["participant_id", "peak_x", "peak_y", "peak_z", "n_voxels", "immediate_slope", "delayed_slope"]])


if __name__ == "__main__":
    main()
