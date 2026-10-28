"""Figure 5C,D: GLM3 condition main effect in the meta-analytic mPFC mask vs curvature gamma, per condition, N = 22."""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import C, Recorder, correlation, lattimore_params  # noqa: E402


def data():
    b = pd.read_csv(C.ROI / "bartra_mpfc.csv").set_index("participant_id")
    g = lattimore_params().loc[b.index]
    return b.join(g)


def main():
    d = data()
    rec = Recorder("main_effect_gamma")
    keys = dict(r="r", ci_low="r_ci_low", ci_high="r_ci_high", p="p", rho="rho", rho_p="rho_p")
    rec.add_many("fig5c", correlation(d.glm3_immediate_main_effect, d.gamma_immediate), keys)
    rec.add_many("fig5d", correlation(d.glm3_delayed_main_effect, d.gamma_delayed), keys)
    rec.save()


if __name__ == "__main__":
    main()
