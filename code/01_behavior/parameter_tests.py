"""Immediate vs delayed Lattimore parameters: paired t-tests and effect sizes (behavioral sample, N = 24).

Differences are immediate minus delayed. Cohen's d is the mean paired difference divided by the standard
deviation of the differences, reported as a magnitude.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import Recorder, lattimore_params, paired  # noqa: E402


def main():
    w = lattimore_params()
    rec = Recorder("parameter_tests")
    for par in ("delta", "gamma"):
        d = paired(w[f"{par}_immediate"], w[f"{par}_delayed"])
        rec.add_many(f"beh_{par}", d, dict(diff="mean", diff_ci_low="ci_low", diff_ci_high="ci_high", t="t", p="p"))
        rec.add(f"beh_{par}_d", abs(d["dz"]))
    rec.add("meth_n_behavior", len(w))
    rec.save()


if __name__ == "__main__":
    main()
