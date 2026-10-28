"""Shared helpers: data loading, statistics and result recording."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config as C  # noqa: E402


# ---------------------------------------------------------------- data loading
def participants(fmri_only=False):
    p = pd.read_csv(C.DATA / "participants.tsv", sep="\t")
    if fmri_only:
        p = p[p.in_fmri_sample == "yes"]
    return p.participant_id.tolist()


def lattimore_params():
    """Wide table: one row per participant, columns delta_immediate, delta_delayed, gamma_immediate, gamma_delayed."""
    t = pd.read_csv(C.BEHAVIOR / "lattimore_params.csv")
    w = t.pivot(index="participant_id", columns="condition", values=["delta", "gamma"])
    w.columns = [f"{a}_{b}" for a, b in w.columns]
    return w


def empirical_weights():
    """Wide table of w_empirical = CE / 20: columns (condition, probability)."""
    t = pd.read_csv(C.BEHAVIOR / "prescan_ce.csv")
    return t.pivot_table(index="participant_id", columns=["condition", "probability"], values="w_empirical")


def delta_w(probs, ids=None):
    """Mean delayed-minus-immediate empirical weight over the given probabilities."""
    w = empirical_weights()
    d = np.mean([w[("delayed", p)] - w[("immediate", p)] for p in probs], axis=0)
    s = pd.Series(d, index=w.index)
    return s if ids is None else s.loc[ids]


def lattimore(p, delta, gamma):
    """Two-parameter probability weighting function of Lattimore, Baker and Witte (1992)."""
    pg = np.power(p, gamma)
    return delta * pg / (delta * pg + np.power(1 - p, gamma))


# ---------------------------------------------------------------- statistics
def one_sample(x):
    x = np.asarray(x, float)
    n = len(x)
    res = stats.ttest_1samp(x, 0)
    sem = stats.sem(x)
    lo, hi = stats.t.interval(0.95, n - 1, loc=x.mean(), scale=sem)
    return dict(n=n, mean=x.mean(), sd=x.std(ddof=1), sem=sem, ci_low=lo, ci_high=hi,
                t=res.statistic, df=n - 1, p=res.pvalue, dz=x.mean() / x.std(ddof=1))


def paired(a, b):
    """Paired t-test of a minus b, with Cohen's d of the paired differences."""
    return one_sample(np.asarray(a, float) - np.asarray(b, float))


def correlation(x, y):
    """Pearson r with Fisher-z 95% CI, and Spearman rho."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    n = len(x)
    r = stats.pearsonr(x, y)
    rho = stats.spearmanr(x, y)
    z, h = np.arctanh(r.statistic), stats.norm.ppf(0.975) / np.sqrt(n - 3)
    return dict(n=n, r=r.statistic, r_ci_low=np.tanh(z - h), r_ci_high=np.tanh(z + h), p=r.pvalue,
                rho=rho.statistic, rho_p=rho.pvalue)


def leave_one_out_r(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    return np.array([stats.pearsonr(np.delete(x, i), np.delete(y, i)).statistic for i in range(len(x))])


# ---------------------------------------------------------------- recording
class Recorder:
    """Collects (stat_id, value) pairs and writes results/stats/<name>.csv."""

    def __init__(self, name):
        self.name, self.rows = name, []

    def add(self, stat_id, value, note=""):
        self.rows.append(dict(stat_id=stat_id, value=value, script=self.name, note=note))

    def add_many(self, prefix, d, keys):
        for k, key in keys.items():
            self.add(f"{prefix}_{k}" if k else prefix, d[key])

    def save(self):
        C.STATS_DIR.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(self.rows).to_csv(C.STATS_DIR / f"{self.name}.csv", index=False)


def save_source(panel, df):
    """Write the plotted data of one figure panel; run_all.py merges these into Source_Data.xlsx."""
    C.SOURCE_DIR.mkdir(parents=True, exist_ok=True)
    df.to_csv(C.SOURCE_DIR / f"{panel}.csv", index=False)


def save_table(name, df):
    out = C.RESULTS / "tables"
    out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / f"{name}.csv", index=False)
