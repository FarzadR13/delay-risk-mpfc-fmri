"""Compare every statistic listed in manuscript_ledger.csv with the reproduced value.

A numeric value matches when the reproduced value, rounded to the printed number of decimals, equals the
printed value. "< x" and "> x" are checked as inequalities, MNI coordinates as strings.
"""
import re
import sys
from decimal import ROUND_HALF_UP, Decimal
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import C  # noqa: E402

HERE = Path(__file__).resolve().parent


def decimals(s):
    m = re.search(r"\.(\d+)", s)
    return len(m.group(1)) if m else 0


def rounded(x, nd):
    return float(Decimal(repr(float(x))).quantize(Decimal(1).scaleb(-nd), rounding=ROUND_HALF_UP))


def compare(printed, value):
    if value is None:
        return "MISSING"
    p = printed.strip()
    ineq = re.match(r"^([<>])\s*(-?[\d.]+)", p)
    if ineq:
        lim = float(ineq.group(2))
        return "MATCH" if (float(value) < lim if ineq.group(1) == "<" else float(value) > lim) else "MISMATCH"
    if p.startswith("["):
        return "MATCH" if p == str(value) else "MISMATCH"
    if re.fullmatch(r"-?[\d.]+", p):
        return "MATCH" if rounded(value, decimals(p)) == float(p) else "MISMATCH"
    if p == "0 (no overlap)":
        return "MATCH" if float(value) == 0 else "MISMATCH"
    if p == "none":
        return "MATCH" if value == "none" else "MISMATCH"
    if p == "significant (p < 0.05)":
        return "MATCH" if float(value) < 0.05 else "MISMATCH"
    if p == "qualitative":
        return "MATCH" if value == "consistent" else "MISMATCH"
    return "UNPARSED"


def main():
    ledger = pd.read_csv(HERE / "manuscript_ledger.csv", dtype=str)
    rep = pd.read_csv(C.RESULTS / "stats_reproduced.csv", dtype={"value": str}).drop_duplicates("stat_id", keep="last")
    vals = dict(zip(rep.stat_id, rep.value))
    rows = []
    for _, r in ledger.iterrows():
        v = vals.get(r.stat_id)
        try:
            vv = float(v)
        except (TypeError, ValueError):
            vv = v
        rows.append(dict(stat_id=r.stat_id, figure_or_section=r.figure_or_section, quantity=r.quantity,
                         manuscript_value=r.value_as_printed, reproduced_value=v, status=compare(r.value_as_printed, vv)))
    out = pd.DataFrame(rows)
    out.to_csv(C.RESULTS / "check_against_manuscript.csv", index=False)
    print(out.status.value_counts().to_string())
    return out


if __name__ == "__main__":
    main()
