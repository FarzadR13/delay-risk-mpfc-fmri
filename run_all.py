"""Regenerate the statistics, figures and source-data sheets of this repository from data/.

    python run_all.py

Outputs: results/stats_reproduced.csv, results/check_against_manuscript.csv, results/figures/,
results/tables/ and source_data/Source_Data.xlsx. The only network access is the one-time download of
third-party atlases in step 0. Whole-brain statistics are read from the SPM tables in data/group_maps
(written by the MATLAB scripts in code/03_second_level).
"""
import runpy
import shutil
import sys
import time

import pandas as pd

import config as C

STEPS = ["00_masks/build_masks.py",
         "01_behavior/fit_lattimore.py", "01_behavior/parameter_tests.py", "01_behavior/weighting_decomposition.py",
         "04_roi_and_brain_behavior/loso_roi.py", "04_roi_and_brain_behavior/heatmap_gradients.py",
         "04_roi_and_brain_behavior/main_effect_gamma.py", "04_roi_and_brain_behavior/wp_bin_profiles.py",
         "04_roi_and_brain_behavior/lowp_brain_behavior.py", "04_roi_and_brain_behavior/whole_brain_summary.py",
         "06_figures/figure4.py", "06_figures/figure5.py", "06_figures/figure6.py", "06_figures/figure7.py"]


def run(script):
    t0 = time.time()
    sys.argv = [script]
    runpy.run_path(str(C.CODE / script), run_name="__main__")
    print(f"  {script} ({time.time() - t0:.1f} s)")


def merge_outputs():
    stats = pd.concat([pd.read_csv(f) for f in sorted(C.STATS_DIR.glob("*.csv"))], ignore_index=True)
    stats.to_csv(C.RESULTS / "stats_reproduced.csv", index=False)
    C.SOURCE_DATA_XLSX.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(C.SOURCE_DATA_XLSX, engine="openpyxl") as xl:
        pd.DataFrame(dict(sheet=[f.stem for f in sorted(C.SOURCE_DIR.glob("*.csv"))])).assign(
            content="data plotted in this figure panel").to_excel(xl, sheet_name="README", index=False)
        for f in sorted(C.SOURCE_DIR.glob("*.csv")):
            pd.read_csv(f).to_excel(xl, sheet_name=f.stem[:31], index=False)


def main():
    t0 = time.time()
    sys.path.insert(0, str(C.CODE))
    for d in (C.STATS_DIR, C.SOURCE_DIR, C.FIGURES, C.RESULTS / "tables"):
        shutil.rmtree(d, ignore_errors=True)
    for s in STEPS:
        run(s)
    merge_outputs()
    runpy.run_path(str(C.CODE / "07_check" / "check_against_manuscript.py"), run_name="__main__")
    print(f"done in {time.time() - t0:.0f} s")

if __name__ == "__main__":
    main()
