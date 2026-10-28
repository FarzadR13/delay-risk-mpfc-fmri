"""Paths and constants shared by every analysis script. All paths are relative to the repository root."""
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
BEHAVIOR = DATA / "behavior"
ROI = DATA / "roi"
GROUP_MAPS = DATA / "group_maps"
OWN_MASKS = DATA / "masks"
CODE = ROOT / "code"
RESULTS = ROOT / "results"
STATS_DIR = RESULTS / "stats"            # one CSV per analysis script
FIGURES = RESULTS / "figures"
SOURCE_DIR = RESULTS / "source_data"     # one CSV per figure panel, merged into Source_Data.xlsx
DERIVED_MASKS = RESULTS / "masks"        # masks rebuilt from third-party atlases (not committed)
CACHE = ROOT / ".cache"                  # downloaded third-party files (not committed)
SOURCE_DATA_XLSX = ROOT / "source_data" / "Source_Data.xlsx"

# Bartra, McGuire and Kable (2013), Figure 3 masks, distributed by the Kable lab
BARTRA_URL = "https://www.kablelab.com/uploads/1/3/0/4/130439175/nifti_roi-masks_fig03.zip"
BARTRA_SHA256 = "d4e41bf2a93109aa5c75b021635e7448d7321b5c91b521e5150dd898f8c4a1eb"
BARTRA_FILE = "fig03A/I_POS_kda10_testStat_mean_cMass_testStat_bin.nii.gz"

PROBABILITIES = [0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95]
DISTANCE_LEVELS = [-30, -20, -10, 0, 10, 20, 30]           # percent of the certainty equivalent
PROBABILITY_BINS = {0.05: 1, 0.10: 1, 0.25: 2, 0.50: 2, 0.75: 2, 0.90: 3, 0.95: 3}
DISTANCE_BINS = {-30: 1, -20: 1, -10: 2, 0: 2, 10: 2, 20: 3, 30: 3}
CONDITIONS = ["immediate", "delayed"]

DF_GROUP = 21                       # one-sample t-tests, N = 22
CLUSTER_P = 0.001                   # cluster-forming threshold
DISPLAY_P = 0.005                   # display threshold for brain figures

# Peaks (MNI mm) used to identify clusters and spheres
PEAK_FIG4A_MPFC = (6, 50, 14)
PEAK_FIG6A = {"mpfc": (-9, 59, 26), "pcc": (-6, -52, 32), "left_ag": (-48, -64, 26)}
LOSO_RADIUS_MM = 8.0

