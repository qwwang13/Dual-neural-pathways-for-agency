# Read data and define analysis settings
from pathlib import Path
import sys
import importlib
import numpy as np
import pandas as pd

ROOT = next(p for p in [Path.cwd(), *Path.cwd().parents]
            if (p / "dual-agency-pathways.Rproj").is_file())
HERE = ROOT / "analysis/06_ieeg"
sys.path.insert(0, str(HERE / "lib"))
import TF_features_ieeg
import ITPC_ieeg
importlib.reload(TF_features_ieeg)
importlib.reload(ITPC_ieeg)

DATA = ROOT / "data/06_ieeg"
OUT = ROOT / "analysis/.cache/derived/06_ieeg"
BEHAVIOR = DATA / "behavior_ieeg.csv"
SIGNAL = DATA / "ieeg_MTGR.csv"

sys.path.insert(0, str(ROOT / "analysis"))
from config import windows
w = windows()
w = w.loc[(w.binding == "outcome") & (w.region == "MTG-R")].copy()
TIME_WINDOW = tuple(w.iloc[0][["primary_lower_ms", "primary_upper_ms"]].astype(float) / 1000)
print("MTG-R window (seconds):", TIME_WINDOW)

SFREQ = 2000.0
FREQS = np.arange(1, 31, dtype=float)
N_CYCLES = FREQS / 2.0
TFR_TIME_WINDOW = (-0.2, 1.0)
BANDS = {"Delta": (1, 4), "Theta": (4, 8),
         "Alpha": (8, 13), "Beta": (13, 30)}

# Extract and save neural features
kwargs = dict(csv_path=SIGNAL, time_window=TIME_WINDOW, sfreq=SFREQ,
              freqs=FREQS, bands=BANDS, n_cycles=N_CYCLES,
              tfr_time_window=TFR_TIME_WINDOW, behavior_csv=BEHAVIOR)
feature_csv, features = TF_features_ieeg.export_trial_features_ieeg(
    **kwargs, out_dir=OUT / "feature")
w.to_csv(OUT / "windows_used.csv", index=False)
counts = features.groupby(["subject", "label", "condition"]).size().reset_index(name="n_trials")
counts

# Compute and save phase angles and ITPC
angles_csv, itpc_csv, angles, itpc_mean = ITPC_ieeg.export_angles_itpc_ieeg(
    **kwargs, out_dir=OUT / "ITPC")
itpc_mean
