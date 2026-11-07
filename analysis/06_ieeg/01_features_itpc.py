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
assert len(w) == 1, "Run 03_functional.Rmd first; Outcome MTG-R must have one window."
TIME_WINDOW = tuple(w.iloc[0][["primary_lower_ms", "primary_upper_ms"]].astype(float) / 1000)
print("MTG-R window (seconds):", TIME_WINDOW)

SFREQ = 2000.0
FREQS = np.arange(1, 31, dtype=float)
N_CYCLES = FREQS / 2.0
TFR_TIME_WINDOW = (-0.2, 1.0)
BANDS = {"Delta": (1, 4), "Theta": (4, 8),
         "Alpha": (8, 13), "Beta": (13, 30)}
CHANNELS = {2: "chan83-chan84", 4: "chan10-chan11"}
CONDITIONS = ["Voluntary_OutcomeOnsetReportSound",
              "Involuntary_OutcomeOnsetReportSound"]

assert TFR_TIME_WINDOW[0] <= TIME_WINDOW[0] < TIME_WINDOW[1] <= TFR_TIME_WINDOW[1]
keys = ["subject", "condition", "trial"]
behavior = pd.read_csv(BEHAVIOR)
assert not behavior.duplicated(keys).any(), "Duplicate behavior trial keys."
behavior = behavior.loc[np.isfinite(behavior.behavior_value)].copy()
neural = pd.read_csv(SIGNAL, usecols=keys + ["label"]).drop_duplicates()
assert set(zip(neural.subject, neural.label)) == set(CHANNELS.items()), "Unexpected patient/channel."
assert set(neural.condition) == set(CONDITIONS), "Unexpected conditions."
assert not neural.duplicated(keys).any(), "More than one channel per patient/trial."
matched = neural.merge(behavior[keys], on=keys, how="inner", validate="one_to_one")
expected_groups = {(s, lab, c) for s, lab in CHANNELS.items() for c in CONDITIONS}
assert set(zip(matched.subject, matched.label, matched.condition)) == expected_groups

# Extract and save neural features
kwargs = dict(csv_path=SIGNAL, time_window=TIME_WINDOW, sfreq=SFREQ,
              freqs=FREQS, bands=BANDS, n_cycles=N_CYCLES,
              tfr_time_window=TFR_TIME_WINDOW, behavior_csv=BEHAVIOR)
feature_csv, features = TF_features_ieeg.export_trial_features_ieeg(
    **kwargs, out_dir=OUT / "feature")
assert not features.duplicated(keys).any()
assert np.isfinite(features.select_dtypes(include="number").to_numpy()).all()
assert set(map(tuple, features[keys].to_numpy())) == set(map(tuple, matched[keys].to_numpy()))
w.to_csv(OUT / "windows_used.csv", index=False)
counts = features.groupby(["subject", "label", "condition"]).size().reset_index(name="n_trials")
counts

# Compute and save phase angles and ITPC
angles_csv, itpc_csv, angles, itpc_mean = ITPC_ieeg.export_angles_itpc_ieeg(
    **kwargs, out_dir=OUT / "ITPC")
assert set(map(tuple, angles[keys].drop_duplicates().to_numpy())) == set(map(tuple, matched[keys].to_numpy()))
itpc_check = features.groupby(["label", "condition"])[
    [f"{b}_bandITPC" for b in BANDS]].first().reset_index().melt(
    id_vars=["label", "condition"], var_name="band", value_name="feature_ITPC")
itpc_check["band"] = itpc_check.band.str.replace("_bandITPC", "", regex=False)
check = itpc_check.merge(itpc_mean, on=["label", "condition", "band"], validate="one_to_one")
assert len(check) == len(itpc_check)
assert np.allclose(check.feature_ITPC, check.itpc_mean)
itpc_mean
