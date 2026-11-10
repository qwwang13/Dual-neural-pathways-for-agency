# Load configuration and shared functions
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import pandas as pd
from config import ROOT, DATA, CACHE, BEHAVIOR, windows

sys.path.insert(0, str(ROOT / "analysis/05_LMM"))
from TF_features import export_trial_features
from ITPC import export_angles_itpc

def export(stage):
    # Read selected windows and matched behavioural trials
    w = windows()
    for name in w.file:
        if not (DATA / "05_LMM" / name).is_file() and not (DATA / "05_LMM" / (name+".gz")).is_file():
            raise FileNotFoundError(DATA / "05_LMM" / name)
    part = "05_LMM/features" if stage == "features" else "05_LMM/ITPC"
    out = CACHE / part; out.mkdir(parents=True, exist_ok=True)
    records = []
    # Extract neural measures and save intermediate data
    for row in w.itertuples():
        source = DATA / "05_LMM" / row.file
        if not source.is_file():
            source = source.with_name(source.name+".gz")
        limits = (row.primary_lower_ms / 1000, row.primary_upper_ms / 1000)
        kwargs = dict(csv_path=source, time_window=limits, out_dir=out, sfreq=1200.,
                      freqs=list(range(1,31)), behavior_csv=BEHAVIOR, exclude_subjects=None)
        print(f"{stage}: {row.binding} {row.region} {limits}")
        if stage == "features":
            path, frame = export_trial_features(**kwargs)
        else:
            path, _, frame, means = export_angles_itpc(**kwargs)
        for cond, g in frame.groupby("condition"):
            records.append(dict(binding=row.binding, roi=row.roi, condition=cond,
                start_ms=row.primary_lower_ms, end_ms=row.primary_upper_ms,
                n_subjects=g.subject.nunique(), n_trials=len(g[["subject","trial"]].drop_duplicates()), file=str(path)))
    w.to_csv(out / "windows_used.csv", index=False)
    summary = pd.DataFrame(records)
    return summary

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=["features", "itpc"])
    export(parser.parse_args().stage)
