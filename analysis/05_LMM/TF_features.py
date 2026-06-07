# TF_features.py
import numpy as np
import pandas as pd
from pathlib import Path
from mne.time_frequency import tfr_array_morlet


def _filter_to_behavior_trials(df, behavior_csv, encoding):
    if behavior_csv is None:
        return df

    key_cols = ["subject", "condition", "trial"]
    beh = pd.read_csv(behavior_csv, encoding=encoding)
    missing = set(key_cols) - set(beh.columns)
    if missing:
        raise ValueError(f"behavior_csv missing columns: {sorted(missing)}")

    if "behavior_value" in beh.columns:
        beh = beh[beh["behavior_value"].notna()].copy()

    key_df = beh[key_cols].drop_duplicates().copy()
    df2 = df.copy()

    tmp_cols = [f"__key_{col}" for col in key_cols]
    for col, tmp in zip(key_cols, tmp_cols):
        key_df[tmp] = key_df[col].astype(str)
        df2[tmp] = df2[col].astype(str)

    filtered = df2.merge(key_df[tmp_cols].drop_duplicates(), on=tmp_cols, how="inner")
    filtered = filtered.drop(columns=tmp_cols)
    if filtered.empty:
        raise ValueError("No neural trials match behavior_csv after filtering")

    return filtered


def compute_trial_features(
    csv_path: str | Path,
    time_window: tuple,
    sfreq: float = 1200.0,
    freqs=None,
    bands=None,
    encoding: str = "utf-8",
    behavior_csv: str | Path | None = None,
    exclude_subjects=None,
):
    if freqs is None:
        freqs = np.arange(1, 31, 1)
    freqs = np.asarray(freqs, dtype=float)
    n_cycles = freqs / 2.0

    if bands is None:
        bands = {
            "Delta": (1, 4),
            "Theta": (4, 8),
            "Alpha": (8, 13),
            "Beta": (13, 30),
        }

    df = pd.read_csv(csv_path, encoding=encoding)

    if "Time" in df.columns and "time" not in df.columns:
        df = df.rename(columns={"Time": "time"})
    if "Value" in df.columns and "value" not in df.columns:
        df = df.rename(columns={"Value": "value"})

    required = {"subject", "trial", "condition", "time", "value"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns: {sorted(missing)}")

    if exclude_subjects:
        exclude_set = {str(x) for x in exclude_subjects}
        df = df[~df["subject"].astype(str).isin(exclude_set)].copy()

    df = _filter_to_behavior_trials(df, behavior_csv, encoding)

    trial_power_rows = []
    itpc_rows = []

    for (subj, cond), g in df.groupby(["subject", "condition"], sort=False):
        trials = sorted(g["trial"].unique())
        if len(trials) == 0:
            continue

        g0 = g[g["trial"] == trials[0]].sort_values("time")
        t = g0["time"].to_numpy()

        epochs = np.empty((len(trials), 1, t.size), dtype=float)
        for i, tr in enumerate(trials):
            gt = g[g["trial"] == tr].sort_values("time")
            if gt.shape[0] != t.size or not np.allclose(gt["time"].to_numpy(), t):
                raise ValueError(f"time axis mismatch: subj={subj}, cond={cond}, trial={tr}")
            epochs[i, 0, :] = gt["value"].to_numpy()

        tmask = (t >= time_window[0]) & (t <= time_window[1])
        if not np.any(tmask):
            raise ValueError(f"time_window has no samples: subj={subj}, cond={cond}")

        amp_per_trial = epochs[:, 0, :][:, tmask].mean(axis=1)

        cfs = tfr_array_morlet(
            epochs,
            sfreq=sfreq,
            freqs=freqs,
            n_cycles=n_cycles,
            output="complex",
            n_jobs=1,
        )

        power = np.abs(cfs) ** 2
        phase = cfs / np.maximum(np.abs(cfs), np.finfo(float).eps)
        itpc = np.abs(phase.mean(axis=0))

        for band_name, (f_lo, f_hi) in bands.items():
            fmask = (freqs >= f_lo) & (freqs < f_hi)
            if not np.any(fmask):
                continue

            band_power_per_trial = power[:, 0, :, :][:, fmask, :][:, :, tmask].mean(axis=(1, 2))
            for tr, bp, amp in zip(trials, band_power_per_trial, amp_per_trial):
                trial_power_rows.append([subj, cond, tr, band_name, float(bp), float(amp)])

            band_itpc = itpc[0, :, :][np.ix_(fmask, tmask)].mean()
            itpc_rows.append([subj, cond, band_name, float(band_itpc)])

    TrialPowerLong = pd.DataFrame(
        trial_power_rows,
        columns=["subject", "condition", "trial", "band", "bandPower", "meanAmplitude"],
    )

    TrialPowerWide = (
        TrialPowerLong.pivot_table(
            index=["subject", "condition", "trial", "meanAmplitude"],
            columns="band",
            values="bandPower",
            aggfunc="mean",
        )
        .rename(columns=lambda b: f"{b}_bandPower")
        .reset_index()
    )

    ITPCSummaryLong = pd.DataFrame(
        itpc_rows,
        columns=["subject", "condition", "band", "bandITPC"],
    )

    ITPCWide = (
        ITPCSummaryLong.pivot_table(
            index=["subject", "condition"],
            columns="band",
            values="bandITPC",
            aggfunc="mean",
        )
        .rename(columns=lambda b: f"{b}_bandITPC")
        .reset_index()
    )

    merged = TrialPowerWide.merge(ITPCWide, on=["subject", "condition"], how="left")
    return merged


def export_trial_features(
    csv_path: str | Path,
    time_window: tuple,
    out_dir: str | Path,
    sfreq: float = 1200.0,
    freqs=None,
    bands=None,
    encoding: str = "utf-8",
    behavior_csv: str | Path | None = None,
    exclude_subjects=None,
):
    csv_path = Path(csv_path)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    merged = compute_trial_features(
        csv_path=csv_path,
        time_window=time_window,
        sfreq=sfreq,
        freqs=freqs,
        bands=bands,
        encoding=encoding,
        behavior_csv=behavior_csv,
        exclude_subjects=exclude_subjects,
    )

    out_csv = out_dir / csv_path.name
    merged.to_csv(out_csv, index=False, encoding="utf-8-sig")
    return out_csv, merged
