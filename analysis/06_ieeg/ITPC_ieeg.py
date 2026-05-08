import numpy as np
import pandas as pd
from pathlib import Path
from mne.time_frequency import tfr_array_morlet


def _wrap_to_2pi(theta):
    return np.mod(theta, 2 * np.pi)


def compute_trial_band_angles_and_subject_itpc_ieeg(
    csv_path: str | Path,
    time_window: tuple,
    sfreq: float = 1200.0,
    freqs=None,
    bands=None,
    encoding: str = "utf-8",
    exclude_subjects=None,
):
    if freqs is None:
        freqs = np.arange(1, 80, 1)
    freqs = np.asarray(freqs, dtype=float)
    n_cycles = freqs / 2.0

    if bands is None:
        bands = {
            "Delta": (1, 4),
            "Theta": (4, 8),
            "Alpha": (8, 13),
            "Beta": (13, 30),
            "Gamma": (30, 80),
        }

    df = pd.read_csv(csv_path, encoding=encoding)

    if "Time" in df.columns and "time" not in df.columns:
        df = df.rename(columns={"Time": "time"})
    if "Value" in df.columns and "value" not in df.columns:
        df = df.rename(columns={"Value": "value"})

    required = {"subject", "label", "trial", "condition", "time", "value"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns: {sorted(missing)}")

    if exclude_subjects:
        df = df[~df["subject"].isin(exclude_subjects)].copy()

    angles_rows = []

    for (subj, lab, cond), g in df.groupby(["subject", "label", "condition"], sort=False):
        trials = sorted(g["trial"].unique())
        if len(trials) == 0:
            continue

        g0 = g[g["trial"] == trials[0]].sort_values("time")
        t = g0["time"].to_numpy()

        tmask = (t >= time_window[0]) & (t <= time_window[1])
        if not np.any(tmask):
            raise ValueError(f"time_window has no samples: subj={subj}, label={lab}, cond={cond}")

        epochs = np.empty((len(trials), 1, t.size), dtype=float)
        for i, tr in enumerate(trials):
            gt = g[g["trial"] == tr].sort_values("time")
            if gt.shape[0] != t.size or not np.allclose(gt["time"].to_numpy(), t):
                raise ValueError(f"time axis mismatch: subj={subj}, label={lab}, cond={cond}, trial={tr}")
            epochs[i, 0, :] = gt["value"].to_numpy()

        cfs = tfr_array_morlet(
            epochs,
            sfreq=sfreq,
            freqs=freqs,
            n_cycles=n_cycles,
            output="complex",
            n_jobs=1,
        )

        amp = np.abs(cfs)
        unit = cfs / np.maximum(amp, np.finfo(float).eps)

        for band_name, (f_lo, f_hi) in bands.items():
            fmask = (freqs >= f_lo) & (freqs <= f_hi)
            if not np.any(fmask):
                continue

            unit_bt = unit[:, 0, :, :][:, fmask, :][:, :, tmask]
            z_trial = unit_bt.mean(axis=(1, 2))
            ang = _wrap_to_2pi(np.angle(z_trial))

            for tr, a in zip(trials, ang):
                angles_rows.append([subj, lab, cond, tr, band_name, float(a)])

    angles_df = pd.DataFrame(
        angles_rows,
        columns=["subject", "label", "condition", "trial", "band", "angle_rad"],
    )

    itpc_sub_rows = []
    for (subj, lab, cond, band), gg in angles_df.groupby(["subject", "label", "condition", "band"], sort=False):
        theta = gg["angle_rad"].to_numpy()
        if theta.size == 0:
            continue
        R = np.exp(1j * theta).mean()
        itpc_sub_rows.append([subj, lab, cond, band, float(np.abs(R))])

    itpc_sub_df = pd.DataFrame(
        itpc_sub_rows,
        columns=["subject", "label", "condition", "band", "itpc_sub"],
    )

    itpc_mean_df = (
        itpc_sub_df.groupby(["condition", "band", "label"], as_index=False)
        .agg(itpc_mean=("itpc_sub", "mean"))
    )

    return angles_df, itpc_mean_df


def export_angles_itpc_ieeg(
    csv_path: str | Path,
    time_window: tuple,
    out_dir: str | Path,
    sfreq: float = 1200.0,
    freqs=None,
    bands=None,
    encoding: str = "utf-8",
    exclude_subjects=None,
):
    csv_path = Path(csv_path)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    angles_df, itpc_mean_df = compute_trial_band_angles_and_subject_itpc_ieeg(
        csv_path=csv_path,
        time_window=time_window,
        sfreq=sfreq,
        freqs=freqs,
        bands=bands,
        encoding=encoding,
        exclude_subjects=exclude_subjects,
    )

    stem = csv_path.stem
    angles_csv = out_dir / f"{stem}_angles.csv"
    itpc_csv = out_dir / f"{stem}_itpc_mean.csv"

    angles_df.to_csv(angles_csv, index=False, encoding="utf-8-sig")
    itpc_mean_df.to_csv(itpc_csv, index=False, encoding="utf-8-sig")

    return angles_csv, itpc_csv, angles_df, itpc_mean_df