import numpy as np
import pandas as pd
from pathlib import Path
from mne.time_frequency import tfr_array_morlet


def _wrap_to_2pi(theta):
    return np.mod(theta, 2 * np.pi)


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


def compute_trial_band_angles_and_subject_itpc(
    csv_path: str | Path,
    time_window: tuple,
    sfreq: float = 1200.0,
    freqs=None,
    bands=None,
    encoding: str = "utf-8",
    exclude_subjects=None,
    behavior_csv: str | Path | None = None,
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

    angles_rows = []
    itpc_sub_rows = []

    for (subj, cond), g in df.groupby(["subject", "condition"], sort=False):
        trials = sorted(g["trial"].unique())
        if len(trials) == 0:
            continue

        g0 = g[g["trial"] == trials[0]].sort_values("time")
        t = g0["time"].to_numpy()

        tmask = (t >= time_window[0]) & (t <= time_window[1])
        if not np.any(tmask):
            raise ValueError(f"time_window has no samples: subj={subj}, cond={cond}")

        epochs = np.empty((len(trials), 1, t.size), dtype=float)
        for i, tr in enumerate(trials):
            gt = g[g["trial"] == tr].sort_values("time")
            if gt.shape[0] != t.size or not np.allclose(gt["time"].to_numpy(), t):
                raise ValueError(f"time axis mismatch: subj={subj}, cond={cond}, trial={tr}")
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
        itpc_tf = np.abs(unit.mean(axis=0))

        for band_name, (f_lo, f_hi) in bands.items():
            fmask = (freqs >= f_lo) & (freqs < f_hi)
            if not np.any(fmask):
                continue

            unit_bt = unit[:, 0, :, :][:, fmask, :][:, :, tmask]
            z_trial = unit_bt.mean(axis=(1, 2))
            ang = _wrap_to_2pi(np.angle(z_trial))
            band_itpc = itpc_tf[0, :, :][np.ix_(fmask, tmask)].mean()
            angle_vector = np.exp(1j * ang).mean()

            itpc_sub_rows.append(
                [
                    subj,
                    cond,
                    band_name,
                    float(band_itpc),
                    float(angle_vector.real),
                    float(angle_vector.imag),
                    float(np.abs(angle_vector)),
                    float(_wrap_to_2pi(np.angle(angle_vector))),
                ]
            )

            for tr, a in zip(trials, ang):
                angles_rows.append([subj, cond, tr, band_name, float(a)])

    angles_df = pd.DataFrame(
        angles_rows, columns=["subject", "condition", "trial", "band", "angle_rad"]
    )

    itpc_sub_df = pd.DataFrame(
        itpc_sub_rows,
        columns=[
            "subject",
            "condition",
            "band",
            "bandITPC",
            "angle_vector_real",
            "angle_vector_imag",
            "angle_vector_mean",
            "angle_mean_rad",
        ],
    )

    itpc_mean_df = (
        itpc_sub_df.groupby(["condition", "band"], as_index=False)
        .agg(
            n_subjects=("subject", "nunique"),
            itpc_mean=("bandITPC", "mean"),
            angle_vector_real=("angle_vector_real", "mean"),
            angle_vector_imag=("angle_vector_imag", "mean"),
        )
    )
    itpc_mean_df["angle_vector_mean"] = np.hypot(
        itpc_mean_df["angle_vector_real"], itpc_mean_df["angle_vector_imag"]
    )
    itpc_mean_df["angle_mean_rad"] = _wrap_to_2pi(
        np.arctan2(
            itpc_mean_df["angle_vector_imag"],
            itpc_mean_df["angle_vector_real"],
        )
    )
    itpc_mean_df = itpc_mean_df[
        [
            "condition",
            "band",
            "n_subjects",
            "itpc_mean",
            "angle_mean_rad",
            "angle_vector_mean",
            "angle_vector_real",
            "angle_vector_imag",
        ]
    ]

    return angles_df, itpc_mean_df


def export_angles_itpc(
    csv_path: str | Path,
    time_window: tuple,
    out_dir: str | Path,
    sfreq: float = 1200.0,
    freqs=None,
    bands=None,
    encoding: str = "utf-8",
    exclude_subjects=None,
    behavior_csv: str | Path | None = None,
):
    csv_path = Path(csv_path)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    angles_df, itpc_mean_df = compute_trial_band_angles_and_subject_itpc(
        csv_path=csv_path,
        time_window=time_window,
        sfreq=sfreq,
        freqs=freqs,
        bands=bands,
        encoding=encoding,
        exclude_subjects=exclude_subjects,
        behavior_csv=behavior_csv,
    )

    stem = csv_path.stem
    angles_csv = out_dir / f"{stem}_angles.csv"
    itpc_csv = out_dir / f"{stem}_itpc_mean.csv"

    angles_df.to_csv(angles_csv, index=False, encoding="utf-8-sig")
    itpc_mean_df.to_csv(itpc_csv, index=False, encoding="utf-8-sig")

    return angles_csv, itpc_csv, angles_df, itpc_mean_df
