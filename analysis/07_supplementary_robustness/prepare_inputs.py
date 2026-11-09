"""Recover participant inputs and the original bootstrap draws from processed data.

This optional preparation step requires NumPy 2.3.5. The statistical analysis
in 07_robustness.R uses the included CSV inputs and requires only base R.
"""

import argparse
import csv
import hashlib
from pathlib import Path

import numpy as np


SPECS = [
    ("action", "PreSMAL", "PreSMA-L", "Theta", 220, 270),
    ("outcome", "PrecentralR", "PreCG-R", "Alpha", 180, 240),
    ("outcome", "MTGR", "MTG-R", "Delta", 780, 830),
    ("outcome", "MTGR", "MTG-R", "Theta", 780, 830),
]
CONDITIONS = {
    "action": [("va_a", "Voluntary"), ("ia_a", "Involuntary")],
    "outcome": [("vs_s", "Voluntary"), ("is_s", "Involuntary")],
}


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path, columns, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(columns)
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path,
                        default=Path(__file__).resolve().parents[2],
                        help="Root of the processed dual-agency-pathways package")
    parser.add_argument("--output-dir", type=Path,
                        default=None)
    args = parser.parse_args()
    root = args.source_root.resolve()
    out = (args.output_dir or root / "data/07_supplementary_robustness").resolve()
    out.mkdir(parents=True, exist_ok=True)
    behavior = root / "data/01_temporal_binding_effect/perceptual_shift_action_outcome.csv"
    feature_root = root / "analysis/.cache/derived/05_LMM/features"
    windows_file = feature_root / "windows_used.csv"
    windows = read_csv(windows_file)
    behavior_rows = read_csv(behavior)
    if sorted(int(r["subject"]) for r in behavior_rows) != list(range(2, 21)):
        raise ValueError("Expected the 19 participants numbered 2 through 20.")
    sources = {behavior, windows_file}
    records = []
    for binding, roi_key, label, band, lower, upper in SPECS:
        window = [r for r in windows if r["binding"] == binding and r["roi"] == roi_key]
        if len(window) != 1 or (float(window[0]["primary_lower_ms"]),
                                float(window[0]["primary_upper_ms"])) != (lower, upper):
            raise ValueError(f"Unexpected analysis window: {roi_key}")
        suffix = "action_ba" if binding == "action" else "outcome_bs"
        source = feature_root / f"{roi_key}_{suffix}.csv"
        sources.add(source)
        rows = read_csv(source)
        trial_keys = [(r["subject"], r["condition"], r["trial"]) for r in rows]
        if len(set(trial_keys)) != len(trial_keys):
            raise ValueError(f"Duplicate trial keys: {source}")
        for condition_key, condition_label in CONDITIONS[binding]:
            values = {}
            for row in rows:
                if row["condition"] != condition_key:
                    continue
                subject = int(row["subject"])
                value = float(row[f"{band}_bandITPC"])
                if not np.isfinite(value) or not 0 <= value <= 1:
                    raise ValueError(f"Invalid ITPC: {source}")
                if subject in values and values[subject] != value:
                    raise ValueError("ITPC must be constant within participant and condition.")
                values[subject] = value
            if sorted(values) != list(range(2, 21)):
                raise ValueError("Every selected feature and condition must contain 19 participants.")
            for subject in range(2, 21):
                records.append([subject, binding, label, band, condition_label,
                                values[subject], lower, upper])
    write_csv(out / "participant_itpc.csv",
              ["subject", "binding", "roi", "band", "condition", "itpc",
               "lower_ms", "upper_ms"], records)
    # These are the same PCG64 draws used in the original participant analysis.
    # Store one-based row positions for direct use in R; row 1 is participant 2.
    draws = np.random.Generator(np.random.PCG64(20261109)).integers(
        0, 19, size=(20_000, 19), dtype=np.int64)
    np.savetxt(out / "bootstrap_indices.csv", draws + 1, fmt="%d", delimiter=",")
    write_csv(out / "source_files.csv", ["source_file", "sha256"],
              [[p.relative_to(root).as_posix(), hashlib.sha256(p.read_bytes()).hexdigest()]
               for p in sorted(sources)])
    # Canonical text hashes remain valid for LF or CRLF Git checkouts.
    local_inputs = [
        ("../01_temporal_binding_effect/perceptual_shift_action_outcome.csv", behavior),
        ("participant_itpc.csv", out / "participant_itpc.csv"),
        ("bootstrap_indices.csv", out / "bootstrap_indices.csv"),
    ]
    write_csv(out / "input_checksums.csv", ["file", "md5"],
              [[name, hashlib.md5(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()]
               for name, path in local_inputs])
    print(f"Prepared 19 behavioural records, {len(records)} ITPC records and 20,000 paired draws.")
    print(f"Data: {out}")


if __name__ == "__main__":
    main()
