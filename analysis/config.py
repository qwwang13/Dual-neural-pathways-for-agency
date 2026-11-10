from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
CACHE = ROOT / "analysis/.cache/derived"
BEHAVIOR = DATA / "05_LMM/behavior_single_trial.csv"
CONDITIONS = {"action": ["va_a", "ia_a"], "outcome": ["vs_s", "is_s"]}

def windows():
    parts = []
    columns = ["region", "primary_lower_ms", "primary_upper_ms"]
    for binding in CONDITIONS:
        one = pd.read_csv(RESULTS / "03_functional_regression" / f"{binding}_binding.csv")
        one = one[columns].drop_duplicates().assign(binding=binding)
        parts.append(one)
    w = pd.concat(parts, ignore_index=True)
    w["primary_window"] = (w.primary_lower_ms.astype(str) + "–" + w.primary_upper_ms.astype(str))
    w = w.sort_values(["binding", "primary_lower_ms", "region"]).reset_index(drop=True)
    w["roi"] = w.region.str.replace("-", "", regex=False)
    w["file"] = w.roi + w.binding.map({"action": "_action_ba.csv", "outcome": "_outcome_bs.csv"})
    if not ((w.primary_lower_ms >= 0) & (w.primary_upper_ms <= 1000)
            & (w.primary_lower_ms < w.primary_upper_ms)).all():
        raise ValueError("Invalid time limits")
    return w
