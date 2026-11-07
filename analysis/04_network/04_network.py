"""Action and outcome PCMCI-GPDC networks."""
# Read data and configure paths
import json
import hashlib
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize, LinearSegmentedColormap
from matplotlib.patches import Ellipse
import network_support as base

ROOT = Path(__file__).resolve().parents[2]
BIN_MS = 10

def _folder(run_name, figures=False):
    return (ROOT / "figures" if figures else ROOT / "analysis/.cache/derived") / "04_network"


def prepare_input(binding, window_ms, rois):
    """Match the original preprocessing: bin within subject, then average 19 subjects."""
    start, end = window_ms
    if not (0 <= start < end <= 1000) or start % BIN_MS or end % BIN_MS:
        raise ValueError(f"{binding}: choose 10-ms boundaries within 0-1000 ms.")
    rois = [r.replace("-", "") for r in rois]
    if len(rois) < 2 or len(set(rois)) != len(rois):
        raise ValueError(f"{binding}: select at least two distinct ROIs.")
    source = ROOT / "data/03_functional_regression" / f"meg_{binding}_difference.csv"
    d = pd.read_csv(source)
    missing = set(rois) - set(d.region)
    if missing:
        raise ValueError(f"{binding}: unavailable ROIs {sorted(missing)}; available: {sorted(d.region.unique())}")
    condition = {"action": "vaa_iaa", "outcome": "vss_iss"}[binding]
    d = d.loc[d.region.isin(rois) & (d.condition == condition) & d.Time.between(0, 1000)].copy()
    if not np.isfinite(d[["Time", "Value", "subject"]].to_numpy()).all():
        raise ValueError(f"{binding}: nonfinite input.")
    if d.duplicated(["region", "subject", "Time"]).any():
        raise ValueError(f"{binding}: duplicate subject/ROI/time rows.")
    edges = np.arange(0, 1000 + BIN_MS, BIN_MS)
    bins = np.searchsorted(edges, d.Time.to_numpy(), side="right") - 1
    d["lower"] = np.minimum(bins, len(edges)-2) * BIN_MS
    d = d.loc[(d.lower >= start) & (d.lower < end)]
    subject_means = d.groupby(["region", "subject", "lower"], as_index=False).Value.mean()
    for region in rois:
        one = subject_means.loc[subject_means.region == region]
        if set(one.subject) != set(range(2, 21)):
            raise ValueError(f"{binding}/{region}: expected all subjects 2-20.")
        if len(one) != 19 * ((end-start)//BIN_MS) or not one.groupby("lower").subject.nunique().eq(19).all():
            raise ValueError(f"{binding}/{region}: incomplete subject/bin grid.")
    matrix = subject_means.groupby(["lower", "region"]).Value.mean().unstack("region")
    matrix = matrix.reindex(index=np.arange(start, end, BIN_MS), columns=rois)
    if not np.isfinite(matrix.to_numpy()).all():
        raise ValueError(f"{binding}: missing bins.")
    return matrix

# Estimate networks and save link results
def calculate(networks, tau_max_values, *, run_name, tau_min=1, pc_alpha=.05, alpha_level=.05):
    if set(networks) != {"action", "outcome"}:
        raise ValueError("Provide one action and one outcome configuration.")
    taus = list(tau_max_values)
    if not taus or len(set(taus)) != len(taus) or any(int(t) != t for t in taus):
        raise ValueError("tau_max_values must contain distinct integers.")
    taus = [int(t) for t in taus]
    if int(tau_min) != tau_min or not 1 <= tau_min <= min(taus):
        raise ValueError("Require 1 <= tau_min <= every tau_max.")
    if not (0 < pc_alpha < 1 and 0 < alpha_level < 1):
        raise ValueError("Significance thresholds must lie between 0 and 1.")
    matrices = {b: prepare_input(b, **networks[b]) for b in ("action", "outcome")}
    for b, matrix in matrices.items():
        if len(matrix) <= 2 * max(taus):
            raise ValueError(f"{b}: {len(matrix)} bins are insufficient for tau_max={max(taus)} (PCMCI cut-off 2*tau_max).")
    out = _folder(run_name)
    out.mkdir(parents=True, exist_ok=True)
    settings = dict(networks={b: dict(window_ms=list(networks[b]["window_ms"]), rois=list(matrices[b].columns))
                              for b in matrices}, tau_min=int(tau_min), tau_max_values=taus,
                    bin_ms=BIN_MS, n_subjects=19, pc_alpha=pc_alpha, alpha_level=alpha_level,
                    significance="analytic", multiple_testing="none")
    settings["settings_key"] = hashlib.sha256(json.dumps(settings,sort_keys=True).encode()).hexdigest()
    summary, all_links = [], []
    for b, matrix in matrices.items():
        regions = list(matrix.columns)
        for maximum in taus:
            print(f"{b}: {networks[b]['window_ms']} ms, tau={tau_min}..{maximum}")
            pcmci = base.PCMCI(dataframe=base.pp.DataFrame(matrix.to_numpy(), var_names=regions),
                               cond_ind_test=base.GPDC(significance="analytic"), verbosity=0)
            result = pcmci.run_pcmci(tau_min=int(tau_min), tau_max=maximum,
                                     pc_alpha=pc_alpha, alpha_level=alpha_level)
            np.savez_compressed(out / f"{b}_tau_max_{maximum}.npz", settings_key=settings["settings_key"], **result)
            links = base.results_to_long_df(result, regions, tau_min=int(tau_min), tau_max=maximum)
            links["binding"] = b
            links["tau_max"] = maximum
            links["lag_ms"] = links.tau * BIN_MS
            links["significant"] = links.p_value <= alpha_level
            all_links.append(links)
            summary.append(dict(binding=b, start_ms=networks[b]["window_ms"][0],
                end_ms=networks[b]["window_ms"][1], tau_min=int(tau_min), tau_max=maximum,
                n_cross_links=int(((links.source != links.target) & links.significant).sum())))
    links = pd.concat(all_links, ignore_index=True)
    result_dir = ROOT / "results/04_network"
    result_dir.mkdir(parents=True, exist_ok=True)
    links.to_csv(result_dir / "all_links.csv", index=False)
    summary = pd.DataFrame(summary)
    (out / "settings.json").write_text(json.dumps(settings, indent=2), encoding="utf-8")
    return summary

# Plot and save network figures
def _positions(n):
    if n in base.NODE_POSITIONS:
        return base.NODE_POSITIONS[n]
    angles = np.pi/2 - np.arange(n) * 2*np.pi/n
    return {"x": .65*np.cos(angles), "y": .55*np.sin(angles)}

def _label(roi):
    names = {"PreSMAL": "PreSMA-L", "PreSMAR": "PreSMA-R", "SMAR": "SMA-R",
             "InsulaR": "Insula-R", "AGR": "AG-R", "DLPFCR": "DLPFC-R",
             "MTGR": "MTG-R", "SPGL": "SPG-L", "PrecuneusR": "Precuneus\nR",
             "PrecentralR": "PreCG-R", "PostcentralL": "PoCG-L",
             "PostcentralR": "PoCG-R"}
    return names.get(roi, roi)

def plot_combined(*, run_name, tau_max_values=None, figsize=(14, 6), dpi=300,
                  arrow_head_size=24, arrow_width=2.8, lag_label_size=11,
                  node_size=.40, node_label_size=10, title_size=12,
                  curve_radius=.16, edge_colors=None, node_colors=None):
    """Yield one saved combined PDF per tau; always use the fitted settings."""
    out = _folder(run_name)
    settings = json.loads((out / "settings.json").read_text(encoding="utf-8"))
    taus = list(settings["tau_max_values"] if tau_max_values is None else tau_max_values)
    if len(set(taus)) != len(taus) or not set(taus) <= set(settings["tau_max_values"]):
        raise ValueError("Choose distinct tau_max values from this saved run.")
    figdir = _folder(run_name, figures=True)
    figdir.mkdir(parents=True, exist_ok=True)
    edge_cmap = LinearSegmentedColormap.from_list("custom_edges", base.EDGE_COLORS if edge_colors is None else edge_colors)
    node_cmap = LinearSegmentedColormap.from_list("custom_nodes", base.NODE_COLORS if node_colors is None else node_colors)
    plt.rcParams.update({"font.family": base.FONT_FAMILY, "pdf.fonttype": 42})
    for maximum in taus:
        fig, axes = plt.subplots(1, 2, figsize=figsize, facecolor="white")
        try:
            fig.subplots_adjust(left=.03, right=.98, top=.75, bottom=.20, wspace=.12)
            for ax, b in zip(axes, ("action", "outcome")):
                cfg = settings["networks"][b]
                regions = cfg["rois"]
                with np.load(out / f"{b}_tau_max_{maximum}.npz") as saved:
                    if saved["settings_key"].item() != settings["settings_key"]:
                        raise ValueError("Incomplete or mixed run; rerun the calculation for this run_name.")
                    base.tp.plot_graph(graph=saved["graph"], val_matrix=saved["val_matrix"],
                        var_names=[_label(r) for r in regions], fig_ax=(fig, ax),
                        node_pos=_positions(len(regions)), node_aspect=1.,
                        node_size=node_size, node_label_size=node_label_size,
                        arrow_linewidth=arrow_width, arrowhead_size=arrow_head_size,
                        curved_radius=curve_radius, link_label_fontsize=lag_label_size,
                        lag_array=np.arange(maximum+1), show_colorbar=False,
                        vmin_edges=0, vmax_edges=1, cmap_edges=edge_cmap,
                        vmin_nodes=0, vmax_nodes=1, cmap_nodes=node_cmap)
                for label in ax.texts:
                    if label.get_text().replace(",", "").replace(" ", "").isdigit():
                        label.set_color("#475569")
                        label.set_path_effects([])
                        label.set_bbox(dict(facecolor="white", edgecolor="none", pad=.6, alpha=.9))
                        label.set_zorder(5)
                for patch in ax.patches:
                    if isinstance(patch, Ellipse):
                        patch.set_edgecolor("white"); patch.set_linewidth(1.5)
                ax.set_xlim(-1, 1); ax.set_ylim(-.8, .8); ax.set_aspect("equal", adjustable="box")
                start, end = cfg["window_ms"]
                ax.set_title(f"{b.capitalize()} binding\n{start:g}-{end:g} ms", loc="left",
                             fontsize=title_size, fontweight="bold", pad=12)
            for rect, cmap, label in [([.14,.08,.30,.025], edge_cmap,"Cross-region GPDC"),
                                      ([.56,.08,.30,.025], node_cmap,"Self-dependence GPDC")]:
                cb = fig.colorbar(ScalarMappable(norm=Normalize(0,1),cmap=cmap),
                                  cax=fig.add_axes(rect), orientation="horizontal")
                cb.set_label(label, fontsize=10); cb.ax.tick_params(labelsize=9)
            fig.suptitle(f"Network comparison | tau_max = {maximum} ({maximum*BIN_MS} ms)",
                         fontsize=16, fontweight="bold", y=.98)
            fig.text(.5,.88,f"Lags {settings['tau_min']}-{maximum} | PC alpha = {settings['pc_alpha']:g} | "
                     f"alpha level = {settings['alpha_level']:g} | Edge labels: lag (1 = {BIN_MS} ms)",
                     ha="center", fontsize=11, color="#475569")
            pdf = figdir / f"combined_tau_max_{maximum}.pdf"
            fig.savefig(pdf, bbox_inches="tight", facecolor="white")
        finally:
            plt.close(fig)
        yield pdf

if __name__ == "__main__":
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument("--plot-only",action="store_true")
    args=parser.parse_args()
    networks={"action":dict(window_ms=(100,600),rois=["PreSMAL","PreSMAR","InsulaR"]),
              "outcome":dict(window_ms=(100,900),rois=["PrecentralR","PostcentralL","MTGR"])}
    if not args.plot_only:
        calculate(networks,[7],run_name="manuscript",tau_min=1,pc_alpha=.05,alpha_level=.05)
    list(plot_combined(run_name="manuscript",tau_max_values=[7]))
