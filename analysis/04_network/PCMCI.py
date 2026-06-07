# pcmci_utils.py
from __future__ import annotations

from pathlib import Path
from typing import Dict, Any, List, Tuple, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import tigramite.data_processing as pp
from tigramite.pcmci import PCMCI
from tigramite.independence_tests.gpdc import GPDC
from tigramite import plotting as tp


def load_pcmci_dataframe(csv_path: str | Path) -> Tuple[pp.DataFrame, List[str], np.ndarray]:
    df = pd.read_csv(csv_path)
    pivot = df.pivot(index="lower", columns="region", values="Value").sort_index()

    if pivot.isna().any().any():
        print("Warning: missing values detected, dropping rows with NaN.")
        pivot = pivot.dropna(axis=0)

    regions = pivot.columns.tolist()
    data = pivot.to_numpy()
    dataframe = pp.DataFrame(data, var_names=regions)
    return dataframe, regions, data


def run_pcmci_gpdc(
    dataframe: pp.DataFrame,
    regions: List[str],
    label: str,
    tau_max: int,
    pc_alpha: float = 0.05,
    alpha_level: float = 0.05,
    label_link: str = "GPDC",
    label_node: str = "auto-GPDC",
    significance: str = "analytic",
    sig_samples: int = 200,
    seed: int = 0,
    save_fig_path: Optional[str | Path] = None,
    show_fig: bool = True,
    verbosity: int = 0,
    print_significant: bool = False,
) -> Dict[str, Any]:
    print(f"\n[{label}] tau_max={tau_max}")

    if significance == "analytic":
        cond_ind_test = GPDC(significance="analytic")
    elif significance == "shuffle_test":
        cond_ind_test = GPDC(significance="shuffle_test", sig_samples=sig_samples, seed=seed)
    else:
        raise ValueError("significance must be 'analytic' or 'shuffle_test'")

    pcmci = PCMCI(dataframe=dataframe, cond_ind_test=cond_ind_test, verbosity=verbosity)
    results = pcmci.run_pcmci(
        tau_min=1,
        tau_max=tau_max,
        pc_alpha=pc_alpha,
        alpha_level=alpha_level,
    )

    if print_significant:
        pcmci.print_significant_links(
            p_matrix=results["p_matrix"],
            val_matrix=results["val_matrix"],
            alpha_level=alpha_level,
        )

    two_nodes = len(regions) == 2
    if two_nodes:
        fig, ax = tp.plot_graph(
            val_matrix=results["val_matrix"],
            graph=results["graph"],
            var_names=regions,
            link_colorbar_label="GPDC",
            node_colorbar_label="auto-GPDC",
            figsize=(6, 3),
            node_pos={"x": [-0.6, 0.6], "y": [0.0, 0.0]},
            node_aspect=1.0,
            node_size=0.2,
            curved_radius=0.05,
            arrow_linewidth=3.0,
            arrowhead_size=10,
            show_colorbar=True,
        )
        ax.set_xlim(-0.7, 0.7)
        ax.set_ylim(-0.5, 0.5)
        ax.set_aspect("equal", adjustable="box")
    else:
        fig, ax = tp.plot_graph(
            val_matrix=results["val_matrix"],
            graph=results["graph"],
            var_names=regions,
            link_colorbar_label=label_link,
            node_colorbar_label=label_node,
            figsize=(6, 6),
        )

    fig.suptitle(f"{label}")

    if save_fig_path is not None:
        save_fig_path = Path(save_fig_path)
        save_fig_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_fig_path, dpi=300, bbox_inches="tight")

    if show_fig:
        plt.show()
    else:
        plt.close(fig)

    return results


def results_to_long_df(
    results: Dict[str, Any],
    regions: List[str],
    tau_min: int = 1,
    tau_max: Optional[int] = None,
) -> pd.DataFrame:
    p = results["p_matrix"]      # shape: [N, N, tau_max+1]
    v = results["val_matrix"]    # shape: [N, N, tau_max+1]
    N = len(regions)
    if tau_max is None:
        tau_max = p.shape[2] - 1

    rows = []
    for i in range(N):
        for j in range(N):
            for tau in range(tau_min, tau_max + 1):
                rows.append(
                    {
                        "source": regions[i],
                        "target": regions[j],
                        "tau": tau,
                        "gpdc_value": float(v[i, j, tau]),
                        "p_value": float(p[i, j, tau]),
                    }
                )
    return pd.DataFrame(rows)
