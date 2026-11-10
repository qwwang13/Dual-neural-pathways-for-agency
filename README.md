# 🧠 Dual Agency Pathways

Neural mechanisms underlying perceptual shifts in action and outcome timing under voluntary and involuntary conditions.

This repository provides the analysis-ready data, computational workflows and figures supporting the manuscript:

> *Dual neural pathways for agency: Distinct spatiotemporal signatures of action binding and outcome binding in the human brain.*

The analyses follow the progression from behavioural temporal binding to regional MEG responses, functional regression, directed connectivity, single-trial modelling and intracranial recordings.

## Repository structure

```text
dual-agency-pathways/
├── analysis/
│   ├── .cache/derived/
│   ├── 01_temporal_binding_effect/
│   ├── 02_meg_results/
│   ├── 03_functional_regression/
│   ├── 04_network/
│   ├── 05_LMM/
│   ├── 06_ieeg/
│   ├── 07_supplementary_robustness/
│   ├── common.R
│   ├── config.py
│   ├── run.R
│   └── README.md
├── data/
├── results/
├── figures/
├── dual-agency-pathways.Rproj
├── renv.lock
└── README.md
```

The numbered modules identify the same analyses across the repository. Each module contains its workflow scripts and helper functions. `data/` holds the analysis inputs, `results/` holds the final numerical tables, and `figures/` holds the combined PDF figures. Generated feature tables, phase angles and network matrices are stored in `analysis/.cache/derived`. This is the repository's intermediate-output directory. For example, `analysis/05_LMM/01_features.py` reads the preprocessed epochs in `data/05_LMM` and saves the nine neural features to `analysis/.cache/derived/05_LMM/features`. The LMMs and supplementary ITPC analysis read those same feature tables. Supplied copies are included and can be regenerated from the existing inputs; see [the feature-data explanation](data/README.md#derived-meg-feature-tables).

## Analytical modules

| Module | Analysis | Main outputs |
|---|---|---|
| `01_temporal_binding_effect` | Quantify baseline-corrected perceptual shifts and compare voluntary and involuntary outcome binding | Signed-rank tests and behavioural plots |
| `02_meg_results` | Paired MEG permutation tests and regional time-course plots | Original Brainstorm script/statistical maps and mean/standard-error plots for 12 regions |
| `03_functional_regression` | Relate participant-level binding contrasts to MEG time courses using sparse functional regression | Selected regional coefficients and time windows |
| `04_network` | Estimate directed interactions using PCMCI–GPDC | Lagged-link estimates and action/outcome network plots |
| `05_LMM` | Extract trial-wise neural features, characterize phase distributions and fit mixed-effects models | Fixed-effect estimates, mean-amplitude panels and four-band ITPC plots |
| `06_ieeg` | Characterize behavioural and neural measures in intracranial MTG recordings | Patient summaries and behavioural/ITPC plots |
| `07_supplementary_robustness` | Bayesian behavioural robustness and participant-level ITPC slope comparisons | JZS Bayes factors and paired-bootstrap confidence intervals |

See [analysis/README.md](analysis/README.md) for the scripts, execution order and analysis settings in each module.

## Running the analyses

Open `dual-agency-pathways.Rproj` in RStudio, or run the commands below from the repository root. Restore the R packages first:

```sh
Rscript analysis/restore.R
```

The full workflow is:

```sh
Rscript analysis/run.R behavior
Rscript analysis/run.R meg
Rscript analysis/run.R functional
python analysis/04_network/04_network.py
python analysis/05_LMM/01_features.py
python analysis/05_LMM/02_itpc.py
Rscript analysis/run.R models
Rscript analysis/run.R coefficients
Rscript analysis/run.R amplitude
python analysis/06_ieeg/01_features_itpc.py
Rscript analysis/run.R ieeg
Rscript --vanilla analysis/run.R robustness
```

Run the phase-plot entry points in MATLAB:

```matlab
run('analysis/05_LMM/plot_itpc.m')
run('analysis/06_ieeg/plot_itpc.m')
```

Functional regression supplies the time windows for MEG feature extraction and the intracranial MTG analysis. Extract the neural features before fitting the LMMs. Scripts within these modules are numbered in execution order.

The R commands execute the report code and save the corresponding tables and figures. To also render an HTML report, append `--html`, for example:

```sh
Rscript analysis/run.R behavior --html
```

### MEG permutation tests and large source-map files

The original Brainstorm permutation export is available as [`a01_permutation_tests.m`](analysis/02_meg_results/a01_permutation_tests.m). Its 76 source-map inputs and two saved statistical maps are organized separately from the Fig. 2 plotting summaries. The large archives are prepared as GitHub Release assets, while their file inventories, pairing list and instructions are tracked in the repository. See [the MEG module](analysis/02_meg_results/README.md) for the Brainstorm setup, [data layout](data/02_meg_results/meg/README.md) and [saved results](results/02_meg_results/README.md). The existing `meg` R command redraws the time-course figure; it does not rerun the permutation calculation.

### Supplementary robustness analyses

Run `Rscript analysis/run.R robustness` for the Bayesian behavioural tests and participant-level ITPC bootstrap. The R Markdown file has two separate sections: Bayesian behavioural analysis and Participant-level ITPC bootstrap. Each section reads its inputs, runs its own analysis and saves its own result table. It reads the original behavioural summary and existing MEG feature tables, and saves two final CSV tables. Python with NumPy generates the paired bootstrap draws from the original fixed seed. See [analysis/README.md](analysis/README.md#supplementary-robustness-analyses) for settings and dependencies.

### Redrawing supplied results

Use the supplied tables and cached features to redraw the figures. Run the `coefficients`, `amplitude` or `ieeg` R step, or either MATLAB phase-plot entry point, to redraw those figures. To redraw the networks, use:

```sh
python analysis/04_network/04_network.py --plot-only
```

The functional-regression report also identifies the setup and plotting chunks needed to redraw its saved results.

## Software environment

| Software | Use | Environment specification |
|---|---|---|
| R 4.4.3 | Behavioural tests, functional regression, LMMs and plotting | `renv.lock`, including glmnet 4.1-10 |
| Python | Time-frequency features, ITPC and PCMCI–GPDC | Module-specific `requirements.txt` files |
| MATLAB | Original Brainstorm MEG permutation export and phase-distribution figures | Brainstorm for `02_meg_results`; plotting entry points in `05_LMM` and `06_ieeg` |

The network module uses Tigramite; neural feature extraction uses MNE. Python dependencies are listed in the [network](analysis/04_network/requirements.txt), [MEG feature](analysis/05_LMM/requirements.txt) and [intracranial feature](analysis/06_ieeg/requirements.txt) requirements files. Pandoc is used for optional HTML rendering.

## Data, results and figures

**Data.** The inputs are preprocessed regional MEG time courses, single-trial epochs, behavioural responses and intracranial bipolar-channel recordings. MEG and behavioural inputs include participants 2–20; intracranial analyses use patients 2 and 4. Dataset fields, condition codes and baseline corrections are documented in [data/README.md](data/README.md).

**Results.** The repository retains the final behavioural tests, functional-regression coefficients, network links, LMM fixed effects, intracranial patient summaries and supplementary Bayes factors/ITPC bootstrap intervals. File contents and column definitions are documented in [results/README.md](results/README.md).

**Figures.** Each figure category is exported as one combined PDF in its numbered module folder. The single-trial module provides mean-amplitude, LMM-coefficient and ITPC figures. Its ITPC figure shows PreSMA-L, PreCG-R and MTG-R, each with Delta, Theta, Alpha and Beta.

## Brain-region labels

Figure labels use the following abbreviations:

| Data identifier | Figure label |
|---|---|
| `Precentral-R` or `PrecentralR` | PreCG-R |
| `Postcentral-R` or `PostcentralR` | PoCG-R |
| `Postcentral-L` or `PostcentralL` | PoCG-L |

