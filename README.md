# 🧠 Dual Agency Pathways

Neural mechanisms underlying perceptual shifts in action and outcome timing under voluntary and involuntary conditions.

This repository provides the analysis-ready data, computational workflows and figures supporting the manuscript:

> *Dual neural pathways for agency: Distinct spatiotemporal signatures of action binding and outcome binding in the human brain.*

The analyses follow the progression from behavioural temporal binding to regional MEG responses, functional regression, directed connectivity, single-trial modelling and intracranial recordings.

## Repository structure

```text
dual-agency-pathways/
├── analysis/
│   ├── 01_temporal_binding_effect/
│   ├── 02_meg_results/
│   ├── 03_functional_regression/
│   ├── 04_network/
│   ├── 05_LMM/
│   ├── 06_ieeg/
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

The numbered modules identify the same analyses across the repository. Each module contains its workflow scripts and helper functions. `data/` holds the analysis inputs, `results/` holds the final numerical tables, and `figures/` holds the combined PDF figures. Intermediate features, phase angles and network matrices are stored in `analysis/.cache/derived`.

## Analytical modules

| Module | Analysis | Main outputs |
|---|---|---|
| `01_temporal_binding_effect` | Quantify baseline-corrected perceptual shifts and compare voluntary and involuntary outcome binding | Signed-rank tests and behavioural plots |
| `02_meg_results` | Display regional MEG amplitude time courses across conditions | Mean and standard-error plots for 12 regions |
| `03_functional_regression` | Relate participant-level binding contrasts to MEG time courses using sparse functional regression | Selected regional coefficients and time windows |
| `04_network` | Estimate directed interactions using PCMCI–GPDC | Lagged-link estimates and action/outcome network plots |
| `05_LMM` | Extract trial-wise neural features, characterize phase distributions and fit mixed-effects models | Fixed-effect estimates, mean-amplitude panels and four-band ITPC plots |
| `06_ieeg` | Characterize behavioural and neural measures in intracranial MTG recordings | Patient summaries and behavioural/ITPC plots |

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
| MATLAB | MEG and intracranial phase-distribution figures | MATLAB entry points in `05_LMM` and `06_ieeg` |

The network module uses Tigramite; neural feature extraction uses MNE. Python dependencies are listed in the [network](analysis/04_network/requirements.txt), [MEG feature](analysis/05_LMM/requirements.txt) and [intracranial feature](analysis/06_ieeg/requirements.txt) requirements files. Pandoc is used for optional HTML rendering.

## Data, results and figures

**Data.** The inputs are preprocessed regional MEG time courses, single-trial epochs, behavioural responses and intracranial bipolar-channel recordings. MEG and behavioural inputs include participants 2–20; intracranial analyses use patients 2 and 4. Dataset fields, condition codes and baseline corrections are documented in [data/README.md](data/README.md).

**Results.** The repository retains the final behavioural tests, functional-regression coefficients, network links, LMM fixed effects and intracranial patient summaries. File contents and column definitions are documented in [results/README.md](results/README.md).

**Figures.** Each figure category is exported as one combined PDF in its numbered module folder. The single-trial module provides mean-amplitude, LMM-coefficient and ITPC figures. Its ITPC figure shows PreSMA-L, PreCG-R and MTG-R, each with Delta, Theta, Alpha and Beta.

## Brain-region labels

Figure labels use the following abbreviations:

| Data identifier | Figure label |
|---|---|
| `Precentral-R` or `PrecentralR` | PreCG-R |
| `Postcentral-R` or `PostcentralR` | PoCG-R |
| `Postcentral-L` or `PostcentralL` | PoCG-L |

