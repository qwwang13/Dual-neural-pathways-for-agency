# Analysis workflows

Run scripts from the repository root. Each module contains its workflow files and helper functions; `common.R` and `config.py` provide shared paths and settings.

## Temporal binding effects

[`01_temporal_binding_effect/01_behavior.Rmd`](01_temporal_binding_effect/01_behavior.Rmd) reads participant-level perceptual shifts, runs normality and exact one-sided signed-rank tests, and generates the combined behavioural figure. Run with `Rscript analysis/run.R behavior`.

Action-binding tests use the greater-than-zero alternative; outcome-binding tests use the less-than-zero alternative. The paired outcome test uses voluntary minus involuntary with a less-than-zero alternative. Reported V is the larger signed-rank sum; both rank sums are retained in the output tables.

## MEG time courses

[`02_meg_results/02_timecourses.Rmd`](02_meg_results/02_timecourses.Rmd) reads the supplied regional means and between-participant standard errors, plots the voluntary and involuntary conditions, and saves the combined 12-region time-course figure. Run with `Rscript analysis/run.R meg`.

## Functional regression

[`03_functional_regression/03_functional.Rmd`](03_functional_regression/03_functional.Rmd) contains three sections: **Read data**, **Run functional regression**, and **Plot and save**. Run with `Rscript analysis/run.R functional`.

| Setting | Value |
|---|---|
| Participants and model domain | All 19 participants; 0–1000 ms |
| Averaging windows | 10, 20 and 40 ms |
| Elastic-net grid | Intercept included; alpha = 0.1, 0.2, …, 1.0; up to 100 lambda values per alpha |
| Continuity thresholds | 1–20 bins at 10 ms; 1–10 at 20 ms; 1–5 at 40 ms |
| Gap bridging | Gaps of at most 20 ms are bridged within each resolution before intersection |
| Common interval | Minimum span of 40 ms |
| Support limits | At most 20% non-zero bins in each resolution; at most five raw 10 ms segments outside the primary interval |

At each continuity threshold, eligible models are ranked by the number of non-zero coefficients, training R-squared, alpha and lambda path order. The three-resolution combination prioritizes the longest qualifying common interval, followed by the support and localization criteria implemented in `functional_binding_functions.R`.

The module saves action and outcome coefficient tables with the selected analysis-window bounds. Fitted model paths are stored in `analysis/.cache/functional`. To redraw saved results, run the setup chunk and the **Plot and save** section.

## Directed connectivity

[`04_network/04_network.py`](04_network/04_network.py) prepares participant-averaged MEG contrast time courses, estimates PCMCI–GPDC networks and saves the combined network figure. Run with `python analysis/04_network/04_network.py`; add `--plot-only` to use the saved network matrices.

Networks use 10 ms bins, lags 1–7, PC alpha = 0.05 and link alpha = 0.05. Action binding uses 100–600 ms; outcome binding uses 100–900 ms. Inputs are the mean of the 19 participant voluntary-minus-involuntary time courses.

## Single-trial features and mixed-effects models

| Entry point in `05_LMM` | Task |
|---|---|
| [`01_features.py`](05_LMM/01_features.py) | Extract mean amplitude, four band powers and four ITPC features |
| [`02_itpc.py`](05_LMM/02_itpc.py) | Export trial phase angles and participant ITPC summaries |
| [`03_models.Rmd`](05_LMM/03_models.Rmd) | Fit voluntary models and the selected involuntary follow-up models |
| [`04_coefficients.Rmd`](05_LMM/04_coefficients.Rmd) | Plot mean-amplitude and ITPC coefficient estimates with 95% confidence intervals |
| [`05_amplitude.Rmd`](05_LMM/05_amplitude.Rmd) | Plot pooled-trial amplitude distributions and paired participant means |
| [`plot_itpc.m`](05_LMM/plot_itpc.m) | Draw the four-band phase distributions in MATLAB |

Feature windows are read from the functional-regression coefficient tables. Trials are matched by participant, condition and trial identifier. The nine neural predictors are z-scored separately within each condition after matching; behavioural responses remain in milliseconds. Models use maximum likelihood and a participant random intercept.

All retained regions are first fitted in the voluntary condition. Regions with at least one neural predictor having p < 0.05 are then fitted in the involuntary condition using the same nine predictors. Coefficient p-values are two-sided.

ITPC is computed across trials within each participant and condition. Figure 4 shows the mean-amplitude and four ITPC coefficients for PreSMA-L, PreCG-R and MTG-R; estimates with p < 0.05 are coloured. Amplitude histograms pool matched trials within condition into 20 shared equal-width bins, and the paired points show participant condition means.

## Intracranial recordings

| Entry point in `06_ieeg` | Task |
|---|---|
| [`01_features_itpc.py`](06_ieeg/01_features_itpc.py) | Extract trial features and phase measures using the selected MTG-R window |
| [`02_ieeg.Rmd`](06_ieeg/02_ieeg.Rmd) | Save patient summaries and plot behavioural distributions with ITPC values |
| [`plot_itpc.m`](06_ieeg/plot_itpc.m) | Draw the four-band phase distributions in MATLAB |

The analysis uses patient 2, channel `chan83-chan84`, and patient 4, channel `chan10-chan11`.

## Input preparation and execution

MEG and behavioural inputs already contain the corresponding participant-specific baseline correction. Each behavioural trial has the participant's baseline-condition average subtracted. Each MEG trial has the participant's baseline-condition mean waveform subtracted at the same region and time point before feature extraction. Voluntary action uses `ba`, involuntary movement uses `bi`, and both outcome conditions use `bs`. Functional-regression and network contrasts are voluntary minus involuntary after baseline correction. See [the data documentation](../data/README.md) for the condition codes.

The R runner executes the report chunks and saves tables and figures. Append `--html` to render an HTML report with Pandoc. The project root is identified by `dual-agency-pathways.Rproj`; `DUAL_AGENCY_ROOT` can specify it explicitly. `DUAL_AGENCY_R_LIBRARY` can point to an existing R package library.
