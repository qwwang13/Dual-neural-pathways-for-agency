# Results

This directory contains the final numerical tables produced by the analysis workflows. Each numbered folder corresponds to the same module under `analysis/`.

## Output files

| Module | File | Contents |
|---|---|---|
| `01_temporal_binding_effect` | `normality.csv` | Shapiro–Wilk tests for the four binding conditions and the paired outcome difference |
| `01_temporal_binding_effect` | `against_zero.csv` | Exact one-sided signed-rank tests of perceptual shifts against zero |
| `01_temporal_binding_effect` | `outcome_paired.csv` | Exact one-sided paired comparison of voluntary and involuntary outcome binding |
| `02_meg_results` | `permutation_results.7z` (Release asset) | Original Brainstorm permutation export and two saved source-level p/t maps |
| `03_functional_regression` | `action_binding.csv` | Selected functional-regression coefficients and analysis windows for action binding |
| `03_functional_regression` | `outcome_binding.csv` | Selected functional-regression coefficients and analysis windows for outcome binding |
| `04_network` | `all_links.csv` | GPDC estimates, p-values and lag information for the tested links |
| `05_LMM` | `fixed_effects.csv` | Voluntary screening-model coefficients and involuntary follow-up coefficients for the selected regions |
| `06_ieeg` | `patient_summary.csv` | Patient-specific trial counts, behavioural summaries and ITPC values |
| `07_supplementary_robustness` | `bayesian_t_tests.csv` | JZS Bayesian one-sample and paired t tests for behavioural shifts |
| `07_supplementary_robustness` | `itpc_bootstrap.csv` | Condition-specific participant ITPC slopes and paired-bootstrap differences/intervals |
| `07_supplementary_robustness` | `RESULTS.md` | Three-decimal summary of both supplementary analyses |

## Behavioural tests

`normality.csv` includes the variable name, participant count, Shapiro–Wilk W and p-value. The signed-rank tables include the comparison, sample size, mean, median, standard deviation, test direction and p-value. `V` is the larger signed-rank sum; `W_positive` and `W_negative` retain both components. Behavioural values are in milliseconds.

## Functional-regression coefficients

| Column | Meaning |
|---|---|
| `region` | Original brain-region identifier |
| `time_window` | Averaging resolution in milliseconds |
| `lower`, `upper` | Bounds of the coefficient bin in milliseconds |
| `coef` | Selected coefficient multiplied by `10 / time_window` |
| `primary_lower_ms`, `primary_upper_ms` | Selected analysis-window bounds used for plotting and downstream feature extraction |

Rows contain the non-zero coefficients of the selected models for the retained regions. The primary window is recorded alongside those coefficients.

## Network links

`source` and `target` identify the directed pair. `gpdc_value` is the estimated GPDC and `p_value` is its empirical p-value. `tau` gives the lag in 10 ms bins; `lag_ms` gives the same lag in milliseconds. `binding` identifies the action or outcome network, and `significant` marks links passing the configured threshold.

## Mixed-effects estimates

`binding`, `roi`, `condition` and `term` identify each coefficient. `estimate`, `std_error`, `df`, `t_value` and `p_value` describe the fixed effect; `conf_low` and `conf_high` give the 95% confidence interval. `signif` records the p-value category.

The table contains voluntary-condition estimates for all retained regions, followed by involuntary-condition estimates for regions with at least one voluntary neural-predictor p < 0.05. Both stages use the same nine predictors, standardized separately within condition. Behaviour remains in milliseconds.

## Intracranial summaries

`subject`, `label` and `condition` identify the patient, bipolar channel and condition. `n_trials` gives the matched trial count; `behavior_mean` and `behavior_sd` summarize perceptual shifts in milliseconds. The four `*_bandITPC` columns give ITPC within each patient and condition.

Intermediate trial features, phase angles, network matrices and window-matching metadata are stored in `analysis/.cache/derived`. Figures are saved under `figures/`. Running an analysis updates its corresponding outputs.

## Supplementary robustness results

In `bayesian_t_tests.csv`, `BF10` compares the two-sided JZS alternative with zero; `BF_positive0` and `BF_negative0` use the corresponding truncated and renormalized priors. `prior_cauchy_scale` is sqrt(2)/2. The paired outcome row is voluntary minus involuntary.

In `itpc_bootstrap.csv`, each condition has a separate participant-level OLS fit with an intercept. `difference_ms_per_unit_ITPC` is the voluntary slope minus the involuntary slope. `ci_low` and `ci_high` are the 95% percentile limits from 20,000 paired participant draws. `bootstrap_p` retains the approximate two-sided zero-crossing tail calculation, with a plus-one correction. The CSV files retain numerical precision; `RESULTS.md` formats continuous statistics to three decimal places.

## MEG permutation results

The original action/outcome maps and their archive members are described in [`02_meg_results/README.md`](02_meg_results/README.md). The saved `pmap` and `tmap` arrays are retained as supplied; Brainstorm statistical display/FDR thresholding is a separate setting. The binary archive is prepared for a GitHub Release rather than an ordinary repository file.
