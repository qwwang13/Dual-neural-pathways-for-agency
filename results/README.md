------------------------------------------------------------------------

# Results

This directory contains statistical outputs corresponding to each analytical module of the manuscript.

All files are provided in `.csv` format.

------------------------------------------------------------------------

## General Column Structure

Common column names include:

-   `Time` -- Time point (in milliseconds) relative to event onset
-   `subject` -- Participant identifier
-   `region` -- Brain region label
-   `trial` -- Trial index
-   `band` -- Frequency band label
-   `condition` -- Experimental condition identifier

Unless otherwise specified, condition labels reflect contrasts relative to their respective baseline conditions.

------------------------------------------------------------------------

## Condition Naming Conventions

The `condition` column encodes experimental manipulations. Unless otherwise specified, condition values represent baseline-corrected signals.

### Action Binding (MEG)

Voluntary:

-   `va_ba`
-   `vaa_ba`
-   `va_a`

Involuntary:

-   `ia_bi`
-   `iaa_bi`
-   `ia_a`

------------------------------------------------------------------------

### Outcome Binding (MEG)

Voluntary:

-   `vs_bs`
-   `vss_bs`
-   `vs_s`

Involuntary:

-   `is_bs`
-   `iss_bs`
-   `is_s`

------------------------------------------------------------------------

### Probability Manipulation (MEG)

75% with outcome:

-   `a75_ay_ba`

75% without outcome:

-   `a75_an_ba`
-   `a75_an`

50% with outcome:

-   `a50_ay_ba`
-   `a50_ay`

50% without outcome:

-   `a50_an_ba`
-   `a50_an`

------------------------------------------------------------------------

### Difference Conditions (Baseline-Corrected Contrasts)

-   Action binding (Voluntary − Involuntary): `vaa_iaa`
-   Outcome binding (Voluntary − Involuntary): `vss_iss`
-   75% without outcome − 50% without outcome: `a75_a50_an`
-   50% with outcome − 50% without outcome: `a50_ay_an`

------------------------------------------------------------------------

### sEEG Conditions

In the sEEG dataset, baseline subtraction is **not applied**.

-   `Voluntary_OutcomeOnsetReportSound`
-   `Involuntary_OutcomeOnsetReportSound`

------------------------------------------------------------------------

# 01_temporal_binding_effect

Statistical evaluation of perceptual shift.

-   `one_sample_tests.csv` One-sample t-tests against zero for each condition.

-   `paired_tests.csv` Paired comparisons between conditions. The `contrast` column specifies the compared conditions.

Files include sample size, test statistics, and p-values.

------------------------------------------------------------------------

# 03_functional_regression

Results of functional regression identifying significant brain regions and associated temporal windows under different condition contrasts.

Files:

-   `action_binding.csv`
-   `outcome_binding.csv`
-   `75_50_without_outcome.csv`
-   `50_with_without_outcome.csv`

Column definitions:

-   `time_window` -- Temporal window width (in ms)
-   `lower` -- Lower bound of the significant time window
-   `upper` -- Upper bound of the significant time window
-   `coef` -- Estimated regression coefficient within the time window

------------------------------------------------------------------------

# 04_network

Directed functional connectivity results.

Column definitions:

-   `source` -- Originating brain region
-   `target` -- Target brain region
-   `tau` -- Temporal lag (order of transfer)
-   `gpdc_value` -- Generalized Partial Directed Coherence statistic
-   `p_value` -- Statistical significance of the connectivity estimate

Files:

-   `action_binding.csv`
-   `outcome_binding.csv`
-   `75_50_without_outcome.csv`

------------------------------------------------------------------------

# 05_LMM

Linear mixed-effects modeling results under different condition contrasts.

Files:

-   `action_binding.csv`
-   `outcome_binding.csv`
-   `75_50_without_outcome.csv`
-   `50_with_without_outcome.csv`

Column definitions:

-   `ROI` -- Brain region identified under the specified condition
-   `term` -- Extracted feature name

For each condition, the following columns are provided:

-   `estimate_condition` -- Estimated fixed-effect coefficient
-   `std_error_condition` -- Standard error
-   `df_error_condition` -- Degrees of freedom
-   `t_value_error_condition` -- Test statistic
-   `p_value_error_condition` -- Associated p-value
-   `signif_error_condition` -- Significance indicator

------------------------------------------------------------------------

# 06_ieeg

Pearson correlation results between sEEG-derived features and behavioral measures.

File:

-   `ieeg_PostcentralL_pearson.csv`

Column definitions:

-   `feature` -- Feature name
-   `n` -- Sample size
-   `r` -- Pearson correlation coefficient
-   `p` -- Associated p-value

------------------------------------------------------------------------
