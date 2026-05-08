------------------------------------------------------------------------

# Data

This directory contains processed and analysis-ready datasets supporting each analytical module of the manuscript.

All files are provided in `.csv` format.

The datasets do **not** contain raw MEG or sEEG recordings. Instead, they comprise preprocessed and derived measures used for modeling and visualization.

------------------------------------------------------------------------

## General Column Structure

Common column names include:

-   `Time`-- Time point (in milliseconds) relative to event onset
-   `subject`-- Participant identifier
-   `region`-- Brain region label
-   `trial`-- Trial index
-   `band`-- Frequency band label
-   `condition`-- Experimental condition identifier

Condition labels reflect experimental contrasts relative to their respective baseline conditions (unless otherwise specified).

Files containing `"behavior"` or `"perceptual_shift"` refer to behavioral measures. All other files correspond to MEG or sEEG-derived neural data.

------------------------------------------------------------------------

## Condition Naming Conventions

The `condition` column encodes experimental manipulations. Unless otherwise specified, condition values represent signals after subtraction of their respective baseline conditions.

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

### Difference Conditions (Baseline-Corrected Contrasts)

These conditions represent contrasts computed after baseline correction:

-   Action binding (Voluntary − Involuntary): `vaa_iaa`
-   Outcome binding (Voluntary − Involuntary): `vss_iss`

------------------------------------------------------------------------

### sEEG Conditions

In the sEEG dataset, baseline subtraction is **not applied**.

Outcome binding conditions:

-   `Voluntary_OutcomeOnsetReportSound`
-   `Involuntary_OutcomeOnsetReportSound`

------------------------------------------------------------------------

## Data Level

-   Sections **01--04** contain subject-level data.
-   Sections **05--06** contain trial-level data.

------------------------------------------------------------------------

# 01_temporal_binding_effect

-   `perceptual_shift_action_outcome.csv` Perceptual shift values under voluntary and involuntary conditions. Column names correspond directly to condition labels.

------------------------------------------------------------------------

# 02_meg_results

This directory contains four CSV files for stimulus-locked MEG amplitude responses:

- `ta_action_va_12.csv`
- `ta_action_ia_12.csv`
- `ta_outcome_vs_12.csv`
- `ta_outcome_is_12.csv`

Each file contains summary time-course data for 12 ROIs, with columns `region`, `Time`, `condition`, `Value`, and `SE`. Here, `Value` denotes the mean amplitude, and `SE` denotes the standard error of the mean.

------------------------------------------------------------------------

# 03_functional_regression

The following files contain MEG amplitude **differences (0--1000 ms)** across conditions:

-   `meg_action_difference.csv`
-   `meg_outcome_difference.csv`

Behavioral difference:

-   `perceptual_shift_difference.csv`

### plot_data/

Contains averaged MEG amplitude differences (averaged across subjects) used for visualization:

-   `meg_mean_action.csv`
-   `meg_mean_outcome.csv`

------------------------------------------------------------------------

# 04_network

Contains mean MEG amplitude differences (averaged across subjects) within predefined time windows.

-   `lower` and `upper` indicate the temporal window boundaries.
-   `Value` represents the mean amplitude difference within that window.

Files include:

-   `action_binding_for_network.csv`
-   `outcome_binding_for_network.csv`

Each file represents a distinct condition contrast.

------------------------------------------------------------------------

# 05_LMM

Trial-level MEG data (−200 ms to 1000 ms).

Files are named as:

```         
RegionName + ExperimentalCondition.csv
```

Each file contains single-trial MEG amplitude values for a specific region and condition. All conditions have been baseline-corrected. Condition differences are **not** computed at this stage.

-   `behavior_single_trial.csv` Behavioral responses at the single-trial level.

### feature/

Contains time-window-specific features extracted from MEG signals:

-   Mean amplitude
-   Band power
-   Band ITPC

Time windows correspond to those identified in `03_functional_regression`.

### plot_data_for_itpc/

-   Files ending with `angles` contain phase angle data.
-   Files ending with `itpc_mean` contain subject-averaged ITPC values.

------------------------------------------------------------------------

# 06_ieeg

Trial-level sEEG data (−200 ms to 1000 ms). An additional column, `label`, denotes the electrode contact label.

Unlike MEG data, these files are **not baseline-corrected**.

-   `ieeg_MTGR.csv`

-   `behavior_ieeg.csv` Single-trial behavioral response times under voluntary and involuntary outcome conditions.

### feature/

Contains extracted features from MTGR:

-   Mean amplitude
-   Band power
-   Band ITPC

Time windows correspond to those identified in `03_functional_regression`.

### plot_data_for_itpc/

-   Files ending with `angles` contain phase angle data.
-   Files ending with `itpc_mean` contain subject-averaged ITPC values.

------------------------------------------------------------------------
