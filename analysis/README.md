# Analysis Modules

This directory contains all computational workflows supporting the manuscript. Each numbered subdirectory corresponds to a dedicated analytical component.

Scripts are organized into:

-   **Workflow files** (e.g., `.Rmd`, `.ipynb`, `.m`)
-   **Utility functions** (e.g., `.R`, `.py`)

File numbering indicates execution order where sequential processing is required.

------------------------------------------------------------------------

## 01_temporal_binding_effect

Behavioral quantification of perceptual shift.

-   `01_temporal_binding_effect.Rmd` Workflow for generating perceptual shift figures.

-   `temporal_binding_effect_utils.R` Utility functions for visualization.

------------------------------------------------------------------------

## 02_meg_results

Visualization of stimulus-locked MEG amplitude responses across conditions.

-   `02_plot_mean_SE_12roi.Rmd` Workflow for generating MEG amplitude figures.

-   `meg_timecourse_plot_utils.R` Utility functions for visualization.

------------------------------------------------------------------------

## 03_functional_regression

Functional linear regression modeling of behavioral outcomes using MEG signals as predictors.

-   `03_1_action_binding_functional_regression.Rmd` Functional regression for action binding.

-   `03_2_outcome_binding_functional_regression.Rmd` Functional regression for outcome binding.

-   `elasticnet_functions.R` Modeling utilities.

-   `fun_plot_functional_regression.R` Plotting utilities.

------------------------------------------------------------------------

## 04_network

Construction and analysis of the functional connectivity network.

-   `04_network.ipynb` Workflow for directed connectivity estimation and network modeling.

-   `PCMCI.py` Core functions for network estimation.

Environment-specific files:

-   `.gitignore`
-   `environment.yaml`
-   `requirements.txt`

------------------------------------------------------------------------

## 05_LMM

Single-trial frequency analysis and linear mixed-effects modeling.

Sequential workflow:

1.  `05_1_features_save.ipynb` Extraction of MEG mean amplitude and time--frequency features. Outputs saved to `data/05_LMM/feature/`. Uses `TF_features.py`.

2.  `05_2_ITPC_save.ipynb` Computation of phase angles for ITPC visualization. Outputs saved to `data/05_LMM/plot_data_for_itpc/`. Uses `ITPC.py`.

3.  `run_plot_itpc_rose_all_05_3.m` Phase visualization workflow. Uses `plot_itpc_rose_cartesian.m`.

4.  `05_4_LMM.Rmd` Linear mixed-effects modeling across conditions. Uses:

    -   `fun_lmm_split.R` (modeling)
    -   `fun_plot_lmm_split.R` (visualization)

Environment-specific files:

-   `.gitignore`
-   `environment.yaml`
-   `requirements.txt`

------------------------------------------------------------------------

## 06_ieeg

Intracranial EEG validation of MEG-derived findings.

Sequential workflow:

1.  `06_1_features_save.ipynb` Extraction of sEEG mean amplitude and time--frequency features. Outputs saved to `data/06_ieeg/feature/`. Uses `TF_features_ieeg.py`.

2.  `06_2_ITPC_save_ieeg.ipynb` Computation of sEEG phase angles. Outputs saved to `data/06_ieeg/plot_data_for_itpc/`. Uses `ITPC_ieeg.py`.

3.  `run_plot_itpc_rose_all_ieeg_06_3.m` Phase visualization workflow. Uses `plot_itpc_rose_cartesian_ieeg.m`.

4.  `06_4_ieeg.Rmd` Statistical validation in MTG-R. Tests whether MEG-derived features replicate in intracranial recordings. Uses `pearson_utils.R`.

Environment-specific files:

-   `.gitignore`
-   `environment.yaml`
-   `requirements.txt`

------------------------------------------------------------------------

# Execution Notes

-   All scripts are executed from the project root directory.
-   Sequential modules (e.g., 05_LMM and 06_ieeg) must be run in numerical order.
-   Outputs are written automatically to the corresponding `data/`, `results/`, and `figures/` directories.

------------------------------------------------------------------------
