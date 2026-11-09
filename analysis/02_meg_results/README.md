# MEG permutation tests and Fig. 2 time courses

## Files

- `a01_permutation_tests.m`: supplied Brainstorm export for the action and outcome paired permutation tests. Input filenames, pairing order and numerical options are unchanged.
- `02_timecourses.Rmd`: plots the supplied condition means and standard errors and exports the combined Fig. 2 PDF.
- `meg_timecourse_plot_utils.R`: plotting functions used by the R report.

## Permutation-test inputs and settings

The input source maps are distributed as the 15 large archives listed in [the data folder](../../data/02_meg_results/meg/README.md). They contain 76 MAT files: 19 participant maps for each of `vab`, `iab`, `vsb` and `isb`. The exact within-participant pairing is preserved in `input_pairs.csv`.

The original export calls Brainstorm `process_test_permutation2p` for each binding component, using a paired Student t statistic, 1,000 randomizations, a two-sided alternative and the 0-1 s time window. Time samples are retained separately and no scout selection is requested. The export does not set a random seed.

## Using the Brainstorm export

Use MATLAB with Brainstorm and a protocol containing the corresponding group studies and cortical surface. The original protocol-relative paths are retained:

```text
Group_analysis/exp1_200hz_221127/  # action source maps
Group_analysis/exp1_200hz_221203/  # outcome source maps
```

After downloading and extracting the input archives, place/register the source maps in the matching Brainstorm studies, preserving the filenames and the paired order in `input_pairs.csv`. The saved statistics reference `@default_subject/tess_cortex_pial_low.mat` and 15,002 cortical vertices. Run the export manually in that protocol when a new permutation calculation is wanted:

```matlab
run('analysis/02_meg_results/a01_permutation_tests.m')
```

Brainstorm writes new results and its processing report to the active protocol. The supplied historical results are described in [results/02_meg_results](../../results/02_meg_results/README.md).

## FDR display and saved maps

The manuscript describes an alpha of 0.05 and FDR control across time. Brainstorm statistical thresholding is configured separately from the permutation-test export. Both supplied MAT results retain the pre-threshold `pmap` and `tmap` arrays with `Correction = 'no'`; the export does not record the Stat-panel correction-dimension selection. Keep the original Stat-panel settings with any exported thresholded maps or screenshots.

See the official [Brainstorm statistics tutorial](https://neuroimage.usc.edu/brainstorm/Tutorials/Statistics) for the Stat panel and the separate Apply statistic threshold process. The repository retains the supplied permutation export and saved maps; no FDR settings have been added to the historical export.

## Drawing Fig. 2

From the repository root, `Rscript analysis/run.R meg` reads the four `ta_*.csv` condition-summary tables and redraws the figure. The highlighted ROI lists are the supplied figure settings in `02_timecourses.Rmd`. This plotting step does not launch the Brainstorm permutation export.

This addition was organized by reading the script, archive listings and saved-result metadata. The MATLAB export, permutation tests and FDR calculations were not executed during organization.
