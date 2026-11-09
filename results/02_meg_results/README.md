# Saved MEG permutation results

`permutation_results.7z` preserves the supplied archive of the original Brainstorm export and its two saved statistical maps. It is a local large-file asset for the repository's [GitHub Releases](https://github.com/qwwang13/Dual-neural-pathways-for-agency/releases), and is excluded from ordinary Git tracking.

| Comparison | Member inside the archive | Contents |
|---|---|---|
| Action: voluntary vs involuntary | `output/presults_no_221128_0151.mat` | Paired permutation-test p and t maps |
| Outcome: voluntary vs involuntary | `output/presults_no_221203_1536.mat` | Paired permutation-test p and t maps |

Each MAT file is 288,321,664 bytes. `pmap` and `tmap` each contain 15,002 source rows and 1,201 time samples over 0-1 s. `Options.TestType` is `ttest_paired`; `Options.TestTail` is `two`. The saved `History` lists match the 19 paired inputs for each condition in the accompanying export.

Both files store `Correction = 'no'`, retaining statistical maps before Brainstorm display thresholding. The manuscript's FDR settings concern the separate statistical display/thresholding stage; the supplied export does not record the Stat-panel dimension selection. See [the analysis README](../../analysis/02_meg_results/README.md).

The archive also contains the original `script/permutation t test.m`. An identical copy is provided under the public filename `analysis/02_meg_results/a01_permutation_tests.m`. The saved statistical values have not been recomputed or modified.
