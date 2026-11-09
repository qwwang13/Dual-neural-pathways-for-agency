# Bayesian robustness and participant-level ITPC analyses

This module reproduces the supplementary analyses reported in the revised manuscript. It contains the processed participant inputs, executable code and numerical results. All 19 MEG participants (IDs 2-20) are retained.

## Run the analyses

Open `07_robustness.R` in RStudio and click **Source**. Alternatively, from the repository root run:

```text
Rscript --vanilla analysis/run.R robustness
```

Only base R is needed. Paths are resolved relative to the repository. The supplied inputs are sufficient to run the analysis without regenerating the neural feature cache. The code has been checked with R 4.4.3.

The script writes two numerical tables and one readable report into `results/07_supplementary_robustness/`. Continuous statistics in the report have three decimal places; CSV outputs retain numerical precision. It verifies the included input checksums and checks that the manuscript's Bayes factors, slope differences and confidence limits reproduce at the reported precision.

## Bayesian robustness analysis

One-sample JZS Bayesian t tests compare each participant-level binding measure with zero. The paired outcome comparison uses each participant's voluntary-minus-involuntary outcome difference. The alternative prior on the standardized effect is Cauchy(0, sqrt(2)/2), conventionally reported as scale 0.707.

`BF10` is two-sided. `BF_positive0` and `BF_negative0` use the positive and negative halves of the same prior, renormalized to integrate to one. Bayes factors are evaluated by numerical integration; the two-sided result is independently checked using both an inverse-gamma mixture and a noncentral-t integral.

The reported outcome comparison is the between-condition outcome-binding difference (BF10 = 4.883). The code implements this specified contrast. Method: Rouder et al. (2009), https://doi.org/10.3758/PBR.16.2.225.

## Participant-level ITPC bootstrap

For each selected region and band, a separate linear regression with an estimated intercept is fitted in each condition:

```text
binding_ms = intercept + slope * ITPC
```

ITPC is kept on its original scale; behavioural binding remains in milliseconds. Behaviour comes from the original participant summaries rather than newly matched trial means. The comparison is the voluntary slope minus the involuntary slope. Each bootstrap draw samples 19 participant rows with replacement and carries both conditions together. The 95% interval uses the 2.5th and 97.5th percentiles (R quantile type 7).

The four manuscript comparisons are PreSMA-L theta (220-270 ms), PreCG-R alpha (180-240 ms), and MTG-R delta and theta (780-830 ms). The fixed set of 20,000 resampling rows is included, preserving the draws used for the manuscript's reported confidence intervals. These draws were generated with NumPy 2.3.5, PCG64 seed 20261109. Positions are one-based: position 1 refers to participant 2, and position 19 to participant 20.

The output also retains the previously calculated approximate bootstrap p value: twice the smaller proportion of differences at or beyond zero, with a plus-one correction and a cap of one. The current manuscript reports the confidence intervals.

## Input definitions

- `data/01_temporal_binding_effect/perceptual_shift_action_outcome.csv`: the existing processed participant behavioural summary, reused by this module. `va_ba` is voluntary action minus action baseline; `ia_bi` is involuntary movement minus its movement baseline; `vs_bs` and `is_bs` are voluntary and involuntary sound-outcome judgements minus sound baseline. Values are in milliseconds.
- `data/07_supplementary_robustness/participant_itpc.csv`: one processed ITPC value per participant, region, selected band and condition (152 rows). Values are recovered from the existing feature tables after verifying that they are constant across that participant's trials. `lower_ms` and `upper_ms` identify the analysis window. `PrecentralR` in the source is displayed as `PreCG-R`.
- `data/07_supplementary_robustness/bootstrap_indices.csv`: the 20,000 paired resampling rows, each containing 19 participant-row positions.
- `data/07_supplementary_robustness/source_files.csv`: relative paths and SHA256 hashes of the processed source files used to prepare the inputs.
- `data/07_supplementary_robustness/input_checksums.csv`: MD5 hashes of the three numerical inputs after CRLF-to-LF normalization, allowing either Git checkout convention.

The included inputs start at the processed participant-value stage. To recreate these files from the processed project rather than use the included copies, the optional preparation script requires Python with NumPy 2.3.5:

```text
python analysis/07_supplementary_robustness/prepare_inputs.py
```

It reads `data/01_temporal_binding_effect/perceptual_shift_action_outcome.csv` and `analysis/.cache/derived/05_LMM/features/` in this repository. Run `07_robustness.R` afterward. No source data are edited.

## Reproduction checks

The original behavioural R analysis and the original participant-level Python analysis were rerun from their processed source inputs. The standalone R outputs agreed across 64 numerical comparisons; the largest absolute difference was 3.2e-12. All 152 extracted ITPC values and all 380,000 bootstrap index entries matched the original analysis.

The repository entry point and direct Source execution were checked in fresh R sessions. A copy containing the supplied analysis inputs was also run without the neural feature cache, using CRLF text inputs as a checkout compatibility check. Both numerical result tables reproduced. These checks were completed on 9 November 2026.
