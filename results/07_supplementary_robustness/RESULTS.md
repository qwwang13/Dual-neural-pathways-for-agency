# Supplementary robustness results

All manuscript Bayes factors, slope differences and percentile confidence intervals reproduce at three decimal places.

## Bayesian t tests

JZS prior: standardized effect ~ Cauchy(0, sqrt(2)/2). BF10 compares the two-sided alternative with zero; BF+0 and BF-0 use the corresponding truncated and renormalized priors.

| Comparison | n | t | df | BF10 | BF+0 | BF-0 |
|---|---:|---:|---:|---:|---:|---:|
| Voluntary action vs zero | 19 | 1.726 | 18 | 0.826 | 1.552 | 0.099 |
| Involuntary action vs zero | 19 | 0.484 | 18 | 0.264 | 0.355 | 0.172 |
| Voluntary outcome vs zero | 19 | -4.924 | 18 | 258.260 | 0.059 | 516.460 |
| Involuntary outcome vs zero | 19 | -4.573 | 18 | 130.875 | 0.061 | 261.689 |
| Outcome: voluntary minus involuntary | 19 | -2.843 | 18 | 4.883 | 0.075 | 9.691 |

## Participant-level ITPC slope differences

Each condition is fitted separately as binding_ms = intercept + slope * ITPC. Differences are voluntary minus involuntary, in ms per unit ITPC. Each of the 20,000 draws resamples participants together with both conditions. Intervals are 95% percentile intervals.

| Region | Band | Voluntary slope | Involuntary slope | Difference | CI low | CI high | Bootstrap p |
|---|---|---:|---:|---:|---:|---:|---:|
| PreSMA-L | Theta | 313.366 | -61.143 | 374.509 | -134.163 | 869.900 | 0.121 |
| PreCG-R | Alpha | -361.534 | 259.470 | -621.004 | -1650.805 | 235.740 | 0.225 |
| MTG-R | Delta | 319.125 | 264.362 | 54.763 | -279.925 | 280.104 | 0.623 |
| MTG-R | Theta | -455.559 | 185.597 | -641.156 | -1505.191 | 581.905 | 0.387 |

The supplementary bootstrap p value is twice the smaller zero-crossing tail proportion, with a plus-one correction and a cap of 1. The manuscript reports the percentile confidence intervals. All four intervals include zero.

## Reproduction details

R: R version 4.4.3 (2025-02-28 ucrt)
R package dependencies: base and stats, supplied with R.
Bootstrap draws: NumPy PCG64, seed 20261109, shape 20000 x 19; included as one-based participant-row positions.
Largest absolute difference between the two Bayes-factor integrations: 5.4e-10.
Observed and sampled regression slopes were also checked with lm(), including its intercept.
Input MD5 checks (LF-normalized text):
- ../01_temporal_binding_effect/perceptual_shift_action_outcome.csv: a3b857f5c4b20dd204fa9f47bb284754
- participant_itpc.csv: 90517f4f1528e8044d608d103d70f5a8
- bootstrap_indices.csv: a548a6ee9828652724fde7d5040375fc

