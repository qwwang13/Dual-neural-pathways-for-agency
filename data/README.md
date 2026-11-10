# Data

The supplied files are preprocessed inputs for the behavioural, MEG and intracranial analyses.

## Dataset organization

| Folder | Contents |
|---|---|
| `01_temporal_binding_effect` | One row per participant, containing baseline-corrected action and outcome perceptual shifts |
| `02_meg_results` | Regional mean/standard-error summaries, plus `meg/` source-map archive inventories and exact permutation input pairs |
| `03_functional_regression` | Participant-level neural and behavioural voluntary-minus-involuntary contrasts, plus mean curves for plotting |
| `05_LMM` | Baseline-corrected regional single-trial MEG epochs (`.csv.gz`) and matched behavioural responses |
| `06_ieeg` | MTG bipolar-channel epochs and behavioural responses for patients 2 and 4 |

## Identifiers and units

Single-trial tables are matched by participant, condition and trial identifier. Behavioural responses are in milliseconds. Neural time columns in the single-trial inputs are in seconds, whereas MEG summary and functional-regression time columns are in milliseconds.

| Field | Meaning |
|---|---|
| `subject` | MEG participant or intracranial patient identifier |
| `condition` | Experimental condition or contrast, as defined below |
| `trial` | Trial identifier within participant and condition |
| `region` | Brain-region identifier |
| `label` | Intracranial bipolar-channel identifier |
| `Time` or `time` | Time relative to the epoch reference event, in the units specified above |
| `Value` or `value` | Neural amplitude in the regional time-course or epoch tables |
| `behavior_value` | Trial-wise perceptual shift in milliseconds |
| `SE` | Between-participant standard error in the MEG summaries |

The MEG single-trial readers accept both `.csv` and `.csv.gz` files. Feature windows are read from the functional-regression coefficient tables.

Network analysis reads the participant difference time courses from `03_functional_regression`.

## Condition and variable names

The tables below describe how the condition names are used in the supplied analyses. `v` denotes voluntary and `i` denotes involuntary conditions; action and sound-outcome measures are used for action and outcome binding, respectively.

### Baseline conditions

| Code | Meaning |
|---|---|
| `ba` | Voluntary action-only baseline |
| `bi` | Involuntary movement-only baseline |
| `bs` | Sound-only baseline |

These codes also appear in filenames such as `PreSMAL_action_ba.csv.gz` and `MTGR_outcome_bs.csv.gz`. Action files contain voluntary trials corrected using `ba` and involuntary trials corrected using `bi`; both outcome conditions use `bs`.

### Single-trial MEG and behavioural labels

The `condition` column in `05_LMM` uses the following labels:

| Label | Analysis condition | Baseline used |
|---|---|---|
| `va_a` | Voluntary action binding | `ba` |
| `ia_a` | Involuntary action binding | `bi` |
| `vs_s` | Voluntary outcome binding | `bs` |
| `is_s` | Involuntary outcome binding | `bs` |

MEG amplitudes and behavioural perceptual shifts in these four conditions are baseline-corrected. Each behavioural trial has the participant's corresponding baseline-condition mean subtracted. Each MEG trial has the participant's baseline-condition mean waveform subtracted at the same region and time point.

### Participant summaries and MEG time courses

Behavioural participant summaries and MEG summaries use different names for the same four baseline-corrected conditions:

| Behavioural column | MEG condition | Meaning |
|---|---|---|
| `va_ba` | `vaa_ba` | Voluntary action binding: voluntary action condition minus `ba` |
| `ia_bi` | `iaa_bi` | Involuntary action binding: involuntary movement condition minus `bi` |
| `vs_bs` | `vss_bs` | Voluntary outcome binding: sound-outcome condition minus `bs` |
| `is_bs` | `iss_bs` | Involuntary outcome binding: sound-outcome condition minus `bs` |

The functional-regression and network contrasts are voluntary minus involuntary, after the corresponding baseline correction:

| Neural condition | Behavioural column | Contrast |
|---|---|---|
| `vaa_iaa` | `vaa_iaa_binding` | Action: `vaa_ba - iaa_bi` for MEG; `va_ba - ia_bi` for behaviour |
| `vss_iss` | `vss_iss_binding` | Outcome: `vss_bs - iss_bs` for MEG; `vs_bs - is_bs` for behaviour |

### Intracranial conditions

- `Voluntary_OutcomeOnsetReportSound`: voluntary condition, epochs referenced to outcome onset, with sound-time reports.
- `Involuntary_OutcomeOnsetReportSound`: involuntary condition, epochs referenced to outcome onset, with sound-time reports.

Intracranial features are extracted directly from the supplied patient epochs.

## Derived MEG feature tables

The preprocessed MEG epochs are stored in `data/05_LMM`. Feature extraction creates a separate output directory, `analysis/.cache/derived/05_LMM/features`. The name `.cache` is used here for generated intermediate outputs. These feature tables are derived from the supplied epochs using `analysis/05_LMM/01_features.py`.

Each regional feature CSV has one row per retained participant-condition-trial and contains `subject`, `condition`, `trial`, `meanAmplitude`, four `*_bandPower` columns and four `*_bandITPC` columns. ITPC is calculated across trials within each participant and condition, so the same ITPC value is repeated across that participant's trial rows. `windows_used.csv` records the brain regions and extraction windows. The tables are saved before condition-specific model standardization.

The LMM module reads these tables to fit the nine-feature models. The supplementary ITPC analysis reads each participant-condition ITPC value once and combines it with the original participant behavioural summary. Supplied feature tables are included in the repository. They can be regenerated by running:

```sh
python analysis/05_LMM/01_features.py
```

## Supplementary analysis inputs

The supplementary analyses reuse `01_temporal_binding_effect/perceptual_shift_action_outcome.csv` and the existing MEG feature tables under `analysis/.cache/derived/05_LMM/features`. Each participant-condition ITPC value is read once from those tables. Bootstrap draws are generated during analysis; no additional supplementary data folder is needed.

## Fig. 2 permutation-test source maps

The `02_meg_results/meg` folder contains the file-pairing list and instructions for 15 large input archives with 76 Brainstorm source-map MAT files. `input_pairs.csv` preserves the original ordering of the 19 paired maps in each comparison. The local binary archives are excluded from ordinary Git tracking and are prepared for separate Release distribution. See [the source-map README](02_meg_results/meg/README.md) for the layout.
