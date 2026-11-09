# MEG source maps for the Fig. 2 permutation tests

## Data files

`1.7z` through `15.7z` contain the 76 source-map MAT files referenced by the supplied Brainstorm script. There are 19 maps per condition: `vab` and `iab` for the action comparison, and `vsb` and `isb` for the outcome comparison. Original archive and MAT filenames are retained.

`archives.csv` lists archive sizes, member counts and SHA256 hashes. The archives total 12,579,638,717 bytes (approximately 12.58 GB). `input_pairs.csv` records the exact two ordered input lists for each paired test, including the archive containing each file. `pair_index` is the row position in the original paired lists, rather than a participant identifier inferred from the timestamp in a filename.

## Large-file distribution

These archives are prepared as separate GitHub Release assets. Upload all 15 archives to the repository's [Releases page](https://github.com/qwwang13/Dual-neural-pathways-for-agency/releases), keeping their current names. The local archives remain in this folder and are excluded from ordinary Git tracking; publishing the code does not upload the Release assets automatically.

After downloading, extract every archive into a local staging folder. These are separate archives with different members, so all 15 are needed. Their combined contents match all 76 script inputs. Use `input_pairs.csv` to place/register the files in the two original Brainstorm study paths described in [the analysis instructions](../../../analysis/02_meg_results/README.md).

The four parent-folder `ta_*.csv` tables contain the regional mean/standard-error summaries used for Fig. 2 plotting. The source-map archives are inputs to the original Brainstorm permutation tests.
