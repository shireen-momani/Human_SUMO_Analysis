# SUMO-Oracle

Predicts SUMOylation sites on lysine residues in human proteins from AlphaFold structures. A Random Forest model scores every lysine in each protein using 20 structural and sequence-based features within a ±2-residue window.



## How it works

For each structure in the input folder, `sumo_site_predictor.py` does the following:

1. **Sequence extraction.** The protein sequence is read from the `.cif` or `.pdb` file. Large proteins (over 2700 residues) are split by AlphaFold into multiple overlapping fragment files, and concatenating their sequences is not straightforward, so for these the full sequence is instead read from the provided FASTA file (`multifrag_prot_idmapping_2026_06_19.fasta`). Fragment filenames are matched against `human_AF_multifragment.txt`, which lists all fragment files for the 215 large proteins in the human proteome.

2. **Sequence based features.** All lysine positions are identified and a 5-mer window centred on each lysine (2 residues either side) is extracted. Selected AAindex1 physicochemical and biochemical descriptors (listed in `model_metadata.json`) are mapped onto positions within the window (13 per-position features) and summed across all 5 positions (4 summed features). Metapredict v3 is then run to compute per-residue intrinsic disorder scores, and the score at each lysine position is merged into the feature matrix.

3. **Structure based features.** Alpha-helix coordinates are read from the DSSP annotation files (`*_dssp.cif`) in the `dssp/` folder. Where overlapping helices occur because of overlapping AlphaFold fragments, only the longest helix covering that region is kept. The minimum distance from each lysine to the nearest helix boundary is added as `dist_nearest_helix`. Absolute solvent accessibility is read from the OpenStructure output CSVs in the `asa/` folder; the files are combined, duplicate rows from overlapping fragments are removed, and the `absolute_asa` value for each lysine is added to the feature matrix.

4. **Prediction.** Missing values are imputed with the training medians (`train_medians.joblib`), then the trained Random Forest (`random_forest.joblib`) assigns each lysine a SUMOylation probability. Two binary calls are derived from it: `sumo_f1_prediction` (threshold 0.452, F1-optimised, higher sensitivity) and `sumo_mcc_prediction` (threshold 0.602, MCC-optimised, higher confidence), where 1 = predicted SUMOylation site and 0 = not predicted. The full feature matrix, probability and both calls are written to a single output CSV.

## Requirements

Python 3.10 to 3.12 (scikit-learn is pinned to 1.3.2 to match the version the model was trained with, and it has no wheels for newer Python versions).

```
pip install -r requirements.txt
```

Metapredict installs its own dependencies (including PyTorch), so the first install can take a few minutes. See https://metapredict.readthedocs.io/en/latest/getting_started.html if you run into trouble.

## Usage

```
python sumo_site_predictor.py -i structures/ -o sumoylation_site_predictions.csv
```


## Output

One row per lysine. The first columns identify the site (`Protein`, `K-pos`, `5-mer`), followed by the feature values, then:

- `sumoylation_probability` raw probability from the Random Forest
- `sumo_f1_prediction` binary call at the F1-optimised threshold (0.452)
- `sumo_mcc_prediction` binary call at the MCC-optimised threshold (0.602)

Use the F1 call when you want sensitivity (more true sites recovered, more false positives) and the MCC call when you want confidence (fewer false positives, some true sites missed).