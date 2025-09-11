# Construct 15-mer Sequences for Motif-x

This folder contains a Python script to generate **15-mer sequences** required for motif-x analysis of SUMOylation sites.

## Description
To extract over-represented sequence motifs, a short peptide window centred on the modified residue is required.  

  - For each SUMOylated lysine, a 15-mer peptide is created by extending 7 residues on each side.  
  - If the lysine is located less than seven residues from the protein N- or C-terminus, the sequence is padded with `"X"` to reach 15 amino acids.  
  - Redundant peptides are removed to keep only unique sequences.  

## How to run

1. Make sure the script and the `example_input/` folder are in the same directory.  
2. Open a terminal and change into that directory:  

```bash
 cd /path/to/folder
```
   
3.	Run the script directly. Output files will be written automatically into the output/ folder, created if it does not already exist:

```bash
python foreground_15-mer.py
```
