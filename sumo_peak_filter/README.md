# Removing SUMO Remnant Fragment Peaks

This folder contains an R script that filters out diagnostic fragment-ion peaks arising from fragmentation of SUMO remnants in MS/MS spectra. 

The goal is to remove highly abundant SUMO remnant ions that add “clutter” to spectra. These ions are not informative for peptide sequencing but can dominate the signal, leading to reduced PSM confidence and less accurate SUMO site localisation.


## Input

- Takes `.mzML` files as input.  
- The provided script currently removes MS2 peaks matching diagnostic-ion m/z values (`DIlist`) for QQTGG/Pyro-QQTGG. 


## How to run

Run the scripts directly. All output files will be saved automatically into the `output/` folder created by each script if it doesn’t already exist:

```bash
Rscript rm_diagnostic_ions.R
```
