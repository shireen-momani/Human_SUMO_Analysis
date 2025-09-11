# Human SUMO Analysis

Scripts and example input files supporting the **human SUMO build manuscript**.
Organised into subfolders, each handling a different step of the analysis:

- **peptide_db/**

Python scripts for generating custom peptide databases. Includes in-silico digestions (Lys-C, Glu-C, Asp-N) allowing up to eight missed cleavages.

- **sumo_peak_filter/**

R script to filter SUMO remnant fragment ion peaks from mzML files. Removes diagnostic ions using the `mzR` package to improve PSM confidence and accurate SUMO site localisation.

- **site_level_filter/**

Python scripts for post mzidFLR processing. Collapse results to the site level, retain sites with FLR < 5%, remove contaminants/decoys, merge parallel searches, and apply Asp-N–specific rules to endogenous datasets.

- **gsb_tiering/**

Python script to classify SUMO sites into quality tiers. Defines **Gold** (≥2 datasets, FLR < 1%), **Silver** (1 dataset, FLR < 1%), and **Bronze** (FLR < 5% but not in Gold/Silver) sets.

- **15mer_motif/**

Python script to generate 15-mer sequences required for motif-x analysis. Used to extract statistically overrepresented sequence motifs around SUMOylated lysines.

- **figures/**

Scripts to generate figures. Includes motif heatmaps, amino acid enrichment plots, structural summaries, and co-modification analyses.
