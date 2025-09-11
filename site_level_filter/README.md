# Site-Based Merging and Filtering

This folder contains Python scripts for collapsing **peptidoform-site mzidFLR results** to the **site level** and applying filtering strategies to reduce false positives.  



## Scripts

- **`single_site_flr.py`**  
	Collapses `mzidFLR` results at the site level and filters by FLR threshold (default 5%).    

- **`merged_site_flr.py`**  
	Merges results from datasets with multiple SUMO footprints or digestion strategies before applying site-level FLR filtering.  

- **`endogenous_site_flr.py`**  
  Applies Asp-N–specific rules for endogenous datasets, filtering C-terminal lysine SUMO sites unless preceded by **D/E**.  


## How to run

1. Make sure the script and the `example_input/` folder are in the same directory.  
2. Open a terminal and change into that directory: 
 
```bash
cd /path/to/folder
```
3. Run the scripts directly. All output files will be saved automatically into the `output/` folder created by each script if it doesn’t already exist:

```bash
python single_site_flr.py
python merged_site_flr.py
python endogenous_site_flr.py
```
