# Preparing Peptide Search Databases

This folder contains Python scripts to perform **in-silico digestion** using different proteases.  

## Background
For datasets requiring up to eight missed cleavages, which exceeds Comet’s *in silico* digestion limit of five missed cleavages, peptide databases were generated. 

- **lysC_digest.py**  
  Simulate Lys-C digestion (cleaves **after K**, unless followed by **P**).
  
- **aspN_digest.py**  
  Simulate Asp-N + N-terminal Glu digestion (cleaves **before D/E**).

- **gluC_digest.py**  
  Simulate Glu-C digestion (cleaves **after D/E**).

## How to run

1. Make sure the script and the `example_input/` folder are in the same directory.  
2. Open a terminal and change into that directory:  
```bash
cd /path/to/folder
```
3. Run the scripts directly. All output files will be saved automatically into the `output/` folder created by each script if it doesn’t already exist:
```bash
python lysC_digest.py
python aspN_digest.py
python gluC_digest.py

```
