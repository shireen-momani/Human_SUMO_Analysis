# Gold–Silver–Bronze Tiering System

This folder contains a Python script that classifies identified SUMO sites based on the frequency of observation across independent reanalyses and the FLR threshold applied.

## Background

SUMO sites identified < 5% FLR were combined and classified into three sets, with different levels of quality:  

- **Gold**: Sites identified in **two or more datasets** with an **FLR < 1%**.  
- **Silver**: Sites found in **one dataset** with an **FLR < 1%**.  
- **Bronze**: Sites not meeting the Gold or Silver criteria, but with an **FLR < 5%**.
This tiering system highlights high-confidence SUMOylation sites while retaining broader coverage at a controlled FLR.

## How to run

1. Make sure the script and the `example_input/` folder are in the same directory.  
2. Open a terminal and change into that directory:  

```bash
cd /path/to/folder
```
   
3.	Run the script directly. Output files will be written automatically into the output/ folder, created if it does not already exist:

```bash
python gsb_tiering.py
```
