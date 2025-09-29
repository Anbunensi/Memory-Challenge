# ImmuCC Merged Signature Generation

This repository contains a script for generating a merged signature from the ImmuCC dataset.

## Files

- `generate_ImmuCC_merged_signatures.R`: An R script that downloads the 9 ImmuCC signature CSVs, merges them (union of genes), fills missing values with 0, optionally normalizes the data, maps gene symbols to Ensembl (GRCm39) using the biomaRt package, and writes two outputs:
  - `ImmuCC_merged_signature_symbols.txt`: Contains the merged signature with gene symbols.
  - `ImmuCC_merged_signature_ensembl_GRCm39.txt`: Contains the merged signature with Ensembl gene IDs.

## How to Run the Script

1. Ensure you have R installed on your system.
2. Install the required packages if you haven't already:
   ```r
   install.packages("biomaRt")
   install.packages("dplyr")
   ```
3. Run the script in R:
   ```r
   source("data/generate_ImmuCC_merged_signatures.R")
   ```
4. Check the output files for the merged signatures.