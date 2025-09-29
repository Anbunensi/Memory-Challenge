# generate_ImmuCC_merged_signatures.R
# Download 9 ImmuCC tissue signature CSVs, merge them (union of genes),
# fill missing with 0, and save two outputs:
# - ImmuCC_merged_signature_symbols.txt
# - ImmuCC_merged_signature_ensembl_GRCm39.txt
#
# Usage: run in R (install packages readr,dplyr,purrr,biomaRt if needed)

library(readr)
library(dplyr)
library(tidyr)
library(purrr)

out_symbol_file <- "data/ImmuCC_merged_signature_symbols.txt"
out_ensembl_file <- "data/ImmuCC_merged_signature_ensembl_GRCm39.txt"
fill_missing_with <- 0
normalize_columns <- FALSE

files_raw <- c(
  BoneMarrow = "https://raw.githubusercontent.com/wuaipinglab/ImmuCC/master/tissue_immucc/SignatureMatrix/BoneMarrow.sig.matrix.csv",
  Kidney = "https://raw.githubusercontent.com/wuaipinglab/ImmuCC/master/tissue_immucc/SignatureMatrix/Kidney.sig.matrix.csv",
  Liver = "https://raw.githubusercontent.com/wuaipinglab/ImmuCC/master/tissue_immucc/SignatureMatrix/Liver.sig.matrix.csv",
  Lung = "https://raw.githubusercontent.com/wuaipinglab/ImmuCC/master/tissue_immucc/SignatureMatrix/Lung.sig.matrix.csv",
  MammaryGland_Pregnancy = "https://raw.githubusercontent.com/wuaipinglab/ImmuCC/master/tissue_immucc/SignatureMatrix/MammaryGland.Pregnancy.sig.matrix.csv",
  Muscle = "https://raw.githubusercontent.com/wuaipinglab/ImmuCC/master/tissue_immucc/SignatureMatrix/Muscle.sig.matrix.csv",
  PeripheralBlood = "https://raw.githubusercontent.com/wuaipinglab/ImmuCC/master/tissue_immucc/SignatureMatrix/PeripheralBlood.sig.matrix.csv",
  SmallIntestine = "https://raw.githubusercontent.com/wuaipinglab/ImmuCC/master/tissue_immucc/SignatureMatrix/SmallIntestine.sig.matrix.csv",
  Spleen = "https://raw.githubusercontent.com/wuaipinglab/ImmuCC/master/tissue_immucc/SignatureMatrix/Spleen.sig.matrix.csv"
)

read_and_prefix <- function(url, prefix) {
  message("Downloading: ", url)
  df <- tryCatch(read_csv(url, show_col_types = FALSE),
                 error = function(e) tryCatch(read_tsv(url, show_col_types = FALSE),
                                              error = function(e2) stop("Cannot read: ", url)))
  colnames(df)[1] <- "Gene"
  df <- df %>% group_by(Gene) %>% summarise(across(everything(), ~ mean(as.numeric(.), na.rm = TRUE))) %>% ungroup()
  value_cols <- setdiff(colnames(df), "Gene")
  new_names <- paste0(prefix, "__", value_cols)
  colnames(df)[match(value_cols, colnames(df))] <- new_names
  return(df)
}

sig_list <- imap(files_raw, ~ read_and_prefix(.x, .y))
merged <- reduce(sig_list, full_join, by = "Gene")
value_cols <- setdiff(colnames(merged), "Gene")
merged <- merged %>% select(Gene, all_of(sort(value_cols)))

if (is.na(fill_missing_with)) {
  merged_filled <- merged
} else {
  merged_filled <- merged %>% mutate(across(-Gene, ~ replace_na(as.numeric(.), fill_missing_with)))
}

if (normalize_columns) {
  vcols <- setdiff(colnames(merged_filled), "Gene")
  merged_filled <- merged_filled %>% mutate(across(all_of(vcols), ~ {
    x <- as.numeric(.); s <- sum(x, na.rm = TRUE)
    if (is.na(s) || s == 0) return(x); x / s
  }))
}

# Ensure data/ exists
if (!dir.exists("data")) dir.create("data")
write.table(merged_filled, out_symbol_file, sep = "\t", quote = FALSE, row.names = FALSE)
message("Saved symbol file: ", out_symbol_file)

# Map to Ensembl (GRCm39) using biomaRt
suppressMessages({
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    message("biomaRt not installed; skipping Ensembl mapping. Install biomaRt to generate Ensembl file.")
  } else {
    library(biomaRt)
    mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    symbols <- merged_filled$Gene
    symbols <- symbols[!is.na(symbols) & symbols != ""]
    map <- getBM(attributes = c("external_gene_name", "ensembl_gene_id"),
                 filters = "external_gene_name", values = symbols, mart = mart)
    if (nrow(map) == 0) {
      message("No mappings found; check gene symbols.")
    } else {
      map_u <- map %>% group_by(external_gene_name) %>% slice(1) %>% ungroup()
      merged_map <- merged_filled %>% inner_join(map_u, by = c("Gene" = "external_gene_name"))
      merged_map <- merged_map %>% select(ensembl_gene_id, everything(), -Gene) %>% rename(Gene = ensembl_gene_id)
      write.table(merged_map, out_ensembl_file, sep = "\t", quote = FALSE, row.names = FALSE)
      message("Saved Ensembl file: ", out_ensembl_file)
    }
  }
})