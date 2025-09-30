#!/usr/bin/env Rscript
# summarize_blast_hits_pigments.R
# --------------------------------
# Classify BLAST subjects as pigment-related (anthocyanin/flavonoid & carotenoid pathways + regulators),
# summarize counts, and write:
#   1) a filtered table of pigment-related hits
#   2) a "best hits" table (top 5 lowest e-values per query_id)
#
# Usage:
#   Rscript summarize_blast_hits_pigments.R \
#     --input path/to/blast_table.tsv \
#     --filtered-out path/to/blast_hits_pigment_related.csv \
#     --best-out path/to/blast_hits_pigment_candidates.tsv
#
# Example:
#   Rscript summarize_blast_hits_pigments.R \
#     --input results/blast_tables/lasius_blast_table.tsv \
#     --filtered-out results/pigment_hits/lasius_blast_hits_pigment_related.csv \
#     --best-out results/pigment_hits/lasius_blast_hits_pigment_candidates.tsv
#
# Input requirements:
#   TSV with columns at least: query_id, subject_id, evalue
#
# Notes:
#   - Only stdout gets summary prints; files get written to provided paths.
#   - Exclusion terms reduce false positives (e.g., ABA/xanthoxin/NCED, Xanthoceras).

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readr)
})

# ---- Small CLI parser -------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
help_msg <- paste0(
  "Usage:\n",
  "  Rscript summarize_blast_hits_pigments.R \\\n",
  "    --input path/to/blast_table.tsv \\\n",
  "    --filtered-out results/pigment_hits/lasius_blast_hits_pigment_related.csv \\\n",
  "    --best-out results/pigment_hits/lasius_blast_hits_pigment_candidates.tsv\n\n",
  "Required flags:\n",
  "  --input         Path to TSV with BLAST hits (query_id, subject_id, evalue required)\n",
  "  --filtered-out  Output CSV of pigment-related hits\n",
  "  --best-out      Output TSV of best pigment-related hits (top 5 per query by e-value)\n"
)

if (length(args) == 0 || any(args %in% c("-h", "--help"))) {
  cat(help_msg)
  quit(status = 0)
}

# Parse flags
get_flag <- function(flag) {
  i <- which(args == flag)
  if (length(i) == 1 && i < length(args)) args[i + 1] else NULL
}

in_path      <- get_flag("--input")
filtered_out <- get_flag("--filtered-out")
best_out     <- get_flag("--best-out")

if (is.null(in_path) || is.null(filtered_out) || is.null(best_out)) {
  cat(help_msg, file = stderr())
  quit(status = 1)
}

# ---- Read input -------------------------------------------------------------
if (!file.exists(in_path)) {
  stop("Input file not found: ", in_path)
}

# Read forgivingly; coerce later as needed
df <- read_tsv(in_path, show_col_types = FALSE, progress = FALSE)

required_cols <- c("query_id", "subject_id", "evalue")
missing <- setdiff(required_cols, names(df))
if (length(missing) > 0) {
  stop("Input missing required columns: ", paste(missing, collapse = ", "))
}

# ---- Regex term lists -------------------------------------------------------
# INCLUDE: enzymes, regulators, transport/decoration, metabolite names, carotenoids
include_terms <- c(
  # Flavonoid/anthocyanin enzymes
  "chalcone synthase\\b", "CHS\\b",
  "chalcone isomerase\\b", "CHI\\b", "CHIL\\b",
  "flavanone 3[- ]hydroxylase\\b", "\\bF3H\\b",
  "flavonoid 3'[,’']?[- ]hydroxylase\\b", "F3H'\\b", "F3'[- ]?H\\b", "\\bTT7\\b",
  "flavonoid 3'[,’']?5'[- ]hydroxylase\\b", "F3'?5'?H\\b",
  "flavonol synthase\\b", "\\bFLS\\b",
  "dihydroflavonol 4[- ]reductase\\b", "\\bDFR\\b",
  "anthocyanidin synthase\\b", "leucoanthocyanidin dioxygenase\\b", "\\bANS\\b", "\\bLDOX\\b",

  # Decoration / transport
  "UDP-?glucosyltransferase\\b", "glycosyltransferase\\b", "\\bUFGT\\b", "\\b3GT\\b", "\\b5GT\\b",
  "O-?methyltransferase\\b", "\\bOMT\\b",
  "glutathione S-?transferase\\b", "\\bGST\\b", "\\bTT19\\b",

  # Metabolite/compound terms
  "anthocyanin(s)?\\b", "anthocyanidin(s)?\\b",
  "pelargonidin\\b", "cyanidin\\b", "delphinidin\\b",
  "flavonoid(s)?\\b", "flavonol(s)?\\b",

  # Regulators (MBW complex)
  "R2R3-?MYB\\b", "\\bPAP1\\b", "MYB(75|90)\\b", "\\bSG6\\b",
  "\\bTT8\\b", "bHLH\\b",
  "\\bTTG1\\b", "WD40\\b",
  "TRANSPARENT TESTA\\b", "\\bTT\\d+\\b",

  # Carotenoid pathway
  "carotenoid(s)?\\b", "xanthophyll(s)?\\b",
  "\\bPSY\\b|phytoene synthase\\b",
  "\\bPDS\\b|phytoene desaturase\\b",
  "\\bZDS\\b|zeta-?carotene desaturase\\b",
  "\\bCRTISO\\b|carotenoid isomerase\\b",
  "lycopene cyclase\\b|\\bLCY[BE]\\b",
  "beta-?carotene hydroxylase\\b|\\bBCH\\b",
  "\\bZEP\\b|zeaxanthin epoxidase\\b",
  "violaxanthin de-epoxidase\\b|\\bVDE\\b"
)

# EXCLUDE: reduce false positives
exclude_terms <- c(
  "Xanthoceras"
)

include_re <- regex(paste0("\\b(", paste(include_terms, collapse = "|"), ")"), ignore_case = TRUE)
exclude_re <- regex(paste(exclude_terms, collapse = "|"), ignore_case = TRUE)

# ---- Classification ---------------------------------------------------------
df <- df %>%
  mutate(
    subject_id = as.character(subject_id),
    pigment_related = str_detect(subject_id, include_re) & !str_detect(subject_id, exclude_re),
    matched_pattern = if_else(
      pigment_related,
      str_extract(subject_id, include_re),
      NA_character_
    )
  )

# ---- Summaries to stdout ----------------------------------------------------
cat("Counts (pigment_related):\n")
print(table(df$pigment_related))

cat("\nNumber of pigment-related hits per query_id (top shown):\n")
pigment_counts <- df %>%
  filter(pigment_related) %>%
  count(query_id, name = "n_pigment_hits", sort = TRUE)
print(utils::head(pigment_counts, 20))

# ---- Write filtered table ---------------------------------------------------
# CSV (so it opens in Excel easily); no row names, no quotes
dir.create(dirname(filtered_out), showWarnings = FALSE, recursive = TRUE)
write.csv(df %>% filter(pigment_related),
          filtered_out, row.names = FALSE, quote = FALSE)

# ---- Best hits per query (lowest e-values) ----------------------------------
# Ensure numeric evalue
df$evalue <- suppressWarnings(as.numeric(df$evalue))

best_hits <- df %>%
  filter(pigment_related) %>%
  group_by(query_id) %>%
  arrange(evalue, .by_group = TRUE) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  left_join(pigment_counts, by = "query_id")

dir.create(dirname(best_out), showWarnings = FALSE, recursive = TRUE)
write.table(best_hits, best_out, sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n[OK] Wrote:\n")
cat(" - Filtered hits: ", normalizePath(filtered_out), "\n", sep = "")
cat(" - Best hits:     ", normalizePath(best_out), "\n", sep = "")
