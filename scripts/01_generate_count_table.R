#!/usr/bin/env Rscript
# Build a per-pegRNA raw count table and reporter-editing table from
# per-read mapping files for one PE HCR-FlowFISH screen.
#
# For each replicate, three FACS bins are sequenced: a bottom-20%-expression
# bin ("b20"), a top-20%-expression bin ("t20"), and a bulk unsorted
# population used to measure library representation and reporter editing
# rate. Reads are already resolved to a single best-matching oligo_id
# (perfect match to a library pegRNA/reporter/barcode combination) upstream
# of this script, by 00_quantify_fastq_counts.py (spacer match in read 1 +
# barcode match in read 2, kept only if both agree).
#
# Usage:
#   Rscript 01_generate_count_table.R <sample_prefix> <n_replicates> <mapped_reads_dir> <reporter_editing_dir> <pegRNA_library.csv> <output_dir>
#
# Expects, for each replicate r in 1..n_replicates and bin in {b20, t20, bulk}:
#   <mapped_reads_dir>/<sample_prefix>_<bin>_r<r>_mapped_reads_perfect_matches.csv  (column: oligo_id)
#   <reporter_editing_dir>/<sample_prefix>_bulk_r<r>_reporter_editing.csv           (columns: oligo_id, total_counts, percent_editing)

suppressMessages(library(dplyr))
suppressMessages(library(data.table))

MIN_BULK_COUNT <- 10

args <- commandArgs(trailingOnly = TRUE)
sample_prefix <- args[1]
n_reps <- as.integer(args[2])
mapped_reads_dir <- args[3]
reporter_editing_dir <- args[4]
library_path <- args[5]
out_dir <- args[6]

count_bin_rep <- function(bin, rep) {
  path <- file.path(mapped_reads_dir, sprintf("%s_%s_r%d_mapped_reads_perfect_matches.csv", sample_prefix, bin, rep))
  col <- sprintf("%s_%s_r%d", sample_prefix, bin, rep)
  fread(path) %>%
    group_by(oligo_id) %>%
    summarise(!!col := n(), .groups = "drop")
}

# raw counts: full_join across bins within each replicate, filtered on bulk
# depth so library dropouts (rather than genuine sorting-driven depletion)
# are excluded before RPM normalization
rep_tables <- lapply(seq_len(n_reps), function(r) {
  bulk_col <- sprintf("%s_bulk_r%d", sample_prefix, r)
  count_bin_rep("b20", r) %>%
    full_join(count_bin_rep("t20", r), by = "oligo_id") %>%
    full_join(count_bin_rep("bulk", r), by = "oligo_id") %>%
    filter(.data[[bulk_col]] >= MIN_BULK_COUNT)
})
raw_counts <- Reduce(function(a, b) full_join(a, b, by = "oligo_id"), rep_tables)

sample_cols <- setdiff(colnames(raw_counts), "oligo_id")
total_counts <- colSums(raw_counts[, sample_cols], na.rm = TRUE)
rpm <- raw_counts %>%
  mutate(across(all_of(sample_cols), ~ .x / total_counts[cur_column()] * 1e6))

# average reporter editing rate per oligo across bulk replicates, filtered
# on read depth per replicate the same way as the count table
editing_reps <- lapply(seq_len(n_reps), function(r) {
  path <- file.path(reporter_editing_dir, sprintf("%s_bulk_r%d_reporter_editing.csv", sample_prefix, r))
  fread(path) %>%
    filter(total_counts >= MIN_BULK_COUNT) %>%
    select(oligo_id, percent_editing) %>%
    rename(!!sprintf("percent_editing_rep%d", r) := percent_editing)
})
editing <- Reduce(function(a, b) full_join(a, b, by = "oligo_id"), editing_reps) %>%
  mutate(average_percent_editing = rowMeans(across(starts_with("percent_editing_rep")), na.rm = TRUE)) %>%
  select(oligo_id, average_percent_editing)

# non-targeting/control oligos install no edit, so their reporter is always
# the reference sequence -- treat as 100% "editing" for downstream scaling
lib <- fread(library_path)
editing <- lib %>%
  filter(class %in% c("non-targeting")) %>%
  transmute(oligo_id = name, average_percent_editing = 100) %>%
  bind_rows(editing %>% filter(oligo_id %in% (lib %>% filter(class != "non-targeting") %>% pull(name))))

# each targeting pegRNA has a matched-control pegRNA installing the REF
# allele at the same site (docs/manuscript_figures/pegRNA_matched_control_schematic.png).
# multiple targeting pegRNAs can share one matched control, so duplicate the
# matched control's counts/editing rate once per targeting pegRNA that
# references it -- this avoids treating the same matched-control reads as
# non-independent replicates of each other in the downstream BEAN model
matched_key <- lib %>%
  filter(class != "matched_control", !is.na(matched_peg_id)) %>%
  select(name, matched_peg_id) %>%
  group_by(matched_peg_id) %>%
  mutate(matched_duplicate_id = if (n() > 1) paste0(matched_peg_id, ".", row_number()) else matched_peg_id) %>%
  ungroup()

duplicated_counts <- matched_key %>%
  filter(matched_duplicate_id != matched_peg_id) %>%
  left_join(raw_counts, by = c("matched_peg_id" = "oligo_id")) %>%
  select(-name, -matched_peg_id) %>%
  rename(oligo_id = matched_duplicate_id)

full_counts <- bind_rows(raw_counts, duplicated_counts) %>%
  mutate(across(all_of(sample_cols), ~ replace(.x, is.na(.x), 0)))

duplicated_editing <- matched_key %>%
  filter(matched_duplicate_id != matched_peg_id) %>%
  left_join(editing, by = c("matched_peg_id" = "oligo_id")) %>%
  select(oligo_id = matched_duplicate_id, average_percent_editing)

full_editing <- bind_rows(editing, duplicated_editing)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
fwrite(full_counts %>% arrange(oligo_id), file.path(out_dir, sprintf("%s_counts.csv", sample_prefix)))
fwrite(rpm %>% arrange(oligo_id), file.path(out_dir, sprintf("%s_RPM.csv", sample_prefix)))
fwrite(full_editing %>% arrange(oligo_id), file.path(out_dir, sprintf("%s_average_editing.csv", sample_prefix)))

cat("pegRNAs passing bulk-depth filter:", nrow(raw_counts), "\n")
cat("pegRNAs after matched-control duplication:", nrow(full_counts), "\n")
