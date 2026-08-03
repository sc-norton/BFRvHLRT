# === PROJECT: BFR v HLRT Proteomics ==============================
# === Script: Normalization

# === PURPOSE: Load raw data and metadata, assemble a DAList for use in
# proteoDA, filter the data for contaminants and missingness, detect and
# remove outlier samples, and normalize for downstream analyses.

# === INPUTS: raw data, meta data, HPA muscle data

# === OUTPUTS: raw sample summary, normalization reports, flagged
# contaminants, pre & post QC reports, filtering effects, outlier
# diagnostics, normalized data (csv), normalized DAList, report_data.rds
# (consumed by 01_normalization_BFR_reports.R)

# --- 0: Setup --------------------------------------------------------
cat("\n>> 0 - Setup\n")
cat("\n Loading required packages...")

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(proteoDA, readxl, readr, dplyr, tidyr, stringr, purrr, ggplot2,
               ggrepel, patchwork, tibble)

cat("\n Packages successfully loaded.\n")

setwd(rprojroot::find_rstudio_root_file())
getwd()

# Configuration
config <- list(
  raw_data   = "00_input/BFR_raw.csv",
  meta_data  = "00_input/BFR_meta.csv",
  hpa_muscle = "00_input/HPA_muscle.csv",
  
  report_dir = "01_normalization/b_reports",
  data_dir   = "01_normalization/c_data",
  
  # Parameters
  miss_pct    = 50,
  mahal_p     = 0.01,
  mad_k       = 3,
  norm_method = "cycloess"
)

dir.create(config$report_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(config$data_dir,   showWarnings = FALSE, recursive = TRUE)

# --- Helper functions (shared with reports script) --------------------

compute_missingness_stats <- function(d) {
  data <- d$data
  list(
    protein_missing = tibble(uniprot_id = rownames(data),
                             pct_missing = rowMeans(is.na(data)) * 100),
    sample_missing  = tibble(sample = colnames(data),
                             n_missing = colSums(is.na(data)))
  )
}

compute_pca <- function(data, log_transform = TRUE) {
  data_for_pca <- data
  for (j in seq_len(ncol(data_for_pca))) {
    nas <- is.na(data_for_pca[, j])
    if (any(nas)) data_for_pca[nas, j] <- median(data_for_pca[, j], na.rm = TRUE)
  }
  if (log_transform) data_for_pca <- log2(data_for_pca + 1)
  prcomp(t(data_for_pca), center = TRUE, scale. = TRUE)
}

log_step <- function(log, step, n_before, n_after) {
  bind_rows(log, tibble(
    step      = step,
    n_before  = n_before,
    n_after   = n_after,
    n_removed = n_before - n_after
  ))
}

cat("\n Base directory: ", getwd(), "\n")

# --- 1: Load & Restructure Data ---------------------------------------
cat("\n>> 1 - Load & Restructure Data\n")

raw  <- read.csv(file.path(config$raw_data))
meta <- read.csv(file.path(config$meta_data))

# Flag cRAP contaminants (searches the full Protein.Group string, before
# the multi-match trimming below potentially drops a cRAP tag past the ";")
raw <- raw %>%
  mutate(is_crap = str_detect(Protein.Group, "cRAP")) %>%
  relocate(is_crap, .after = N.Proteotypic.Sequences)

cat("\n", sum(raw$is_crap), "protein group(s) flagged as cRAP contaminants.\n")

# Remove multiple matches
raw <- raw %>%
  mutate(Protein.Group = trimws(sub(";.*", "", Protein.Group)),
         Protein.Names = trimws(sub(";.*", "", Protein.Names)),
         Genes         = trimws(sub(";.*", "", Genes)))

# Rename annotation columns
raw <- raw %>%
  rename(uniprot_id = Protein.Group,
         protein     = Protein.Names,
         gene        = Genes,
         description = First.Protein.Description,
         n_seq       = N.Sequences,
         n_pro_seq   = N.Proteotypic.Sequences)

# Clean sample column names: leading X -> S, CL_MR_ -> S, "." -> "_"
sample_col_idx <- 8:138
names(raw)[sample_col_idx] <- names(raw)[sample_col_idx] %>%
  sub("^X", "S", .) %>%
  sub("^CL_MR_", "S", .) %>%
  gsub("\\.", "_", .)

# --- 1b: Flag reruns, then finalize sample column names ----------------
cat("\n>> 1b - Flag reruns\n")

sample_cols <- names(raw)[sample_col_idx]

# Rerun marker = one or more trailing "r"/"R" characters at the end of the
# (already cleaned) column name.
#   rerun_level 0 -> original measurement, not a rerun
#   rerun_level 1 -> rerun
#   rerun_level 2 -> rerun-of-a-rerun (e.g. a name ending "...T2rR")
rerun_suffix <- str_extract(sample_cols, "[rR]+$")
rerun_level  <- ifelse(is.na(rerun_suffix), 0, nchar(rerun_suffix))
base_sample  <- sub("[rR]+$", "", sample_cols)

# Rename raw's sample columns to the finalized base name (rerun marker dropped)
names(raw)[sample_col_idx] <- base_sample

rerun_info <- tibble(
  sample   = base_sample,
  rerun_level = rerun_level,
  is_rerun    = rerun_level > 0
)

cat("\n", sum(rerun_info$is_rerun), "rerun sample column(s) identified.\n")

# --- 1c: Update metadata -------------------------------------------------
cat("\n>> 1c - Update metadata\n")

# "sample" holds the full S-prefixed sample identifier throughout (matches
# the raw data's cleaned column names and the reporting helper functions).
meta <- meta %>%
  mutate(sample = paste0("S", sample)) %>%
  left_join(rerun_info, by = "sample")

# QC: every metadata row should have matched a raw data column
if (any(is.na(meta$is_rerun))) {
  warning("Some metadata samples did not match a raw data column: ",
          paste(meta$sample[is.na(meta$is_rerun)], collapse = ", "))
}

# QC: every raw sample column should have a matching metadata row
unmatched_raw <- setdiff(names(raw)[sample_col_idx], meta$sample)
if (length(unmatched_raw) > 0) {
  warning("Some raw data columns have no matching metadata row: ",
          paste(unmatched_raw, collapse = ", "))
}

annot_cols <- c("uniprot_id", "protein", "gene", "description", "n_seq",
                "n_pro_seq", "is_crap")

intensity   <- raw[, setdiff(names(raw), annot_cols)]
intensity.b <- intensity
rownames(intensity.b) <- make.unique(raw$gene)
annotation  <- raw[, annot_cols]

cat(sprintf("\n Raw data: %d proteins x %d samples\n", nrow(raw), ncol(intensity)))

rownames(meta) <- meta$sample
meta <- meta[order(meta$sample), ]
intensity   <- intensity[, meta$sample]
intensity.b <- intensity.b[, meta$sample]

stopifnot("Sample mismatch" = setequal(colnames(intensity), meta$sample))
cat("   Validation passed:", ncol(intensity), "samples aligned.\n")

sample_summary <- as.data.frame(meta %>% count(group))
write.csv(sample_summary, file.path(config$data_dir, "00_raw_sample_summary.csv"),
          row.names = FALSE)
write.csv(raw,        file.path(config$data_dir, "01_raw_data_annotated.csv"), row.names = FALSE)
write.csv(annotation, file.path(config$data_dir, "02_annotation_data.csv"),    row.names = FALSE)

# --- 2a: Contaminant Assessment ----------------------------------
cat("\n>> 2 - Contaminant Assessment\n")

# --- A: Human Protein Atlas - Skeletal Muscle
HPAm <- read.csv(file.path(config$hpa_muscle))

hpa_lookup <- HPAm %>%
  dplyr::select(Gene, Ensembl, Uniprot) %>%
  distinct(Uniprot, .keep_all = TRUE)

annotation <- annotation %>%
  left_join(hpa_lookup, by = c("uniprot_id" = "Uniprot")) %>%
  mutate(in_hpa_muscle = uniprot_id %in% hpa_lookup$Uniprot)

n_flagged_hpa <- sum(!annotation$in_hpa_muscle)
cat(sprintf("\n Proteins flagged as contaminants (HPA): %d\n", n_flagged_hpa))

# --- B: Blood Protein Blacklist ---
contaminant_immunoglobulins <- annotation %>%
  filter(grepl("^IGH|^IGK|^IGL", gene)) %>%
  pull(gene)

contaminant_keratins <- annotation %>%
  filter(grepl("^KRT", gene)) %>%
  pull(gene)

blood_blacklist <- c(
  # Albumin & transport
  "ALB", "TTR", "APOA1", "APOA2", "APOA4", "APOB", "APOC1", "APOC2", "APOC3",
  "APOC4", "APOD", "APOE", "APOH", "APOL1", "APOL2", "APOM", "CLU", "TF", "GC",
  # Coagulation
  "FGA", "FGB", "FGG", "F2", "F5", "F7", "F10", "F13A1", "F13B", "SERPINC1",
  "SERPIND1", "PROS1", "PROC", "VWF", "PLG", "KNG1",
  # Complement System
  "C1QA", "C1QB", "C1QC", "C3", "C4A", "C4B", "C4BPA", "C4BPB", "C5", "C6",
  "C7", "C8A", "C8B", "C8G", "C9", "CFH", "CFB", "CFI", "SERPING1",
  # True Hemoglobins (explicitly named to avoid breaking growth factors)
  "HBA1", "HBA2", "HBB", "HBD", "HBG1", "HBG2", "HBM", "HBZ", "HBE1", "HBQ1",
  # Acute phase / Inter-alpha-trypsin inhibitors
  "CRP", "SAA1", "SAA2", "ORM1", "ORM2", "HP", "HPX", "A1BG", "A2M", "ITIH1",
  "ITIH2", "ITIH3", "ITIH4", "SERPINA1",
  # Platelet / Neutrophil markers
  "PPBP", "PF4", "THBS1", "SELP", "LYZ", "DEFA3", "JCHAIN",
  # Others
  "AHSG", "AHSP",
  # Dynamic vectors
  contaminant_immunoglobulins,
  contaminant_keratins
) %>% unique()

annotation <- annotation %>%
  mutate(in_blood_blacklist = gene %in% blood_blacklist)

flagged_hpa <- annotation %>%
  filter(!in_hpa_muscle) %>%
  mutate(flag_source = "HPA")

flagged_blood <- annotation %>%
  filter(in_blood_blacklist) %>%
  mutate(flag_source = "Blood blacklist")

cat(sprintf("  Blood blacklist: %d proteins to remove\n", sum(annotation$in_blood_blacklist)))

# --- 2b: Contamination Index ----------------------------
cat("\n>> 2b - Contamination index\n")

blood_markers  <- c("ALB", "HBB", "HBA1")
muscle_markers <- c("MB", "CKM")
all_markers    <- c(blood_markers, muscle_markers)

marker_rows <- which(rownames(intensity.b) %in% all_markers)

if (length(marker_rows) > 0) {
  blood_intensity  <- colSums(intensity.b[rownames(intensity.b) %in% blood_markers, ],  na.rm = TRUE)
  muscle_intensity <- colSums(intensity.b[rownames(intensity.b) %in% muscle_markers, ], na.rm = TRUE)
  total_intensity  <- colSums(intensity.b, na.rm = TRUE)
  blood_pct        <- (blood_intensity / total_intensity) * 100
  muscle_pct       <- (muscle_intensity / total_intensity) * 100
  B_M_ratio        <- (blood_pct / muscle_pct)
  
  marker_pct_df <- as.data.frame(t(intensity.b[marker_rows, ])) %>%
    mutate(across(everything(), ~ . / total_intensity * 100)) %>%
    rename_with(~ paste0(.x, "_pct")) %>%
    mutate(sample  = colnames(intensity.b),
           blood_pct  = blood_pct,
           muscle_pct = muscle_pct,
           B_M_ratio  = B_M_ratio)
  
  cat(sprintf("  Blood %% of total intensity — min: %.2f%%, median: %.2f%%, max: %.2f%%\n",
              min(blood_pct), median(blood_pct), max(blood_pct)))
  cat(sprintf("  Muscle %% of total intensity — min: %.2f%%, median: %.2f%%, max: %.2f%%\n",
              min(muscle_pct), median(muscle_pct), max(muscle_pct)))
} else {
  stop("Markers not found in rownames of intensity.b — check gene names.")
}

meta <- meta %>% left_join(marker_pct_df, by = "sample")
rownames(meta) <- meta$sample

# --- 2c: Blood Score (correlation) ----------------------------
blood_score <- colMeans(intensity.b[rownames(intensity.b) %in% blood_markers, ], na.rm = TRUE)

cor_with_blood <- apply(intensity.b, 1, function(x) cor(x, blood_score, use = "pairwise"))
cor_df <- data.frame(gene = rownames(intensity.b), cor = cor_with_blood)

top_genes <- cor_df %>%
  arrange(desc(abs(cor))) %>%
  dplyr::slice(1:50) %>%
  pull(gene)

# --- 2d: Top 5 Abundance proteins per sample ---------------------
top5_per_sample <- apply(intensity.b, 2, function(x) {
  idx <- order(x, decreasing = TRUE, na.last = TRUE)[1:5]
  tibble(gene      = rownames(intensity.b)[idx],
         intensity = x[idx],
         rank      = 1:5)
}) %>%
  bind_rows(.id = "sample")

# --- 3: Assemble DAList -----------------------------------------
cat("\n>> 3 - Assemble DAList \n")

cat(" Checking for duplicate uniprot ids...\n")
dup_ids <- annotation$uniprot_id[duplicated(annotation$uniprot_id)]

if (length(dup_ids) > 0) {
  cat(sprintf("   Found %d duplicate uniprot_ids — deduplicating by highest mean intensity.\n",
              length(dup_ids)))
  
  intensity_num <- as.data.frame(lapply(intensity, as.numeric))
  annotation$row_mean <- rowMeans(intensity_num, na.rm = TRUE)
  
  keep_idx <- annotation %>%
    mutate(row_idx = row_number()) %>%
    group_by(uniprot_id) %>%
    slice_max(row_mean, n = 1, with_ties = FALSE) %>%
    pull(row_idx)
  
  annotation <- annotation[keep_idx, ]
  intensity  <- intensity[keep_idx, ]
  annotation$row_mean <- NULL
  
  cat(sprintf("   After deduplication: %d unique proteins\n", nrow(annotation)))
} else {
  cat("   No duplicate uniprot_ids found.\n")
}

dal <- DAList(
  data       = intensity,
  metadata   = meta,
  annotation = annotation
)

# Named list of one -> keeps this script's downstream code (and the
# reports script) identical in shape to a multi-comparison pipeline,
# in case a second BFR comparison is added later.
dals          <- list(BFR = dal)
grouping_cols <- c(BFR = "group")

cat(sprintf(" BFR DAList: %d proteins x %d samples\n", nrow(dal$data), ncol(dal$data)))

# --- 3b: Preserve raw (pre-filter) data for missingness/PCA reporting -----
cat("\n>> 3b - Snapshot raw data for reporting\n")
dals_raw <- dals

missingness_raw <- lapply(dals_raw, compute_missingness_stats)
pca_raw <- lapply(dals_raw, function(d) compute_pca(d$data, log_transform = TRUE))

# --- 4: Filter Data -------------------------------------
cat("\n>> 4 - Filter Data \n")

f.dals <- lapply(names(dals), function(nm) {
  d   <- dals[[nm]]
  grp <- grouping_cols[[nm]]
  
  log     <- tibble(step = character(), n_before = integer(),
                    n_after = integer(), n_removed = integer())
  n_start <- nrow(d$data)
  
  # Step 1: Convert 0s to NA
  cat(sprintf(" [%s] Step 1: Convert 0s to NA\n", nm))
  d <- zero_to_missing(d)
  
  # Step 2: HPA muscle filter
  cat(sprintf(" [%s] Step 2: HPA muscle contaminant removal\n", nm))
  n_before <- nrow(d$data)
  d        <- filter_proteins_by_annotation(d, in_hpa_muscle)
  n_after  <- nrow(d$data)
  log      <- log_step(log, "HPA muscle filter", n_before, n_after)
  cat(sprintf("   %d -> %d (%d removed)\n", n_before, n_after, n_before - n_after))
  
  # Step 3: Blood blacklist filter
  cat(sprintf(" [%s] Step 3: Blood blacklist filter\n", nm))
  n_before <- nrow(d$data)
  d        <- filter_proteins_by_annotation(d, !in_blood_blacklist)
  n_after  <- nrow(d$data)
  log      <- log_step(log, "Blood blacklist filter", n_before, n_after)
  cat(sprintf("   %d -> %d (%d removed)\n", n_before, n_after, n_before - n_after))
  
  # Step 4: Missingness filter (within 'group': BFR vs HLRT)
  cat(sprintf(" [%s] Step 4: Missingness filter\n", nm))
  n_before <- nrow(d$data)
  d        <- filter_proteins_by_proportion(d, min_prop = config$miss_pct / 100,
                                            grouping_column = grp)
  n_after  <- nrow(d$data)
  log      <- log_step(log, "Missingness filter", n_before, n_after)
  cat(sprintf("   %d -> %d (%d removed)\n", n_before, n_after, n_before - n_after))
  
  log <- bind_rows(
    tibble(step = "Raw input", n_before = NA_integer_,
           n_after = n_start, n_removed = NA_integer_),
    log
  ) %>% mutate(pct_retained = round(n_after / n_start * 100, 1))
  
  list(dal = d, log = log)
})

names(f.dals) <- names(dals)
dals        <- lapply(f.dals, `[[`, "dal")
filter_logs <- lapply(f.dals, `[[`, "log")

cat("\n   Filtering summary:\n")
lapply(names(filter_logs), function(nm) {
  cat(sprintf("\n  -- %s --\n", nm))
  print(filter_logs[[nm]])
})

# --- 5: Outlier Detection & Removal -----------------------------------------
cat("\n>> 5 - Outlier Detection & Removal\n")

# NOTE: BFR metadata has no timepoint/pid pairing (unlike the CR muscle
# SURV comparison), so this uses simple (unpaired) missingness rather than
# a paired delta-missingness flag. Consensus threshold is n_flags >= 2 of
# {missingness, PCA Mahalanobis, MAD intensity}.
run_outlier_detection <- function(d, nm) {
  cat(sprintf("\n [%s] Running outlier detection...\n", nm))
  
  pct_missing    <- colMeans(is.na(d$data)) * 100
  miss_threshold <- quantile(pct_missing, 0.75) + 1.5 * IQR(pct_missing)
  
  missing <- tibble(
    sample   = names(pct_missing),
    pct_missing = pct_missing,
    miss_flag   = pct_missing > miss_threshold
  )
  
  # PCA
  data_for_pca <- d$data
  for (j in seq_len(ncol(data_for_pca))) {
    nas <- is.na(data_for_pca[, j])
    if (any(nas)) data_for_pca[nas, j] <- median(data_for_pca[, j], na.rm = TRUE)
  }
  data_log2 <- log2(data_for_pca + 1)
  pca_res   <- prcomp(t(data_log2), center = TRUE, scale. = TRUE)
  
  pc_scores  <- pca_res$x[, 1:3]
  mahal_dist <- mahalanobis(pc_scores, colMeans(pc_scores), cov(pc_scores))
  
  pca_flags <- tibble(
    sample  = colnames(d$data),
    mahal_dist = mahal_dist,
    pca_flag   = mahal_dist > qchisq(1 - config$mahal_p, df = 3)
  )
  
  # MAD
  sample_medians <- apply(log2(d$data + 1), 2, median, na.rm = TRUE)
  global_median  <- median(sample_medians)
  mad_val        <- mad(sample_medians)
  
  mad_flags <- tibble(
    sample     = names(sample_medians),
    sample_median = sample_medians,
    mad_deviation = abs(sample_medians - global_median),
    mad_flag      = abs(sample_medians - global_median) > config$mad_k * mad_val
  )
  
  outlier_diag <- missing %>%
    dplyr::select(sample, pct_missing, miss_flag) %>%
    left_join(pca_flags, by = "sample") %>%
    left_join(mad_flags,  by = "sample") %>%
    mutate(n_flags = miss_flag + pca_flag + mad_flag,
           consensus_outlier = n_flags >= 2)
  
  n_outliers <- sum(outlier_diag$consensus_outlier)
  cat(sprintf("   Outliers detected: %d\n", n_outliers))
  
  list(
    outlier_diag = outlier_diag,
    pca_res      = pca_res,
    thresholds   = list(miss_threshold = miss_threshold)
  )
}

run_outlier_removal <- function(d, outlier_diag, nm) {
  n_outliers <- sum(outlier_diag$consensus_outlier)
  
  if (n_outliers > 0) {
    outlier_ids <- outlier_diag %>% filter(consensus_outlier) %>% pull(sample)
    cat(sprintf("   [%s] Removing %d outlier sample(s): %s\n",
                nm, length(outlier_ids), paste(outlier_ids, collapse = ", ")))
    keep_idx   <- !colnames(d$data) %in% outlier_ids
    d$data     <- d$data[, keep_idx, drop = FALSE]
    d$metadata <- d$metadata[d$metadata$sample %in% colnames(d$data), ]
    cat(sprintf("   After outlier removal: %d samples remain\n", ncol(d$data)))
  } else {
    cat(sprintf("   [%s] No outlier samples removed.\n", nm))
  }
  
  d
}

outlier_diag <- lapply(names(dals),
                       function(nm) run_outlier_detection(dals[[nm]], nm)) |>
  setNames(names(dals))

# Save data before outlier removal for quality control
dals_pre_outlier <- dals

dals <- lapply(names(dals),
               function(nm) run_outlier_removal(dals[[nm]], outlier_diag[[nm]]$outlier_diag, nm)) |>
  setNames(names(dals))

# --- 6: Normalize Data ----------------------------
cat("\n>> 6 - Normalize Data\n")

lapply(names(dals), function(nm) {
  out_dir <- config$data_dir
  saveRDS(dals[[nm]], file.path(out_dir, sprintf("00_DAList_pre_norm_%s.rds", nm)))
})

lapply(names(dals), function(nm) {
  report_dir <- config$report_dir
  write_norm_report(dals[[nm]], grouping_column = grouping_cols[nm],
                    output_dir = report_dir,
                    filename   = sprintf("02_norm_comparison_%s.pdf", nm),
                    overwrite  = TRUE)
  write_qc_report(dals[[nm]], color_column = grouping_cols[nm], label_column = "sample",
                  output_dir = report_dir,
                  filename   = sprintf("01_qc_report_pre_norm_%s.pdf", nm),
                  overwrite  = TRUE)
})

dals <- lapply(dals, normalize_data, norm_method = config$norm_method)

cat("\n Diagnostic: post-normalization data range (to check if already log-scale)\n")
lapply(names(dals), function(nm) {
  rng <- range(dals[[nm]]$data, na.rm = TRUE)
  cat(sprintf("  [%s] range: %.2f to %.2f\n", nm, rng[1], rng[2]))
})

lapply(names(dals), function(nm) {
  report_dir <- config$report_dir
  write_qc_report(dals[[nm]], color_column = grouping_cols[nm], label_column = "sample",
                  output_dir = report_dir,
                  filename   = sprintf("03_qc_report_post_norm_%s.pdf", nm),
                  overwrite  = TRUE)
})

lapply(names(dals), function(nm) {
  cat(sprintf("  [%s] Final: %d proteins x %d samples\n", nm,
              nrow(dals[[nm]]$data), ncol(dals[[nm]]$data)))
})

pca_postnorm <- lapply(dals, function(d) compute_pca(d$data, log_transform = FALSE)) # set TRUE if diagnostic shows linear scale

# --- 7: Export Normalized Data ----------------------
cat("\n>> 7 - Export Normalized Data\n")

invisible(lapply(names(dals), function(nm) {
  d        <- dals[[nm]]
  data_dir <- config$data_dir
  export   <- bind_cols(
    d$annotation %>% dplyr::select(uniprot_id, gene, protein, description),
    d$data
  )
  
  write.csv(export, file.path(data_dir, sprintf("01_normalized_data_%s.csv", nm)))
  saveRDS(d,        file.path(data_dir, sprintf("02_DAList_normalized_%s.RDS", nm)))
  
  cat(sprintf("  [%s] Exported: %d proteins x %d samples\n", nm, nrow(d$data), ncol(d$data)))
}))

# --- 8: Store Data for Reports --------------
report_data <- list(
  config           = config,
  flagged_hpa      = flagged_hpa,
  flagged_blood    = flagged_blood,
  marker_pct_df    = marker_pct_df,
  outlier_diag     = outlier_diag,
  filter_logs      = filter_logs,
  dals_pre_outlier = dals_pre_outlier,
  sample_summary   = sample_summary,
  dals             = dals,
  grouping_cols    = grouping_cols,
  dals_raw         = dals_raw,
  missingness_raw  = missingness_raw,
  pca_raw          = pca_raw,
  pca_postnorm     = pca_postnorm
)

saveRDS(report_data, file.path(config$data_dir, "00_report_data.rds"))

writeLines(capture.output(sessionInfo()), file.path(config$data_dir, "sessionInfo.txt"))

cat("\n>> Normalization complete. Run 01_normalization_BFR_reports.R next.\n")