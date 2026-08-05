# GSEA Pathway Consolidation — one canonical pass over all comparisons x databases
#
# PURPOSE
# Every downstream figure that touches GSEA results (bubble plot, heatmap-
# Sankey, ORA butterfly, reversal NES plots, etc.) needs the same three
# things: (1) the full tested-pathway table across every database, (2) the
# significance flag, and (3) which pathways survive redundancy filtering.
# This script runs that logic ONCE and writes a single RDS (+ companion
# xlsx) that every figure script loads, instead of each figure script
# recomputing (and potentially drifting on) the same filtering logic.
#
# WHAT LIVES HERE vs WHAT STAYS IN FIGURE SCRIPTS
# - HERE: loading raw GSEA RDS, building the full tested-pathway table,
#   the significance flag, and the pooled greedy redundancy filter
#   (Jaccard/overlap coefficient on @geneSets). This is "what counts as a
#   real, non-redundant signal" -- an analysis decision, not a plotting one.
# - FIGURE SCRIPTS: TOP_N selection, per-contrast vs pooled ranking,
#   pathway ordering/clustering for display, and all ggplot calls. This is
#   "what subset do I show and how" -- a display decision that legitimately
#   differs per figure.
#
# NOTE ON STRUCTURE: the BFR pipeline is a single dataset with one family
# of comparisons (Baseline, Training_BFR, Training_HLRT, Post_training,
# Interaction) -- there's no tissue or arm (SURV/CRE) split here, so this
# is organized around a single "group" dimension (named "BFR" below) for
# consistency with the rest of the pipeline's list-based structure, and to
# make it a one-line change (add another config$gsea_dirs entry +
# comp_defs entry) if a second comparison family is ever added.

# --- 0: Setup ------------------------------------------------------------
setwd(rprojroot::find_rstudio_root_file())
getwd()
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(readxl, readr, dplyr, tidyr, tibble, stringr, purrr, writexl)

# --- Configuration ---------------------------------------------------------
# gsea_dirs: group -> directory of per-comparison gsea_results_*.rds files.
# Add another entry here (with a matching comp_defs entry) if a second
# comparison family is ever run -- no other code in this script needs to
# change.
config <- list(
  gsea_dirs = list(
    BFR = "03_analysis/c_data/GSEA"
  ),
  out_dir = "03_analysis/c_data/GSEA/consolidated"
)
dir.create(config$out_dir, showWarnings = FALSE, recursive = TRUE)

# Comparisons pooled together per group -- redundancy is assessed on the
# UNION of significant pathways across all comps in the group, so a shared
# pathway is filtered consistently regardless of which comp(s) flagged it
# as significant.
comp_defs <- list(
  BFR = c("Baseline", "Training_BFR", "Training_HLRT",
          "Post_training", "Interaction")
)

comp_labels <- c(
  Baseline      = "BFR_PRE vs HLRT_PRE",
  Training_BFR  = "BFR_POST vs BFR_PRE",
  Training_HLRT = "HLRT_POST vs HLRT_PRE",
  Post_training = "BFR_POST vs HLRT_POST",
  Interaction   = "Interaction"
)

db_slots <- c(
  gsea_bp_clean  = "GO:BP",
  gsea_cc_clean  = "GO:CC",
  gsea_mf_clean  = "GO:MF",
  gsea_bp_slim   = "GO Slim",
  gsea_kegg      = "KEGG",
  gsea_reactome  = "Reactome",
  gsea_hall      = "Hallmark"
)
slot_for_database <- setNames(names(db_slots), unname(db_slots))

adjp_thresh <- 0.05

# --- Redundancy filter parameters -------------------------------------------
# Centralized here so every figure inherits the same definition of
# "redundant" instead of each figure script picking its own threshold.
JACCARD_THRESH <- 0.5
OVERLAP_THRESH <- 0.5
MAX_N_FOR_REDUNDANCY <- Inf  # cap per group/database if O(n^2) gets slow

# --- 1: Load ----------------------------------------------------------------
# One RDS per comparison, each holding a named list of gseaResult objects
# keyed by db_slots names, i.e. gsea_list[[comparison]][[slot]] -> gseaResult
load_gsea_rds <- function(dir) {
  files <- list.files(dir, pattern = "^gsea_results_.*\\.rds$", full.names = TRUE)
  nms   <- gsub("^gsea_results_|\\.rds$", "",
                list.files(dir, pattern = "^gsea_results_.*\\.rds$"))
  setNames(lapply(files, readRDS), nms)
}

# group -> loaded gsea_list. Only groups with a configured directory in
# config$gsea_dirs are loaded.
gsea_lists <- list()
for (group_name in names(config$gsea_dirs)) {
  gsea_lists[[group_name]] <- load_gsea_rds(config$gsea_dirs[[group_name]])
}

# --- Diagnostic check --------------------------------------------------------
# Run before trusting anything downstream -- if names don't match
# comp_defs/db_slots, fix the mismatch here, not by patching symptoms below.
for (group_name in names(gsea_lists)) {
  gl <- gsea_lists[[group_name]]
  cat("=== ", group_name, " ===\n", sep = "")
  cat("Comparisons found:", paste(names(gl), collapse = ", "), "\n")
  cat("Slots in first comparison:", paste(names(gl[[1]]), collapse = ", "), "\n")
  
  expected_comps <- comp_defs[[group_name]]
  if (is.null(expected_comps)) {
    stop("No comp_defs entry for group '", group_name,
         "' -- add one before this group can be consolidated.")
  }
  if (!all(expected_comps %in% names(gl))) {
    stop(group_name, " is missing comparisons: ",
         paste(setdiff(expected_comps, names(gl)), collapse = ", "))
  }
  if (!all(names(db_slots) %in% names(gl[[1]]))) {
    stop(group_name, " is missing db_slots entries: ",
         paste(setdiff(names(db_slots), names(gl[[1]])), collapse = ", "))
  }
}

# --- 2: Data wrangling -------------------------------------------------------
# Full tested-pathway table across every database in db_slots (not just
# whatever a given figure plots), with core_enrichment retained for the
# redundancy-filter fallback below.
make_gsea_full_data <- function(gsea_list, group_name, comps) {
  bind_rows(lapply(comps, function(comp) {
    bind_rows(lapply(names(db_slots), function(slot) {
      res <- gsea_list[[comp]][[slot]]
      if (is.null(res)) return(NULL)
      df <- as.data.frame(res)
      if (nrow(df) == 0) return(NULL)
      if (!"core_enrichment" %in% colnames(df)) df$core_enrichment <- NA_character_
      df %>%
        dplyr::select(ID, Description, NES, pvalue, p.adjust, setSize, core_enrichment) %>%
        dplyr::mutate(
          comparison = comp,
          database   = db_slots[[slot]],
          group      = group_name
        )
    }))
  }))
}

gsea_supp_full <- purrr::map_dfr(names(gsea_lists), function(group_name) {
  make_gsea_full_data(gsea_lists[[group_name]], group_name, comp_defs[[group_name]])
}) %>%
  mutate(
    comparison_label  = comp_labels[as.character(comparison)],
    database          = factor(database, levels = unname(db_slots)),
    Description_clean = Description %>%
      str_remove("^HALLMARK_") %>%
      str_replace_all("_", " ") %>%
      str_to_sentence(),
    significant = p.adjust < adjp_thresh
  )

cat("\n=== Combined data check ===\n")
cat("Rows:", nrow(gsea_supp_full), "\n")
print(table(gsea_supp_full$database, gsea_supp_full$group))

gsea_supp_raw <- gsea_supp_full %>%
  dplyr::select(group, database,
                comparison = comparison_label,
                ID, Description = Description_clean,
                NES, pvalue, p.adjust, setSize, significant) %>%
  arrange(group, database, comparison, p.adjust)

gsea_supp_summary <- gsea_supp_raw %>%
  group_by(group, database) %>%
  summarise(
    n_tested      = n(),
    n_significant = sum(significant, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n=== Summary sheet check (pre-redundancy) ===\n")
print(gsea_supp_summary)

# --- 3: Pooled redundancy filtering (Jaccard + overlap coefficient) ---------
# Greedy significance-ordered FILTER, not clustering: pathways ranked by
# significance, then a pathway is dropped if it's redundant (Jaccard OR
# overlap coefficient >= threshold) with an already-kept, more-significant
# pathway. No transitive grouping -- A and C are never both dropped/merged
# just because each is separately similar to some B.
#
# Filtering is done ONCE on the pooled (union) set of significant pathways
# across all comps in a group, per group x database -- NOT independently
# per contrast (see rationale in header comment above).
#
# Similarity is computed on FULL pathway gene-set membership (gseaResult's
# @geneSets slot), not core_enrichment/leading-edge genes, since leading
# edge is contrast-specific. core_enrichment is used only as a fallback if
# @geneSets isn't available, and flagged as such in the output.
parse_gene_set <- function(x) {
  if (is.na(x) || x == "") return(character(0))
  str_split(x, "/")[[1]]
}

jaccard_index <- function(a, b) {
  if (length(a) == 0 || length(b) == 0) return(0)
  length(intersect(a, b)) / length(union(a, b))
}

overlap_coefficient <- function(a, b) {
  if (length(a) == 0 || length(b) == 0) return(0)
  length(intersect(a, b)) / min(length(a), length(b))
}

get_reference_gene_sets <- function(gsea_list, comp_priority, slot) {
  for (comp in comp_priority) {
    res <- gsea_list[[comp]][[slot]]
    if (is.null(res)) next
    gs <- tryCatch(res@geneSets, error = function(e) NULL)
    if (!is.null(gs) && length(gs) > 0) return(gs)
  }
  NULL
}

empty_redundancy_result <- tibble(
  ID = character(), Description_clean = character(), min_padj = double(),
  group = character(), database = character(),
  is_representative = logical(), redundant_with = character(), gene_source = character()
)

build_redundancy_mapping <- function(group_name, database_name) {
  slot <- slot_for_database[[database_name]]
  gene_sets_ref <- get_reference_gene_sets(
    gsea_lists[[group_name]], comp_defs[[group_name]], slot
  )
  used_fallback <- is.null(gene_sets_ref)
  
  sig_df <- gsea_supp_full %>%
    filter(group == group_name, database == database_name, significant) %>%
    group_by(ID, Description_clean) %>%
    summarise(
      min_padj            = min(p.adjust, na.rm = TRUE),
      core_enrichment_any = dplyr::first(na.omit(core_enrichment)),
      .groups = "drop"
    ) %>%
    arrange(min_padj)
  
  n <- nrow(sig_df)
  if (n == 0) return(empty_redundancy_result)
  
  if (is.finite(MAX_N_FOR_REDUNDANCY) && n > MAX_N_FOR_REDUNDANCY) {
    cat("Note:", group_name, "/", database_name, "--", n,
        "significant pathways exceeds MAX_N_FOR_REDUNDANCY; filtering only the top",
        MAX_N_FOR_REDUNDANCY, "by p.adjust.\n")
    sig_df <- sig_df %>% slice_head(n = MAX_N_FOR_REDUNDANCY)
    n <- nrow(sig_df)
  }
  
  gene_source <- if (used_fallback) "core_enrichment" else "full_gene_set"
  if (used_fallback) {
    cat("Note:", group_name, "/", database_name,
        "-- no @geneSets found on the stored gseaResult; falling back to",
        "core_enrichment (leading edge), which is contrast-specific and a",
        "less exact measure of pathway overlap.\n")
    gene_lists <- lapply(sig_df$core_enrichment_any, parse_gene_set)
  } else {
    gene_lists <- gene_sets_ref[sig_df$ID]
  }
  
  if (n == 1) {
    return(sig_df %>%
             dplyr::select(-core_enrichment_any) %>%
             mutate(group = group_name, database = database_name,
                    is_representative = TRUE, redundant_with = NA_character_,
                    gene_source = gene_source))
  }
  
  keep           <- rep(TRUE, n)
  redundant_with <- rep(NA_character_, n)
  for (i in 2:n) {
    for (j in 1:(i - 1)) {
      if (!keep[j]) next
      a <- gene_lists[[i]]; b <- gene_lists[[j]]
      if (is.null(a) || is.null(b) || length(a) == 0 || length(b) == 0) next
      if (jaccard_index(a, b) >= JACCARD_THRESH || overlap_coefficient(a, b) >= OVERLAP_THRESH) {
        keep[i] <- FALSE
        redundant_with[i] <- sig_df$ID[j]
        break
      }
    }
  }
  
  sig_df %>%
    dplyr::select(-core_enrichment_any) %>%
    mutate(group = group_name, database = database_name,
           is_representative = keep, redundant_with = redundant_with,
           gene_source = gene_source)
}

redundancy_grid <- gsea_supp_full %>%
  distinct(group, database) %>%
  mutate(database = as.character(database)) %>%
  arrange(group, database)

pathway_redundancy_mapping <- purrr::pmap_dfr(
  redundancy_grid,
  ~ build_redundancy_mapping(group_name = ..1, database_name = ..2)
) %>%
  relocate(group, database, .before = ID) %>%
  arrange(group, database, min_padj)

pathway_redundancy_summary <- redundancy_grid %>%
  left_join(
    pathway_redundancy_mapping %>%
      group_by(group, database) %>%
      summarise(n_sig_unique = n(), n_representative = sum(is_representative), .groups = "drop"),
    by = c("group", "database")
  ) %>%
  mutate(
    n_sig_unique     = coalesce(n_sig_unique, 0L),
    n_representative = coalesce(n_representative, 0L)
  )

cat("\n=== Pooled redundancy filtering (Jaccard >=", JACCARD_THRESH,
    "or Overlap Coefficient >=", OVERLAP_THRESH, ") ===\n")
print(pathway_redundancy_summary)

# representative_lookup: THE lookup every figure script should use to
# restrict candidate pathways to redundancy-filter survivors before doing
# any figure-specific TOP_N/ranking.
representative_lookup <- pathway_redundancy_mapping %>%
  filter(is_representative) %>%
  distinct(group, database, ID)

# Fold unique/representative counts into the same summary as n_tested/n_significant
gsea_supp_summary <- gsea_supp_summary %>%
  mutate(database = as.character(database)) %>%
  left_join(pathway_redundancy_summary, by = c("group", "database"))

cat("\n=== Summary sheet check (with redundancy columns) ===\n")
print(gsea_supp_summary)

# --- 4: Save -----------------------------------------------------------------
# ONE canonical object. Every figure script loads this instead of
# recomputing sections 1-3 above.
gsea_pathway_data <- list(
  gsea_supp_full             = gsea_supp_full,
  gsea_supp_raw              = gsea_supp_raw,
  gsea_supp_summary          = gsea_supp_summary,
  pathway_redundancy_mapping = pathway_redundancy_mapping,
  representative_lookup      = representative_lookup,
  pathway_redundancy_summary = pathway_redundancy_summary,
  params = list(
    adjp_thresh           = adjp_thresh,
    jaccard_thresh         = JACCARD_THRESH,
    overlap_thresh         = OVERLAP_THRESH,
    max_n_for_redundancy   = MAX_N_FOR_REDUNDANCY,
    comp_defs              = comp_defs,
    comp_labels            = comp_labels,
    db_slots               = db_slots,
    generated_at           = Sys.time()
  )
)

saveRDS(gsea_pathway_data, file.path(config$out_dir, "gsea_pathway_data.rds"))

writexl::write_xlsx(
  list(
    Summary            = gsea_supp_summary,
    Redundancy_Mapping = pathway_redundancy_mapping,
    Results            = gsea_supp_raw
  ),
  path = file.path(config$out_dir, "gsea_pathway_data_supplement.xlsx")
)

cat("\nSaved consolidated GSEA pathway data to:\n  ",
    file.path(config$out_dir, "gsea_pathway_data.rds"), "\n  ",
    file.path(config$out_dir, "gsea_pathway_data_supplement.xlsx"), "\n", sep = "")