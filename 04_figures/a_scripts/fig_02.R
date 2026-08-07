# Fig 02 — BFR
#
# NOTE: GSEA data comes from the GSEA consolidation script's saved RDS
# rather than being loaded and recomputed from raw per-comparison RDS
# files here. Section 5's "Significant Pathways (n)" panel counts only
# redundancy-filter survivors (see GSEA_COUNT_MODE below), consistent with
# what the bubble plot / other GSEA figures show.
#
# ADAPTED FROM THE CR MUSCLE/PLASMA VERSION: that script faceted every
# panel by tissue (Muscle/Plasma) and used 3 SURV comparisons. BFR is a
# single dataset with no tissue split, so every panel here is single-panel,
# and the comparison set is BFR's 5 contrasts (Baseline, Training_BFR,
# Training_HLRT, Post_training, Interaction) matching the GSEA/consolidation
# scripts.
#
# ASSUMPTION FLAGGED FOR REVIEW: section 6 (CV plot) needs a metadata
# column that crosses arm (BFR/HLRT) with timepoint (PRE/POST). The
# normalization script only established a "group" column with 2 levels
# (BFR/HLRT); if your metadata encodes timepoint under a different column
# name (or already has a single 4-level combined column), update
# GROUP_COL/TIMEPOINT_COL in section 6 accordingly.

# --- 0: Setup ------------------------
setwd(rprojroot::find_rstudio_root_file())
getwd()
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(readxl, readr, dplyr, tidyr, tibble, stringr, purrr, ggplot2,
               ggrepel, patchwork, cowplot, scales, enrichplot, grid)

# --- Configuration -----------------------------------------------------------
config <- list(
  limma_data = "03_analysis/c_data/limma/limma_results.xlsx",
  gsea_data  = "03_analysis/c_data/GSEA/consolidated/gsea_pathway_data.rds",
  norm_data  = "01_normalization/c_data/02_DAList_normalized_BFR.RDS",
  fig_dir    = "04_figures/fig_02",
  data_dir   = "04_figures/data/fig_02"
)
dir.create(config$fig_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(config$data_dir, showWarnings = FALSE, recursive = TRUE)

# --- Color palette + theme ----------------------------------------------------
source("04_figures/a_scripts/00_theme.R")

# --- Shared constants --------------------------------------------------------
thresh      <- 0.05
adjp_thresh <- 0.05
pi_thresh   <- 0.05

# "representative": count only redundancy-filter survivors (excludes
#   near-duplicate pathways within a comparison x database, per the pooled
#   greedy filter in the GSEA consolidation script). Matches what the
#   bubble plot and other GSEA figures show.
# "raw": count every significant pathway regardless of redundancy, i.e.
#   the original behavior of this panel. Kept as a one-line toggle in case
#   you want the two numbers side by side for a methods comparison.
GSEA_COUNT_MODE <- "representative"

gsea_group <- "BFR"  # matches the group column added by the consolidation script

gsea_comps <- c("Baseline", "Training_BFR", "Training_HLRT",
                "Post_training", "Interaction")

comp_labels <- c(
  Baseline      = "Baseline",
  Training_BFR  = "Training_BFR",
  Training_HLRT = "Training_HLRT",
  Post_training = "Post_training",
  Interaction   = "Interaction"
)

SET_LABELS <- comp_labels

db_slots <- c(
  gsea_bp_clean = "GO:BP",
  gsea_cc_clean = "GO:CC",
  gsea_mf_clean = "GO:MF",
  gsea_bp_slim   = "GO Slim",
  gsea_kegg      = "KEGG",
  gsea_reactome  = "Reactome",
  gsea_hall      = "Hallmark"
)

db_colors <- c(
  "GO:BP"    = "#1b9e77",
  "GO:CC"    = "#e6ab02",
  "GO:MF"    = "#7570b3",
  "GO Slim"  = "#66a61e",
  "KEGG"     = "#d95f02",
  "Reactome" = "#a6761d",
  "Hallmark" = "#e7298a"
)

DIR_COLORS <- c(Up = colors$red, Down = colors$blue, Mixed = colors$purple)

# --- 1: Load Data ------------------------------------------------------------

# NOTE: assumes the limma workbook's 5 sheets are in this order (matching
# the GSEA script's load order). Adjust if your workbook differs.
limma_sheets <- c(
  Baseline      = 1,
  Training_BFR  = 2,
  Training_HLRT = 3,
  Post_training = 4,
  Interaction   = 5
)

limma_results <- setNames(
  lapply(limma_sheets, function(s) read_excel(config$limma_data, sheet = s)),
  names(limma_sheets)
)

# GSEA: read from the consolidated pathway data instead of raw per-
# comparison RDS files. gsea_supp_full carries every tested pathway across
# all db_slots databases (already tagged with group/comparison/database/
# significant); representative_lookup flags which (group, database, ID)
# combos survived the pooled redundancy filter.
if (!file.exists(config$gsea_data)) {
  stop("Consolidated GSEA data not found at ", config$gsea_data,
       " -- run the GSEA consolidation script first.")
}
gsea_pathway_data <- readRDS(config$gsea_data)

gsea_supp_full        <- gsea_pathway_data$gsea_supp_full        %>% filter(group == gsea_group)
representative_lookup <- gsea_pathway_data$representative_lookup %>% filter(group == gsea_group)

norm_dal <- readRDS(config$norm_data)

# --- 2: Data wrangling -------------------------------------------------------
make_da_long <- function(limma_results) {
  bind_rows(lapply(names(limma_results), function(comp) {
    df <- limma_results[[comp]]
    data.frame(
      comparison = comp,
      pval_n     = sum(df$P.Value   < thresh,    na.rm = TRUE),
      adjp_n     = sum(df$adj.P.Val < thresh,    na.rm = TRUE),
      pival_n    = sum(df$Pi.Val    < pi_thresh, na.rm = TRUE)
    )
  })) |>
    pivot_longer(cols = c(pval_n, adjp_n, pival_n),
                 names_to = "method", values_to = "n") |>
    mutate(
      method     = factor(method,
                          levels = c("pval_n", "adjp_n", "pival_n"),
                          labels = c("p-value", "FDR", "\u03a0-value")),
      comparison = factor(comparison, levels = names(limma_results))
    )
}

da_long <- make_da_long(limma_results)

# --- 3: DA counts plot -------------------------------------------------------
n_methods   <- length(levels(da_long$method))
dodge_width <- 0.7
bar_width   <- dodge_width / n_methods

text_data_da <- da_long |>
  mutate(
    rank_in_full = as.integer(method),
    x_dodge      = as.numeric(comparison) +
      (rank_in_full - (n_methods + 1) / 2) * bar_width,
    offset       = max(n) * 0.03
  )

p_da <- ggplot(da_long, aes(x = comparison, y = n)) +
  geom_col(aes(fill = method),
           position = position_dodge(width = dodge_width), width = 0.6) +
  geom_text(
    data        = filter(text_data_da, n > 0),
    aes(x = x_dodge, y = n + offset, label = n),
    vjust       = 0, size = 3, fontface = "bold",
    inherit.aes = FALSE
  ) +
  geom_text(
    data        = filter(text_data_da, n == 0),
    aes(x = x_dodge, y = offset, label = "0"),
    vjust       = 0, size = 3, fontface = "bold",
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = c(colors$blue_light, colors$blue, colors$blue_dark)) +
  scale_x_discrete(labels = comp_labels) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "DAPs (n)", title = "DAPs per Contrast", fill = NULL) +
  theme_cr(base_size = 10, title_hjust = 0) +
  theme(
    legend.position      = c(0.97, 0.97),
    legend.justification = c(1, 1),
    legend.background    = element_rect(fill = NA, linewidth = 0),
    axis.text.x          = element_text(angle = 20, hjust = 1)
  )

# --- 4: UpSet plot -----------------------------------------------------------
# Generalized from the original (hardcoded to 3 named columns:
# Baseline_SURVvCTL/Training_SURVvCTL/Training_SURV) to work for any number
# of contrasts, since BFR has 5. Also folds in what used to be the separate
# "get_upset_components" helper from make_combined_upset -- with a single
# dataset there's no second tissue to stitch side-by-side, so that
# stitching machinery (make_combined_upset) is gone entirely.
make_upset_plot <- function(limma_results, pi_thresh, SET_LABELS, DIR_COLORS, title = NULL) {
  CONTRASTS <- names(limma_results)
  
  sig_sets <- list()
  dir_map  <- list()
  
  for (ctr in CONTRASTS) {
    df      <- limma_results[[ctr]]
    is_sig  <- !is.na(df$Pi.Val) & df$Pi.Val < pi_thresh
    sig_sets[[ctr]] <- df$uniprot[is_sig]
    dir_map[[ctr]]  <- setNames(
      ifelse(df$logFC[is_sig] > 0, "Up", "Down"),
      df$uniprot[is_sig]
    )
  }
  
  # Dynamic combination grid (2^n_contrasts - 1 non-empty combos), replaces
  # the original's hardcoded expand.grid(Baseline_SURVvCTL=..., ...).
  comb_grid <- expand.grid(rep(list(c(TRUE, FALSE)), length(CONTRASTS))) |>
    setNames(CONTRASTS) |>
    filter(rowSums(across(everything())) > 0)
  
  intersection_data <- bind_rows(lapply(seq_len(nrow(comb_grid)), function(i) {
    in_sets  <- CONTRASTS[which(as.logical(comb_grid[i, ]))]
    out_sets <- CONTRASTS[which(!as.logical(comb_grid[i, ]))]
    members  <- Reduce(intersect, sig_sets[in_sets])
    if (length(out_sets) > 0)
      members <- setdiff(members, unique(na.omit(unlist(sig_sets[out_sets]))))
    if (length(members) == 0) return(NULL)
    dirs <- sapply(members, function(prot) {
      d <- sapply(in_sets, \(ctr) dir_map[[ctr]][prot])
      d <- d[!is.na(d)]
      if (all(d == "Up")) "Up" else if (all(d == "Down")) "Down" else "Mixed"
    })
    tibble(
      in_sets = list(in_sets), degree = length(in_sets),
      up      = sum(dirs == "Up"),
      down    = sum(dirs == "Down"),
      mixed   = sum(dirs == "Mixed"),
      total   = length(members)
    )
  })) |>
    filter(total > 0) |>
    arrange(desc(total)) |>
    mutate(x = row_number())
  
  bar_long <- intersection_data |>
    mutate(mixed_only = up == 0 & down == 0 & mixed > 0) |>
    dplyr::select(x, up, down, mixed_only, total) |>
    pivot_longer(c(up, down), names_to = "direction", values_to = "count") |>
    mutate(
      direction = factor(str_to_title(direction), levels = c("Down", "Up")),
      count     = ifelse(mixed_only, 0L, count)
    )
  
  mixed_bar <- intersection_data |>
    filter(up == 0 & down == 0 & mixed > 0) |>
    dplyr::select(x, total)
  
  dot_df <- intersection_data |>
    dplyr::select(x, in_sets) |>
    mutate(set = map(in_sets, \(s) tibble(set = CONTRASTS, active = CONTRASTS %in% s))) |>
    unnest(set) |>
    mutate(
      set_label = factor(SET_LABELS[set], levels = rev(unname(SET_LABELS))),
      ynum      = as.numeric(set_label)
    )
  
  seg_df <- dot_df |>
    filter(active) |>
    group_by(x) |>
    filter(n() > 1) |>
    summarise(ymin = min(ynum), ymax = max(ynum), .groups = "drop")
  
  y_max           <- ceiling(max(c(bar_long$count, mixed_bar$total), na.rm = TRUE) * 1.15 / 10) * 10
  max_label_chars <- max(nchar(unname(SET_LABELS)))
  left_pad        <- max_label_chars * 5.5
  n_dirs          <- 2
  bar_width       <- 0.7 / n_dirs
  
  text_data_bars <- bar_long |>
    filter(count > 0) |>
    group_by(x) |>
    mutate(
      rank_in_full = as.integer(direction),
      x_dodge      = x + (rank_in_full - (n_dirs + 1) / 2) * bar_width,
      offset       = y_max * 0.03
    ) |>
    ungroup()
  
  p_bars <- ggplot(bar_long, aes(x = x, y = count, fill = direction)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6, linewidth = 0.3) +
    geom_col(
      data        = mixed_bar,
      aes(x = x, y = total, fill = "Mixed"),
      width       = 0.6,
      inherit.aes = FALSE
    ) +
    geom_text(
      data     = mixed_bar,
      aes(x = x, y = total + y_max * 0.03, label = total),
      vjust    = 0, size = 3, fontface = "bold",
      inherit.aes = FALSE
    ) +
    geom_text(
      data        = text_data_bars,
      aes(x = x_dodge, y = count + offset, label = count),
      vjust = 0, size = 3, inherit.aes = FALSE
    ) +
    scale_fill_manual(values = DIR_COLORS) +
    scale_x_continuous(breaks = intersection_data$x,
                       expand = expansion(add = 0.5)) +
    scale_y_continuous(limits = c(0, y_max),
                       expand = expansion(mult = c(0, 0)),
                       breaks = scales::breaks_pretty(n = 4)) +
    labs(y = "Intersection size", x = NULL, fill = NULL, title = title) +
    theme_cr(base_size = 10, title_hjust = 0) +
    theme(
      axis.text.x          = element_blank(),
      axis.ticks.x         = element_blank(),
      panel.grid.major.x   = element_blank(),
      panel.grid.minor     = element_blank(),
      legend.position      = c(0.97, 0.97),
      legend.justification = c(1, 1),
      plot.margin          = margin(4, 4, 0, left_pad)
    )
  
  n_sets <- length(CONTRASTS)
  
  p_dots <- ggplot() +
    annotate("rect",
             xmin = 0.5, xmax = max(intersection_data$x) + 0.5,
             ymin = seq_len(n_sets) - 0.45, ymax = seq_len(n_sets) + 0.45,
             fill = rep(c("grey96", "grey90"), length.out = n_sets)) +
    {if (nrow(seg_df) > 0)
      geom_segment(data = seg_df,
                   aes(x = x, xend = x, y = ymin, yend = ymax),
                   linewidth = 0.7, color = "grey25")} +
    geom_point(data = filter(dot_df, !active),
               aes(x = x, y = ynum), color = "grey82", size = 2.5) +
    geom_point(data = filter(dot_df, active),
               aes(x = x, y = ynum), color = "grey15", size = 2.5) +
    scale_x_continuous(expand = expansion(add = 0.5)) +
    scale_y_continuous(breaks = seq_len(n_sets),
                       labels = levels(dot_df$set_label),
                       expand = expansion(add = 0.4)) +
    labs(x = NULL, y = NULL) +
    theme_cr(base_size = 10, title_hjust = 0) +
    theme(
      axis.text.x  = element_blank(),
      axis.ticks   = element_blank(),
      panel.grid   = element_blank(),
      panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.3),
      axis.text.y  = element_text(size = 9, face = "bold"),
      plot.margin  = margin(0, 4, 4, 4)
    )
  
  (p_bars / p_dots) + plot_layout(heights = c(0.75, 0.25))
}

p_upset <- make_upset_plot(
  limma_results = limma_results,
  pi_thresh     = pi_thresh,
  SET_LABELS    = SET_LABELS,
  DIR_COLORS    = DIR_COLORS,
  title         = "Contrast Overlap (UpSet)"
)

# --- 5: GSEA counts plot -----------------------------------------------------
# Sourced from gsea_supp_full (consolidated). Counting logic is controlled
# by GSEA_COUNT_MODE:
#   "representative" -- a pathway counts only if it's both significant in
#     THIS comparison AND flagged representative for its (group, database)
#     by the pooled redundancy filter (representative_lookup). The
#     redundancy flag itself is pooled across all 5 gsea_comps (see the
#     consolidation script), but whether it counts toward THIS bar is
#     still gated on significance in this specific comparison -- so a
#     pathway significant in two comparisons still contributes to both
#     bars, it just won't be double-counted against a near-duplicate
#     pathway within the same bar.
#   "raw" -- original behavior: every significant pathway counts,
#     redundant or not.
make_gsea_counts <- function() {
  base <- gsea_supp_full %>%
    filter(comparison %in% gsea_comps, significant)
  
  if (GSEA_COUNT_MODE == "representative") {
    base <- base %>%
      semi_join(representative_lookup, by = c("group", "database", "ID"))
  } else if (GSEA_COUNT_MODE != "raw") {
    stop("GSEA_COUNT_MODE must be 'representative' or 'raw'")
  }
  
  counts <- base %>%
    mutate(database = as.character(database)) %>%
    group_by(comparison, database) %>%
    summarise(n = n(), .groups = "drop") %>%
    # ensure every comparison x database combo appears (0 counts included),
    # same as the old per-slot loop always emitting a row
    tidyr::complete(
      comparison = gsea_comps,
      database   = unname(db_slots),
      fill       = list(n = 0L)
    ) |>
    mutate(
      comparison = factor(comparison, levels = gsea_comps),
      database   = factor(database,   levels = unname(db_slots))
    )
  
  n_databases <- length(levels(counts$database))
  bar_width   <- 0.7 / n_databases
  
  text_data_gsea <- counts |>
    group_by(comparison) |>
    mutate(
      rank_in_full = as.integer(database),
      x_dodge      = as.numeric(comparison) +
        (rank_in_full - (n_databases + 1) / 2) * bar_width
    ) |>
    ungroup() |>
    mutate(offset = max(n) * 0.03)
  
  ggplot(counts, aes(x = comparison, y = n)) +
    geom_col(aes(fill = database),
             position  = position_dodge(width = 0.7), width = 0.6,
             color     = "black", linewidth = 0.3) +
    geom_text(
      data        = filter(text_data_gsea, n > 0),
      aes(x = x_dodge, y = n + offset, label = n),
      vjust       = 0, size = 3, fontface = "bold",
      inherit.aes = FALSE
    ) +
    geom_text(
      data        = filter(text_data_gsea, n == 0),
      aes(x = x_dodge, y = offset, label = "0"),
      vjust       = 0, size = 3, fontface = "bold",
      inherit.aes = FALSE
    ) +
    scale_fill_manual(values = db_colors) +
    scale_x_discrete(labels = comp_labels) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(x = NULL, y = "Significant Pathways (n)",
         title = "Pathway Enrichment", fill = NULL,
         caption = if (GSEA_COUNT_MODE == "representative")
           "Redundant pathways (Jaccard/overlap-coefficient filtered) excluded"
         else NULL) +
    theme_cr(base_size = 10, title_hjust = 0) +
    theme(
      panel.grid.major.x   = element_blank(),
      panel.grid.minor     = element_blank(),
      legend.key.size      = unit(0.4, "cm"),
      legend.text          = element_text(size = 7),
      legend.position      = c(0.97, 0.9),
      legend.justification = c(1, 1),
      legend.background    = element_rect(fill = NA, linewidth = 0.3),
      axis.text.x          = element_text(angle = 20, hjust = 1),
      plot.caption         = element_text(size = 6, face = "plain",
                                          hjust = 0, color = colors$gray_dark)
    )
}

p_gsea <- make_gsea_counts()

# --- 6: CV Plot --------------------------------------------------------------
# REVIEW THIS SECTION: adapted from CR's cancer/supp_time grouping to a
# generic "arm x timepoint" grouping for BFR. GROUP_COL is confirmed (the
# normalization script filters missingness "within 'group': BFR vs
# HLRT"); TIMEPOINT_COL is assumed -- update the name below (and the
# BFR_GROUPS labels/order) to match your actual metadata columns if
# different.
GROUP_COL     <- "group"      # BFR / HLRT
TIMEPOINT_COL <- "time"  # PRE / POST -- confirm/rename to match your metadata

build_cv_data <- function(data, annot, meta, group_col = GROUP_COL, timepoint_col = TIMEPOINT_COL) {
  group_map <- meta %>%
    mutate(group_time = paste0(.data[[group_col]], "_", .data[[timepoint_col]])) %>%
    dplyr::select(sample, group = group_time) %>%
    filter(!is.na(group))
  
  rn   <- rownames(data)
  data <- as.data.frame(base::apply(data, 2, as.numeric))
  rownames(data) <- rn
  
  data %>%
    rownames_to_column("uniprot_id") %>%
    left_join(annot %>% dplyr::select(uniprot_id, gene), by = "uniprot_id") %>%
    filter(!is.na(gene)) %>%
    dplyr::select(gene, any_of(group_map$sample)) %>%
    pivot_longer(-gene, names_to = "sample", values_to = "intensity") %>%
    left_join(group_map, by = "sample", relationship = "many-to-many") %>%
    filter(!is.na(group), !is.na(intensity)) %>%
    group_by(gene, group) %>%
    summarise(
      mean_int = mean(intensity, na.rm = TRUE),
      sd_int   = sd(intensity,   na.rm = TRUE),
      n        = sum(!is.na(intensity)),
      cv       = (sd_int / mean_int) * 100,
      .groups  = "drop"
    ) %>%
    filter(mean_int > 0, n >= 3)
}

BFR_GROUPS <- c("HLRT_PRE", "HLRT_POST", "BFR_PRE", "BFR_POST")

GROUP_COLOURS <- c(
  HLRT_PRE  = colors$blue,
  HLRT_POST = colors$blue_dark,
  BFR_PRE   = colors$red,
  BFR_POST  = colors$red_dark
)

cv_df <- build_cv_data(
  data  = norm_dal$data,
  annot = norm_dal$annotation,
  meta  = norm_dal$metadata
) 

# --- 6b: Bootstrap median CV with 95% CI ------------------------------------
bootstrap_median_cv <- function(cv_df, groups, n_boot = 2000, seed = 42) {
  set.seed(seed)
  
  bind_rows(lapply(groups, function(grp) {
    vals <- cv_df$cv[cv_df$group == grp]
    vals <- vals[!is.na(vals)]
    
    boot_medians <- replicate(n_boot, median(sample(vals, replace = TRUE)))
    
    tibble(
      group  = grp,
      median = median(vals),
      ci_lo  = quantile(boot_medians, 0.025),
      ci_hi  = quantile(boot_medians, 0.975),
      n      = length(vals)
    )
  }))
}

boot_cv <- bootstrap_median_cv(cv_df, BFR_GROUPS)

cat("\nBFR — bootstrap median CV (95% CI):\n")
print(boot_cv)

# --- 6c: Brown-Forsythe global test (for methods reporting) ------------------
if (!requireNamespace("car", quietly = TRUE)) install.packages("car")
library(car)

bf_cv <- car::leveneTest(
  cv ~ factor(group, levels = BFR_GROUPS),
  data   = cv_df,
  center = median
)

cat("\nBFR — Brown-Forsythe test (equal CV spread across groups):\n")
print(bf_cv)

fstat <- round(bf_cv$`F value`[1], 2)
pval  <- bf_cv$`Pr(>F)`[1]
psym  <- ifelse(pval < 0.001, "p < 0.001", sprintf("p = %.3f", pval))
cat(sprintf("Brown-Forsythe: F = %s, %s\n", fstat, psym))

# --- 6d: Violin fill palette -------------------------------------------------
make_shades <- function(base_color, n, lightest = 0.3, darkest = 0.85) {
  colorRampPalette(c("white", base_color))(100)[
    round(seq(lightest, darkest, length.out = n) * 100)
  ]
}

# Two shades per arm (PRE lighter, POST darker)
hlrt_shades <- make_shades(colors$blue, n = 2)
bfr_shades  <- make_shades(colors$red,  n = 2)

violin_fills <- c(
  HLRT_PRE  = hlrt_shades[1],
  HLRT_POST = hlrt_shades[2],
  BFR_PRE   = bfr_shades[1],
  BFR_POST  = bfr_shades[2]
)

# --- 6e: CV violin plot ------------------------------------------------------
make_cv_violin <- function(cv_df, boot_cv) {
  cv_df   <- cv_df   %>% mutate(group = factor(group, levels = BFR_GROUPS))
  boot_cv <- boot_cv %>% mutate(group = factor(group, levels = BFR_GROUPS))
  
  violin_tops <- cv_df %>%
    group_by(group) %>%
    summarise(y_top = quantile(cv, 0.95, na.rm = TRUE), .groups = "drop")
  
  label_df <- boot_cv %>% left_join(violin_tops, by = "group")
  
  ggplot(cv_df, aes(x = group, y = cv)) +
    geom_violin(aes(fill = group), alpha = 0.75, linewidth = 0.4, trim = FALSE) +
    geom_boxplot(aes(fill = group), width = 0.12, alpha = 0.95,
                 linewidth = 0.4, outlier.shape = NA, colour = colors$gray_dark) +
    geom_point(
      data        = boot_cv,
      aes(x = group, y = median),
      inherit.aes = FALSE,
      size        = 2.5,
      shape       = 18,
      color       = colors$black
    ) +
    geom_errorbar(
      data        = boot_cv,
      aes(x = group, ymin = ci_lo, ymax = ci_hi),
      inherit.aes = FALSE,
      width       = 0.12,
      linewidth   = 0.8,
      color       = colors$black
    ) +
    geom_text(
      data        = label_df,
      aes(x = group, y = y_top, label = sprintf("%.1f%%", median)),
      inherit.aes = FALSE,
      vjust       = -0.5,
      size        = 3,
      fontface    = "bold",
      color       = colors$black
    ) +
    scale_fill_manual(values = violin_fills, guide = "none") +
    scale_x_discrete(labels = c(HLRT_PRE = "HLRT PRE", HLRT_POST = "HLRT POST",
                                BFR_PRE = "BFR PRE", BFR_POST = "BFR POST")) +
    scale_y_continuous(limits = c(0, NA),
                       expand = expansion(mult = c(0, 0.15))) +
    labs(
      x       = NULL,
      y       = "CV (%)",
      title   = "Per-Group CV \u2014 Normalized",
      caption = "\u25c6 Median CV (95% bootstrap CI, n = 2,000 resamples)"
    ) +
    theme_cr(base_size = 10, title_hjust = 0) +
    theme(
      axis.text.x  = element_text(face = "bold", size = 9, angle = 20, hjust = 1),
      plot.caption = element_text(size = 7, face = "plain",
                                  hjust = 0, color = colors$gray_dark)
    )
}

p_cv <- make_cv_violin(cv_df, boot_cv)

# --- 7: Assemble & Save ------------------------------------------------------
f2 <- plot_grid(
  p_da, p_gsea, p_upset, p_cv,
  labels = "AUTO",
  nrow   = 2,
  ncol   = 2
)

ggsave(file.path(config$fig_dir, "fig_02.pdf"),
       f2, width = 16, height = 12, units = "in", device = cairo_pdf)
ggsave(file.path(config$fig_dir, "fig_02.png"),
       f2, width = 16, height = 12, units = "in", dpi = 300)

ggsave(file.path(config$fig_dir, "fig_02_A.pdf"), p_da,    width = 6,  height = 6, device = cairo_pdf)
ggsave(file.path(config$fig_dir, "fig_02_B.pdf"), p_gsea,  width = 6,  height = 6, device = cairo_pdf)
ggsave(file.path(config$fig_dir, "fig_02_C.pdf"), p_upset, width = 10, height = 8, device = cairo_pdf)
ggsave(file.path(config$fig_dir, "fig_02_D.pdf"), p_cv,    width = 5,  height = 5, device = cairo_pdf)

ggsave(file.path(config$fig_dir, "fig_02_A.png"), p_da,    width = 6,  height = 6,  dpi = 300)
ggsave(file.path(config$fig_dir, "fig_02_B.png"), p_gsea,  width = 6,  height = 6,  dpi = 300)
ggsave(file.path(config$fig_dir, "fig_02_C.png"), p_upset, width = 10, height = 8,  dpi = 300)
ggsave(file.path(config$fig_dir, "fig_02_D.png"), p_cv,    width = 5,  height = 5,  dpi = 300)