# Protein-Level Intensity Bar Plots: Training_BFR / Training_HLRT, Up/Down
#
# For each of Training_BFR and Training_HLRT, and separately for UP- and
# DOWN-regulated proteins (per DE.pi in the limma results): a bar plot of
# mean normalized intensity per protein, with PRE/POST as differently
# colored dodged bars, proteins ordered along the x-axis by significance
# (most significant first). One plot per (contrast x direction) = 4 total.
#
# Intensity comes from the normalized DAList, not the limma summary
# stats -- each contrast only uses that contrast's OWN arm's samples
# (Training_BFR -> BFR arm samples only, both PRE and POST; Training_HLRT
# -> HLRT arm samples only), matching what that contrast actually tests.

# --- 0: Setup ------------------------
setwd(rprojroot::find_rstudio_root_file())
getwd()
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(readxl, dplyr, tidyr, tibble, purrr, ggplot2)

# --- Configuration -----------------------------------------------------------
config <- list(
  limma_data = "03_analysis/c_data/limma/limma_results.xlsx",
  norm_data  = "01_normalization/c_data/02_DAList_normalized_BFR.RDS",
  fig_dir    = "04_figures/fig_protein_intensity",
  data_dir   = "04_figures/data/fig_protein_intensity"
)
dir.create(config$fig_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(config$data_dir, showWarnings = FALSE, recursive = TRUE)

source("04_figures/a_scripts/00_theme.R")

pi_thresh        <- 0.05
TOP_N_PROTEINS   <- 20  # proteins shown per plot, most significant first
SHOW_ERROR_BARS  <- TRUE  # SEM error bars on each bar

# Same metadata-column assumptions used in fig_02_BFR.R's CV section and
# the sample report (01_normalization_BFR_reports.R). Confirm these match
# your actual metadata via the diagnostic printout in section 1 below.
GROUP_COL        <- "group"       # BFR / HLRT
TIMEPOINT_COL    <- "time"   # PRE / POST
TIMEPOINT_VALUES <- c(pre = "PRE", post = "POST")

# Bar colors (edit these -- any colors$... entry from 00_theme.R)
PRE_COLOR  <- colors$blue_light
POST_COLOR <- colors$blue_dark

# Colors for the "significant in both contrasts" plot (4 groups: arm x
# timepoint). Edit these -- any colors$... entry from 00_theme.R.
BOTH_GROUP_COLORS <- c(
  BFR_PRE   = colors$red_light,
  BFR_POST  = colors$red_dark,
  HLRT_PRE  = colors$blue_light,
  HLRT_POST = colors$blue_dark
)
TOP_N_BOTH <- 20  # proteins shown in the "significant in both" plots

contrasts_to_plot <- c(Training_BFR = "BFR", Training_HLRT = "HLRT")  # contrast -> arm

comp_labels <- c(
  Training_BFR  = "BFR_POST vs BFR_PRE",
  Training_HLRT = "HLRT_POST vs HLRT_PRE"
)

# --- 1: Load data --------------------------------------------------------
limma_results <- list(
  Training_BFR  = read_excel(config$limma_data, sheet = 2),
  Training_HLRT = read_excel(config$limma_data, sheet = 3)
)

norm_dal <- readRDS(config$norm_data)

cat("\n=== Metadata column check ===\n")
cat("Metadata columns available:", paste(names(norm_dal$metadata), collapse = ", "), "\n")
if (all(c(GROUP_COL, TIMEPOINT_COL) %in% names(norm_dal$metadata))) {
  cat("Group x timepoint crosstab:\n")
  print(table(norm_dal$metadata[[GROUP_COL]], norm_dal$metadata[[TIMEPOINT_COL]], useNA = "ifany"))
} else {
  stop("GROUP_COL ('", GROUP_COL, "') and/or TIMEPOINT_COL ('", TIMEPOINT_COL,
       "') not found in metadata. Actual columns: ", paste(names(norm_dal$metadata), collapse = ", "))
}

# --- 2: Intensity extraction ----------------------------------------------
# Matches significant proteins (by gene symbol) back to the normalized
# DAList's annotation/data, restricted to one arm's samples.
get_protein_intensity_long <- function(genes, arm) {
  row_idx <- match(genes, norm_dal$annotation$gene)
  valid   <- !is.na(row_idx)
  
  if (any(!valid)) {
    warning(sprintf("%d gene(s) not found in normalized data annotation and will be dropped: %s",
                    sum(!valid), paste(genes[!valid], collapse = ", ")))
  }
  
  row_idx       <- row_idx[valid]
  genes_matched <- genes[valid]
  
  meta_arm <- norm_dal$metadata %>% filter(.data[[GROUP_COL]] == arm)
  
  # Hard check: if this arm's samples don't actually contain BOTH
  # timepoints, every gene collapses to a single group_by(gene, timepoint)
  # row downstream -- i.e. exactly the "only one bar per protein" symptom,
  # rendered silently instead of erroring. Fail loudly here instead, with
  # the actual values found, so a TIMEPOINT_COL/TIMEPOINT_VALUES mismatch
  # is immediately diagnosable rather than looking like a plotting bug.
  actual_timepoints   <- unique(meta_arm[[TIMEPOINT_COL]])
  expected_timepoints <- unname(TIMEPOINT_VALUES)
  if (!all(expected_timepoints %in% actual_timepoints)) {
    stop(sprintf(
      "Arm '%s': expected both '%s' and '%s' in metadata column '%s', but only found: %s.\n",
      arm, TIMEPOINT_VALUES[["pre"]], TIMEPOINT_VALUES[["post"]], TIMEPOINT_COL,
      paste(actual_timepoints, collapse = ", ")
    ), "Fix TIMEPOINT_COL and/or TIMEPOINT_VALUES at the top of the script to match ",
    "the crosstab printed in section 1, then re-run.")
  }
  
  intensity_sub <- norm_dal$data[row_idx, meta_arm$sample, drop = FALSE]
  rownames(intensity_sub) <- genes_matched
  
  as.data.frame(intensity_sub) %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "sample", values_to = "intensity") %>%
    left_join(meta_arm %>% dplyr::select(sample, timepoint = all_of(TIMEPOINT_COL)), by = "sample")
}

summarise_intensity <- function(long_df) {
  long_df %>%
    filter(!is.na(intensity), !is.na(timepoint)) %>%
    group_by(gene, timepoint) %>%
    summarise(
      mean_intensity = mean(intensity, na.rm = TRUE),
      sem            = sd(intensity, na.rm = TRUE) / sqrt(sum(!is.na(intensity))),
      n              = sum(!is.na(intensity)),
      .groups = "drop"
    )
}

# --- 2b: Intensity extraction (both arms) ---------------------------------
# Same idea as get_protein_intensity_long(), but pulls ALL samples (both
# arms, both timepoints) instead of restricting to one arm -- needed for
# the 4-bar "significant in both contrasts" plot below.
get_protein_intensity_long_both_arms <- function(genes) {
  row_idx <- match(genes, norm_dal$annotation$gene)
  valid   <- !is.na(row_idx)
  
  if (any(!valid)) {
    warning(sprintf("%d gene(s) not found in normalized data annotation and will be dropped: %s",
                    sum(!valid), paste(genes[!valid], collapse = ", ")))
  }
  
  row_idx       <- row_idx[valid]
  genes_matched <- genes[valid]
  
  meta_all <- norm_dal$metadata
  
  # Same hard-check principle as the single-arm version, but across both
  # arms and both timepoints at once -- a 4-bar-per-protein plot needs
  # all four (group x timepoint) combinations to actually exist.
  actual_groups     <- unique(meta_all[[GROUP_COL]])
  actual_timepoints <- unique(meta_all[[TIMEPOINT_COL]])
  if (!all(unname(contrasts_to_plot) %in% actual_groups)) {
    stop("Expected both '", contrasts_to_plot[["Training_BFR"]], "' and '",
         contrasts_to_plot[["Training_HLRT"]], "' in metadata column '", GROUP_COL,
         "', but found: ", paste(actual_groups, collapse = ", "))
  }
  if (!all(unname(TIMEPOINT_VALUES) %in% actual_timepoints)) {
    stop("Expected both '", TIMEPOINT_VALUES[["pre"]], "' and '", TIMEPOINT_VALUES[["post"]],
         "' in metadata column '", TIMEPOINT_COL, "', but found: ",
         paste(actual_timepoints, collapse = ", "))
  }
  
  intensity_sub <- norm_dal$data[row_idx, meta_all$sample, drop = FALSE]
  rownames(intensity_sub) <- genes_matched
  
  as.data.frame(intensity_sub) %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "sample", values_to = "intensity") %>%
    left_join(
      meta_all %>% dplyr::select(sample, group = all_of(GROUP_COL), timepoint = all_of(TIMEPOINT_COL)),
      by = "sample"
    ) %>%
    mutate(group_timepoint = paste0(group, "_", timepoint))
}

summarise_intensity_both <- function(long_df) {
  long_df %>%
    filter(!is.na(intensity), !is.na(group_timepoint)) %>%
    group_by(gene, group_timepoint) %>%
    summarise(
      mean_intensity = mean(intensity, na.rm = TRUE),
      sem            = sd(intensity, na.rm = TRUE) / sqrt(sum(!is.na(intensity))),
      n              = sum(!is.na(intensity)),
      .groups = "drop"
    )
}

# --- 3: Plot builder ------------------------------------------------------
make_intensity_bar_plot <- function(summary_df, gene_order, title) {
  summary_df <- summary_df %>%
    mutate(
      gene      = factor(gene, levels = gene_order),
      timepoint = factor(timepoint, levels = c(TIMEPOINT_VALUES[["pre"]], TIMEPOINT_VALUES[["post"]]))
    )
  
  p <- ggplot(summary_df, aes(x = gene, y = mean_intensity, fill = timepoint)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.7, color = "black", linewidth = 0.2)
  
  if (SHOW_ERROR_BARS) {
    p <- p +
      geom_errorbar(
        aes(ymin = mean_intensity - sem, ymax = mean_intensity + sem),
        position = position_dodge(width = 0.75), width = 0.25, linewidth = 0.4
      )
  }
  
  p +
    scale_fill_manual(values = c(PRE = PRE_COLOR, POST = POST_COLOR), name = "Timepoint") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(x = NULL, y = "Normalized Intensity", title = title) +
    theme_cr(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
      legend.position = "right"
    )
}

# --- 3b: Plot builder (both contrasts, 4 bars per protein) -----------------
make_intensity_bar_plot_both <- function(summary_df, gene_order, title) {
  summary_df <- summary_df %>%
    mutate(
      gene            = factor(gene, levels = gene_order),
      group_timepoint = factor(group_timepoint, levels = names(BOTH_GROUP_COLORS))
    )
  
  p <- ggplot(summary_df, aes(x = gene, y = mean_intensity, fill = group_timepoint)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75, color = "black", linewidth = 0.2)
  
  if (SHOW_ERROR_BARS) {
    p <- p +
      geom_errorbar(
        aes(ymin = mean_intensity - sem, ymax = mean_intensity + sem),
        position = position_dodge(width = 0.8), width = 0.2, linewidth = 0.4
      )
  }
  
  p +
    scale_fill_manual(values = BOTH_GROUP_COLORS, name = "Group", drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(x = NULL, y = "Normalized Intensity", title = title) +
    theme_cr(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
      legend.position = "right"
    )
}

# --- 4: Build all 4 plots --------------------------------------------------
for (cmp in names(contrasts_to_plot)) {
  arm <- contrasts_to_plot[[cmp]]
  df  <- limma_results[[cmp]]
  
  for (direction in c("UP", "DOWN")) {
    sig <- df %>%
      filter(DE.pi == direction) %>%
      arrange(Pi.Val) %>%
      slice_head(n = TOP_N_PROTEINS)
    
    if (nrow(sig) == 0) {
      cat(sprintf("[%s / %s] No significant proteins -- skipping.\n", cmp, direction))
      next
    }
    
    gene_order <- sig$gene  # already ordered by significance (most sig first)
    
    long_df <- get_protein_intensity_long(gene_order, arm)
    summary_df <- summarise_intensity(long_df)
    
    plot_title <- sprintf("%s \u2014 %s (Top %d, %s-regulated)",
                          cmp, comp_labels[[cmp]], length(gene_order), direction)
    
    p <- make_intensity_bar_plot(summary_df, gene_order, plot_title)
    
    file_stub <- sprintf("intensity_%s_%s", cmp, tolower(direction))
    ggsave(file.path(config$fig_dir, paste0(file_stub, ".pdf")),
           p, width = max(8, length(gene_order) * 0.4), height = 6, device = cairo_pdf)
    ggsave(file.path(config$fig_dir, paste0(file_stub, ".png")),
           p, width = max(8, length(gene_order) * 0.4), height = 6, dpi = 300, bg = "white")
    
    write.csv(summary_df, file.path(config$data_dir, paste0(file_stub, "_data.csv")), row.names = FALSE)
    
    cat(sprintf("Saved: %s.pdf/.png (%d proteins)\n", file_stub, length(gene_order)))
  }
}

# --- 5: Proteins significant in BOTH contrasts (concordant direction) -----
# "Significant in both" here means the SAME direction in both contrasts
# (concordant), consistent with the Concordant/Discordant framing already
# used elsewhere in this pipeline (the pathway/protein concordance
# scripts) -- a protein flagged UP in Training_BFR and DOWN in
# Training_HLRT wouldn't make sense to show as a single "up-regulated"
# bar set, so those are excluded here rather than included ambiguously.
df_bfr_both  <- limma_results$Training_BFR  %>% dplyr::select(gene, Pi.Val_BFR = Pi.Val, DE.pi_BFR = DE.pi)
df_hlrt_both <- limma_results$Training_HLRT %>% dplyr::select(gene, Pi.Val_HLRT = Pi.Val, DE.pi_HLRT = DE.pi)

df_both <- inner_join(df_bfr_both, df_hlrt_both, by = "gene") %>%
  mutate(max_pi = pmax(Pi.Val_BFR, Pi.Val_HLRT))  # ranked so BOTH p-values must be small, not just one

for (direction in c("UP", "DOWN")) {
  sig_both <- df_both %>%
    filter(DE.pi_BFR == direction, DE.pi_HLRT == direction) %>%
    arrange(max_pi) %>%
    slice_head(n = TOP_N_BOTH)
  
  if (nrow(sig_both) == 0) {
    cat(sprintf("[Both contrasts / %s] No proteins significant (same direction) in both -- skipping.\n", direction))
    next
  }
  
  gene_order <- sig_both$gene
  
  long_df    <- get_protein_intensity_long_both_arms(gene_order)
  summary_df <- summarise_intensity_both(long_df)
  
  plot_title <- sprintf("Concordant %s-regulated \u2014 Significant in Both Contrasts (Top %d)",
                        direction, length(gene_order))
  
  p <- make_intensity_bar_plot_both(summary_df, gene_order, plot_title)
  
  file_stub <- sprintf("intensity_both_concordant_%s", tolower(direction))
  ggsave(file.path(config$fig_dir, paste0(file_stub, ".pdf")),
         p, width = max(10, length(gene_order) * 0.5), height = 6, device = cairo_pdf)
  ggsave(file.path(config$fig_dir, paste0(file_stub, ".png")),
         p, width = max(10, length(gene_order) * 0.5), height = 6, dpi = 300, bg = "white")
  
  write.csv(summary_df, file.path(config$data_dir, paste0(file_stub, "_data.csv")), row.names = FALSE)
  
  cat(sprintf("Saved: %s.pdf/.png (%d proteins)\n", file_stub, length(gene_order)))
}