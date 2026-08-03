# === PROJECT: BFR v HLRT Proteomics ==============================
# === Script: Normalization Reports

# === PURPOSE: Generate QC reports, contaminant summaries, filter logs,
# outlier diagnostics, CV plots, and protein/sample assessments from
# the saved normalization report data.

# === INPUTS:  00_report_data.rds (from 01_normalization_BFR.R)
# === OUTPUTS: plots, CSVs, and PDFs in report_dir / data_dir

# --- 0: Setup ------------------------------------------------------------
setwd(rprojroot::find_rstudio_root_file())

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(proteoDA, readxl, readr, dplyr, tidyr, stringr, purrr, ggplot2,
               ggrepel, patchwork, tibble, ggpubr, Hmisc, scales, vegan)

config <- list(
  report_data = "01_normalization/c_data/00_report_data.rds",
  report_dir  = "01_normalization/b_reports",
  data_dir    = "01_normalization/c_data"
)

dir.create(config$report_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(config$data_dir,   showWarnings = FALSE, recursive = TRUE)

# Color palette (matches the CR muscle report palette for visual consistency)
colors <- list(
  red          = "#C0392B", red_light    = "#E6A19A", red_dark     = "#7B241C",
  blue         = "#2471A3", blue_light   = "#A9CCE3", blue_dark    = "#154360",
  green        = "#4E8B6F", green_light  = "#A9C9B8", green_dark   = "#2E5E4E",
  pink         = "#D98C9C", pink_light   = "#F2C6CF", pink_dark    = "#A65C6B",
  purple       = "#7B5EA7", purple_light = "#C5B4DD", purple_dark  = "#4A3A73",
  yellow       = "#D4A017", yellow_light = "#F1D77A", yellow_dark  = "#9A6F0A",
  orange       = "#C4693A", orange_light = "#E5A07B", orange_dark  = "#8A4321",
  aqua         = "#5DAFA3", aqua_light   = "#A9D6CF", aqua_dark    = "#3A7F76",
  gray_light   = "#F5F5F5", gray         = "#9E9E9E", gray_dark    = "#4A4A4A",
  black        = "#000000", white        = "#FFFFFF"
)

panel_border_theme <- theme(
  panel.border = element_rect(colour = colors$gray_dark, fill = NA, linewidth = 0.6)
)

# --- 1: Load Data --------------------------------------------------------
report_data <- readRDS(config$report_data)

flagged_hpa     <- report_data$flagged_hpa
flagged_blood   <- report_data$flagged_blood
marker_pct_df   <- report_data$marker_pct_df
sample_summary  <- report_data$sample_summary
outlier_diag    <- report_data$outlier_diag
filter_logs     <- report_data$filter_logs
dals            <- report_data$dals
annotation      <- dals$BFR$annotation
meta            <- dals$BFR$metadata
grouping_cols   <- report_data$grouping_cols
dals_raw        <- report_data$dals_raw
missingness_raw <- report_data$missingness_raw
pca_raw         <- report_data$pca_raw
pca_postnorm    <- report_data$pca_postnorm

dal_prenorm_BFR <- readRDS(file.path(config$data_dir, "00_DAList_pre_norm_BFR.rds"))

# --- 1b: Removed Sample Sets (for bolding flagged labels) ----------------

removed_samples <- map(names(dals), function(nm) {
  raw_ids   <- dals_raw[[nm]]$metadata$sample
  final_ids <- dals[[nm]]$metadata$sample
  setdiff(raw_ids, final_ids)
}) %>% set_names(names(dals))

outlier_samples <- map(names(outlier_diag), function(nm) {
  outlier_diag[[nm]]$outlier_diag %>%
    filter(consensus_outlier) %>%
    pull(sample)
}) %>% set_names(names(outlier_diag))

removed_samples_all <- unique(unlist(removed_samples))

# Bold flagged sample IDs using base R plotmath expressions (renders bold +
# rotated natively; ggtext::element_markdown() can't do markdown at angle != 0)
bold_flagged_labels <- function(x, flagged) {
  lapply(x, function(lbl) {
    if (lbl %in% flagged) bquote(bold(.(lbl))) else lbl
  })
}

# --- 2: Flagged Contaminants ---------------------------------------------

flagged_contaminants_all <- bind_rows(flagged_hpa, flagged_blood) %>%
  group_by(gene) %>%
  summarise(
    protein             = dplyr::first(protein),
    gene                = dplyr::first(gene),
    description         = dplyr::first(description),
    in_hpa_muscle       = dplyr::first(in_hpa_muscle),
    in_blood_blacklist  = dplyr::first(in_blood_blacklist),
    .groups = "drop"
  )

write.csv(flagged_hpa,   file.path(config$report_dir, "01_flagged_contaminants_hpa.csv"), row.names = FALSE)
write.csv(flagged_blood, file.path(config$report_dir, "01_flagged_contaminants_blood.csv"), row.names = FALSE)
write.csv(flagged_contaminants_all, file.path(config$report_dir, "01_flagged_contaminants_all.csv"), row.names = FALSE)
cat("\n Saved flagged contaminant CSVs to:", config$report_dir, "\n")

# --- 3a: Blood vs Muscle Intensity Plot -----------------------------------

marker_levels <- c("ALB", "HBB", "HBA1", "MB", "CKM")
marker_colors <- c(
  ALB  = colors$red_dark,
  HBB  = colors$red_light,
  HBA1 = colors$red,
  MB   = colors$aqua,
  CKM  = colors$aqua_dark
)

make_bm_plot <- function(marker_df, sample_ids, label, flagged_samples = character(0)) {
  df <- marker_df %>% filter(sample %in% sample_ids)
  
  label_df <- df %>%
    dplyr::select(sample, ALB_pct, HBB_pct, HBA1_pct, MB_pct, CKM_pct, B_M_ratio) %>%
    mutate(
      blood_total = ALB_pct + HBB_pct + HBA1_pct,
      total_pct   = blood_total + MB_pct + CKM_pct
    ) %>%
    arrange(desc(blood_total)) %>%
    mutate(sample = factor(sample, levels = unique(sample)))
  
  sample_order <- levels(label_df$sample)
  
  df %>%
    dplyr::select(sample, ALB_pct, HBB_pct, HBA1_pct, MB_pct, CKM_pct) %>%
    pivot_longer(-sample, names_to = "marker", values_to = "pct") %>%
    mutate(
      marker    = gsub("_pct", "", marker),
      marker    = factor(marker, levels = marker_levels),
      sample = factor(sample, levels = sample_order)
    ) %>%
    ggplot(aes(x = sample, y = pct, fill = marker)) +
    geom_col(width = 0.8) +
    scale_fill_manual(name = "Marker", values = marker_colors) +
    scale_x_discrete(labels = function(x) bold_flagged_labels(x, flagged_samples)) +
    geom_text(
      data        = label_df,
      aes(x = sample, y = total_pct, label = round(B_M_ratio, 2)),
      vjust       = -0.4, size = 2.5, inherit.aes = FALSE
    ) +
    labs(x = "Sample", y = "% of total intensity",
         title = sprintf("Blood vs muscle marker %% of total intensity (%s)", label)) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
          legend.position = "bottom") +
    panel_border_theme
}

p_bm <- make_bm_plot(marker_pct_df, marker_pct_df$sample, "all samples",
                     flagged_samples = removed_samples_all)

ggsave(file.path(config$report_dir, "01_blood_v_muscle_intensity_plot_raw.pdf"),
       p_bm, width = 14, height = 9)

write.csv(marker_pct_df, file.path(config$data_dir, "01_blood_contamination.csv"), row.names = FALSE)
cat(" Saved: 01_blood_v_muscle_intensity_plot_raw.pdf\n")

# --- 4: Filter Logs & Filtered Protein Lists -----------------------------

walk(names(filter_logs), function(nm) {
  log     <- filter_logs[[nm]]
  d       <- dals[[nm]]
  out_dir <- config$data_dir
  rep_dir <- config$report_dir
  
  write_csv(log, file.path(rep_dir, sprintf("02_filtering_effects_%s.csv", nm)))
  
  kept_uniprots    <- rownames(d$data)
  removed_uniprots <- setdiff(annotation$uniprot_id, kept_uniprots)
  
  filtered_proteins <- annotation %>%
    filter(uniprot_id %in% removed_uniprots) %>%
    dplyr::select(uniprot_id, gene, description, in_hpa_muscle, in_blood_blacklist) %>%
    mutate(reason = if_else(!in_hpa_muscle | in_blood_blacklist, "Contaminant", "Missingness"))
  
  write_csv(filtered_proteins, file.path(rep_dir, sprintf("03_filtered_proteins_%s.csv", nm)))
  
  cat(sprintf(" Saved: 02_filtering_effects_%s.csv\n", nm))
  cat(sprintf(" Saved: 03_filtered_proteins_%s.csv\n", nm))
})

# --- 5: Filter Plots -----------------------------------------------------

make_filter_plot <- function(log, label) {
  log_plot <- log %>% filter(!is.na(n_after))
  ggplot(mutate(log_plot, step = factor(step, levels = step)),
         aes(x = step, y = n_after)) +
    geom_col(fill = colors$blue, width = 0.5) +
    geom_text(aes(label = n_after), vjust = -0.3, size = 3) +
    labs(x = NULL, y = "Number of Proteins", title = "Proteins Retained Through Filtering") +
    theme_minimal(base_size = 13) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1)) +
    panel_border_theme
}

walk(names(filter_logs), function(nm) {
  rep_dir <- config$report_dir
  p <- make_filter_plot(filter_logs[[nm]], nm)
  ggsave(file.path(rep_dir, sprintf("04_filtering_plot_%s.pdf", nm)), p, width = 8, height = 5)
  cat(sprintf(" Saved: 04_filtering_plot_%s.pdf\n", nm))
})

make_protein_missingness_plot <- function(df, label) {
  ggplot(df, aes(x = pct_missing)) +
    geom_histogram(binwidth = 5, boundary = 0, fill = colors$blue, color = "white") +
    labs(x = "% Missing", y = "Number of Proteins", title = "Per-Protein Missingness (Raw)") +
    theme_minimal(base_size = 13) +
    panel_border_theme
}

make_sample_missingness_plot <- function(df, label, flagged_samples = character(0)) {
  ggplot(df, aes(x = reorder(sample, n_missing), y = n_missing)) +
    geom_col(fill = colors$orange, width = 0.7) +
    scale_x_discrete(labels = function(x) bold_flagged_labels(x, flagged_samples)) +
    labs(x = "Sample", y = "Number of Missing Proteins",
         title = "Per-Sample Missingness (Raw, Pre-Normalization)") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6)) +
    panel_border_theme
}

make_pca_plot <- function(pca_res, meta, group_col, label, title_suffix, group_colors = NULL,
                          run_permanova = TRUE, permanova_permutations = 999) {
  var_explained <- round(summary(pca_res)$importance[2, 1:2] * 100, 1)
  pc_df <- as.data.frame(pca_res$x[, 1:2]) %>%
    rownames_to_column("sample") %>%
    left_join(meta %>% dplyr::select(sample, group = all_of(group_col)), by = "sample")
  
  p <- ggplot(pc_df, aes(x = PC1, y = PC2, color = group)) +
    stat_ellipse(aes(fill = group), type = "norm", level = 0.95,
                 geom = "polygon", alpha = 0.12, color = NA, show.legend = FALSE) +
    geom_point(size = 3) +
    labs(x = sprintf("PC1 (%.1f%%)", var_explained[1]),
         y = sprintf("PC2 (%.1f%%)", var_explained[2]),
         title = title_suffix, color = NULL) +
    theme_minimal(base_size = 11) +
    panel_border_theme
  
  if (!is.null(group_colors)) {
    p <- p +
      scale_color_manual(values = group_colors) +
      scale_fill_manual(values = group_colors, guide = "none")
  }
  
  if (run_permanova) {
    perm_df <- as.data.frame(pca_res$x) %>%
      rownames_to_column("sample") %>%
      left_join(meta %>% dplyr::select(sample, group = all_of(group_col)), by = "sample") %>%
      filter(!is.na(group))
    
    permanova_label <- tryCatch({
      dist_mat <- vegan::vegdist(perm_df %>% dplyr::select(-sample, -group), method = "euclidean")
      set.seed(42)
      res <- vegan::adonis2(dist_mat ~ group, data = perm_df, permutations = permanova_permutations)
      r2_val  <- round(res$R2[1], 3)
      p_val   <- res$`Pr(>F)`[1]
      p_label <- if (p_val < 0.001) "p < 0.001" else sprintf("p = %.3f", p_val)
      sprintf("PERMANOVA\nR\u00b2 = %.3f\n%s", r2_val, p_label)
    }, error = function(e) {
      message(sprintf("PERMANOVA failed for %s (%s): %s", label, title_suffix, e$message))
      NULL
    })
    
    if (!is.null(permanova_label)) {
      p <- p +
        annotate("label", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.1,
                 label = permanova_label, size = 3, lineheight = 0.9,
                 color = colors$gray_dark, fill = colors$white, label.size = 0.2)
    }
  }
  
  p
}

make_bm_pca_combined <- function(pca_res, meta, marker_df, group_col, label, stage,
                                 group_colors = NULL) {
  pc1_df <- as.data.frame(pca_res$x[, 1, drop = FALSE]) %>%
    rownames_to_column("sample")
  
  df <- pc1_df %>%
    left_join(marker_df %>% dplyr::select(sample, B_M_ratio), by = "sample") %>%
    left_join(meta %>% dplyr::select(sample, group = all_of(group_col)), by = "sample") %>%
    filter(!is.na(B_M_ratio), !is.na(group))
  
  fmt_p <- function(p) ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p))
  
  overall_test  <- cor.test(df$PC1, df$B_M_ratio)
  overall_label <- sprintf("All: r = %.2f, %s", overall_test$estimate, fmt_p(overall_test$p.value))
  
  group_stats <- df %>%
    group_by(group) %>%
    filter(n() >= 3) %>%
    group_modify(~ {
      test <- cor.test(.x$PC1, .x$B_M_ratio)
      tibble(r = unname(test$estimate), p = test$p.value)
    }) %>%
    ungroup() %>%
    mutate(label_txt = sprintf("%s: r = %.2f, %s", group, r, fmt_p(p)))
  
  annotation_text <- paste(c(overall_label, group_stats$label_txt), collapse = "\n")
  
  p <- ggplot(df, aes(x = PC1, y = B_M_ratio)) +
    geom_point(aes(color = group), size = 2.2, alpha = 0.85) +
    geom_smooth(aes(color = group, fill = group), method = "lm", se = TRUE,
                linewidth = 0.7, alpha = 0.12) +
    geom_smooth(method = "lm", color = colors$black, fill = colors$black,
                se = TRUE, linewidth = 1.6, alpha = 0.15) +
    labs(x = "PC1", y = "B:M Ratio",
         title = sprintf("B:M Ratio vs PC1 (%s)", stage), color = NULL) +
    theme_minimal(base_size = 11) +
    panel_border_theme +
    annotate("label", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.05,
             label = annotation_text, size = 2.8, lineheight = 0.95,
             color = colors$gray_dark, fill = colors$white, label.size = 0.2)
  
  if (!is.null(group_colors)) {
    p <- p +
      scale_color_manual(values = group_colors) +
      scale_fill_manual(values = group_colors, guide = "none")
  }
  
  p
}

# --- 6: Contaminant Classification Plot ----------------------------------

flagged_contaminants <- flagged_contaminants_all %>%
  mutate(
    contaminant_category = case_when(
      grepl("KERATIN|KRT",        description, ignore.case = TRUE) ~ "Keratins",
      grepl("SERUM|ALBUMIN",      description, ignore.case = TRUE) ~ "Serum Proteins",
      grepl("IMMUNOGLOBULIN|IGG", description, ignore.case = TRUE) ~ "Immunoglobulins",
      grepl("TRYP",               description, ignore.case = TRUE) ~ "Trypsins",
      grepl("HEMOGLOBIN",         description, ignore.case = TRUE) ~ "Hemoglobins",
      grepl("LIPOPROTEIN",        description, ignore.case = TRUE) ~ "Lipoproteins",
      TRUE                                                         ~ "Other"
    )
  )

contaminant_summary <- flagged_contaminants %>%
  count(contaminant_category) %>%
  mutate(percentage = round(n / sum(n) * 100, 1)) %>%
  arrange(contaminant_category)

p_contaminants <- ggplot(contaminant_summary, aes(x = 2, y = n, fill = contaminant_category)) +
  geom_col(color = "white") +
  coord_polar(theta = "y") +
  xlim(c(0.5, 2.5)) +
  scale_fill_brewer(
    palette = "Set2",
    labels  = paste0(contaminant_summary$contaminant_category, " (", contaminant_summary$percentage, "%)")
  ) +
  theme_void() +
  theme(legend.position = "right", legend.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  panel_border_theme +
  labs(title = "Contaminant Classification by Type", fill = NULL)

ggsave(file.path(config$report_dir, "05_contaminant_classification.pdf"), p_contaminants, width = 8, height = 12)
ggsave(file.path(config$report_dir, "05_contaminant_classification.png"), p_contaminants, width = 8, height = 12)
cat(" Saved: 05_contaminant_classification.pdf/.png\n")

# --- 7: Outlier Diagnostics ----------------------------------------------

walk(names(outlier_diag), function(nm) {
  out_dir <- config$data_dir
  write_csv(outlier_diag[[nm]]$outlier_diag,
            file.path(out_dir, sprintf("06_outlier_diagnostics_%s.csv", nm)))
})

walk(names(outlier_diag), function(nm) {
  res <- outlier_diag[[nm]]
  df  <- res$outlier_diag
  thr <- res$thresholds
  pca_res <- res$pca_res
  
  miss_threshold <- thr$miss_threshold
  sample_medians <- df$sample_median
  global_median  <- median(sample_medians, na.rm = TRUE)
  mad_val        <- mad(sample_medians, na.rm = TRUE)
  
  # A: Missingness
  p1 <- ggplot(df, aes(x = reorder(sample, pct_missing), y = pct_missing)) +
    geom_col(aes(fill = consensus_outlier)) +
    geom_hline(yintercept = miss_threshold, linetype = "dashed", color = "red", alpha = 0.5) +
    scale_fill_manual(values = c("FALSE" = "gray40", "TRUE" = "red")) +
    labs(x = "Sample", y = "% Missing", title = "A: Sample Missingness") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6), legend.position = "none") +
    panel_border_theme
  
  # B: PCA Mahalanobis
  var_explained <- round(summary(pca_res)$importance[2, 1:2] * 100, 1)
  pc_df <- as.data.frame(pca_res$x[, 1:2]) %>%
    rownames_to_column("sample") %>%
    left_join(df %>% dplyr::select(sample, consensus_outlier), by = "sample")
  
  p2 <- ggplot(pc_df, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = consensus_outlier), size = 3) +
    scale_color_manual(values = c("FALSE" = "gray40", "TRUE" = "red")) +
    labs(x = sprintf("PC1 (%.1f%%)", var_explained[1]),
         y = sprintf("PC2 (%.1f%%)", var_explained[2]),
         title = "B: PCA Outliers (Mahalanobis)") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none") +
    panel_border_theme
  
  # C: MAD intensity
  p3 <- ggplot(df, aes(x = reorder(sample, sample_median), y = sample_median)) +
    geom_point(aes(color = consensus_outlier), size = 2.5) +
    geom_hline(yintercept = global_median, color = "black") +
    geom_hline(yintercept = global_median + c(-3, 3) * mad_val, linetype = "dashed", color = "red", alpha = 0.5) +
    scale_color_manual(values = c("FALSE" = "gray40", "TRUE" = "red")) +
    labs(x = "Sample", y = "Median log2 intensity", title = "C: MAD Intensity Outliers") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6), legend.position = "none") +
    panel_border_theme
  
  rep_dir <- config$report_dir
  ggsave(file.path(rep_dir, sprintf("06_outlier_diagnostic_plots_%s.pdf", nm)),
         p1 + p2 + p3 + plot_layout(ncol = 1), width = 10, height = 12)
  cat(sprintf(" Saved: 06_outlier_diagnostic_plots_%s.pdf\n", nm))
})

# --- 8: CV Assessment ----------------------------------------------------
# NOTE: BFR has a single grouping factor ("group": BFR vs HLRT), unlike the
# CR muscle script's cancer/supp/timepoint splits, so this section is
# simplified to two groups pulled dynamically from the metadata.

GROUPS <- sort(unique(meta$group))

build_cv_data <- function(data, annot, meta, groups) {
  group_map <- meta %>% dplyr::select(sample, group) %>% filter(group %in% groups)
  
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

make_cv_scatter <- function(cv_wide, x_group, y_group, n_label = 20, delta_col = NULL, label = NULL) {
  missing <- setdiff(c(x_group, y_group), colnames(cv_wide))
  if (length(missing) > 0) stop("Groups not found in cv_wide: ", paste(missing, collapse = ", "))
  
  if (is.null(delta_col)) {
    delta_col <- paste0("cv_delta_", x_group, "v", y_group)
    cv_wide   <- cv_wide %>% mutate(!!delta_col := .data[[y_group]] - .data[[x_group]])
  }
  
  df       <- cv_wide %>% filter(!is.na(.data[[x_group]]), !is.na(.data[[y_group]]))
  cor_test <- cor.test(df[[x_group]], df[[y_group]])
  r_val    <- unname(cor_test$estimate)
  p_val    <- cor_test$p.value
  p_label  <- if (p_val < 0.001) "p < 0.001" else sprintf("p = %.3f", p_val)
  delta_lim <- max(abs(df[[delta_col]]), na.rm = TRUE)
  df_label  <- df %>% slice_max(abs(.data[[delta_col]]), n = n_label)
  
  ggplot(df, aes(x = .data[[x_group]], y = .data[[y_group]], colour = .data[[delta_col]])) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey40", linewidth = 0.6) +
    geom_point(size = 1.8, alpha = 0.8) +
    geom_smooth(method = "lm", colour = colors$black, fill = "grey80", alpha = 0.15, linewidth = 0.8) +
    geom_text_repel(
      data = df_label, aes(label = gene), colour = "black", fontface = "bold", size = 2.8,
      max.overlaps = Inf, box.padding = 0.4, point.padding = 0.2,
      segment.color = "grey60", segment.size = 0.3, show.legend = FALSE, seed = 42
    ) +
    scale_colour_gradientn(
      colours = c(colors$purple_dark, colors$purple, colors$orange, colors$yellow,
                  colors$orange, colors$purple, colors$purple_dark),
      values  = scales::rescale(c(-delta_lim, -delta_lim * 0.5, -delta_lim * 0.15, 0,
                                  delta_lim * 0.15, delta_lim * 0.5, delta_lim)),
      limits  = c(-delta_lim, delta_lim),
      name    = sprintf("\u0394CV\n(%s \u2212 %s)", y_group, x_group)
    ) +
    annotate("label", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.1,
             label = sprintf("r = %.2f\n%s", r_val, p_label), size = 3.2, lineheight = 0.95,
             fontface = "bold", colour = colors$gray_dark, fill = colors$white, label.size = 0.2) +
    labs(x = bquote(CV[.(x_group)]~"(%)"), y = bquote(CV[.(y_group)]~"(%)"),
         title = sprintf("Per-protein CV: %s vs %s", x_group, y_group), subtitle = label) +
    theme_minimal() +
    panel_border_theme
}

GROUP_COLOURS <- setNames(c(colors$purple, colors$orange), GROUPS)

compute_cv_summary_stats <- function(df, groups, n_boot = 2000) {
  df <- df %>% filter(group %in% groups, !is.na(cv))
  
  pooled_skew <- with(df, mean((cv - mean(cv))^3) / sd(cv)^3)
  use_median  <- abs(pooled_skew) > 1
  
  cat(sprintf("  CV summary stat: pooled skewness = %.2f -> using %s \u00b1 95%% CI\n",
              pooled_skew, ifelse(use_median, "Median (bootstrap)", "Mean (normal-theory)")))
  
  set.seed(42)
  df %>%
    group_by(group) %>%
    group_modify(~ {
      x <- .x$cv
      n <- length(x)
      if (use_median) {
        boot_medians <- replicate(n_boot, median(sample(x, n, replace = TRUE)))
        tibble(n = n, stat_type = "Median", estimate = median(x),
               ci_lower = quantile(boot_medians, 0.025, names = FALSE),
               ci_upper = quantile(boot_medians, 0.975, names = FALSE),
               y_max = max(x))
      } else {
        ci <- Hmisc::smean.cl.normal(x, conf.int = 0.95)
        tibble(n = n, stat_type = "Mean", estimate = ci[["Mean"]],
               ci_lower = ci[["Lower"]], ci_upper = ci[["Upper"]], y_max = max(x))
      }
    }) %>%
    ungroup()
}

make_cv_violin <- function(cv_df, groups, label, n_boot = 2000) {
  df <- cv_df %>% filter(group %in% groups) %>% mutate(group = factor(group, levels = groups))
  stats_df <- compute_cv_summary_stats(df, groups, n_boot = n_boot) %>%
    mutate(group = factor(group, levels = groups))
  
  y_range  <- diff(range(df$cv, na.rm = TRUE))
  stats_df <- stats_df %>%
    mutate(label_y = y_max + y_range * 0.08,
           label_txt = sprintf("%.1f [%.1f\u2013%.1f]", estimate, ci_lower, ci_upper))
  
  ggplot(df, aes(x = group, y = cv, fill = group, color = group)) +
    geom_violin(alpha = 0.5, linewidth = 0.5, trim = FALSE) +
    geom_boxplot(width = 0.15, alpha = 0.8, linewidth = 0.5, outlier.shape = NA, colour = colors$gray_dark) +
    geom_label(data = stats_df, aes(x = group, y = label_y, label = label_txt), inherit.aes = FALSE,
               size = 2.8, lineheight = 0.9, fontface = "bold",
               color = colors$gray_dark, fill = colors$white, label.size = 0.2) +
    scale_fill_manual(values = GROUP_COLOURS[groups], guide = "none") +
    scale_colour_manual(values = GROUP_COLOURS[groups], guide = "none") +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0.02, 0.2))) +
    labs(x = NULL, y = "CV (%)", title = "Per-protein CV by group") +
    theme_minimal() +
    panel_border_theme
}

intensity_raw <- report_data$dals_pre_outlier$BFR$data

cv_df_raw <- build_cv_data(dals_raw$BFR$data, dals_raw$BFR$annotation, dals_raw$BFR$metadata, GROUPS)
cv_df_bfr <- build_cv_data(dals$BFR$data,     dals$BFR$annotation,     dals$BFR$metadata,     GROUPS)
cv_df_bfr_prenorm <- build_cv_data(dal_prenorm_BFR$data, dal_prenorm_BFR$annotation, dal_prenorm_BFR$metadata, GROUPS)

build_cv_wide <- function(cv_df, groups) {
  cv_df %>%
    filter(group %in% groups) %>%
    dplyr::select(gene, group, cv) %>%
    pivot_wider(names_from = group, values_from = cv) %>%
    mutate(cv_delta = .data[[groups[2]]] - .data[[groups[1]]])
}

cv_wide_bfr     <- build_cv_wide(cv_df_bfr, GROUPS)
cv_wide_bfr_raw <- build_cv_wide(cv_df_raw, GROUPS)

cv_violin_specs <- list(
  list(data = cv_df_bfr_prenorm, groups = GROUPS, label = "pre-normalization", file = "07_cv_violin_BFR_prenorm.pdf"),
  list(data = cv_df_bfr,         groups = GROUPS, label = "normalized",       file = "07_cv_violin_BFR_normalized.pdf")
)

walk(cv_violin_specs, function(s) {
  p <- make_cv_violin(s$data, s$groups, s$label)
  ggsave(file.path(config$report_dir, s$file), p, width = 6, height = 6, device = cairo_pdf)
  cat(sprintf(" Saved: %s\n", s$file))
})

cv_scatter_specs <- list(
  list(data = cv_wide_bfr,     label = "normalized", n = 30),
  list(data = cv_wide_bfr_raw, label = "raw",        n = 50)
)

walk(cv_scatter_specs, function(s) {
  plot_name <- sprintf("08_cv_scatter_%s_v_%s_%s.pdf", GROUPS[1], GROUPS[2], s$label)
  ggsave(
    file.path(config$report_dir, plot_name),
    make_cv_scatter(s$data, GROUPS[1], GROUPS[2], delta_col = "cv_delta", label = s$label, n_label = s$n),
    width = 7, height = 6, device = cairo_pdf
  )
  cat(sprintf(" Saved: %s\n", plot_name))
})

# --- 8b: Combined CV Report (scatter above violin) -------------------------

make_cv_report <- function(cv_wide, cv_violin_df, groups, stage_label, n_label) {
  scatter_subtitle <- sub("^[A-Za-z0-9]+,\\s*", "", stage_label)
  p_scatter <- make_cv_scatter(cv_wide, groups[1], groups[2], delta_col = "cv_delta",
                               label = scatter_subtitle, n_label = n_label)
  p_violin  <- make_cv_violin(cv_violin_df, groups, stage_label)
  p_scatter / p_violin + plot_layout(heights = c(1, 1))
}

cv_report_normalized <- make_cv_report(cv_wide_bfr, cv_df_bfr, GROUPS,
                                       stage_label = "BFR, normalized", n_label = 30)
cv_report_prenorm <- make_cv_report(cv_wide_bfr_raw, cv_df_bfr_prenorm, GROUPS,
                                    stage_label = "BFR, pre-normalization", n_label = 30)

ggsave(file.path(config$report_dir, "14_cv_report_BFR_normalized.pdf"), cv_report_normalized, width = 8, height = 12, device = cairo_pdf)
ggsave(file.path(config$report_dir, "14_cv_report_BFR_prenorm.pdf"),    cv_report_prenorm,    width = 8, height = 12, device = cairo_pdf)

cat(" Saved: 14_cv_report_BFR_normalized.pdf\n")
cat(" Saved: 14_cv_report_BFR_prenorm.pdf\n")

# --- 9: Protein Assessment -----------------------------------------------

protein_summary <- map_dfr(names(dals), function(nm) {
  d        <- dals[[nm]]
  miss_pct <- rowMeans(is.na(d$data))
  tibble(dataset = nm, Total = nrow(d$data), Complete = sum(miss_pct == 0), Incomplete = sum(miss_pct > 0))
})

write.csv(protein_summary, file.path(config$data_dir, "11_protein_report_all.csv"), row.names = FALSE)
cat(" Saved: 11_protein_report_all.csv\n")
print(protein_summary)

walk(names(dals), function(nm) {
  out_dir  <- config$data_dir
  d        <- dals[[nm]]
  miss_pct <- rowMeans(is.na(d$data))
  report   <- tibble(Total = nrow(d$data), Complete = sum(miss_pct == 0), Incomplete = sum(miss_pct > 0))
  write.csv(report, file.path(out_dir, sprintf("11_protein_report_%s.csv", nm)), row.names = FALSE)
  cat(sprintf("  [%s] Total: %d | Complete: %d | Incomplete: %d\n", nm, report$Total, report$Complete, report$Incomplete))
})

# --- 10: Sample Assessment -----------------------------------------------

compute_sample_report <- function(meta, nm) {
  tibble(
    dataset       = nm,
    total_samples = nrow(meta),
    n_rerun       = sum(meta$is_rerun, na.rm = TRUE),
    !!!setNames(
      as.list(table(meta$group)[GROUPS]),
      paste0("n_", GROUPS)
    )
  )
}

sample_report_all <- map_dfr(names(dals), function(nm) compute_sample_report(dals[[nm]]$metadata, nm))

write.csv(sample_report_all, file.path(config$report_dir, "12_sample_report_all.csv"), row.names = FALSE)
cat(" Saved: 12_sample_report_all.csv\n")
print(sample_report_all)

walk(names(dals), function(nm) {
  rep_dir <- config$report_dir
  write.csv(compute_sample_report(dals[[nm]]$metadata, nm),
            file.path(rep_dir, sprintf("12_sample_report_%s.csv", nm)), row.names = FALSE)
})

writeLines(capture.output(sessionInfo()), file.path(config$data_dir, "sessionInfo_reports.txt"))

cat("\n>> Reporting complete.\n")

# --- 11: Combined Normalization Summary (per dataset) ---------------------
cat("\n>> 11 - Combined normalization summary\n")

pca_group_colors <- GROUP_COLOURS

walk(names(dals), function(nm) {
  rep_dir <- config$report_dir
  
  p_filter     <- make_filter_plot(filter_logs[[nm]], nm)
  p_protmiss   <- make_protein_missingness_plot(missingness_raw[[nm]]$protein_missing, nm)
  p_samplemiss <- make_sample_missingness_plot(missingness_raw[[nm]]$sample_missing, nm,
                                               flagged_samples = removed_samples[[nm]])
  p_pca_raw    <- make_pca_plot(pca_raw[[nm]], dals_raw[[nm]]$metadata, grouping_cols[[nm]], nm,
                                "PCA — Raw", group_colors = pca_group_colors)
  p_pca_norm   <- make_pca_plot(pca_postnorm[[nm]], dals[[nm]]$metadata, grouping_cols[[nm]], nm,
                                "PCA — Post-Normalization", group_colors = pca_group_colors)
  
  ggsave(file.path(rep_dir, sprintf("09_protein_missingness_%s.pdf", nm)), p_protmiss, width = 7, height = 5)
  ggsave(file.path(rep_dir, sprintf("09_pca_raw_%s.pdf", nm)), p_pca_raw, width = 6, height = 5)
  ggsave(file.path(rep_dir, sprintf("09_pca_postnorm_%s.pdf", nm)), p_pca_norm, width = 6, height = 5)
  ggsave(file.path(rep_dir, sprintf("09_sample_missingness_%s.pdf", nm)), p_samplemiss, width = 10, height = 5)
  
  combined <- (p_filter + p_protmiss) /
    (p_pca_raw + p_pca_norm) /
    p_samplemiss +
    plot_layout(heights = c(1, 1, 1))
  
  ggsave(file.path(rep_dir, sprintf("10_normalization_summary_%s.pdf", nm)), combined, width = 14, height = 20)
  
  cat(sprintf(" Saved: 10_normalization_summary_%s.pdf (+ individual panels)\n", nm))
})

# --- 12: Contamination Report ---------------------------------------------
cat("\n>> 12 - Contamination report\n")

walk(names(dals), function(nm) {
  rep_dir <- config$report_dir
  grp_col <- grouping_cols[[nm]]
  
  p_bm_nm <- make_bm_plot(marker_pct_df, dals_raw[[nm]]$metadata$sample, nm,
                          flagged_samples = removed_samples[[nm]])
  
  p_pca_raw_combined  <- make_bm_pca_combined(pca_raw[[nm]], dals_raw[[nm]]$metadata, marker_pct_df,
                                              group_col = grp_col, label = nm, stage = "Raw",
                                              group_colors = pca_group_colors)
  p_pca_norm_combined <- make_bm_pca_combined(pca_postnorm[[nm]], dals[[nm]]$metadata, marker_pct_df,
                                              group_col = grp_col, label = nm, stage = "Post-Normalization",
                                              group_colors = pca_group_colors)
  
  ggsave(file.path(rep_dir, sprintf("13_bm_pca_raw_combined_%s.pdf", nm)), p_pca_raw_combined, width = 7, height = 6)
  ggsave(file.path(rep_dir, sprintf("13_bm_pca_postnorm_combined_%s.pdf", nm)), p_pca_norm_combined, width = 7, height = 6)
  
  combined_contam <- p_bm_nm /
    (p_pca_raw_combined + p_pca_norm_combined) +
    plot_layout(heights = c(1.4, 1))
  
  ggsave(file.path(rep_dir, sprintf("13_contamination_report_%s.pdf", nm)), combined_contam, width = 14, height = 14)
  
  cat(sprintf(" Saved: 13_contamination_report_%s.pdf (+ individual panels)\n", nm))
})