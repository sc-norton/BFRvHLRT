# Fig 03 — Volcano plots (BFR)
#
# ADAPTED FROM THE CR MUSCLE/PLASMA VERSION: that script faceted by
# tissue (rows) x comparison (columns) for 3 SURV comparisons. BFR is a
# single dataset with no tissue split, so this is a single-row facet over
# BFR's 5 comparisons instead.

# --- 0: Setup ------------------------
setwd(rprojroot::find_rstudio_root_file())
getwd()
# Load packages
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(readxl, readr, dplyr, tidyr, tibble, stringr, purrr, ggplot2,
               ggrepel, ggvenn, patchwork, ggtext, cowplot, ggh4x, ggalluvial,
               scales, gtsummary, limma, Hmisc, nlme, ggpubr, tidytext, ggsankey,
               enrichplot)

# Configuration
config <- list(
  # Data files
  limma_data = "03_analysis/c_data/limma/limma_results.xlsx",
  # Directories
  fig_dir  = "04_figures/fig_03",
  data_dir = "04_figures/data/fig_03"
)

dir.create(config$fig_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(config$data_dir, showWarnings = FALSE, recursive = TRUE)

# --- Color palette + theme ----------------------------------------------------
source("04_figures/a_scripts/00_theme.R")

make_strip <- function(fills, type = "x") {
  n <- length(fills)
  if (type == "x") {
    strip_themed(
      background_x = elem_list_rect(fill = fills, linewidth = rep(0, n)),
      text_x       = elem_list_text(face = rep("bold", n))
    )
  } else if (type == "xy") {
    strip_themed(
      background_x = elem_list_rect(fill = fills$x, linewidth = rep(0, length(fills$x))),
      background_y = elem_list_rect(fill = fills$y, linewidth = rep(0, length(fills$y))),
      text_x       = elem_list_text(face = rep("bold", length(fills$x))),
      text_y       = elem_list_text(face = rep("bold", length(fills$y)))
    )
  }
}

# --- 0e: Comparison Labels ---------------------------------------
comp_labels <- c(
  Baseline      = "Baseline",
  Training_BFR  = "Training_BFR",
  Training_HLRT = "Training_HLRT",
  Post_training = "Post_training",
  Interaction   = "Interaction"
)

# --- 1: Load Data --------------------------
# NOTE: assumes the limma workbook's 5 sheets are in this order (matching
# the GSEA/fig_02 scripts' load order). Adjust if your workbook differs.
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

# --- 2: Restructure Data -------------------
# --- 2a: Volcano data ---
# No tissue dimension to loop over -- just tag each sheet with its
# comparison name and bind.
df_vol <- bind_rows(lapply(names(limma_results), function(nm) {
  limma_results[[nm]] %>% mutate(comparison = nm)
}))

# --- 3: Volcano plot function ---------------
make_volcano_plot <- function(df_vol, selected_comparisons,
                              n_labels = 5, pval_thresh = 0.05, lfc_thresh = 0.6) {
  df <- df_vol %>%
    filter(comparison %in% selected_comparisons) %>%
    mutate(comparison = factor(comparison, levels = selected_comparisons))

  df_label <- df %>%
    mutate(nudge_x = ifelse(DE.pi == "UP", 2.5, -2.5), nudge_y = 0.5) %>%
    filter(DE.pi != "NDE") %>%
    group_by(comparison, DE.pi) %>%
    slice_max(abs(logFC), n = n_labels) %>%
    ungroup()

  # Count UP and DOWN per facet
  df_counts <- df %>%
    filter(DE.pi != "NDE") %>%
    group_by(comparison, DE.pi) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(
      x     = ifelse(DE.pi == "UP",  Inf, -Inf),
      y     = Inf,
      hjust = ifelse(DE.pi == "UP",  1.2, -0.2),
      label = as.character(n),
      color = ifelse(DE.pi == "UP", colors$red, colors$blue)
    )

  strip <- make_strip(
    fills = rep(colors$light_gray, length(selected_comparisons)),
    type  = "x"
  )

  ggplot(df, aes(logFC, -log10(P.Value), color = DE.pi)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, color = colors$gray) +
    geom_point(size = 1.2, alpha = 0.8) +
    geom_text_repel(data = df_label, aes(label = gene),
                    color = "black", fontface = "bold", size = 3,
                    show.legend = FALSE, max.overlaps = Inf,
                    box.padding = 0.4, point.padding = 0.3,
                    force = 2, force_pull = 0.01,
                    min.segment.length = 0.05, segment.color = "black",
                    arrow = arrow(length = unit(0.01, "npc")), seed = 42) +
    geom_text(
      data        = df_counts,
      aes(x = x, y = y, label = label, hjust = hjust),
      vjust       = 1.5,
      fontface    = "bold",
      size        = 3,
      color       = df_counts$color,
      inherit.aes = FALSE
    ) +
    scale_color_manual(values   = c(UP = colors$red, DOWN = colors$blue),
                       na.value = colors$gray) +
    guides(color = guide_legend(override.aes = list(size = 4))) +
    coord_cartesian(xlim = c(-4, 4)) +
    facet_wrap2(~ comparison, nrow = 1, strip = strip,
               labeller = labeller(comparison = comp_labels)) +
    labs(x = bquote(bold("log2(Fold Change)")), y = bquote(bold("-log10(p-value)")),
         color = "Direction") +
    theme_cr(base_size = 10)
}

# Build volcano
volcano_all <- make_volcano_plot(
  df_vol,
  c("Baseline", "Training_BFR", "Training_HLRT", "Post_training", "Interaction")
)

ggsave(file.path(config$fig_dir, "panel_3A.pdf"), volcano_all, width = 16, height = 4, device = "pdf")
ggsave(file.path(config$fig_dir, "panel_3A.png"), volcano_all, width = 16, height = 4, device = "png", dpi = 300)
