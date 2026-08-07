# Fig 03 — Volcano plots (BFR)
#
# ADAPTED FROM THE CR MUSCLE/PLASMA VERSION: that script faceted by
# tissue (rows) x comparison (columns) for 3 SURV comparisons.
#
# RESTRUCTURED (this version): rather than one ggplot faceted across all
# 5 comparisons, each comparison is now built as its own standalone
# volcano plot (identical styling/logic to what each facet produced
# before -- same points, labels, UP/DOWN counts, color scale), and the 5
# are assembled into a lettered A-E grid via patchwork. Facet-specific
# machinery (make_strip, facet_wrap2) is gone since there's no longer a
# single faceted plot to theme strips for -- each panel gets its own
# title instead of a facet strip label.

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

# --- 0f: Per-contrast background shading (edit these to change panel colors) ---
# One solid background tint per panel, keyed by comparison name -- makes
# it easy to tell the 5 panels apart at a glance in the combined figure.
# Swap any of these for a different colors$... entry from 00_theme.R.
PANEL_BG_COLORS <- c(
  Baseline      = colors$gray,
  Training_BFR  = colors$red_light,
  Training_HLRT = colors$blue_light,
  Post_training = colors$green_light,
  Interaction   = colors$purple_light
)
PANEL_BG_ALPHA <- 0.15  # opacity -- higher = more visible tint

# --- 3: Volcano plot function (single comparison) ---------------
# Builds ONE panel for ONE comparison -- same internals as the old
# faceted version's per-facet content, just without facet_wrap2. The
# comparison name becomes a plot title instead of a facet strip label.
make_volcano_plot <- function(df_vol, comparison_name, title, y_limit,
                              n_labels = 5, pval_thresh = 0.05, lfc_thresh = 0.6) {
  df <- df_vol %>% filter(comparison == comparison_name)
  
  df_label <- df %>%
    filter(DE.pi != "NDE") %>%
    group_by(DE.pi) %>%
    slice_max(abs(logFC), n = n_labels) %>%
    ungroup()
  
  df_counts <- df %>%
    filter(DE.pi != "NDE") %>%
    group_by(DE.pi) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(
      x     = ifelse(DE.pi == "UP",  Inf, -Inf),
      y     = Inf,
      hjust = ifelse(DE.pi == "UP",  1.2, -0.2),
      label = as.character(n),
      color = ifelse(DE.pi == "UP", colors$red, colors$blue)
    )
  
  ggplot(df, aes(logFC, -log10(P.Value), color = DE.pi)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
             fill = PANEL_BG_COLORS[[comparison_name]], alpha = PANEL_BG_ALPHA) +
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
                       limits   = c("UP", "DOWN"),
                       na.value = colors$gray, name = "Direction", drop = FALSE) +
    guides(color = guide_legend(override.aes = list(size = 4))) +
    coord_cartesian(xlim = c(-2.5, 2.5), ylim = c(0, y_limit)) +
    labs(x = bquote(bold("log2(Fold Change)")), y = bquote(bold("-log10(p-value)")),
         title = title) +
    theme_cr(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 11, hjust = 0.5))
}

# --- 4: Build 5 panels and assemble A-E grid ---------------------------
comparisons_ordered <- c("Baseline", "Training_BFR", "Training_HLRT",
                         "Post_training", "Interaction")

# Shared y-axis limit across ALL 5 contrasts (not per-panel), so every
# panel uses the same -log10(p) scale and bar heights/point positions are
# directly comparable at a glance. A small multiplier leaves headroom
# above the tallest point so it isn't flush against the panel edge.
y_limit <- max(-log10(df_vol$P.Value), na.rm = TRUE) * 1.05

volcano_panels <- lapply(comparisons_ordered, function(cmp) {
  make_volcano_plot(df_vol, cmp, title = comp_labels[[cmp]], y_limit = y_limit)
})

# guides = "collect" is supposed to merge identical legends across panels
# automatically, but in practice it can fail to recognize two legends as
# identical if the underlying DE.pi values aren't perfectly consistent
# across contrasts (stray capitalization, an extra category present in
# only one sheet, etc.) -- rather than debug that per-dataset quirk, this
# builds ONE legend from a clean synthetic UP/DOWN dataframe (guaranteed
# correct, independent of whatever inconsistency is in the real data) and
# attaches it manually instead of relying on automatic collection.
#
# Custom layout: A/B/C across the top row, D/E centered on the bottom row
# with blank ("#") cells flanking them. wrap_plots() matches the plot
# list to design letters in order (volcano_panels is already in
# comparisons_ordered order: Baseline=A, Training_BFR=B, Training_HLRT=C,
# Post_training=D, Interaction=E).
layout_design <- "
AAAABBBBCCCC
##DDDDEEEE##
"

volcano_grid <- wrap_plots(volcano_panels) +
  plot_layout(design = layout_design) +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "none")

legend_source <- ggplot(
  tibble(Direction = factor(c("UP", "DOWN"), levels = c("UP", "DOWN"))),
  aes(x = 1, y = 1, color = Direction)
) +
  geom_point(size = 4) +
  scale_color_manual(values = c(UP = colors$red, DOWN = colors$blue), name = "Direction") +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme_cr(base_size = 10) +
  theme(legend.position = "right")

shared_legend <- cowplot::get_legend(legend_source)

volcano_all <- cowplot::plot_grid(
  volcano_grid, shared_legend,
  ncol = 2, rel_widths = c(1, 0.08)
)

ggsave(file.path(config$fig_dir, "panel_3A.pdf"), volcano_all, width = 16, height = 9, device = "pdf")
ggsave(file.path(config$fig_dir, "panel_3A.png"), volcano_all, width = 16, height = 9, device = "png", dpi = 300)