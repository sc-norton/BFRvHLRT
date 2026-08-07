# Fig GSEA Bubble — Hallmark + GO:BP, pathway x contrast (BFR)
#
# ADAPTED FROM THE CR MUSCLE/PLASMA VERSION, with one structural change
# beyond tissue removal: that version loaded raw per-comparison GSEA RDS
# files and rebuilt the full pathway table + redundancy filter itself
# (its sections 1-3c). That's exactly the logic 02_gsea_consolidation_BFR.R
# now centralizes into gsea_pathway_data.rds -- the same file fig_02_BFR.R
# already reads from. This version does the same: load the consolidated
# RDS, and skip straight to bubble-plot-specific selection/plotting logic.
#
# Also: single dataset (group = "BFR"), no tissue split, and restricted to
# Training_BFR / Training_HLRT only (not all 5 BFR contrasts) so every
# plot in this script compares just those two training responses.

# --- 0: Setup ------------------------
setwd(rprojroot::find_rstudio_root_file())
getwd()
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(readxl, readr, dplyr, tidyr, tibble, stringr, purrr, ggplot2,
               patchwork, cowplot, scales, enrichplot, writexl)

# --- Configuration -----------------------------------------------------------
config <- list(
  gsea_data = "03_analysis/c_data/GSEA/consolidated/gsea_pathway_data.rds",
  fig_dir   = "04_figures/fig_gsea_bubble",
  data_dir  = "04_figures/data/fig_gsea_bubble"
)
dir.create(config$fig_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(config$data_dir, showWarnings = FALSE, recursive = TRUE)

# --- Color palette + theme ------------------------------------------------
source("04_figures/a_scripts/00_theme.R")

# --- Shared constants --------------------------------------------------------
adjp_thresh <- 0.05

gsea_group <- "BFR"  # matches the group column added by the consolidation script

# Restricted to the two training-response contrasts -- every plot in this
# script's x-axis (and pathway selection, since selection is ranked from
# gsea_supp_full below, which is filtered to just these two) compares
# only Training_BFR vs Training_HLRT.
gsea_comps <- c("Training_BFR", "Training_HLRT")

comp_labels <- c(
  Training_BFR  = "Training_BFR",
  Training_HLRT = "Training_HLRT"
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

# Panels for THIS figure -- Hallmark + GO:BP only.
bubble_slots <- c(gsea_hall = "Hallmark", gsea_bp_clean = "GO:BP")

# NES fill uses the same Up/Down red-blue convention as DIR_COLORS in fig_02.
NES_LOW  <- colors$blue   # negative NES
NES_HIGH <- colors$red    # positive NES

# --- 1: Load consolidated GSEA data -------------------------------------
# Everything sections 1-3c did in the original version (loading raw RDS,
# building the full pathway table, significance flagging, the pooled
# greedy Jaccard/overlap redundancy filter) already happened once in
# 02_gsea_consolidation_BFR.R. Nothing here recomputes any of it.
if (!file.exists(config$gsea_data)) {
  stop("Consolidated GSEA data not found at ", config$gsea_data,
       " -- run the GSEA consolidation script first.")
}
gsea_pathway_data <- readRDS(config$gsea_data)

gsea_supp_full         <- gsea_pathway_data$gsea_supp_full         %>% filter(group == gsea_group)
representative_lookup  <- gsea_pathway_data$representative_lookup  %>% filter(group == gsea_group)

# Restrict to Training_BFR / Training_HLRT ONLY, right at the source --
# every selection/ranking function below reads from gsea_supp_full, so
# filtering once here means pathway selection itself is based purely on
# these two contrasts' evidence rather than being influenced by
# Baseline/Post_training/Interaction significance. Explicit factor levels
# guarantee consistent x-axis ordering (Training_BFR then Training_HLRT)
# across every plot in this script.
gsea_supp_full <- gsea_supp_full %>%
  filter(comparison %in% gsea_comps) %>%
  mutate(comparison = factor(comparison, levels = gsea_comps))

cat("=== Loaded consolidated GSEA data (group = '", gsea_group, "') ===\n", sep = "")
cat("Rows:", nrow(gsea_supp_full), "\n")
print(table(gsea_supp_full$database, gsea_supp_full$comparison))

stopifnot(
  "gsea_supp_full comparisons must match gsea_comps" =
    all(gsea_comps %in% unique(as.character(gsea_supp_full$comparison))),
  "gsea_supp_full must contain all bubble_slots databases" =
    all(unname(bubble_slots) %in% unique(as.character(gsea_supp_full$database)))
)

# --- 2: Select pathways to plot (Hallmark + GO:BP, per-database TOP_N) ------
# adjp_thresh (shared constant, 0.05)  : "significant in >=1 contrast" gate
# representative_lookup (from consolidation) : drops non-representative
#                                               pathways from redundant
#                                               clusters BEFORE ranking
# TOP_N_HALLMARK / TOP_N_GOBP          : cap per database, ranked by best p.adjust
TOP_N_HALLMARK <- 30
TOP_N_GOBP     <- 30
TOP_N <- c(Hallmark = TOP_N_HALLMARK, `GO:BP` = TOP_N_GOBP)

select_pathways <- function(database_name, top_n) {
  df <- gsea_supp_full %>% filter(database == database_name)
  
  ranked <- df %>%
    group_by(ID, Description_clean) %>%
    summarise(min_padj = min(p.adjust, na.rm = TRUE), .groups = "drop") %>%
    filter(min_padj < adjp_thresh) %>%
    arrange(min_padj)
  
  if (nrow(ranked) == 0) {
    warning("No pathways pass adjp_thresh = ", adjp_thresh, " for ", database_name,
            " -- falling back to top ", top_n, " by p.adjust regardless of significance.")
    ranked <- df %>%
      group_by(ID, Description_clean) %>%
      summarise(min_padj = min(p.adjust, na.rm = TRUE), .groups = "drop") %>%
      arrange(min_padj)
  } else {
    # Restrict candidates to survivors of the pooled greedy redundancy
    # filter BEFORE capping to top_n.
    rep_ids <- representative_lookup %>% filter(database == database_name) %>% pull(ID)
    if (length(rep_ids) > 0) {
      n_before <- nrow(ranked)
      ranked <- ranked %>% filter(ID %in% rep_ids)
      cat("  [", database_name, "] redundancy filter: ",
          n_before, " -> ", nrow(ranked), " candidate pathways\n", sep = "")
    }
  }
  
  ranked %>% slice_head(n = top_n) %>% pull(Description_clean)
}

build_bubble_selection <- function() {
  keep <- bind_rows(lapply(names(TOP_N), function(db) {
    tibble(database = db, Description_clean = select_pathways(db, TOP_N[[db]]))
  }))
  
  gsea_supp_full %>%
    filter(database %in% names(TOP_N)) %>%
    semi_join(keep, by = c("database", "Description_clean"))
}

# --- 3: Order pathways within each database (hierarchical clustering on NES) -
order_pathways <- function(df, key_col = "Description_clean") {
  nes_mat <- df %>%
    dplyr::select(all_of(key_col), comparison, NES) %>%
    pivot_wider(names_from = comparison, values_from = NES) %>%
    column_to_rownames(key_col) %>%
    as.matrix()
  
  nes_mat[is.na(nes_mat)] <- 0
  
  if (nrow(nes_mat) > 2) {
    hc <- hclust(dist(nes_mat), method = "average")
    rownames(nes_mat)[hc$order]
  } else {
    rownames(nes_mat)
  }
}

# --- 4: Plot ---------------------------------------------------------------
make_gsea_bubble_plot <- function() {
  sel <- build_bubble_selection()
  
  ordered_levels <- unlist(lapply(names(TOP_N), function(db) {
    order_pathways(filter(sel, database == db))
  }))
  
  sel <- sel %>%
    mutate(
      Description_clean = factor(Description_clean, levels = ordered_levels),
      significant = p.adjust < adjp_thresh
    )
  
  n_pathways <- length(ordered_levels)
  
  p <- ggplot(sel, aes(x = comparison, y = Description_clean)) +
    geom_point(aes(size = -log10(p.adjust), fill = NES,
                   color = significant, stroke = significant),
               shape = 21) +
    scale_fill_gradient2(
      low = NES_LOW, mid = colors$white, high = NES_HIGH, midpoint = 0,
      name = "NES"
    ) +
    scale_color_manual(
      values = c(`TRUE` = colors$black, `FALSE` = colors$gray_light),
      labels = c(`TRUE` = paste0("FDR < ", adjp_thresh),
                 `FALSE` = paste0("FDR \u2265 ", adjp_thresh)), 
      name = "Significance"
    ) +
    scale_discrete_manual(
      aesthetics = "stroke",
      values = c(`TRUE` = 1, `FALSE` = 0.3),
      guide = "none"
    ) +
    scale_size_continuous(name = expression(-log[10]~"(FDR)"), range = c(1, 7)) +
    scale_x_discrete(labels = comp_labels) +
    facet_grid(rows = vars(database), scales = "free_y", space = "free_y") +
    labs(title = "GSEA \u2014 BFR vs HLRT (Hallmark + GO:BP)", x = NULL, y = NULL) +
    theme_cr(base_size = 9, title_hjust = 0) +
    theme(
      axis.text.y      = element_text(size = 7, face = "plain"),
      axis.text.x      = element_text(size = 8, face = "bold"),
      legend.position   = "right",
      panel.grid.major  = element_line(color = colors$gray_light),
      panel.grid.minor  = element_blank(),
      strip.background  = element_rect(fill = colors$gray_light, color = NA),
      strip.text.y      = element_text(angle = 0)
    )
  
  list(plot = p, n_pathways = n_pathways)
}

fig_bfr <- make_gsea_bubble_plot()

plot_height <- max(4, 0.22 * fig_bfr$n_pathways + 1.5)

ggsave(file.path(config$fig_dir, "gsea_bubble_BFR.pdf"),
       fig_bfr$plot, width = 12, height = plot_height, device = cairo_pdf)
ggsave(file.path(config$fig_dir, "gsea_bubble_BFR.png"),
       fig_bfr$plot, width = 12, height = plot_height, dpi = 300, bg = "white")

cat("\nSaved to:", config$fig_dir, "\n")

# --- 5: All-database top-N bubble plot -----------------------------------------
# Pools ALL databases (db_slots, not just bubble_slots) and selects a
# single GLOBAL top N by p.adjust across every database combined --
# database becomes purely a facet label rather than each database getting
# its own quota. A database with unusually strong/abundant signal (e.g.
# GO:BP, just by virtue of having many more granular terms than Hallmark)
# can end up supplying most or all of the 30 shown; that's intended, not a
# bug -- use the per-database TOP_N figure above if you want every
# database guaranteed some representation regardless of relative signal.
#
# Redundancy filtering is still applied WITHIN each database before the
# pooled ranking -- comparing gene-set overlap between, say, a Hallmark
# term and a KEGG term isn't meaningful, so cross-database redundancy is
# out of scope; only the final top-N selection is pooled across databases.
TOP_N_ALL_DB <- 30

# "padj" -- rank significant survivors by p.adjust (most statistically
#           confident first). Safer default.
# "NES"  -- rank significant survivors by |NES| among pathways already
#           cleared for significance. NOT a way to bypass the significance
#           gate -- that still applies either way, only the sort order
#           among survivors changes.
RANK_METRIC <- "padj"

select_pathways_global <- function(top_n) {
  ranked <- gsea_supp_full %>%
    group_by(database, ID, Description_clean) %>%
    summarise(
      min_padj    = min(p.adjust, na.rm = TRUE),
      max_abs_nes = max(abs(NES), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    # Restrict to redundancy-filter survivors within each database (falls
    # through unfiltered for any database with 0 significant pathways).
    dplyr::group_by(database) %>%
    dplyr::group_modify(~ {
      db_name <- as.character(.y$database)
      rep_ids <- representative_lookup %>% filter(database == db_name) %>% pull(ID)
      if (length(rep_ids) == 0) return(.x)
      dplyr::filter(.x, ID %in% rep_ids)
    }) %>%
    dplyr::ungroup() %>%
    filter(min_padj < adjp_thresh)
  
  ranked <- if (RANK_METRIC == "NES") {
    ranked %>% arrange(desc(max_abs_nes))
  } else {
    ranked %>% arrange(min_padj)
  }
  
  if (nrow(ranked) == 0) {
    warning("No pathways pass adjp_thresh = ", adjp_thresh, " in any database.")
    return(ranked %>% dplyr::select(database, Description_clean))
  }
  
  ranked %>% slice_head(n = top_n) %>% dplyr::select(database, Description_clean)
}

build_bubble_selection_global <- function() {
  keep <- select_pathways_global(TOP_N_ALL_DB)
  gsea_supp_full %>% semi_join(keep, by = c("database", "Description_clean"))
}

make_gsea_bubble_plot_alldb <- function() {
  sel <- build_bubble_selection_global()
  
  present_dbs <- unname(db_slots)[unname(db_slots) %in% unique(as.character(sel$database))]
  
  # Pooling across databases means two DIFFERENT pathways can share the
  # same cleaned display text purely by coincidence -- pathway_key
  # (database::ID) is unique by construction; Description_clean is kept
  # only as the printed label via scale_y_discrete(labels = ...) below.
  sel <- sel %>% mutate(pathway_key = paste(database, ID, sep = "::"))
  
  ordered_levels <- unlist(lapply(present_dbs, function(db) {
    order_pathways(filter(sel, database == db), key_col = "pathway_key")
  }))
  
  label_lookup <- setNames(sel$Description_clean, sel$pathway_key)
  
  sel <- sel %>%
    mutate(
      database    = factor(as.character(database), levels = present_dbs),
      pathway_key = factor(pathway_key, levels = ordered_levels),
      significant = p.adjust < adjp_thresh
    )
  
  n_pathways <- length(ordered_levels)
  
  p <- ggplot(sel, aes(x = comparison, y = pathway_key)) +
    geom_point(aes(size = -log10(p.adjust), fill = NES,
                   color = significant, stroke = significant),
               shape = 21) +
    scale_fill_gradient2(
      low = NES_LOW, mid = colors$white, high = NES_HIGH, midpoint = 0,
      name = "NES"
    ) +
    scale_color_manual(
      values = c(`TRUE` = colors$black, `FALSE` = colors$gray_light),
      labels = c(`TRUE` = paste0("FDR < ", adjp_thresh),
                 `FALSE` = paste0("FDR \u2265 ", adjp_thresh)),
      name = "Significance"
    ) +
    scale_discrete_manual(
      aesthetics = "stroke",
      values = c(`TRUE` = 1, `FALSE` = 0.3),
      guide = "none"
    ) +
    scale_size_continuous(name = expression(-log[10]~"(FDR)"), range = c(1, 7)) +
    scale_x_discrete(labels = comp_labels) +
    scale_y_discrete(labels = label_lookup) +
    facet_grid(rows = vars(database), scales = "free_y", space = "free_y") +
    labs(title = paste0("GSEA \u2014 BFR vs HLRT (Top ", TOP_N_ALL_DB, ", pooled across all databases)"),
         x = NULL, y = NULL) +
    theme_cr(base_size = 9, title_hjust = 0) +
    theme(
      axis.text.y      = element_text(size = 7, face = "plain"),
      axis.text.x      = element_text(size = 8, face = "bold"),
      legend.position  = "right",
      panel.grid.major = element_line(color = colors$gray_light),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = colors$gray_light, color = NA),
      strip.text.y     = element_text(angle = 0)
    )
  
  list(plot = p, n_pathways = n_pathways)
}

fig_bfr_alldb <- make_gsea_bubble_plot_alldb()
plot_height_alldb <- max(4, 0.22 * fig_bfr_alldb$n_pathways + 1.5)

if (fig_bfr_alldb$n_pathways > 0) {
  ggsave(file.path(config$fig_dir, "gsea_bubble_BFR_top30_alldb.pdf"),
         fig_bfr_alldb$plot, width = 12, height = plot_height_alldb, device = cairo_pdf)
  ggsave(file.path(config$fig_dir, "gsea_bubble_BFR_top30_alldb.png"),
         fig_bfr_alldb$plot, width = 12, height = plot_height_alldb, dpi = 300, bg = "white")
  cat("\nSaved all-database top-", TOP_N_ALL_DB, " figure to: ", config$fig_dir, "\n", sep = "")
} else {
  cat("\nAll-database top-", TOP_N_ALL_DB, " plot skipped -- 0 pathways passed adjp_thresh.\n", sep = "")
}

# --- 6: Per-Contrast Top-N Bubble Plots (pooled across all databases) --------
# Same redundancy-filtered candidate pool as section 5, but selection is
# driven by EACH CONTRAST'S OWN p.adjust rather than pooling significance
# across gsea_comps -- i.e. "top 30 pathways diagnostic for Training_BFR
# specifically" rather than "top 30 pathways significant anywhere." Loops
# over gsea_comps, which is now just Training_BFR/Training_HLRT, so this
# produces 2 plots instead of 5.
#
# The resulting plot still shows both comparison columns (not just the
# selecting contrast) so the selected pathways' behavior in the other
# contrast stays visible.
TOP_N_PER_CONTRAST <- 30

select_pathways_per_contrast <- function(contrast_name, top_n) {
  df <- gsea_supp_full %>%
    filter(comparison == contrast_name) %>%
    dplyr::group_by(database) %>%
    dplyr::group_modify(~ {
      db_name <- as.character(.y$database)
      rep_ids <- representative_lookup %>% filter(database == db_name) %>% pull(ID)
      if (length(rep_ids) == 0) return(.x)
      dplyr::filter(.x, ID %in% rep_ids)
    }) %>%
    dplyr::ungroup() %>%
    filter(p.adjust < adjp_thresh) %>%
    arrange(p.adjust)
  
  if (nrow(df) == 0) {
    warning("No pathways pass adjp_thresh = ", adjp_thresh, " in any database for ", contrast_name)
    return(df %>% dplyr::select(database, Description_clean))
  }
  
  df %>% slice_head(n = top_n) %>% dplyr::select(database, Description_clean)
}

build_bubble_selection_per_contrast <- function(contrast_name) {
  keep <- select_pathways_per_contrast(contrast_name, TOP_N_PER_CONTRAST)
  gsea_supp_full %>% semi_join(keep, by = c("database", "Description_clean"))
}

make_gsea_bubble_plot_percontrast <- function(contrast_name) {
  sel <- build_bubble_selection_per_contrast(contrast_name)
  
  present_dbs <- unname(db_slots)[unname(db_slots) %in% unique(as.character(sel$database))]
  
  sel <- sel %>% mutate(pathway_key = paste(database, ID, sep = "::"))
  
  ordered_levels <- unlist(lapply(present_dbs, function(db) {
    order_pathways(filter(sel, database == db), key_col = "pathway_key")
  }))
  
  label_lookup <- setNames(sel$Description_clean, sel$pathway_key)
  
  sel <- sel %>%
    mutate(
      database    = factor(as.character(database), levels = present_dbs),
      pathway_key = factor(pathway_key, levels = ordered_levels),
      significant = p.adjust < adjp_thresh
    )
  
  n_pathways <- length(ordered_levels)
  
  p <- ggplot(sel, aes(x = comparison, y = pathway_key)) +
    geom_point(aes(size = -log10(p.adjust), fill = NES,
                   color = significant, stroke = significant),
               shape = 21) +
    scale_fill_gradient2(
      low = NES_LOW, mid = colors$white, high = NES_HIGH, midpoint = 0,
      name = "NES"
    ) +
    scale_color_manual(
      values = c(`TRUE` = colors$black, `FALSE` = colors$gray_light),
      labels = c(`TRUE` = paste0("FDR < ", adjp_thresh),
                 `FALSE` = paste0("FDR \u2265 ", adjp_thresh)),
      name = "Significance"
    ) +
    scale_discrete_manual(
      aesthetics = "stroke",
      values = c(`TRUE` = 1, `FALSE` = 0.3),
      guide = "none"
    ) +
    scale_size_continuous(name = expression(-log[10]~"(FDR)"), range = c(1, 7)) +
    scale_x_discrete(labels = comp_labels) +
    scale_y_discrete(labels = label_lookup) +
    facet_grid(rows = vars(database), scales = "free_y", space = "free_y") +
    labs(title = paste0("GSEA \u2014 Top ", TOP_N_PER_CONTRAST,
                        " for ", comp_labels[[contrast_name]]),
         subtitle = "Deduplicated, pooled across 7 databases",
         x = NULL, y = NULL) +
    theme_cr(base_size = 9, title_hjust = 0) +
    theme(
      axis.text.y      = element_text(size = 7, face = "plain"),
      axis.text.x      = element_text(size = 8, face = "bold"),
      legend.position  = "right",
      panel.grid.major = element_line(color = colors$gray_light),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = colors$gray_light, color = NA),
      strip.text.y     = element_text(angle = 0)
    )
  
  list(plot = p, n_pathways = n_pathways)
}

# --- 6b: Build, save, and export companion data -------------------------------
percontrast_selection_data <- purrr::map_dfr(gsea_comps, function(contrast_name) {
  df <- build_bubble_selection_per_contrast(contrast_name)
  if (nrow(df) == 0) return(NULL)
  df %>%
    mutate(selected_for = comp_labels[[contrast_name]]) %>%
    dplyr::select(selected_for, database, ID,
                  Description = Description_clean,
                  comparison = comparison_label,
                  NES, pvalue, p.adjust, setSize, significant)
})

writexl::write_xlsx(
  list(BFR = percontrast_selection_data),
  path = file.path(config$data_dir, sprintf("gsea_bubble_percontrast_top%d_selection.xlsx",
                                            TOP_N_PER_CONTRAST))
)
cat("\nPer-contrast selection data saved to:",
    file.path(config$data_dir, sprintf("gsea_bubble_percontrast_top%d_selection.xlsx",
                                       TOP_N_PER_CONTRAST)), "\n")

walk(gsea_comps, function(contrast_name) {
  res <- make_gsea_bubble_plot_percontrast(contrast_name)
  
  if (res$n_pathways == 0) {
    cat("\n", contrast_name, " per-contrast top-", TOP_N_PER_CONTRAST,
        " plot skipped -- 0 pathways passed adjp_thresh.\n", sep = "")
    return(invisible(NULL))
  }
  
  plot_height <- max(4, 0.22 * res$n_pathways + 1.5)
  file_stub <- sprintf("gsea_bubble_BFR_%s_top%d_percontrast",
                       contrast_name, TOP_N_PER_CONTRAST)
  
  ggsave(file.path(config$fig_dir, paste0(file_stub, ".pdf")),
         res$plot, width = 7, height = plot_height, device = cairo_pdf)
  ggsave(file.path(config$fig_dir, paste0(file_stub, ".png")),
         res$plot, width = 7, height = plot_height, dpi = 300, bg = "white")
  
  cat(" Saved: ", file_stub, ".pdf/.png (", res$n_pathways, " pathways)\n", sep = "")
})
