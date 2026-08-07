# =============================================================================
# fig_05_panel_A_protein_BFR.R
# Panel A: BFR vs HLRT training-response concordance scatter + CONCORDANT/
# DIVERGENT ORA flanking bars.
# Depends on fig_05_style.R being sourced first (for fig_theme_scaled,
# scale_text, BASE_*, PANEL_MD, pal_shared, build_quad_label_df,
# geom_quadrant_label, strip_for_composite, colors -- unchanged from the
# original, not redefined here).
#
# REFRAMED FROM THE CR MUSCLE VERSION: the original plotted "reversal" --
# whether training pushed cancer patients' proteins back toward what an
# untrained control group looked like at baseline (x = PRE vs CTL,
# y = POST vs PRE within SURV). BFR has no untrained control; both arms
# train. Per your choice, the two axes here are instead the two arms'
# OWN training responses: x = Training_BFR (BFR_POST vs BFR_PRE),
# y = Training_HLRT (HLRT_POST vs HLRT_PRE). The question this answers is
# "do the two training modalities move a given protein the same way, or
# different ways" -- concordant vs. divergent, rather than reversed vs.
# exacerbated relative to a baseline gap.
#
# The quadrant geometry is identical to the original (a diagonal split:
# same-sign quadrants vs. opposite-sign quadrants), just relabeled:
#   TR (both up) + BL (both down)   -> Concordant   (was: Exacerbated, red)
#   TL (BFR down/HLRT up) + BR (BFR up/HLRT down) -> Divergent (was: Reversed, blue)
# =============================================================================

pacman::p_load(dplyr, tidyr, tibble, ggplot2, ggrepel, patchwork)

# Multiplier applied on top of BASE_PATHWAY, just for the ORA flanking-bar
# labels in Panel A -- bump this if you want them bigger/smaller without
# touching BASE_PATHWAY globally.
PATHWAY_LABEL_SCALE <- 1.45

ora_label_map <- c(
  "L-amino acid catabolic process"                                  = "L-AA Catabolism",
  "branched-chain amino acid catabolic process"                     = "BCAA Catabolism",
  "quinone metabolic process"                                       = "Quinone Metabolism",
  "pyruvate family amino acid catabolic process"                    = "Pyruvate AA Catabolism",
  "spermatid differentiation"                                       = "Spermatid Diff.",
  "calcium ion transmembrane transport"                             = "Calcium Transport",
  "actin cytoskeleton organization"                                 = "Actin Cytoskeleton Org.",
  "cell adhesion"                                                   = "Cell Adhesion",
  "vesicle-mediated transport"                                      = "Vesicle Transport",
  "cell motility"                                                   = "Cell Motility",
  "cell-substrate adhesion"                                         = "Cell-Substrate\nAdhesion",
  "regulation of actin cytoskeleton organization"                   = "Regulation of\nActin Cytoskeleton",
  "response to toxic substance"                                     = "Toxic Substance Response",
  "glutathione metabolic process"                                   = "Glutathione Metabolism",
  "Ras protein signal transduction"                                 = "Ras Signaling",
  "aldehyde metabolic process"                                      = "Aldehyde Metabolism",
  "regulation of proteolysis involved in protein catabolic process" = "Regulation of\nProteolysis",
  "proteolysis"                                                     = "Proteolysis",
  "proton transmembrane transport"                                  = "Proton Transport",
  "energy derivation by oxidation of organic compounds"             = "Oxidative Energy Deriv.",
  "mitochondrial electron transport, NADH to ubiquinone"            = "Mitochondrial\nElectron Trans.",
  "nucleoside phosphate biosynthetic process"                       = "Nucleoside-P\nBiosynthesis",
  "translation"                                                     = "Translation",
  "NADH dehydrogenase complex assembly"                             = "NADH\nDehydrogenase Assembly"
)

apply_labels <- function(df, map) {
  df |> mutate(label = if_else(Description %in% names(map), map[Description], Description))
}

# Kept independent of the pooled-Jaccard filter used for GSEA pathways --
# CONCORDANT/DIVERGENT ORA is protein-level ORA, not one of the pooled
# GSEA contrasts.
jaccard_filter_ora <- function(res, threshold = 0.5,
                               padj_col = "p.adjust", id_col = "ID", gene_col = "geneID") {
  jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))
  df <- as.data.frame(res)
  df <- df[grepl("^GO:", df[[id_col]]), ]
  df <- df[order(df[[padj_col]]), ]
  sets <- strsplit(df[[gene_col]], "/")
  n    <- nrow(df); keep <- rep(TRUE, n)
  for (i in seq_len(n - 1)) {
    if (!keep[i]) next
    for (j in seq(i + 1, n)) {
      if (!keep[j]) next
      if (jaccard(sets[[i]], sets[[j]]) >= threshold) keep[j] <- FALSE
    }
  }
  df[keep, ]
}

make_ora_right <- function(df, fill_low, fill_high, panel_width_mm, show_x_title = FALSE,
                           outside_labels = character(0)) {
  df <- df |>
    mutate(neg_log10_p = -log10(p.adjust),
           y_order = reorder(label, neg_log10_p),
           sig_label = case_when(p.adjust < 0.001 ~ "***", p.adjust < 0.01 ~ "**",
                                 p.adjust < 0.05 ~ "*", TRUE ~ "ns"),
           is_outside = label %in% outside_labels)

  df_in  <- filter(df, !is_outside)
  df_out <- filter(df, is_outside)

  ggplot(df, aes(x = neg_log10_p, y = y_order)) +
    geom_col(aes(fill = neg_log10_p), width = 0.7) +
    # in-bar labels: white text, anchored near the base of the bar
    geom_text(data = df_in, aes(label = label), x = 0.05, hjust = 0,
              size = scale_text(BASE_PATHWAY * PATHWAY_LABEL_SCALE, panel_width_mm),
              colour = "white", fontface = "bold") +
    # flagged labels: dark text, anchored just past the bar's tip
    geom_text(data = df_out, aes(x = neg_log10_p, label = label), hjust = -0.05,
              size = scale_text(BASE_PATHWAY * PATHWAY_LABEL_SCALE, panel_width_mm),
              colour = "grey15", fontface = "bold") +
    geom_text(aes(x = neg_log10_p, label = sig_label), hjust = -0.2, vjust = 0.75,
              size = scale_text(BASE_COUNT, panel_width_mm), colour = "black", fontface = "bold") +
    scale_fill_gradient(low = fill_low, high = fill_high, guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, if (nrow(df_out) > 0) 0.32 else 0.18)), position = "bottom") +
    scale_y_discrete(position = "left") +
    coord_cartesian(clip = "off") +
    labs(x = if (show_x_title) expression(-log[10](p[adj])) else NULL, y = NULL) +
    fig_theme_scaled(panel_width_mm) +
    theme(
      axis.text.y = element_blank(), axis.ticks.y = element_blank(),
      panel.grid = element_blank(), legend.position = "none",
      plot.margin = margin(2, 10, 2, 0)
    )
}

make_ora_left <- function(df, fill_low, fill_high, panel_width_mm, show_x_title = FALSE,
                          outside_labels = character(0)) {
  df <- df |>
    mutate(neg_log10_p = -log10(p.adjust), neg_log10_p_neg = -neg_log10_p,
           y_order = reorder(label, neg_log10_p),
           sig_label = case_when(p.adjust < 0.001 ~ "***", p.adjust < 0.01 ~ "**",
                                 p.adjust < 0.05 ~ "*", TRUE ~ "ns"),
           is_outside = label %in% outside_labels)

  df_in  <- filter(df, !is_outside)
  df_out <- filter(df, is_outside)

  ggplot(df, aes(x = neg_log10_p_neg, y = y_order)) +
    geom_col(aes(fill = neg_log10_p), width = 0.7) +
    # in-bar labels: white text, anchored near the base of the bar
    geom_text(data = df_in, aes(x = neg_log10_p_neg * 0.97, label = label), hjust = 0,
              size = scale_text(BASE_PATHWAY * PATHWAY_LABEL_SCALE, panel_width_mm),
              colour = "white", fontface = "bold") +
    # flagged labels: dark text, anchored just past the bar's tip (further left)
    geom_text(data = df_out, aes(x = neg_log10_p_neg, label = label), hjust = 1,
              size = scale_text(BASE_PATHWAY * PATHWAY_LABEL_SCALE, panel_width_mm),
              colour = "grey15", fontface = "bold") +
    geom_text(aes(x = neg_log10_p_neg, label = sig_label), hjust = 1.4, vjust = 0.75,
              size = scale_text(BASE_COUNT, panel_width_mm), colour = "black", fontface = "bold") +
    scale_fill_gradient(low = fill_low, high = fill_high, guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(if (nrow(df_out) > 0) 0.32 else 0.18, 0)),
                       labels = function(x) abs(x), position = "bottom") +
    scale_y_discrete(position = "right") +
    coord_cartesian(clip = "off") +
    labs(x = if (show_x_title) expression(-log[10](p[adj])) else NULL, y = NULL) +
    fig_theme_scaled(panel_width_mm) +
    theme(
      axis.text.y = element_blank(), axis.ticks.y = element_blank(),
      panel.grid = element_blank(), legend.position = "none",
      plot.margin = margin(2, 0, 2, 10)
    )
}

#' Build Panel A: BFR vs HLRT training-response concordance scatter + ORA flanking bars
#'
#' @param limma_training_bfr  Raw Training_BFR sheet (BFR_POST vs BFR_PRE)
#' @param limma_training_hlrt Raw Training_HLRT sheet (HLRT_POST vs HLRT_PRE)
#' @param gsea_data Named list from load_gsea_rds(), must contain
#'   $CONCORDANT / $DIVERGENT (see NOTE below -- these don't exist yet and
#'   need a corresponding ORA run added to the GSEA script).
#' @param panel_width_mm Total rendered width of Panel A in the final composite (mm)
#' @param pal Named color vector for the 4 sig categories. MUST have names
#'   "Both", "BFR only", "HLRT only", "Neither" -- if you're reusing
#'   pal_shared from the SURV version, its names are currently
#'   "PRE v CTL only"/"POST v PRE only" and need updating to match.
build_panel_A <- function(limma_training_bfr, limma_training_hlrt, gsea_data,
                          panel_width_mm = PANEL_MD, pal = pal_shared,
                          ora_outside = list()) {

  df_bfr <- limma_training_bfr |>
    dplyr::select(uniprot_id, gene, logFC_BFR = logFC, Pi.Val_BFR = Pi.Val, adj.P.Val_BFR = adj.P.Val)
  df_hlrt <- limma_training_hlrt |>
    dplyr::select(uniprot_id, gene, logFC_HLRT = logFC, Pi.Val_HLRT = Pi.Val, adj.P.Val_HLRT = adj.P.Val)

  df_train <- inner_join(df_bfr, df_hlrt, by = c("uniprot_id", "gene")) |>
    mutate(sig = case_when(
      Pi.Val_BFR < 0.05 & Pi.Val_HLRT < 0.05 ~ "Both",
      Pi.Val_BFR < 0.05                      ~ "BFR only",
      Pi.Val_HLRT < 0.05                     ~ "HLRT only",
      TRUE                                   ~ "Neither"
    ))

  lim <- max(abs(c(df_train$logFC_BFR, df_train$logFC_HLRT)), na.rm = TRUE)
  r_spear <- suppressWarnings(cor.test(df_train$logFC_BFR, df_train$logFC_HLRT, method = "spearman", exact = FALSE))

  protein_label_df <- df_train |>
    filter(sig != "Neither", abs(logFC_BFR) > 1 | abs(logFC_HLRT) > 1) |>
    mutate(
      Pi.Val_HLRT = pmax(Pi.Val_HLRT, 1e-300), Pi.Val_BFR = pmax(Pi.Val_BFR, 1e-300),
      direction  = if_else(abs(logFC_BFR) >= abs(logFC_HLRT),
                           if_else(logFC_BFR > 0, "Up", "Down"),
                           if_else(logFC_HLRT > 0, "Up", "Down")),
      score = pmax(-log10(Pi.Val_HLRT), -log10(Pi.Val_BFR))
    ) |>
    group_by(sig, direction) |> slice_max(score, n = 5, with_ties = FALSE) |> ungroup()

  # Diagonal split, same geometry as the original: same-sign quadrants
  # (TR/BL, "Concordant") vs. opposite-sign quadrants (TL/BR, "Divergent").
  quad_counts <- df_train |>
    mutate(q = case_when(
      logFC_BFR < 0 & logFC_HLRT > 0 ~ "TL", logFC_BFR > 0 & logFC_HLRT > 0 ~ "TR",
      logFC_BFR > 0 & logFC_HLRT < 0 ~ "BR", logFC_BFR < 0 & logFC_HLRT < 0 ~ "BL"
    )) |>
    filter(!is.na(q)) |> group_by(q) |>
    summarise(n_total = n(), n_sig = sum(sig != "Neither", na.rm = TRUE), .groups = "drop")

  quad_ann <- build_quad_label_df(
    counts = quad_counts, quad_ids = c("TL", "TR", "BR", "BL"),
    titles = c("Divergent (HLRT \u2191 / BFR \u2193)", "Concordant \u2191",
              "Divergent (BFR \u2191 / HLRT \u2193)", "Concordant \u2193"),
    label_fn = function(title, row) sprintf("%s\n%d / %d", title, row$n_sig, row$n_total),
    blue_quads = c("TL", "BR"), red_quads = c("TR", "BL")
  )

  th <- fig_theme_scaled(panel_width_mm)

  p_scatter <- ggplot(df_train, aes(x = logFC_BFR, y = logFC_HLRT, color = sig, fill = sig)) +
    annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf, fill = colors$blue, alpha = 0.2) +
    annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf, fill = colors$red, alpha = 0.2) +
    annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0, fill = colors$red, alpha = 0.2) +
    annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0, fill = colors$blue, alpha = 0.2) +
    geom_hline(yintercept = 0, linewidth = 0.6, color = "grey30", linetype = "dashed") +
    geom_vline(xintercept = 0, linewidth = 0.6, color = "grey30", linetype = "dashed") +
    geom_point(data = filter(df_train, sig == "Neither"), size = 1.2, alpha = 0.5, shape = 16) +
    geom_point(data = filter(df_train, sig != "Neither"), size = 2.2, alpha = 0.85, shape = 16) +
    scale_color_manual(values = pal, name = NULL) +
    scale_fill_manual(values = pal, guide = "none") +
    geom_quadrant_label(quad_ann, size = scale_text(BASE_QUADRANT, panel_width_mm)) +
    geom_label_repel(
      data = protein_label_df, aes(label = gene, fill = sig),
      color = colors$white, size = scale_text(BASE_GENE, panel_width_mm),
      max.overlaps = 20, box.padding = 0.3, label.padding = unit(0.15, "lines"),
      label.size = NA, segment.color = "grey60", segment.size = 0.3, show.legend = FALSE,
      xlim = c(-lim * 0.85, lim * 0.85), ylim = c(-lim * 0.85, lim * 0.85)
    ) +
    annotate("text", x = -lim * 0.80, y = 0.01, label = expression(log[2]*FC~"(BFR: POST vs PRE)"),
             hjust = 0.5, vjust = -0.5, size = scale_text(BASE_STAT, panel_width_mm), fontface = "bold", color = "grey20") +
    annotate("text", x = -0.01, y = -lim * 0.80, label = expression(log[2]*FC~"(HLRT: POST vs PRE)"),
             hjust = 0.5, vjust = 1.3, angle = 90, size = scale_text(BASE_STAT, panel_width_mm), fontface = "bold", color = "grey20") +
    coord_cartesian(xlim = c(-lim, lim), ylim = c(-lim, lim), clip = "off") +
    th +
    theme(legend.position = "bottom", axis.title = element_blank(),
          axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank())

  axis_breaks <- pretty(c(-lim, lim))
  axis_breaks <- axis_breaks[axis_breaks != 0 & abs(axis_breaks) < lim]
  p_scatter <- p_scatter +
    annotate("segment", x = axis_breaks, xend = axis_breaks, y = -lim * 0.015, yend = lim * 0.015,
             color = "grey30", linewidth = 0.4) +
    annotate("text", x = axis_breaks, y = -lim * 0.03, label = axis_breaks,
             vjust = 1, size = scale_text(BASE_STAT * 0.85, panel_width_mm), color = "grey20") +
    annotate("segment", y = axis_breaks, yend = axis_breaks, x = -lim * 0.015, xend = lim * 0.015,
             color = "grey30", linewidth = 0.4) +
    annotate("text", y = axis_breaks, x = -lim * 0.03, label = axis_breaks,
             hjust = 1, size = scale_text(BASE_STAT * 0.85, panel_width_mm), color = "grey20")

  # --- ORA flanking bars -------------------------------------------------
  # NOTE: gsea_data$CONCORDANT / $DIVERGENT do NOT exist in your GSEA
  # script yet -- this is new relative to the SURV version, which had
  # analogous (also not fully wired up) REVERSAL/EXACERBATION calls. You
  # need to add, in the GSEA script, something like:
  #
  #   concordant_df <- df_train %>% filter(sig != "Neither",
  #                        (logFC_BFR > 0) == (logFC_HLRT > 0)) %>%
  #                      mutate(logFC = logFC_HLRT)   # up/down split key
  #   divergent_df  <- df_train %>% filter(sig != "Neither",
  #                        (logFC_BFR > 0) != (logFC_HLRT > 0)) %>%
  #                      mutate(logFC = logFC_HLRT)
  #   run_ora_pipeline(concordant_df, "CONCORDANT", bg_df = bg_full)
  #   run_ora_pipeline(divergent_df,  "DIVERGENT",  bg_df = bg_full)
  #
  # (using logFC_HLRT as the up/down split key, matching the y-axis, the
  # same convention the original REVERSAL/EXACERBATION split used for RET)
  ora_conc_up   <- jaccard_filter_ora(gsea_data$CONCORDANT$ora_bp_up@result)   |> arrange(p.adjust) |> slice_head(n = 6) |> apply_labels(ora_label_map)
  ora_conc_down <- jaccard_filter_ora(gsea_data$CONCORDANT$ora_bp_down@result) |> arrange(p.adjust) |> slice_head(n = 6) |> apply_labels(ora_label_map)
  ora_div_up    <- jaccard_filter_ora(gsea_data$DIVERGENT$ora_bp_up@result)    |> arrange(p.adjust) |> slice_head(n = 6) |> apply_labels(ora_label_map)
  ora_div_down  <- jaccard_filter_ora(gsea_data$DIVERGENT$ora_bp_down@result)  |> arrange(p.adjust) |> slice_head(n = 6) |> apply_labels(ora_label_map)

  n_enrich <- sum(gsea_data$CONCORDANT$ora_bp_up@result$p.adjust < 0.05, na.rm = TRUE) +
              sum(gsea_data$CONCORDANT$ora_bp_down@result$p.adjust < 0.05, na.rm = TRUE) +
              sum(gsea_data$DIVERGENT$ora_bp_up@result$p.adjust < 0.05, na.rm = TRUE) +
              sum(gsea_data$DIVERGENT$ora_bp_down@result$p.adjust < 0.05, na.rm = TRUE)

  # ORA flank panels sit at ~1.1/4.2 = ~0.26x the total composite's rendered
  # width (see plot_layout(widths) below: 1.1 : 2 : 1.1), so their own text
  # scale target is narrower than the scatter's.
  flank_width_mm <- panel_width_mm * (1.1 / 4.2)

  get_outside <- function(nm) if (is.null(ora_outside[[nm]])) character(0) else ora_outside[[nm]]

  # Spatial mapping: TL sits upper-left of the scatter -> left column, top
  # slot; BL -> left column, bottom slot; TR -> right column, top slot;
  # BR -> right column, bottom slot. Fill gradient follows the same
  # blue=Divergent / red=Concordant split used for the quadrant rects above.
  left_col  <- make_ora_left(ora_div_up,    colors$blue_light, colors$blue_dark, flank_width_mm, FALSE,
                             outside_labels = get_outside("div_up")) /
    make_ora_left(ora_conc_down, colors$red_light,  colors$red_dark,  flank_width_mm, TRUE,
                  outside_labels = get_outside("conc_down"))
  right_col <- make_ora_right(ora_conc_up,   colors$red_light,  colors$red_dark,  flank_width_mm, FALSE,
                              outside_labels = get_outside("conc_up")) /
    make_ora_right(ora_div_down, colors$blue_light, colors$blue_dark, flank_width_mm, TRUE,
                   outside_labels = get_outside("div_down"))

  fig_compound <- (left_col | p_scatter | right_col) +
    plot_layout(widths = c(1.1, 2, 1.1)) &
    theme(plot.margin = margin(2, 0, 2, 0))

  ttl <- "BFR vs HLRT Training Concordance"
  sub <- sprintf("N = %d | %d sig (\u03a0) | %d ORA terms enriched (FDR) | \u03c1 = %.2f",
                 nrow(df_train), sum(df_train$sig != "Neither"), n_enrich, r_spear$estimate)

  list(
    plot           = fig_compound,                       # for standalone save (own legend)
    plot_composite = strip_for_composite(fig_compound),   # no legend/title, for the 380mm grid
    title = ttl, subtitle = sub,
    df_train = df_train
  )
}
