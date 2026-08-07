# GSEA Differential Pathway Assessment: Training_BFR vs Training_HLRT
#
# Classifies every redundancy-filtered pathway (pooled across all 7 GSEA
# databases) into one of:
#   - "BFR only"    : significant in Training_BFR, not in Training_HLRT
#   - "HLRT only"   : significant in Training_HLRT, not in Training_BFR
#   - "Discordant"  : significant in BOTH, but opposite NES direction
#   - "Concordant"  : significant in BOTH, same NES direction
#   - "Neither"     : not significant in either
#
# "Different between the two contrasts" = BFR only + HLRT only + Discordant.
# Concordant pathways are explicitly EXCLUDED from that set -- they're
# significant in both arms and move the same direction, i.e. exactly the
# pathways that are NOT different between the two training modalities.
#
# This is the pathway-level version of the protein-level concordance
# scatter (logfc_scatter_BFR_v_HLRT.R) -- same underlying question, same
# BFR-only/HLRT-only/discordant framing, applied to GSEA results instead
# of individual protein logFCs.

# --- 0: Setup ------------------------
setwd(rprojroot::find_rstudio_root_file())
getwd()
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(dplyr, tidyr, tibble, stringr, purrr, readr, writexl, ggplot2)

# --- Configuration -----------------------------------------------------------
config <- list(
  gsea_data = "03_analysis/c_data/GSEA/consolidated/gsea_pathway_data.rds",
  fig_dir   = "04_figures/fig_gsea_differential",
  data_dir  = "04_figures/data/fig_gsea_differential"
)
dir.create(config$fig_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(config$data_dir, showWarnings = FALSE, recursive = TRUE)

source("04_figures/a_scripts/00_theme.R")

adjp_thresh <- 0.05
gsea_group  <- "BFR"

# Same muted qualitative palette used in fig_gsea_updown_BFR.R, kept
# consistent across figures that color by database.
db_colors <- c(
  "GO:BP"    = "#1b9e77",
  "GO:CC"    = "#e6ab02",
  "GO:MF"    = "#7570b3",
  "GO Slim"  = "#66a61e",
  "KEGG"     = "#d95f02",
  "Reactome" = "#a6761d",
  "Hallmark" = "#e7298a"
)

status_colors <- c(
  "BFR only"   = "#d95f02",
  "HLRT only"  = "#7570b3",
  "Discordant" = "#e7298a",
  "Concordant" = "#66a61e",
  "Neither"    = "grey80"
)

# --- 1: Load consolidated GSEA data -------------------------------------
if (!file.exists(config$gsea_data)) {
  stop("Consolidated GSEA data not found at ", config$gsea_data,
       " -- run the GSEA consolidation script first.")
}
gsea_pathway_data <- readRDS(config$gsea_data)

gsea_supp_full        <- gsea_pathway_data$gsea_supp_full        %>% filter(group == gsea_group)
representative_lookup <- gsea_pathway_data$representative_lookup %>% filter(group == gsea_group)

# --- 2: Build BFR vs HLRT comparison table -----------------------------------
# NOTE: GO:CC and GO:MF were run with pvalueCutoff = 0.05 in the GSEA
# script (unlike the other 5 databases, which use pvalueCutoff = 1), so a
# GO:CC/GO:MF pathway that's non-significant in one contrast may be
# genuinely absent from that contrast's result table rather than just
# "tested and not significant." Practically this means "X only" for
# GO:CC/GO:MF pathways can mean either "significant in X and not in Y" OR
# "significant in X and not even returned for Y" -- both get treated the
# same way here (as non-significant via the NA-safe check below), which is
# the right call for this classification, just worth knowing if you dig
# into a specific GO:CC/GO:MF pathway's NA cells in the exported table.
bfr_df <- gsea_supp_full %>%
  filter(comparison == "Training_BFR") %>%
  transmute(pathway_key = paste(database, ID, sep = "::"),
            database, ID, Description_clean,
            NES_BFR = NES, padj_BFR = p.adjust)

hlrt_df <- gsea_supp_full %>%
  filter(comparison == "Training_HLRT") %>%
  transmute(pathway_key = paste(database, ID, sep = "::"),
            database, ID, Description_clean,
            NES_HLRT = NES, padj_HLRT = p.adjust)

df <- full_join(bfr_df, hlrt_df, by = "pathway_key", suffix = c("_bfr", "_hlrt")) %>%
  mutate(
    database          = coalesce(database_bfr, database_hlrt),
    ID                = coalesce(ID_bfr, ID_hlrt),
    Description_clean = coalesce(Description_clean_bfr, Description_clean_hlrt)
  ) %>%
  dplyr::select(pathway_key, database, ID, Description_clean,
                NES_BFR, padj_BFR, NES_HLRT, padj_HLRT)

# Restrict to redundancy-filter survivors (same representative_lookup used
# throughout the rest of the pipeline), so near-duplicate GO terms don't
# inflate either side of this comparison.
rep_keys <- representative_lookup %>%
  transmute(pathway_key = paste(database, ID, sep = "::")) %>%
  pull(pathway_key)

df <- df %>% filter(pathway_key %in% rep_keys)

# --- 3: Classify --------------------------------------------------------
df <- df %>%
  mutate(
    sig_bfr    = !is.na(padj_BFR)  & padj_BFR  < adjp_thresh,
    sig_hlrt   = !is.na(padj_HLRT) & padj_HLRT < adjp_thresh,
    delta_NES  = NES_BFR - NES_HLRT,
    status = case_when(
      sig_bfr  & sig_hlrt & sign(NES_BFR) == sign(NES_HLRT) ~ "Concordant",
      sig_bfr  & sig_hlrt & sign(NES_BFR) != sign(NES_HLRT) ~ "Discordant",
      sig_bfr  & !sig_hlrt                                  ~ "BFR only",
      !sig_bfr & sig_hlrt                                   ~ "HLRT only",
      TRUE                                                  ~ "Neither"
    )
  )

# The actual deliverable: pathways that differ between the two training
# modalities. Ranked by |delta_NES| (effect-size divergence) rather than
# significance, consistent with the NES-based ranking used elsewhere.
differential_df <- df %>%
  filter(status %in% c("BFR only", "HLRT only", "Discordant")) %>%
  arrange(desc(abs(delta_NES)))

cat(sprintf("\n=== Classification summary (group = '%s') ===\n", gsea_group))
print(table(df$status))
cat(sprintf("\n%d pathways classified as 'different' (BFR only + HLRT only + Discordant)\n",
            nrow(differential_df)))

# --- 4: Export ----------------------------------------------------------
export_cols <- c("status", "database", "ID", "Description_clean",
                 "NES_BFR", "padj_BFR", "NES_HLRT", "padj_HLRT", "delta_NES")

write_xlsx(
  list(
    Differential = differential_df %>% dplyr::select(all_of(export_cols)),
    All_Pathways = df %>% dplyr::select(all_of(export_cols)) %>% arrange(status, desc(abs(delta_NES))),
    Summary      = df %>% count(status, database) %>% pivot_wider(names_from = database, values_from = n, values_fill = 0)
  ),
  path = file.path(config$data_dir, "gsea_differential_BFR_v_HLRT.xlsx")
)
cat("\nSaved:", file.path(config$data_dir, "gsea_differential_BFR_v_HLRT.xlsx"), "\n")

cat("\nTop 15 most divergent pathways:\n")
print(differential_df %>% dplyr::select(status, database, Description_clean, NES_BFR, NES_HLRT, delta_NES) %>%
        head(15))

# --- 5: Summary count plot ------------------------------------------------
count_summary <- df %>%
  filter(status != "Neither") %>%
  count(status, database) %>%
  mutate(
    status   = factor(status, levels = c("BFR only", "HLRT only", "Discordant", "Concordant")),
    database = factor(database, levels = names(db_colors))
  )

p_counts <- ggplot(count_summary, aes(x = status, y = n, fill = database)) +
  geom_col(position = "stack", color = "white", linewidth = 0.3) +
  geom_text(
    data = count_summary %>% group_by(status) %>% summarise(total = sum(n), .groups = "drop"),
    aes(x = status, y = total, label = total), inherit.aes = FALSE,
    vjust = -0.4, fontface = "bold", size = 3.5
  ) +
  scale_fill_manual(values = db_colors, name = "Database") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    x = NULL, y = "Number of pathways",
    title = "GSEA Pathway Classification: Training_BFR vs Training_HLRT",
    subtitle = "\"Different\" = BFR only + HLRT only + Discordant (excludes Concordant)"
  ) +
  theme_cr(base_size = 10) +
  theme(legend.position = "right")

ggsave(file.path(config$fig_dir, "gsea_differential_counts.pdf"),
       p_counts, width = 8, height = 5.5, device = cairo_pdf)
ggsave(file.path(config$fig_dir, "gsea_differential_counts.png"),
       p_counts, width = 8, height = 5.5, dpi = 300, bg = "white")

# --- 6: NES-vs-NES scatter, highlighting differential pathways --------------
# Pathway-level companion to logfc_scatter_BFR_v_HLRT.R -- same axes/
# framing (x = BFR training response, y = HLRT training response), one
# point per pathway instead of per protein. Only plotted for pathways with
# an NES on both sides (a "BFR only"/"HLRT only" pathway can still have a
# non-significant NES on the missing side, which is what lets it appear
# here at all; a pathway entirely absent from one side's result table --
# possible for GO:CC/GO:MF per the note in section 2 -- can't be placed
# and is dropped from this plot only, not from the exported tables above).
scatter_df <- df %>% filter(!is.na(NES_BFR), !is.na(NES_HLRT))

lim <- max(abs(c(scatter_df$NES_BFR, scatter_df$NES_HLRT)), na.rm = TRUE) * 1.05

p_scatter <- ggplot(scatter_df, aes(x = NES_BFR, y = NES_HLRT, color = status)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4, color = colors$gray_dark) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, color = colors$gray_dark) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", linewidth = 0.5, color = colors$gray_dark) +
  geom_point(data = filter(scatter_df, status %in% c("Neither", "Concordant")), size = 1.3, alpha = 0.4) +
  geom_point(data = filter(scatter_df, status %in% c("BFR only", "HLRT only", "Discordant")),
             size = 2.2, alpha = 0.85) +
  scale_color_manual(values = status_colors, name = "Status") +
  coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  labs(
    x = "NES (Training_BFR)", y = "NES (Training_HLRT)",
    title = "Pathway-Level Training Response: BFR vs HLRT"
  ) +
  theme_cr(base_size = 11) +
  theme(legend.position = "right")

ggsave(file.path(config$fig_dir, "gsea_differential_scatter.pdf"),
       p_scatter, width = 8, height = 7, device = cairo_pdf)
ggsave(file.path(config$fig_dir, "gsea_differential_scatter.png"),
       p_scatter, width = 8, height = 7, dpi = 300, bg = "white")

cat("\nSaved figures to:", config$fig_dir, "\n")