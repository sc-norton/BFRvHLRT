# GSEA / ORA Script (BFR)

# --- 0: Setup --------------------------
# Set working directory
setwd(rprojroot::find_rstudio_root_file())
getwd()

# Load packages
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(proteoDA, readxl, readr, dplyr, tidyr, stringr, purrr, limma,
               writexl, clusterProfiler, org.Hs.eg.db, enrichplot, enrichR,
               msigdbr, ReactomePA, biomaRt, GSEABase, GO.db, AnnotationDbi)

# Configuration
config <- list(
  # Directories
  report_dir = "03_analysis/b_reports/GSEA",
  data_dir = "03_analysis/c_data/GSEA",
  # Data files
  BFR_limma = "03_analysis/c_data/limma/limma_results.xlsx"
)

dir.create(config$report_dir, showWarnings = F, recursive = T)
dir.create(config$data_dir, showWarnings = F, recursive = T)

# Set organism for GSEA
organism   <- org.Hs.eg.db

# 1. Connect to the human Ensembl database
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

# 2. Extract the mapping (e.g., matching Ensembl Gene IDs to GO Slim annotations)
human_goslim <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "goslim_goa_accession", "goslim_goa_description"),
  mart = mart
)

slim_bp <- human_goslim %>%
  dplyr::select(goslim_goa_accession, goslim_goa_description) %>%
  dplyr::distinct()

slim_bp_ids <- slim_bp %>%
  filter(goslim_goa_accession != "") %>%
  pull(goslim_goa_accession)

# TERM2GENE matching for GO:SLIM GSEA run
slim_term2gene <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys      = slim_bp_ids,
  columns   = "ENSEMBL",
  keytype   = "GOALL"
) %>%
  filter(ONTOLOGYALL == "BP") %>%
  dplyr::select(GOALL, ENSEMBL) %>%
  dplyr::rename(term = GOALL, gene = ENSEMBL) %>%
  distinct()

slim_term2name <- human_goslim %>%
  filter(goslim_goa_accession != "") %>%
  dplyr::select(goslim_goa_accession, goslim_goa_description) %>%
  dplyr::rename(term = goslim_goa_accession, name = goslim_goa_description) %>%
  distinct()

# --- 1:  Load limma Data for pathway analysis ----
limma_results <- list(
  Baseline          = read_excel(file.path(config$BFR_limma),
                                 sheet = 1),
  Training_BFR      = read_excel(file.path(config$BFR_limma),
                                 sheet = 2),
  Training_HLRT     = read_excel(file.path(config$BFR_limma),
                                 sheet = 3),
  Post_training     = read_excel(file.path(config$BFR_limma),
                                 sheet = 4),
  Interaction  = read_excel(file.path(config$BFR_limma),
                            sheet = 5)
)

# --- 6:  Pathway / Enrichment Analysis ----
# --- GSEA Pipeline Function ---

run_gsea_pipeline <- function(df, label) {
  
  cat(sprintf("\n>> Running GSEA pipeline for: %s\n", label))
  
  set.seed(26)
  
  # --- Rank metric ---
  df$rank                   <- df$t # limma t-statistic
  original_gene_list        <- df$rank
  names(original_gene_list) <- df$Ensembl
  gene_list                 <- na.omit(original_gene_list)
  gene_list                 <- sort(gene_list, decreasing = TRUE)
  
  # --- Hallmark GSEA ---
  hallmark <- msigdbr(collection = "H") %>% dplyr::select(gs_name, ensembl_gene)
  h2       <- msigdbr(collection = "H") %>% dplyr::select(gs_name, gs_description)
  
  gsea_hall <- GSEA(
    geneList  = gene_list,
    TERM2GENE = hallmark,
    TERM2NAME = h2,
    eps = 0,
    minGSSize = 15,
    maxGSSize = 500,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    seed = 26,
    verbose = FALSE
  )
  
  # --- Reactome GSEA ---
  react <- msigdbr(collection = "C2", subcollection = "CP:REACTOME") %>% dplyr::select(gs_name, ensembl_gene)
  react2 <- msigdbr(collection = "C2", subcollection = "CP:REACTOME") %>% dplyr::select(gs_name, gs_description)
  
  gsea_reactome <- GSEA(
    geneList  = gene_list,
    TERM2GENE = react,
    TERM2NAME = react2,
    eps = 0,
    minGSSize = 15,
    maxGSSize = 500,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    seed = 26,
    verbose = FALSE
  )
  
  # --- GO GSEA (BP) ---
  gsea_bp <- gseGO(
    geneList      = gene_list,
    ont           = "BP",
    keyType       = "ENSEMBL",
    seed          = 26,
    eps           = 0,
    minGSSize     = 15,
    maxGSSize     = 500,
    pvalueCutoff  = 1,
    verbose       = FALSE,
    OrgDb         = organism,
    pAdjustMethod = "BH"
  )
  
  gsea_bp_clean        <- gsea_bp
  gsea_bp_clean@result <- gsea_bp@result %>% filter(!is.na(NES), !is.na(pvalue))
  # Simplify on significant terms only
  gsea_bp_sig          <- gsea_bp_clean
  gsea_bp_sig@result   <- gsea_bp_clean@result %>% filter(p.adjust < 0.05)
  gsea_bp_simple       <- simplify(gsea_bp_clean, cutoff = 0.7, by = "p.adjust", select_fun = min)
  
  # --- GO GSEA (CC) ---
  # NOTE: pvalueCutoff changed from 0.05 -> 1 to match every other database
  # (BP, MF, Slim, KEGG, Reactome, Hallmark). With the cutoff at 0.05, this
  # result object only ever contained already-significant terms, which
  # breaks the "full tested-pathway table across every database" the
  # consolidation script's n_tested/n_significant counts assume -- CC (and
  # MF below) were silently pre-truncated relative to every other database.
  gsea_cc <- gseGO(
    geneList      = gene_list,
    ont           = "CC",
    keyType       = "ENSEMBL",
    seed          = 26,
    eps           = 0,
    minGSSize     = 15,
    maxGSSize     = 500,
    pvalueCutoff  = 1,
    verbose       = FALSE,
    OrgDb         = organism,
    pAdjustMethod = "BH"
  )
  gsea_cc_clean        <- gsea_cc
  gsea_cc_clean@result <- gsea_cc@result %>% filter(!is.na(NES), !is.na(pvalue))
  # Simplify on significant terms only
  gsea_cc_sig          <- gsea_cc_clean
  gsea_cc_sig@result   <- gsea_cc_clean@result %>% filter(p.adjust < 0.05)
  gsea_cc_simple       <- simplify(gsea_cc_clean, cutoff = 0.7, by = "p.adjust", select_fun = min)
  
  # --- GO GSEA (MF) ---
  # NOTE: same pvalueCutoff fix as CC above (0.05 -> 1).
  gsea_mf <- gseGO(
    geneList      = gene_list,
    ont           = "MF",
    keyType       = "ENSEMBL",
    seed          = 26,
    eps           = 0,
    minGSSize     = 15,
    maxGSSize     = 500,
    pvalueCutoff  = 1,
    verbose       = FALSE,
    OrgDb         = organism,
    pAdjustMethod = "BH"
  )
  gsea_mf_clean        <- gsea_mf
  gsea_mf_clean@result <- gsea_mf@result %>% filter(!is.na(NES), !is.na(pvalue))
  # Simplify on significant terms only
  gsea_mf_sig          <- gsea_mf_clean
  gsea_mf_sig@result   <- gsea_mf_clean@result %>% filter(p.adjust < 0.05)
  gsea_mf_simple       <- simplify(gsea_mf_clean, cutoff = 0.7, by = "p.adjust", select_fun = min)
  
  # --- GO Slim GSEA (BP) ---
  gsea_bp_slim <- GSEA(
    geneList      = gene_list,
    TERM2GENE     = slim_term2gene,
    TERM2NAME     = slim_term2name,
    eps           = 0,
    minGSSize     = 15,
    maxGSSize     = 500,
    pvalueCutoff  = 1,
    pAdjustMethod = "BH",
    seed          = 26,
    verbose       = FALSE
  )
  
  # --- KEGG GSEA ---
  ids       <- bitr(names(original_gene_list), fromType = "ENSEMBL",
                    toType = "ENTREZID", OrgDb = organism)
  dedup_ids <- ids[!duplicated(ids$ENSEMBL), ]
  df2       <- df[df$Ensembl %in% dedup_ids$ENSEMBL, ]
  df2$Y     <- dedup_ids$ENTREZID[match(df2$Ensembl, dedup_ids$ENSEMBL)]
  
  entrez_gene_list        <- df2$t
  names(entrez_gene_list) <- df2$Y
  entrez_gene_list        <- na.omit(entrez_gene_list)
  entrez_gene_list        <- sort(entrez_gene_list, decreasing = TRUE)
  
  kegg <- gseKEGG(
    geneList      = entrez_gene_list,
    organism      = "hsa",
    eps           = 0,
    minGSSize     = 15,
    maxGSSize     = 500,
    pvalueCutoff  = 1,
    seed          = 26,
    pAdjustMethod = "BH",
    keyType       = "ncbi-geneid"
  )
  
  # --- GO ORA ---
  sig_genes  <- df %>% filter(Pi.Val < 0.05)
  
  sig_genes_up   <- sig_genes %>% filter(logFC > 0)
  sig_genes_down <- sig_genes %>% filter(logFC < 0)
  
  sig_entrez_up   <- bitr(sig_genes_up$Ensembl,   fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = organism)
  sig_entrez_down <- bitr(sig_genes_down$Ensembl, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = organism)
  
  bg_entrez  <- bitr(df$Ensembl, fromType = "ENSEMBL",
                     toType = "ENTREZID", OrgDb = organism)
  
  run_go_ora <- function(entrez_ids, ont_type, direction) {
    tryCatch(
      enrichGO(
        gene          = entrez_ids$ENTREZID,
        universe      = bg_entrez$ENTREZID, # remove to have full GO background
        OrgDb         = organism,
        ont           = ont_type,
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.2,
        readable      = TRUE
      ),
      error = function(e) {
        cat(sprintf("   GO ORA (%s, %s) failed for %s: %s\n", ont_type, direction, label, e$message))
        NULL
      }
    )
  }
  
  ora_bp_up   <- run_go_ora(sig_entrez_up,   "BP", "UP")
  ora_cc_up   <- run_go_ora(sig_entrez_up,   "CC", "UP")
  ora_mf_up   <- run_go_ora(sig_entrez_up,   "MF", "UP")
  ora_bp_down <- run_go_ora(sig_entrez_down, "BP", "DOWN")
  ora_cc_down <- run_go_ora(sig_entrez_down, "CC", "DOWN")
  ora_mf_down <- run_go_ora(sig_entrez_down, "MF", "DOWN")
  
  # --- Save .rds ---
  results <- list(
    label           = label,
    gsea_bp_clean   = gsea_bp_clean,
    gsea_bp_simple  = gsea_bp_simple,
    gsea_bp_slim    = gsea_bp_slim,
    gsea_cc_clean   = gsea_cc_clean,
    gsea_cc_simple  = gsea_cc_simple,
    gsea_mf_clean   = gsea_mf_clean,
    gsea_mf_simple  = gsea_mf_simple,
    gsea_kegg       = kegg,
    gsea_hall       = gsea_hall,
    gsea_reactome   = gsea_reactome,
    ora_bp_up       = ora_bp_up,
    ora_cc_up       = ora_cc_up,
    ora_mf_up       = ora_mf_up,
    ora_bp_down     = ora_bp_down,
    ora_cc_down     = ora_cc_down,
    ora_mf_down     = ora_mf_down
  )
  
  gsea_dir <- file.path(config$data_dir)
  if (!dir.exists(gsea_dir)) dir.create(gsea_dir, recursive = TRUE)
  saveRDS(results, file.path(gsea_dir, sprintf("gsea_results_%s.rds", label)))
  cat(sprintf("   Saved .rds for: %s\n", label))
  
  invisible(results)
}


# --- ORA ONLY Pipeline Function ---
run_ora_pipeline <- function(df, label, bg_df) {
  
  cat(sprintf("\n>> Running ORA pipeline for: %s\n", label))
  
  # --- Background: all detected proteins in the full results ---
  bg_entrez <- bitr(bg_df$Ensembl, fromType = "ENSEMBL",
                    toType = "ENTREZID", OrgDb = organism)
  
  # --- Gene lists: split by direction if mixed, or use all if already directional ---
  sig_entrez_up   <- bitr(df %>% filter(logFC > 0) %>% pull(Ensembl),
                          fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = organism)
  sig_entrez_down <- bitr(df %>% filter(logFC < 0) %>% pull(Ensembl),
                          fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = organism)
  
  run_go_ora <- function(entrez_ids, ont_type, direction) {
    if (nrow(entrez_ids) == 0) return(NULL)
    tryCatch(
      enrichGO(
        gene          = entrez_ids$ENTREZID,
        universe      = bg_entrez$ENTREZID,
        OrgDb         = organism,
        ont           = ont_type,
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.2,
        readable      = TRUE
      ),
      error = function(e) {
        cat(sprintf("   GO ORA (%s, %s) failed for %s: %s\n", ont_type, direction, label, e$message))
        NULL
      }
    )
  }
  
  ora_bp_up   <- run_go_ora(sig_entrez_up,   "BP", "UP")
  ora_cc_up   <- run_go_ora(sig_entrez_up,   "CC", "UP")
  ora_mf_up   <- run_go_ora(sig_entrez_up,   "MF", "UP")
  ora_bp_down <- run_go_ora(sig_entrez_down, "BP", "DOWN")
  ora_cc_down <- run_go_ora(sig_entrez_down, "CC", "DOWN")
  ora_mf_down <- run_go_ora(sig_entrez_down, "MF", "DOWN")
  
  results <- list(
    label       = label,
    ora_bp_up   = ora_bp_up,
    ora_cc_up   = ora_cc_up,
    ora_mf_up   = ora_mf_up,
    ora_bp_down = ora_bp_down,
    ora_cc_down = ora_cc_down,
    ora_mf_down = ora_mf_down
  )
  
  gsea_dir <- file.path(config$data_dir)
  if (!dir.exists(gsea_dir)) dir.create(gsea_dir, recursive = TRUE)
  saveRDS(results, file.path(gsea_dir, sprintf("gsea_results_%s.rds", label)))
  cat(sprintf("   Saved .rds for: %s\n", label))
  
  invisible(results)
}

# --- 6a: Run GSEA pipeline for a specific dataset ----
run_gsea_pipeline(limma_results$Baseline, "Baseline")
run_gsea_pipeline(limma_results$Training_BFR, "Training_BFR")
run_gsea_pipeline(limma_results$Training_HLRT, "Training_HLRT")
run_gsea_pipeline(limma_results$Post_training, "Post_training")
run_gsea_pipeline(limma_results$Interaction, "Interaction")

# --- 6b: CONCORDANT / DIVERGENT ORA (Training_BFR vs Training_HLRT) --------
# Added to support fig_05_panel_A_protein_BFR.R, which plots a
# Training_BFR-vs-Training_HLRT protein concordance scatter and needs
# ORA results for the two flanking-quadrant categories it defines:
#   - CONCORDANT: significant in both contrasts, same direction (TR+BL
#     quadrants of that scatter)
#   - DIVERGENT: significant in both contrasts, opposite direction (TL+BR
#     quadrants)
# This replaces the previous REVERSAL/EXACERBATION block, which was
# leftover from the CR pipeline (referenced diff_only_list_str, rev_df,
# ex_df, and limma_results$Baseline_SURVvCTL -- none of which exist here
# -- and would have errored if run).
#
# Full detected proteome as background. Any contrast's limma sheet works
# for this (same tested universe); Training_BFR is used arbitrarily.
bg_full <- limma_results$Training_BFR

df_bfr_ora <- limma_results$Training_BFR %>%
  dplyr::select(Ensembl, gene, logFC_BFR = logFC, Pi.Val_BFR = Pi.Val)
df_hlrt_ora <- limma_results$Training_HLRT %>%
  dplyr::select(Ensembl, gene, logFC_HLRT = logFC, Pi.Val_HLRT = Pi.Val)

df_train_ora <- inner_join(df_bfr_ora, df_hlrt_ora, by = c("Ensembl", "gene")) %>%
  mutate(sig = case_when(
    Pi.Val_BFR < 0.05 & Pi.Val_HLRT < 0.05 ~ "Both",
    Pi.Val_BFR < 0.05                      ~ "BFR only",
    Pi.Val_HLRT < 0.05                     ~ "HLRT only",
    TRUE                                   ~ "Neither"
  ))

# logFC = logFC_HLRT is the up/down split key passed into run_ora_pipeline
# (matches the y-axis convention used in fig_05_panel_A_protein_BFR.R's
# scatter, where HLRT is the y-axis contrast).
concordant_df <- df_train_ora %>%
  filter(sig == "Both", (logFC_BFR > 0) == (logFC_HLRT > 0)) %>%
  mutate(logFC = logFC_HLRT)

divergent_df <- df_train_ora %>%
  filter(sig == "Both", (logFC_BFR > 0) != (logFC_HLRT > 0)) %>%
  mutate(logFC = logFC_HLRT)

cat(sprintf("\nCONCORDANT proteins: %d | DIVERGENT proteins: %d\n",
            nrow(concordant_df), nrow(divergent_df)))

run_ora_pipeline(concordant_df, "CONCORDANT", bg_df = bg_full)
run_ora_pipeline(divergent_df,  "DIVERGENT",  bg_df = bg_full)