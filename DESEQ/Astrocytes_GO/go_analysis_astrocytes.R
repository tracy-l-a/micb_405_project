# =============================================================================
# GO enrichment analysis — astrocytes, LPS vs saline
# MICB 405 Group 10 (Srinivasan et al. 2016, mouse, FACS-sorted astrocytes)
#
# Input : DESEQ/Astrocytes_DESeq2_results.csv from DESEQ/Deseq.R
#         (results(..., lfcThreshold = 1, altHypothesis = "greaterAbs"))
#         Also accepts an .rds of a DESeq2 results object.
# Output: tables (CSV) + figures (PDF) in OUT_DIR
#
# Two complementary analyses:
#   1. ORA  (over-representation, enrichGO) — up- and down-regulated genes
#           tested separately, against a background of all genes DESeq2 tested
#   2. GSEA (gseGO) — all tested genes ranked by the Wald statistic; no cutoff
#
# Note on the lfcThreshold test used in Deseq.R:
#   - 'padj' already tests |log2FC| > 1, so every gene with padj < 0.05 also
#     has |log2FC| > 1. LFC_CUT below is a no-op for this input (kept so the
#     script still works on results run without lfcThreshold).
#   - 'stat' is still the ordinary Wald statistic (log2FC / lfcSE), but
#     'pvalue' is from the threshold test. GSEA therefore ranks on 'stat';
#     do not switch it to a p-value-based ranking for this input.
# =============================================================================

# ---- 0. Packages ------------------------------------------------------------
# One-time install (uncomment if needed):
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("clusterProfiler", "org.Mm.eg.db", "enrichplot", "DOSE"))
# install.packages(c("ggplot2", "dplyr", "readr"))

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(dplyr)
  library(readr)
})

# ---- 1. Settings — edit these -----------------------------------------------
# Paths are relative to the project root (where micb_405_project.Rproj lives).
INPUT_FILE <- "DESEQ/Astrocytes_DESeq2_results.csv"
OUT_DIR    <- "DESEQ/Astrocytes_GO"
GENE_COL   <- "gene_id"   # unversioned Ensembl IDs (ENSMUSG...) written by Deseq.R
PADJ_CUT   <- 0.05     # DE gene cutoff for ORA
LFC_CUT    <- 1        # |log2FC| cutoff for ORA; already enforced by padj here (see header)
ONTOLOGIES <- c("BP", "MF", "CC")   # BP is usually the one you'll report
SEED       <- 405

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
set.seed(SEED)

# ---- 2. Load DESeq2 results -------------------------------------------------
if (grepl("\\.rds$", INPUT_FILE, ignore.case = TRUE)) {
  res <- as.data.frame(readRDS(INPUT_FILE))
  res$gene_id <- rownames(res)
} else {
  res <- read.csv(INPUT_FILE, check.names = FALSE, stringsAsFactors = FALSE)
  if (is.null(GENE_COL)) {
    # write.csv(as.data.frame(res)) puts rownames in an unnamed first column
    GENE_COL <- names(res)[1]
  }
  if (!GENE_COL %in% names(res)) stop("Column '", GENE_COL, "' not found in ", INPUT_FILE)
  res$gene_id <- as.character(res[[GENE_COL]])
}

stopifnot(all(c("log2FoldChange", "pvalue", "padj") %in% names(res)))

# Strip Ensembl version suffixes (ENSMUSG00000000001.4 -> ENSMUSG00000000001)
res$gene_id <- sub("\\.\\d+$", "", res$gene_id)

# Auto-detect ID type
detect_keytype <- function(ids) {
  ids <- na.omit(ids)
  if (mean(grepl("^ENSMUSG", ids)) > 0.5) return("ENSEMBL")
  if (mean(grepl("^\\d+$",   ids)) > 0.5) return("ENTREZID")
  "SYMBOL"
}
KEYTYPE  <- detect_keytype(res$gene_id)
READABLE <- KEYTYPE != "SYMBOL"   # convert IDs to symbols in output tables
message("Detected gene ID type: ", KEYTYPE)

# Gene symbols for labelling plots
if (KEYTYPE == "SYMBOL") {
  res$symbol <- res$gene_id
} else if ("gene_name" %in% names(res)) {
  res$symbol <- res$gene_name   # symbols from the GTF, already in Deseq.R output
} else {
  res$symbol <- unname(mapIds(org.Mm.eg.db, keys = res$gene_id,
                              column = "SYMBOL", keytype = KEYTYPE,
                              multiVals = "first"))
}

# ---- 3. Define gene sets ----------------------------------------------------
# Background = every gene DESeq2 actually tested (padj not NA after
# independent filtering). Using the whole genome as background inflates
# enrichment, so this matters.
tested <- res %>%
  filter(!is.na(padj), !is.na(gene_id)) %>%
  distinct(gene_id, .keep_all = TRUE)

universe <- tested$gene_id
up   <- tested %>% filter(padj < PADJ_CUT, log2FoldChange >  LFC_CUT) %>% pull(gene_id)
down <- tested %>% filter(padj < PADJ_CUT, log2FoldChange < -LFC_CUT) %>% pull(gene_id)

message(sprintf("Background: %d genes | Up in LPS: %d | Down in LPS: %d",
                length(universe), length(up), length(down)))

# Named log2FC vector (by symbol) for cnetplot colouring
fc_by_symbol <- setNames(tested$log2FoldChange, tested$symbol)
fc_by_symbol <- fc_by_symbol[!is.na(names(fc_by_symbol)) & !duplicated(names(fc_by_symbol))]

# ---- Helpers ----------------------------------------------------------------
n_terms <- function(x) if (is.null(x)) 0 else nrow(as.data.frame(x))

save_plot <- function(p, file, n = 15, width = 8) {
  ggsave(file.path(OUT_DIR, file), p, width = width,
         height = max(4, 0.3 * n + 2), limitsize = FALSE)
}

run_ora <- function(genes, ont) {
  if (length(genes) < 5) return(NULL)
  ego <- enrichGO(gene          = genes,
                  universe      = universe,
                  OrgDb         = org.Mm.eg.db,
                  keyType       = KEYTYPE,
                  ont           = ont,
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.2,
                  minGSSize     = 10,
                  maxGSSize     = 500,
                  readable      = READABLE)
  if (n_terms(ego) == 0) return(NULL)
  # Collapse redundant GO terms (semantic similarity > 0.7).
  # Namespaced because igraph (and others) also export simplify().
  clusterProfiler::simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min)
}

# Term counts for the summary table (records zeros too)
term_counts <- c()

# ---- 4. ORA: up- and down-regulated genes separately ------------------------
ora <- list()
for (ont in ONTOLOGIES) {
  for (direction in c("up", "down")) {
    genes <- if (direction == "up") up else down
    key   <- paste("ORA", direction, ont, sep = "_")
    message(key)
    ego <- run_ora(genes, ont)
    term_counts[key] <- n_terms(ego)

    if (n_terms(ego) == 0) {
      message("  no significant terms")
      next
    }
    ora[[key]] <- ego
    write_csv(as.data.frame(ego), file.path(OUT_DIR, paste0(key, ".csv")))

    title <- sprintf("Astrocytes: GO %s, %s in LPS", ont,
                     ifelse(direction == "up", "up-regulated", "down-regulated"))
    save_plot(dotplot(ego, showCategory = 15) + ggtitle(title),
              paste0(key, "_dotplot.pdf"))
  }
}

# Extra BP figures: term network + gene–concept network
for (direction in c("up", "down")) {
  ego <- ora[[paste("ORA", direction, "BP", sep = "_")]]
  if (n_terms(ego) < 2) next

  ego_sim <- pairwise_termsim(ego)
  save_plot(emapplot(ego_sim, showCategory = 30) +
              ggtitle(sprintf("Astrocytes: GO BP term map (%s in LPS)", direction)),
            sprintf("ORA_%s_BP_emapplot.pdf", direction), n = 30, width = 10)

  # 'foldChange =' is accepted by every enrichplot version (deprecated with a
  # warning in 1.20–1.24, current again from 1.26); 'color.params' is not.
  p_cnet <- suppressWarnings(
    cnetplot(ego, showCategory = 5, foldChange = fc_by_symbol)
  )
  save_plot(p_cnet + ggtitle(sprintf("Astrocytes: top GO BP terms and genes (%s)", direction)),
            sprintf("ORA_%s_BP_cnetplot.pdf", direction), n = 30, width = 11)
}

# Up vs down side by side (BP) — handy as a single summary figure
cc_input <- list(`Up in LPS` = up, `Down in LPS` = down)
cc_input <- cc_input[lengths(cc_input) >= 5]
if (length(cc_input) > 0) {
  ck <- compareCluster(geneClusters  = cc_input,
                       fun           = "enrichGO",
                       universe      = universe,
                       OrgDb         = org.Mm.eg.db,
                       keyType       = KEYTYPE,
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       readable      = READABLE)
  if (n_terms(ck) > 0) {
    ck <- clusterProfiler::simplify(ck, cutoff = 0.7, by = "p.adjust", select_fun = min)
    write_csv(as.data.frame(ck), file.path(OUT_DIR, "ORA_compare_up_down_BP.csv"))
    save_plot(dotplot(ck, showCategory = 10) +
                ggtitle("Astrocytes: GO BP, LPS vs saline"),
              "ORA_compare_up_down_BP_dotplot.pdf", n = 20, width = 9)
  }
}

# ---- 5. GSEA: all tested genes, ranked -------------------------------------
# Rank by the DESeq2 Wald statistic if present (combines effect size and
# confidence). If the table came from lfcShrink() there is no 'stat' column,
# so fall back to sign(log2FC) * -log10(pvalue).
rank_df <- tested %>% filter(!is.na(pvalue))
if ("stat" %in% names(rank_df) && !all(is.na(rank_df$stat))) {
  rank_df$score <- rank_df$stat
  message("GSEA ranking metric: Wald stat")
} else {
  rank_df$score <- sign(rank_df$log2FoldChange) * -log10(pmax(rank_df$pvalue, 1e-300))
  message("GSEA ranking metric: sign(log2FC) * -log10(pvalue)")
}
gene_list <- rank_df %>%
  filter(is.finite(score)) %>%
  arrange(desc(score)) %>%
  { setNames(.$score, .$gene_id) }

gsea <- list()
for (ont in ONTOLOGIES) {
  message("GSEA: ", ont)
  gse <- gseGO(geneList     = gene_list,
               OrgDb        = org.Mm.eg.db,
               keyType      = KEYTYPE,
               ont          = ont,
               minGSSize    = 10,
               maxGSSize    = 500,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH",
               eps          = 0,
               seed         = TRUE,
               verbose      = FALSE)
  if (n_terms(gse) == 0) {
    term_counts[paste0("GSEA_", ont)] <- 0
    message("  no significant terms")
    next
  }
  gse <- clusterProfiler::simplify(gse, cutoff = 0.7, by = "p.adjust", select_fun = min)
  term_counts[paste0("GSEA_", ont)] <- n_terms(gse)
  if (READABLE) gse <- setReadable(gse, OrgDb = org.Mm.eg.db, keyType = KEYTYPE)
  gsea[[ont]] <- gse

  write_csv(as.data.frame(gse), file.path(OUT_DIR, sprintf("GSEA_%s.csv", ont)))

  # Activated (NES > 0, higher in LPS) vs suppressed (NES < 0)
  p <- dotplot(gse, showCategory = 10, split = ".sign") +
    facet_grid(. ~ .sign) +
    ggtitle(sprintf("Astrocytes: GSEA GO %s (LPS vs saline)", ont))
  save_plot(p, sprintf("GSEA_%s_dotplot.pdf", ont), n = 20, width = 11)
}

# Running-score plots for the top BP gene sets
if (!is.null(gsea$BP)) {
  top_ids <- head(as.data.frame(gsea$BP) %>% arrange(p.adjust) %>% pull(ID), 3)
  p <- gseaplot2(gsea$BP, geneSetID = top_ids, pvalue_table = TRUE,
                 title = "Astrocytes: top GO BP gene sets")
  ggsave(file.path(OUT_DIR, "GSEA_BP_top3_running_score.pdf"), p,
         width = 10, height = 7)
}

# ---- 6. Summary + reproducibility -------------------------------------------
summary_tbl <- tibble(analysis = names(term_counts),
                      n_terms  = unname(term_counts))
print(summary_tbl)
write_csv(summary_tbl, file.path(OUT_DIR, "summary_n_terms.csv"))

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
message("Done. Results in: ", normalizePath(OUT_DIR))
