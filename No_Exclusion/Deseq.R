
# ============================================================
# 1. LOAD PACKAGES AND SET ANALYSIS THRESHOLDS
# ============================================================

library(tidyverse)
library(DESeq2)
library(pheatmap)
library(ggplot2)
library(png)
library(grid)
library(ggrepel)

# Needed for the inflammatory-gene (GO:0006954) bold labels.
# Checked here so the script stops before the slow DESeq2 steps.
# org.Mm.eg.db is used via its namespace (not attached) so that
# AnnotationDbi::select() does not mask dplyr::select().
missing_bioc <- c("AnnotationDbi", "org.Mm.eg.db")[
  !vapply(c("AnnotationDbi", "org.Mm.eg.db"),
          requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_bioc) > 0) {
  stop(
    "Missing package(s): ", paste(missing_bioc, collapse = ", "),
    "\nInstall with:\n  BiocManager::install(c(\"",
    paste(missing_bioc, collapse = "\", \""), "\"))"
  )
}

LFC_THRESHOLD <- 1
PADJ_THRESHOLD <- 0.05

# Samples removed from the analysis (state the QC reason in the Methods)
EXCLUDED_SAMPLES <- c("SRR2938964")

# Treatment colours used in every plot
COND_COLOURS <- c("LPS" = "#D55E00", "Saline" = "#0072B2")

# Number of top DEGs per direction shown in each figure.
# Volcano labels (5 up + 5 down) are a subset of the heatmap genes
# (10 up + 10 down), because both use the same ranking (select_top_DEGs).
N_LABEL_PER_DIRECTION <- 5
N_HEATMAP_PER_DIRECTION <- 10


# ============================================================
# 2. LOAD GENE COUNT FILE
# ============================================================

counts <- read.delim(
  "gene_counts.tsv",
  comment.char = "#",
  header = TRUE,
  check.names = FALSE
)


# ============================================================
# 3. REMOVE UNNECESSARY COLUMNS
# ============================================================

# Remove Chr, Start, End, Strand
counts <- counts[, -c(2:5)]

# Rename first two columns
colnames(counts)[1:2] <- c("gene", "length")

# Extract SRR accession numbers
colnames(counts)[3:ncol(counts)] <- sub(
  ".*(SRR[0-9]+).*",
  "\\1",
  colnames(counts)[3:ncol(counts)]
)


# ============================================================
# 4. DEFINE SAMPLE GROUPS
# ============================================================

microglia_LPS <- c(
  "SRR2938958", "SRR2938964", "SRR2938974",
  "SRR2938976", "SRR2938981"
)

microglia_Saline <- c(
  "SRR2938960", "SRR2938961", "SRR2938967",
  "SRR2938971", "SRR2938980"
)

astrocytes_LPS <- c(
  "SRR2938963", "SRR2938973",
  "SRR2938977", "SRR2938985"
)

astrocytes_Saline <- c(
  "SRR2938962", "SRR2938968",
  "SRR2938969", "SRR2938972"
)

# Remove excluded samples from every group
microglia_LPS     <- setdiff(microglia_LPS, EXCLUDED_SAMPLES)
microglia_Saline  <- setdiff(microglia_Saline, EXCLUDED_SAMPLES)
astrocytes_LPS    <- setdiff(astrocytes_LPS, EXCLUDED_SAMPLES)
astrocytes_Saline <- setdiff(astrocytes_Saline, EXCLUDED_SAMPLES)

cat("Excluded samples:", paste(EXCLUDED_SAMPLES, collapse = ", "), "\n")
cat("Microglia: LPS n =", length(microglia_LPS),
    "| Saline n =", length(microglia_Saline), "\n")
cat("Astrocytes: LPS n =", length(astrocytes_LPS),
    "| Saline n =", length(astrocytes_Saline), "\n")


# ============================================================
# 5. FILTER COUNT DATA
# ============================================================

counts <- counts[, setdiff(colnames(counts), EXCLUDED_SAMPLES)]

sample_cols <- colnames(counts)[3:ncol(counts)]

counts[sample_cols] <- lapply(
  counts[sample_cols],
  as.numeric
)

counts <- counts[
  complete.cases(counts[, sample_cols]),
]

counts <- counts[
  rowSums(counts[, sample_cols]) > 0,
]

cat("Dimensions after count filtering:\n")
print(dim(counts))

cat("Remaining NA values:", sum(is.na(counts[, sample_cols])), "\n")
cat("All-zero genes:", sum(rowSums(counts[, sample_cols]) == 0), "\n")


# ============================================================
# 6. CHECK SAMPLE NAMES
# ============================================================

all_samples <- c(
  microglia_LPS,
  microglia_Saline,
  astrocytes_LPS,
  astrocytes_Saline
)

missing_samples <- setdiff(all_samples, colnames(counts))

if (length(missing_samples) > 0) {
  stop(
    paste(
      "Missing sample columns:",
      paste(missing_samples, collapse = ", ")
    )
  )
}

stopifnot(!anyDuplicated(sample_cols))


# ============================================================
# 7. ANNOTATE GENES USING GENCODE M39
# ============================================================

gtf <- read.delim(
  "gencode.vM39.primary_assembly.annotation.gtf",
  header = FALSE,
  sep = "\t",
  comment.char = "#",
  quote = "",
  stringsAsFactors = FALSE
)

colnames(gtf) <- c(
  "chromosome", "source", "feature", "start", "end",
  "score", "strand", "frame", "attributes"
)

gene_annotation <- gtf %>%
  filter(feature == "gene") %>%
  mutate(
    gene_id = stringr::str_match(
      attributes,
      'gene_id "([^"]+)"'
    )[, 2],
    gene_name = stringr::str_match(
      attributes,
      'gene_name "([^"]+)"'
    )[, 2]
  ) %>%
  select(gene_id, gene_name) %>%
  filter(!is.na(gene_id)) %>%
  mutate(gene_id = sub("\\.[0-9]+$", "", gene_id)) %>%
  distinct(gene_id, .keep_all = TRUE)

counts <- counts %>%
  mutate(gene_id = sub("\\.[0-9]+$", "", gene))

matched_genes <- sum(counts$gene_id %in% gene_annotation$gene_id)

cat(
  "Genes matched:", matched_genes, "out of", nrow(counts), "\n"
)

cat(
  "Percentage matched:",
  round(100 * matched_genes / nrow(counts), 2),
  "%\n"
)

counts <- counts %>%
  left_join(gene_annotation, by = "gene_id") %>%
  mutate(
    gene_name = if_else(
      is.na(gene_name) | gene_name == "",
      gene,
      gene_name
    )
  )


# ============================================================
# 8. CREATE FOUR SEPARATE DATASETS
# ============================================================

microglia_LPS_counts <- counts %>%
  select(gene, gene_id, gene_name, length, all_of(microglia_LPS))

microglia_Saline_counts <- counts %>%
  select(gene, gene_id, gene_name, length, all_of(microglia_Saline))

astrocytes_LPS_counts <- counts %>%
  select(gene, gene_id, gene_name, length, all_of(astrocytes_LPS))

astrocytes_Saline_counts <- counts %>%
  select(gene, gene_id, gene_name, length, all_of(astrocytes_Saline))


# ============================================================
# 9. EXPORT FOUR COUNT CSV FILES
# ============================================================

write.csv(microglia_LPS_counts, "Microglia_LPS.csv", row.names = FALSE)
write.csv(microglia_Saline_counts, "Microglia_Saline.csv", row.names = FALSE)
write.csv(astrocytes_LPS_counts, "Astrocytes_LPS.csv", row.names = FALSE)
write.csv(astrocytes_Saline_counts, "Astrocytes_Saline.csv", row.names = FALSE)


# ============================================================
# 10. CREATE FUNCTION TO RUN DESEQ2
# ============================================================

run_deseq <- function(LPS, Saline, cell_type) {
  
  genes <- LPS$gene
  
  stopifnot(identical(LPS$gene, Saline$gene))
  
  LPS_counts <- LPS %>%
    select(-gene, -gene_id, -gene_name, -length)
  
  Saline_counts <- Saline %>%
    select(-gene, -gene_id, -gene_name, -length)
  
  count_matrix <- as.matrix(cbind(LPS_counts, Saline_counts))
  
  stopifnot(
    all(is.finite(count_matrix)),
    all(count_matrix >= 0),
    all(count_matrix == round(count_matrix))
  )
  
  storage.mode(count_matrix) <- "integer"
  rownames(count_matrix) <- genes
  
  sample_info <- data.frame(
    condition = factor(
      c(
        rep("LPS", ncol(LPS_counts)),
        rep("Saline", ncol(Saline_counts))
      ),
      levels = c("Saline", "LPS")
    )
  )
  
  rownames(sample_info) <- colnames(count_matrix)
  
  dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = sample_info,
    design = ~ condition
  )
  
  # Filter genes with fewer than 10 counts in at least 2 samples
  keep <- rowSums(counts(dds) >= 10) >= 2
  dds <- dds[keep, ]
  
  # Run DESeq2
  dds <- DESeq(dds)
  
  # IMPORTANT:
  # Test whether the absolute log2 fold change is greater than 1.
  # The adjusted p-values correspond to this threshold-based test.
  res <- results(
    dds,
    contrast = c("condition", "LPS", "Saline"),
    lfcThreshold = LFC_THRESHOLD,
    altHypothesis = "greaterAbs",
    alpha = PADJ_THRESHOLD
  )
  
  res_df <- as.data.frame(res) %>%
    rownames_to_column("gene") %>%
    mutate(gene_id = sub("\\.[0-9]+$", "", gene)) %>%
    left_join(gene_annotation, by = "gene_id") %>%
    mutate(
      gene_name = if_else(
        is.na(gene_name) | gene_name == "",
        gene,
        gene_name
      )
    ) %>%
    arrange(padj)
  
  # Export all DESeq2 results
  write.csv(
    res_df,
    paste0(cell_type, "_DESeq2_results.csv"),
    row.names = FALSE
  )
  
  # Extract and export normalized counts
  normalized_counts <- counts(dds, normalized = TRUE) %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    mutate(gene_id = sub("\\.[0-9]+$", "", gene)) %>%
    left_join(gene_annotation, by = "gene_id") %>%
    mutate(
      gene_name = if_else(
        is.na(gene_name) | gene_name == "",
        gene,
        gene_name
      )
    )
  
  write.csv(
    normalized_counts,
    paste0(cell_type, "_normalized_counts.csv"),
    row.names = FALSE
  )
  
  return(
    list(
      dds = dds,
      results = res_df,
      normalized_counts = normalized_counts
    )
  )
}


# ============================================================
# 11. RUN DESEQ2 FOR MICROGLIA
# ============================================================

microglia_results <- run_deseq(
  LPS = microglia_LPS_counts,
  Saline = microglia_Saline_counts,
  cell_type = "Microglia"
)


# ============================================================
# 12. RUN DESEQ2 FOR ASTROCYTES
# ============================================================

astrocytes_results <- run_deseq(
  LPS = astrocytes_LPS_counts,
  Saline = astrocytes_Saline_counts,
  cell_type = "Astrocytes"
)


# ============================================================
# 13. EXTRACT DESEQ2 OBJECTS AND RESULTS
# ============================================================

microglia_dds <- microglia_results$dds
astrocytes_dds <- astrocytes_results$dds

microglia_dat <- microglia_results$results
astrocytes_dat <- astrocytes_results$results


# ============================================================
# 14. IDENTIFY SIGNIFICANT DEGS
# ============================================================

# Apply both adjusted p-value and absolute fold-change criteria.
microglia_DEGs <- microglia_dat %>%
  filter(
    !is.na(padj),
    !is.na(log2FoldChange),
    padj < PADJ_THRESHOLD,
    abs(log2FoldChange) >= LFC_THRESHOLD
  )

astrocytes_DEGs <- astrocytes_dat %>%
  filter(
    !is.na(padj),
    !is.na(log2FoldChange),
    padj < PADJ_THRESHOLD,
    abs(log2FoldChange) >= LFC_THRESHOLD
  )

write.csv(
  microglia_DEGs,
  "Microglia_significant_DEGs.csv",
  row.names = FALSE
)

write.csv(
  astrocytes_DEGs,
  "Astrocytes_significant_DEGs.csv",
  row.names = FALSE
)

cat("\nNumber of significant DEGs:\n")
cat("Microglia:", nrow(microglia_DEGs), "\n")
cat("Astrocytes:", nrow(astrocytes_DEGs), "\n")


# ============================================================
# 15. COUNT UPREGULATED AND DOWNREGULATED GENES
# ============================================================

microglia_DEGs %>%
  summarise(
    Upregulated = sum(log2FoldChange >= LFC_THRESHOLD),
    Downregulated = sum(log2FoldChange <= -LFC_THRESHOLD)
  ) %>%
  print()

astrocytes_DEGs %>%
  summarise(
    Upregulated = sum(log2FoldChange >= LFC_THRESHOLD),
    Downregulated = sum(log2FoldChange <= -LFC_THRESHOLD)
  ) %>%
  print()


# ============================================================
# 15b. SHARED TOP-DEG SELECTION (used by volcano, top-5 bar plot
#      and heatmap so all three figures show the same genes)
# ============================================================

# Rank significant DEGs by adjusted p-value (from the |log2FC| > 1
# threshold test, so it reflects both effect size and reliability),
# breaking ties by larger |log2FC|.
# Ranking by raw log2FC alone favours low-count genes whose fold
# changes are unstable (e.g. a gene with baseMean ~28 and log2FC -23).
select_top_DEGs <- function(dat, n_per_direction) {

  sig <- dat %>%
    filter(
      !is.na(padj),
      !is.na(log2FoldChange),
      padj < PADJ_THRESHOLD,
      abs(log2FoldChange) >= LFC_THRESHOLD
    ) %>%
    mutate(
      direction = if_else(log2FoldChange > 0, "UP", "DOWN"),
      display_name = if_else(
        is.na(gene_name) | gene_name == "",
        gene,
        gene_name
      )
    )

  top_up <- sig %>%
    filter(direction == "UP") %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    slice_head(n = n_per_direction)

  top_down <- sig %>%
    filter(direction == "DOWN") %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    slice_head(n = n_per_direction)

  top <- bind_rows(top_up, top_down)

  dup <- duplicated(top$display_name) |
    duplicated(top$display_name, fromLast = TRUE)

  top$display_name[dup] <- paste0(
    top$display_name[dup], " (", top$gene_id[dup], ")"
  )

  top
}


# ============================================================
# 15c. INFLAMMATORY GENE SET (for bold labels only)
# ============================================================

# Defined a priori from GO:0006954 "inflammatory response"
# (including all child terms) in org.Mm.eg.db.
# This list is ONLY used to bold labels. It does not change which
# genes are selected or plotted.
# org.Mm.eg.db is called via its namespace and is NOT attached,
# so AnnotationDbi::select() does not mask dplyr::select().
INFLAMMATORY_GO_ID <- "GO:0006954"

inflammatory_map <- suppressMessages(
  AnnotationDbi::select(
    org.Mm.eg.db::org.Mm.eg.db,
    keys = INFLAMMATORY_GO_ID,
    keytype = "GOALL",
    columns = c("ENSEMBL", "SYMBOL")
  )
)

inflammatory_ensembl <- unique(na.omit(inflammatory_map$ENSEMBL))
inflammatory_symbols <- unique(na.omit(inflammatory_map$SYMBOL))

cat(
  "Inflammatory response genes (", INFLAMMATORY_GO_ID, "):",
  length(inflammatory_symbols), "symbols\n"
)

is_inflammatory <- function(gene_id, gene_name) {
  gene_id %in% inflammatory_ensembl | gene_name %in% inflammatory_symbols
}


# ============================================================
# 16. CREATE OUTPUT FOLDERS
# ============================================================

dir.create("Microglia_plots", showWarnings = FALSE)
dir.create("Astrocytes_plots", showWarnings = FALSE)


# ============================================================
# 17. PCA PLOT
# ============================================================

plot_PCA <- function(dds, cell_type, output_folder) {
  
  rld <- rlog(dds, blind = FALSE)
  
  pcaData <- plotPCA(
    rld,
    intgroup = "condition",
    returnData = TRUE
  )
  
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  p <- ggplot(
    pcaData,
    aes(x = PC1, y = PC2, colour = condition)
  ) +
    geom_point(size = 4) +
    scale_colour_manual(values = COND_COLOURS) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    theme_test() +
    theme(
      axis.text = element_text(colour = "black", size = 11),
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      panel.border = element_rect(linewidth = 1.2),
      plot.title = element_text(face = "bold", size = 13)
    ) +
    labs(
      title = paste0(
        "Global Gene Expression Profiles of ",
        cell_type,
        ": LPS vs Saline"
      ),
      subtitle = "Principal component analysis of biological replicates",
      colour = "Treatment"
    )
  
  print(p)
  
  ggsave(
    paste0(output_folder, "/", cell_type, "_PCA.png"),
    plot = p,
    width = 7,
    height = 6,
    dpi = 300
  )
}


# ============================================================
# 18. VOLCANO PLOT
# ============================================================

plot_volcano <- function(dat, cell_type, output_folder) {
  
  labelled_dat <- dat %>%
    filter(
      !is.na(padj),
      !is.na(log2FoldChange)
    ) %>%
    mutate(
      up_down = case_when(
        padj < PADJ_THRESHOLD &
          log2FoldChange >= LFC_THRESHOLD ~ "UP",
        
        padj < PADJ_THRESHOLD &
          log2FoldChange <= -LFC_THRESHOLD ~ "DOWN",
        
        TRUE ~ "NONE"
      ),
      neg_log10_padj = -log10(
        pmax(padj, .Machine$double.xmin)
      ),
      gene_label = if_else(
        is.na(gene_name) | gene_name == "",
        gene,
        gene_name
      )
    )
  
  # Same genes/ranking as the top-5 bar plot and the heatmap
  top_genes <- select_top_DEGs(dat, N_LABEL_PER_DIRECTION)

  label_dat <- labelled_dat %>%
    filter(gene %in% top_genes$gene) %>%
    left_join(
      top_genes %>% select(gene, display_name),
      by = "gene"
    ) %>%
    mutate(
      gene_label = display_name,
      label_face = if_else(
        is_inflammatory(gene_id, gene_name),
        "bold.italic",
        "italic"
      )
    )

  x_range <- range(labelled_dat$log2FoldChange, na.rm = TRUE)
  y_range <- range(labelled_dat$neg_log10_padj, na.rm = TRUE)
  
  x_padding <- max(diff(x_range) * 0.15, 1)
  y_padding <- max(diff(y_range) * 0.15, 1)
  
  p <- ggplot(
    labelled_dat,
    aes(
      x = log2FoldChange,
      y = neg_log10_padj,
      colour = up_down
    )
  ) +
    geom_point(alpha = 0.7, size = 1.6) +
    
    ggrepel::geom_label_repel(
      data = label_dat,
      aes(label = gene_label, fontface = label_face),
      inherit.aes = TRUE,
      colour = "black",
      fill = "white",
      size = 3,
      label.padding = unit(0.15, "lines"),
      label.r = unit(0.1, "lines"),
      box.padding = 0.8,
      point.padding = 0.6,
      force = 3,
      force_pull = 0.5,
      direction = "both",
      min.segment.length = 0,
      segment.colour = "grey50",
      segment.linewidth = 0.4,
      max.overlaps = Inf,
      max.time = 5,
      show.legend = FALSE
    ) +
    
    # Adjusted p-value threshold
    geom_hline(
      yintercept = -log10(PADJ_THRESHOLD),
      linetype = "dashed",
      linewidth = 0.5,
      colour = "grey40"
    ) +
    
    # Absolute log2 fold-change thresholds
    geom_vline(
      xintercept = c(-LFC_THRESHOLD, LFC_THRESHOLD),
      linetype = "dashed",
      linewidth = 0.5,
      colour = "grey40"
    ) +
    
    coord_cartesian(
      xlim = c(
        x_range[1] - x_padding,
        x_range[2] + x_padding
      ),
      ylim = c(
        max(0, y_range[1] - y_padding),
        y_range[2] + y_padding
      ),
      clip = "off"
    ) +
    
    scale_colour_manual(
      values = c(
        "UP" = "#D73027",
        "DOWN" = "#4575B4",
        "NONE" = "grey75"
      ),
      breaks = c("UP", "DOWN", "NONE"),
      labels = c("Upregulated", "Downregulated", "Not significant")
    ) +
    
    theme_classic() +
    
    labs(
      title = paste(
        "Differential Gene Expression in",
        cell_type,
        "Following LPS Treatment"
      ),
      subtitle = paste0(
        "DEGs: padj < ", PADJ_THRESHOLD,
        " and |log2FC| >= ", LFC_THRESHOLD,
        "; top ", N_LABEL_PER_DIRECTION,
        " up/down by adjusted p-value labelled"
      ),
      caption = paste0(
        "Bold labels: genes annotated to ", INFLAMMATORY_GO_ID,
        " (inflammatory response)"
      ),
      x = "Log2 Fold Change (LPS vs Saline)",
      y = "-Log10 Adjusted P-value",
      colour = "Expression Direction"
    ) +
    
    theme(
      axis.title = element_text(face = "bold", size = 11),
      axis.text = element_text(colour = "black", size = 10),
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9),
      legend.title = element_text(face = "bold"),
      plot.margin = margin(15, 25, 15, 15)
    )
  
  print(p)
  
  ggsave(
    paste0(output_folder, "/", cell_type, "_Volcano.png"),
    plot = p,
    width = 10,
    height = 8,
    dpi = 300,
    bg = "white"
  )
}


# ============================================================
# 19. DEG BAR PLOT
# ============================================================

plot_DEG_bar <- function(dat, cell_type, output_folder) {
  
  labelled_dat <- dat %>%
    filter(
      !is.na(padj),
      !is.na(log2FoldChange),
      padj < PADJ_THRESHOLD,
      abs(log2FoldChange) >= LFC_THRESHOLD
    ) %>%
    mutate(
      abs_l2fc = abs(log2FoldChange),
      up_down = if_else(
        log2FoldChange >= LFC_THRESHOLD,
        "UP",
        "DOWN"
      ),
      l2fc = case_when(
        abs_l2fc >= 4 ~ ">=4",
        abs_l2fc >= 3 ~ "3-<4",
        abs_l2fc >= 2 ~ "2-<3",
        TRUE ~ "1-<2"
      )
    )
  
  if (nrow(labelled_dat) == 0) {
    message(paste("No significant DEGs found for", cell_type))
    return(invisible(NULL))
  }
  
  labelled_dat$l2fc <- factor(
    labelled_dat$l2fc,
    levels = c("1-<2", "2-<3", "3-<4", ">=4")
  )
  
  p <- ggplot(
    labelled_dat,
    aes(x = l2fc, fill = up_down)
  ) +
    geom_bar(colour = "black") +
    scale_fill_manual(
      values = c("UP" = "#D73027", "DOWN" = "#4575B4"),
      labels = c("UP" = "Upregulated", "DOWN" = "Downregulated")
    ) +
    theme_classic() +
    labs(
      title = paste(
        "Distribution of Significant Gene Expression Changes in",
        cell_type
      ),
      subtitle = paste0(
        "DEGs: padj < ", PADJ_THRESHOLD,
        " and |log2FC| >= ", LFC_THRESHOLD
      ),
      x = "Absolute Log2 Fold Change",
      y = "Number of Significant Genes",
      fill = "Expression Direction"
    ) +
    theme(
      axis.title = element_text(face = "bold"),
      axis.text = element_text(colour = "black"),
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9)
    )
  
  print(p)
  
  ggsave(
    paste0(output_folder, "/", cell_type, "_DEG_barplot.png"),
    plot = p,
    width = 7,
    height = 6,
    dpi = 300
  )
}


# ============================================================
# 20. TOP 5 UPREGULATED AND DOWNREGULATED GENES
# ============================================================

plot_top5 <- function(dat, cell_type, output_folder) {
  
  # Same genes/ranking as the volcano labels
  top_genes <- select_top_DEGs(dat, N_LABEL_PER_DIRECTION) %>%
    mutate(
      direction = if_else(
        log2FoldChange >= LFC_THRESHOLD,
        "Upregulated",
        "Downregulated"
      )
    )
  
  if (nrow(top_genes) == 0) {
    message(paste("No significant DEGs found for", cell_type))
    return(invisible(NULL))
  }
  
  top_genes <- top_genes %>%
    arrange(log2FoldChange) %>%
    mutate(display_name = factor(display_name, levels = display_name))
  
  p <- ggplot(
    top_genes,
    aes(x = display_name, y = log2FoldChange, fill = direction)
  ) +
    geom_hline(yintercept = 0, linewidth = 0.5) +
    geom_col(colour = "black", width = 0.8) +
    geom_errorbar(
      aes(
        ymin = log2FoldChange - lfcSE,
        ymax = log2FoldChange + lfcSE
      ),
      width = 0.3
    ) +
    coord_flip() +
    scale_fill_manual(
      values = c(
        "Upregulated" = "#D73027",
        "Downregulated" = "#4575B4"
      )
    ) +
    theme_classic() +
    labs(
      title = paste("Top 5 Upregulated and Downregulated Genes in", cell_type),
      subtitle = paste0(
        "DEGs: padj < ", PADJ_THRESHOLD,
        " and |log2FC| >= ", LFC_THRESHOLD,
        "; ranked by adjusted p-value (same genes as volcano labels)"
      ),
      x = "Gene name",
      y = "Log2 Fold Change (LPS vs Saline)",
      fill = "Expression Direction"
    ) +
    theme(
      axis.text = element_text(colour = "black"),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9)
    )
  
  print(p)
  
  ggsave(
    filename = paste0(output_folder, "/", cell_type, "_Top5_DEGs.png"),
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )
}


# ============================================================
# 21. NORMALIZED COUNTS PLOT
# ============================================================

plot_normalized_counts <- function(dds, dat, cell_type, output_folder) {
  
  normalized_counts <- counts(dds, normalized = TRUE) %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    pivot_longer(
      cols = -gene,
      names_to = "Sample",
      values_to = "count"
    )
  
  sample_info <- as.data.frame(colData(dds)) %>%
    rownames_to_column("Sample")
  
  normalized_counts <- normalized_counts %>%
    left_join(sample_info, by = "Sample")
  
  # Select the most statistically significant DEG meeting both criteria
  gene_of_interest <- dat %>%
    filter(
      !is.na(padj),
      !is.na(log2FoldChange),
      padj < PADJ_THRESHOLD,
      abs(log2FoldChange) >= LFC_THRESHOLD
    ) %>%
    arrange(padj) %>%
    slice_head(n = 1) %>%
    pull(gene)
  
  if (length(gene_of_interest) == 0) {
    message(paste("No significant DEGs found for", cell_type))
    return(invisible(NULL))
  }
  
  gene_symbol <- dat %>%
    filter(gene == gene_of_interest) %>%
    pull(gene_name)
  
  gene_symbol <- gene_symbol[1]
  
  if (length(gene_symbol) == 0 ||
      is.na(gene_symbol) ||
      gene_symbol == "") {
    gene_symbol <- gene_of_interest
  }
  
  gene_data <- normalized_counts %>%
    filter(gene == gene_of_interest)
  
  p <- ggplot(
    gene_data,
    aes(
      x = condition,
      y = count,
      colour = condition
    )
  ) +
    geom_jitter(width = 0.2, size = 3) +
    scale_colour_manual(values = COND_COLOURS) +
    scale_y_continuous(trans = "log1p") +
    theme_classic() +
    labs(
      title = paste(
        "Normalized Expression of",
        gene_symbol,
        "in",
        cell_type
      ),
      subtitle = paste0(
        "Selected DEG: padj < ", PADJ_THRESHOLD,
        " and |log2FC| >= ", LFC_THRESHOLD
      ),
      x = "Treatment",
      y = "Normalized Counts (log1p)",
      colour = "Treatment"
    ) +
    theme(
      axis.title = element_text(face = "bold"),
      axis.text = element_text(colour = "black"),
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9)
    )
  
  print(p)
  
  ggsave(
    paste0(output_folder, "/", cell_type, "_Normalized_Counts.png"),
    plot = p,
    width = 7,
    height = 6,
    dpi = 300
  )
}


# ============================================================
# 22. TOP 20 DEG HEATMAP (10 up + 10 down)
# ============================================================

# Previously this plotted the 20 genes with the highest VST variance
# across samples. That is an unsupervised selection, unrelated to the
# DESeq2 ranking used for the volcano plot, so the two figures showed
# different genes (and 2 of the 20 "variable" genes were not DEGs).
# Now the heatmap shows the top DEGs from the SAME results table and
# SAME ranking as the volcano labels; the 5 + 5 labelled volcano genes
# are a subset of these 10 + 10.

plot_heatmap <- function(dds, dat, cell_type, output_folder) {
  
  top_genes <- select_top_DEGs(dat, N_HEATMAP_PER_DIRECTION)
  
  if (nrow(top_genes) < 2) {
    message(paste("Not enough DEGs for heatmap in", cell_type))
    return(invisible(NULL))
  }
  
  vsd <- vst(dds, blind = FALSE)
  vsd_matrix <- assay(vsd)[top_genes$gene, , drop = FALSE]
  
  # Row z-scores of VST expression
  Z <- t(scale(t(vsd_matrix)))
  
  keep_rows <- apply(Z, 1, function(x) all(is.finite(x)))
  Z <- Z[keep_rows, , drop = FALSE]
  top_genes <- top_genes[keep_rows, ]
  
  rownames(Z) <- top_genes$display_name
  
  # Gene symbols in italics; inflammatory-response genes (GO:0006954)
  # in bold italics. Labels only, gene selection is unchanged.
  inflam <- is_inflammatory(top_genes$gene_id, top_genes$gene_name)
  
  row_labels <- as.expression(
    mapply(
      function(lbl, bold) {
        if (bold) bquote(bolditalic(.(lbl))) else bquote(italic(.(lbl)))
      },
      top_genes$display_name,
      inflam,
      SIMPLIFY = FALSE,
      USE.NAMES = FALSE
    )
  )
  
  annotation_col <- as.data.frame(
    colData(dds)[, "condition", drop = FALSE]
  )
  colnames(annotation_col) <- "Treatment"
  
  annotation_colours <- list(Treatment = COND_COLOURS)
  
  n_up <- sum(top_genes$direction == "UP")
  n_down <- sum(top_genes$direction == "DOWN")
  
  pheatmap(
    Z,
    main = paste0("Top ", nrow(Z), " DEGs in ", cell_type, ": LPS vs Saline"),
    labels_row = row_labels,
    annotation_col = annotation_col,
    annotation_colors = annotation_colours,
    scale = "none",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    fontsize = 10,
    filename = paste0(
      output_folder,
      "/",
      cell_type,
      "_Top20_Heatmap.png"
    ),
    width = 10,
    height = 9
  )
}


# ============================================================
# 23. MICROGLIA VISUALIZATIONS
# ============================================================

plot_PCA(microglia_dds, "Microglia", "Microglia_plots")

plot_volcano(microglia_dat, "Microglia", "Microglia_plots")

plot_DEG_bar(microglia_dat, "Microglia", "Microglia_plots")

plot_top5(microglia_dat, "Microglia", "Microglia_plots")

plot_normalized_counts(
  microglia_dds,
  microglia_dat,
  "Microglia",
  "Microglia_plots"
)

plot_heatmap(microglia_dds, microglia_dat, "Microglia", "Microglia_plots")


# ============================================================
# 24. ASTROCYTES VISUALIZATIONS
# ============================================================

plot_PCA(astrocytes_dds, "Astrocytes", "Astrocytes_plots")

plot_volcano(astrocytes_dat, "Astrocytes", "Astrocytes_plots")

plot_DEG_bar(astrocytes_dat, "Astrocytes", "Astrocytes_plots")

plot_top5(astrocytes_dat, "Astrocytes", "Astrocytes_plots")

plot_normalized_counts(
  astrocytes_dds,
  astrocytes_dat,
  "Astrocytes",
  "Astrocytes_plots"
)

plot_heatmap(astrocytes_dds, astrocytes_dat, "Astrocytes", "Astrocytes_plots")


# ============================================================
# 25. COLLECT ALL PNG FILES
# ============================================================

all_plots <- c(
  list.files(
    "Microglia_plots",
    pattern = "\\.png$",
    full.names = TRUE
  ),
  list.files(
    "Astrocytes_plots",
    pattern = "\\.png$",
    full.names = TRUE
  )
)

cat("\nNumber of plots:", length(all_plots), "\n")
print(basename(all_plots))


# ============================================================
# 26. SAVE ALL PLOTS INTO ONE ZIP FILE
# ============================================================

if (length(all_plots) > 0) {
  zip("All_RNAseq_Plots.zip", files = all_plots)
}


# ============================================================
# 27. COMBINE ALL PLOTS INTO ONE PDF
# ============================================================

if (length(all_plots) > 0) {
  
  pdf("All_RNAseq_Plots.pdf", width = 10, height = 8)
  
  for (plot_file in all_plots) {
    
    img <- png::readPNG(plot_file)
    
    grid.newpage()
    
    grid.text(
      tools::file_path_sans_ext(basename(plot_file)),
      y = unit(0.97, "npc"),
      gp = gpar(fontsize = 12, fontface = "bold")
    )
    
    grid.raster(
      img,
      y = unit(0.46, "npc"),
      width = unit(0.95, "npc"),
      height = unit(0.85, "npc"),
      interpolate = TRUE
    )
  }
  
  dev.off()
}


# ============================================================
# 28. CONFIRM OUTPUT
# ============================================================

cat("\n========================================\n")
cat("RNA-seq ANALYSIS COMPLETED\n")
cat("========================================\n")

cat("DEG criteria:\n")
cat("Adjusted p-value <", PADJ_THRESHOLD, "\n")
cat("|Log2 Fold Change| >=", LFC_THRESHOLD, "\n\n")

cat("Microglia significant DEGs:", nrow(microglia_DEGs), "\n")
cat("Astrocytes significant DEGs:", nrow(astrocytes_DEGs), "\n\n")

cat("Number of plots saved:", length(all_plots), "\n")
cat("Created: All_RNAseq_Plots.zip\n")
cat("Created: All_RNAseq_Plots.pdf\n")
