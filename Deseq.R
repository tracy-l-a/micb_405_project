# ============================================================
# 1. LOAD PACKAGES
# ============================================================

library(tidyverse)
library(DESeq2)
library(pheatmap)
library(ggplot2)
library(png)
library(grid)
library(gridExtra)


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


# ============================================================
# 5. FILTER COUNT DATA
# ============================================================

# Identify sample columns
sample_cols <- colnames(counts)[3:ncol(counts)]

# Convert sample counts to numeric
counts[sample_cols] <- lapply(
  counts[sample_cols],
  as.numeric
)

# Remove genes containing NA values
counts <- counts[
  complete.cases(counts[, sample_cols]),
]

# Remove genes with zero counts across all samples
counts <- counts[
  rowSums(counts[, sample_cols]) > 0,
]

# Check remaining genes
dim(counts)

# Check for NA values
sum(is.na(counts[, sample_cols]))

# Check for all-zero genes
sum(rowSums(counts[, sample_cols]) == 0)


# ============================================================
# 6. CHECK SAMPLE NAMES
# ============================================================

all_samples <- c(
  microglia_LPS,
  microglia_Saline,
  astrocytes_LPS,
  astrocytes_Saline
)

missing_samples <- setdiff(
  all_samples,
  colnames(counts)
)

if (length(missing_samples) > 0) {
  stop(
    paste(
      "Missing sample columns:",
      paste(missing_samples, collapse = ", ")
    )
  )
}

# Check for duplicate sample names
stopifnot(
  !anyDuplicated(colnames(counts)[3:ncol(counts)])
)


# ============================================================
# 7. ANNOTATE GENES USING GENCODE M39
# ============================================================

# Load GENCODE M39 GTF file
gtf <- read.delim(
  "gencode.vM39.primary_assembly.annotation.gtf",
  header = FALSE,
  sep = "\t",
  comment.char = "#",
  quote = "",
  stringsAsFactors = FALSE
)

# Assign column names
colnames(gtf) <- c(
  "chromosome",
  "source",
  "feature",
  "start",
  "end",
  "score",
  "strand",
  "frame",
  "attributes"
)

# Extract gene-level annotations
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
  mutate(
    gene_id = sub("\\.[0-9]+$", "", gene_id)
  ) %>%
  distinct(gene_id, .keep_all = TRUE)

# Add version-free Ensembl IDs to count data
counts <- counts %>%
  mutate(
    gene_id = sub("\\.[0-9]+$", "", gene)
  )

# Check gene annotation matching
matched_genes <- sum(
  counts$gene_id %in% gene_annotation$gene_id
)

cat(
  "Genes matched:",
  matched_genes,
  "out of",
  nrow(counts),
  "\n"
)

cat(
  "Percentage matched:",
  round(
    100 * matched_genes / nrow(counts),
    2
  ),
  "%\n"
)

# Add gene symbols to count data
counts <- counts %>%
  left_join(
    gene_annotation,
    by = "gene_id"
  ) %>%
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
  select(
    gene,
    gene_id,
    gene_name,
    length,
    all_of(microglia_LPS)
  )

microglia_Saline_counts <- counts %>%
  select(
    gene,
    gene_id,
    gene_name,
    length,
    all_of(microglia_Saline)
  )

astrocytes_LPS_counts <- counts %>%
  select(
    gene,
    gene_id,
    gene_name,
    length,
    all_of(astrocytes_LPS)
  )

astrocytes_Saline_counts <- counts %>%
  select(
    gene,
    gene_id,
    gene_name,
    length,
    all_of(astrocytes_Saline)
  )


# ============================================================
# 9. EXPORT FOUR CSV FILES
# ============================================================

write.csv(
  microglia_LPS_counts,
  "Microglia_LPS.csv",
  row.names = FALSE
)

write.csv(
  microglia_Saline_counts,
  "Microglia_Saline.csv",
  row.names = FALSE
)

write.csv(
  astrocytes_LPS_counts,
  "Astrocytes_LPS.csv",
  row.names = FALSE
)

write.csv(
  astrocytes_Saline_counts,
  "Astrocytes_Saline.csv",
  row.names = FALSE
)

# ============================================================
# 10. CREATE FUNCTION TO RUN DESEQ2
# ============================================================

run_deseq <- function(LPS, Saline, cell_type) {
  
  # Extract gene IDs
  genes <- LPS$gene
  
  # Check that gene IDs match between conditions
  stopifnot(
    identical(LPS$gene, Saline$gene)
  )
  
  # Extract sample columns by excluding annotation columns
  LPS_counts <- LPS %>%
    select(-gene, -gene_id, -gene_name, -length)
  
  Saline_counts <- Saline %>%
    select(-gene, -gene_id, -gene_name, -length)
  
  # Combine LPS and Saline counts
  count_matrix <- as.matrix(
    cbind(LPS_counts, Saline_counts)
  )
  
  # Check counts
  stopifnot(
    all(is.finite(count_matrix)),
    all(count_matrix >= 0),
    all(count_matrix == round(count_matrix))
  )
  
  # Convert to integer
  storage.mode(count_matrix) <- "integer"
  
  # Assign gene IDs as row names
  rownames(count_matrix) <- genes
  
  # Create sample information
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
  
  # Create DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = sample_info,
    design = ~ condition
  )
  
  # Filter low-count genes
  keep <- rowSums(counts(dds) >= 10) >= 2
  
  dds <- dds[keep, ]
  
  # Run DESeq2
  dds <- DESeq(dds)
  
  # Extract differential expression results
  res <- results(
    dds,
    contrast = c("condition", "LPS", "Saline"),
    alpha = 0.05
  )
  
  # Convert results to data frame
  res_df <- as.data.frame(res) %>%
    rownames_to_column("gene") %>%
    mutate(
      gene_id = sub("\\.[0-9]+$", "", gene)
    ) %>%
    left_join(
      gene_annotation,
      by = "gene_id"
    ) %>%
    mutate(
      gene_name = if_else(
        is.na(gene_name) | gene_name == "",
        gene,
        gene_name
      )
    ) %>%
    arrange(padj)
  
  # Export DESeq2 results
  write.csv(
    res_df,
    paste0(cell_type, "_DESeq2_results.csv"),
    row.names = FALSE
  )
  
  # Extract normalized counts
  normalized_counts <- counts(
    dds,
    normalized = TRUE
  ) %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    mutate(
      gene_id = sub("\\.[0-9]+$", "", gene)
    ) %>%
    left_join(
      gene_annotation,
      by = "gene_id"
    ) %>%
    mutate(
      gene_name = if_else(
        is.na(gene_name) | gene_name == "",
        gene,
        gene_name
      )
    )
  
  # Export normalized counts
  write.csv(
    normalized_counts,
    paste0(cell_type, "_normalized_counts.csv"),
    row.names = FALSE
  )
  
  # Return results
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

microglia_DEGs <- microglia_dat %>%
  filter(
    !is.na(padj),
    padj < 0.05
  )

astrocytes_DEGs <- astrocytes_dat %>%
  filter(
    !is.na(padj),
    padj < 0.05
  )

# Export significant DEGs
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


# ============================================================
# 15. COUNT UPREGULATED AND DOWNREGULATED GENES
# ============================================================

microglia_DEGs %>%
  summarise(
    Upregulated = sum(log2FoldChange > 0),
    Downregulated = sum(log2FoldChange < 0)
  )

astrocytes_DEGs %>%
  summarise(
    Upregulated = sum(log2FoldChange > 0),
    Downregulated = sum(log2FoldChange < 0)
  )

# ============================================================
# 16. CREATE OUTPUT FOLDERS
# ============================================================

dir.create(
  "Microglia_plots",
  showWarnings = FALSE
)

dir.create(
  "Astrocytes_plots",
  showWarnings = FALSE
)

plot_PCA <- function(dds, cell_type, output_folder) {
  
  rld <- rlog(
    dds,
    blind = FALSE
  )
  
  pcaData <- plotPCA(
    rld,
    intgroup = "condition",
    returnData = TRUE
  )
  
  percentVar <- round(
    100 * attr(pcaData, "percentVar")
  )
  
  p <- ggplot(
    pcaData,
    aes(
      x = PC1,
      y = PC2,
      colour = condition
    )
  ) +
    geom_point(size = 4) +
    xlab(
      paste0(
        "PC1: ",
        percentVar[1],
        "% variance"
      )
    ) +
    ylab(
      paste0(
        "PC2: ",
        percentVar[2],
        "% variance"
      )
    ) +
    theme_test() +
    theme(
      axis.text = element_text(
        colour = "black",
        size = 11
      ),
      axis.title = element_text(
        face = "bold"
      ),
      legend.title = element_text(
        face = "bold"
      ),
      panel.border = element_rect(
        linewidth = 1.2
      ),
      plot.title = element_text(
        face = "bold",
        size = 13
      )
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
    paste0(
      output_folder,
      "/",
      cell_type,
      "_PCA.png"
    ),
    plot = p,
    width = 7,
    height = 6,
    dpi = 300
  )
}






plot_volcano <- function(dat, cell_type, output_folder) {
  
  # Prepare data
  labelled_dat <- dat %>%
    filter(
      !is.na(padj),
      !is.na(log2FoldChange)
    ) %>%
    mutate(
      up_down = case_when(
        padj < 0.05 & log2FoldChange > 1 ~ "UP",
        padj < 0.05 & log2FoldChange < -1 ~ "DOWN",
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
  
  # Select top 5 upregulated genes
  top_up <- labelled_dat %>%
    filter(up_down == "UP") %>%
    arrange(desc(log2FoldChange), padj) %>%
    slice_head(n = 5)
  
  # Select top 5 downregulated genes
  top_down <- labelled_dat %>%
    filter(up_down == "DOWN") %>%
    arrange(log2FoldChange, padj) %>%
    slice_head(n = 5)
  
  # Combine selected genes
  label_dat <- bind_rows(
    top_up,
    top_down
  )
  
  # Distinguish duplicate gene symbols
  duplicate_labels <- duplicated(label_dat$gene_label) |
    duplicated(label_dat$gene_label, fromLast = TRUE)
  
  label_dat$gene_label[duplicate_labels] <- paste0(
    label_dat$gene_label[duplicate_labels],
    " (",
    label_dat$gene_id[duplicate_labels],
    ")"
  )
  
  # Calculate plot boundaries
  x_range <- range(labelled_dat$log2FoldChange, na.rm = TRUE)
  y_range <- range(labelled_dat$neg_log10_padj, na.rm = TRUE)
  
  x_padding <- diff(x_range) * 0.15
  y_padding <- diff(y_range) * 0.15
  
  # Avoid zero-width limits
  if (x_padding == 0) x_padding <- 1
  if (y_padding == 0) y_padding <- 1
  
  # Create volcano plot
  p <- ggplot(
    labelled_dat,
    aes(
      x = log2FoldChange,
      y = neg_log10_padj,
      colour = up_down
    )
  ) +
    
    # Plot all genes
    geom_point(
      alpha = 0.7,
      size = 1.6
    ) +
    
    # Label top 5 upregulated and downregulated genes
    ggrepel::geom_label_repel(
      data = label_dat,
      aes(label = gene_label),
      inherit.aes = TRUE,
      colour = "black",
      fill = "white",
      size = 3,
      fontface = "italic",
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
    
    # Significance threshold
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed",
      linewidth = 0.5,
      colour = "grey40"
    ) +
    
    # Fold-change thresholds
    geom_vline(
      xintercept = c(-1, 1),
      linetype = "dashed",
      linewidth = 0.5,
      colour = "grey40"
    ) +
    
    # Give labels additional room around plot boundaries
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
    
    theme_classic() +
    
    labs(
      title = paste(
        "Differential Gene Expression in",
        cell_type,
        "Following LPS Treatment"
      ),
      subtitle = "Top 5 upregulated and downregulated significant genes are labelled",
      x = "Log2 Fold Change (LPS vs Saline)",
      y = "-Log10 Adjusted P-value",
      colour = "Expression Direction"
    ) +
    
    theme(
      axis.title = element_text(
        face = "bold",
        size = 11
      ),
      axis.text = element_text(
        colour = "black",
        size = 10
      ),
      plot.title = element_text(
        face = "bold",
        size = 13
      ),
      plot.subtitle = element_text(
        size = 9
      ),
      legend.title = element_text(
        face = "bold"
      ),
      plot.margin = margin(
        15, 25, 15, 15
      )
    )
  
  print(p)
  
  # Save high-resolution plot
  ggsave(
    paste0(
      output_folder,
      "/",
      cell_type,
      "_Volcano.png"
    ),
    plot = p,
    width = 10,
    height = 8,
    dpi = 300,
    bg = "white"
  )
}



plot_DEG_bar <- function(dat, cell_type, output_folder) {
  
  labelled_dat <- dat %>%
    filter(
      !is.na(padj),
      !is.na(log2FoldChange),
      padj < 0.05
    ) %>%
    mutate(
      abs_l2fc = abs(log2FoldChange),
      
      up_down = if_else(
        log2FoldChange > 0,
        "UP",
        "DOWN"
      ),
      
      l2fc = case_when(
        abs_l2fc >= 4 ~ ">=4",
        abs_l2fc >= 3 ~ "3-<4",
        abs_l2fc >= 2 ~ "2-<3",
        abs_l2fc >= 1 ~ "1-<2",
        TRUE ~ "<1"
      )
    )
  
  if (nrow(labelled_dat) == 0) {
    message(
      paste("No significant DEGs found for", cell_type)
    )
    return(invisible(NULL))
  }
  
  labelled_dat$l2fc <- factor(
    labelled_dat$l2fc,
    levels = c(
      "<1",
      "1-<2",
      "2-<3",
      "3-<4",
      ">=4"
    )
  )
  
  p <- ggplot(
    labelled_dat,
    aes(
      x = l2fc,
      fill = up_down
    )
  ) +
    geom_bar(
      colour = "black"
    ) +
    theme_classic() +
    labs(
      title = paste(
        "Distribution of Significant Gene Expression Changes in",
        cell_type
      ),
      subtitle = "Significant DEGs (adjusted p-value < 0.05), grouped by absolute log2 fold change",
      x = "Absolute Log2 Fold Change",
      y = "Number of Significant Genes",
      fill = "Expression Direction"
    ) +
    theme(
      axis.title = element_text(
        face = "bold"
      ),
      axis.text = element_text(
        colour = "black"
      ),
      plot.title = element_text(
        face = "bold",
        size = 13
      ),
      plot.subtitle = element_text(
        size = 9
      )
    )
  
  print(p)
  
  ggsave(
    paste0(
      output_folder,
      "/",
      cell_type,
      "_DEG_barplot.png"
    ),
    plot = p,
    width = 7,
    height = 6,
    dpi = 300
  )
}

plot_top5 <- function(dat, cell_type, output_folder) {
  # Keep significant genes with valid fold-change estimates
  sig_dat <- dat %>%
    filter(
      !is.na(padj),
      !is.na(log2FoldChange),
      !is.na(lfcSE),
      padj < 0.05
    ) %>%
    mutate(
      # Use the annotated gene symbol; fall back to Ensembl ID if unavailable
      display_name = if_else(
        is.na(gene_name) | gene_name == "",
        gene,
        gene_name
      )
    )
  
  # Select the five most upregulated and five most downregulated genes
  top5_up <- sig_dat %>%
    filter(log2FoldChange > 0) %>%
    arrange(desc(log2FoldChange)) %>%
    slice_head(n = 5)
  
  top5_down <- sig_dat %>%
    filter(log2FoldChange < 0) %>%
    arrange(log2FoldChange) %>%
    slice_head(n = 5)
  
  top_genes <- bind_rows(top5_up, top5_down) %>%
    mutate(
      direction = if_else(log2FoldChange > 0, "Upregulated", "Downregulated")
    )
  
  if (nrow(top_genes) == 0) {
    message(paste("No significant DEGs found for", cell_type))
    return(invisible(NULL))
  }
  
  # If multiple Ensembl genes share a symbol, append the Ensembl ID
  duplicate_names <- duplicated(top_genes$display_name) |
    duplicated(top_genes$display_name, fromLast = TRUE)
  
  top_genes$display_name[duplicate_names] <- paste0(
    top_genes$display_name[duplicate_names],
    " (",
    top_genes$gene_id[duplicate_names],
    ")"
  )
  
  # Order bars by log2 fold change
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
    theme_classic() +
    labs(
      title = paste("Top 5 Upregulated and Downregulated Genes in", cell_type),
      subtitle = "Top significant DEGs ranked by log2 fold change (adjusted p-value < 0.05)",
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

plot_normalized_counts <- function(
    dds,
    dat,
    cell_type,
    output_folder
) {
  
  # Extract normalized counts
  normalized_counts <- counts(
    dds,
    normalized = TRUE
  ) %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    pivot_longer(
      cols = -gene,
      names_to = "Sample",
      values_to = "count"
    )
  
  # Extract sample information
  sample_info <- as.data.frame(
    colData(dds)
  ) %>%
    rownames_to_column("Sample")
  
  # Add treatment information
  normalized_counts <- normalized_counts %>%
    left_join(
      sample_info,
      by = "Sample"
    )
  
  # Select the most significant DEG
  gene_of_interest <- dat %>%
    filter(
      !is.na(padj),
      !is.na(log2FoldChange),
      padj < 0.05
    ) %>%
    arrange(padj) %>%
    slice_head(n = 1) %>%
    pull(gene)
  
  if (length(gene_of_interest) == 0) {
    message(
      paste("No significant DEGs found for", cell_type)
    )
    return(invisible(NULL))
  }
  
  # Retrieve actual gene symbol without using first(), which can be
  # masked by another loaded package
  gene_symbol <- dat %>%
    filter(gene == gene_of_interest) %>%
    pull(gene_name)
  
  gene_symbol <- gene_symbol[1]
  
  if (length(gene_symbol) == 0 || is.na(gene_symbol) || gene_symbol == "") {
    gene_symbol <- gene_of_interest
  }
  
  # Extract selected gene expression
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
    geom_jitter(
      width = 0.2,
      size = 3
    ) +
    scale_y_continuous(
      trans = "log1p"
    ) +
    theme_classic() +
    labs(
      title = paste(
        "Normalized Expression of",
        gene_symbol,
        "in",
        cell_type
      ),
      subtitle = "Comparison between LPS-treated and Saline-treated biological replicates",
      x = "Treatment",
      y = "Normalized Counts (log1p)",
      colour = "Treatment"
    ) +
    theme(
      axis.title = element_text(
        face = "bold"
      ),
      axis.text = element_text(
        colour = "black"
      ),
      plot.title = element_text(
        face = "bold",
        size = 13
      ),
      plot.subtitle = element_text(
        size = 9
      )
    )
  
  print(p)
  
  ggsave(
    paste0(
      output_folder,
      "/",
      cell_type,
      "_Normalized_Counts.png"
    ),
    plot = p,
    width = 7,
    height = 6,
    dpi = 300
  )
}

plot_heatmap <- function(dds, cell_type, output_folder) {
  
  # Variance stabilizing transformation
  vsd <- vst(
    dds,
    blind = FALSE
  )
  
  # Extract transformed counts
  vsd_matrix <- assay(vsd)
  
  # Calculate gene-wise Z-scores
  Z <- t(
    scale(
      t(vsd_matrix)
    )
  )
  
  # Remove genes with non-finite Z-scores
  Z_no_NA <- Z[
    apply(
      Z,
      1,
      function(x) all(is.finite(x))
    ),
    ,
    drop = FALSE
  ]
  
  if (nrow(Z_no_NA) == 0) {
    message(
      paste("No suitable genes for heatmap in", cell_type)
    )
    return(invisible(NULL))
  }
  
  # Calculate gene variance
  gene_variance <- apply(
    vsd_matrix[
      rownames(Z_no_NA),
      ,
      drop = FALSE
    ],
    1,
    var
  )
  
  # Select top 20 variable genes
  topVarGenes <- names(
    sort(
      gene_variance,
      decreasing = TRUE
    )
  ) %>%
    head(20)
  
  # Extract Z-scores
  Z_topGenes <- Z_no_NA[
    topVarGenes,
    ,
    drop = FALSE
  ]
  
  # Match Ensembl IDs to gene symbols
  heatmap_gene_names <- gene_annotation$gene_name[
    match(
      sub(
        "\\.[0-9]+$",
        "",
        rownames(Z_topGenes)
      ),
      gene_annotation$gene_id
    )
  ]
  
  # Use Ensembl ID if gene symbol is unavailable
  missing_names <- is.na(heatmap_gene_names) |
    heatmap_gene_names == ""
  
  heatmap_gene_names[missing_names] <- rownames(
    Z_topGenes
  )[missing_names]
  
  # Distinguish duplicate gene symbols
  duplicate_names <- duplicated(heatmap_gene_names) |
    duplicated(
      heatmap_gene_names,
      fromLast = TRUE
    )
  
  heatmap_gene_names[duplicate_names] <- paste0(
    heatmap_gene_names[duplicate_names],
    " (",
    rownames(Z_topGenes)[duplicate_names],
    ")"
  )
  
  # Replace heatmap row names
  rownames(Z_topGenes) <- heatmap_gene_names
  
  # Create sample annotation
  annotation_col <- as.data.frame(
    colData(dds)[, "condition", drop = FALSE]
  )
  
  # Generate and save heatmap
  pheatmap(
    Z_topGenes,
    main = paste(
      "Expression Patterns of the Top 20 Most Variable Genes in",
      cell_type
    ),
    annotation_col = annotation_col,
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

plot_PCA(
  microglia_dds,
  "Microglia",
  "Microglia_plots"
)

plot_volcano(
  microglia_dat,
  "Microglia",
  "Microglia_plots"
)

plot_DEG_bar(
  microglia_dat,
  "Microglia",
  "Microglia_plots"
)

plot_top5(
  microglia_dat,
  "Microglia",
  "Microglia_plots"
)

plot_normalized_counts(
  microglia_dds,
  microglia_dat,
  "Microglia",
  "Microglia_plots"
)

plot_heatmap(
  microglia_dds,
  "Microglia",
  "Microglia_plots"
)


# ============================================================
# 24. ASTROCYTES VISUALIZATIONS
# ============================================================

plot_PCA(
  astrocytes_dds,
  "Astrocytes",
  "Astrocytes_plots"
)

plot_volcano(
  astrocytes_dat,
  "Astrocytes",
  "Astrocytes_plots"
)

plot_DEG_bar(
  astrocytes_dat,
  "Astrocytes",
  "Astrocytes_plots"
)

plot_top5(
  astrocytes_dat,
  "Astrocytes",
  "Astrocytes_plots"
)

plot_normalized_counts(
  astrocytes_dds,
  astrocytes_dat,
  "Astrocytes",
  "Astrocytes_plots"
)

plot_heatmap(
  astrocytes_dds,
  "Astrocytes",
  "Astrocytes_plots"
)

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

# Check number of plots
print(length(all_plots))

# Display filenames
print(basename(all_plots))


# ============================================================
# 26. SAVE ALL PLOTS INTO ONE ZIP FILE
# ============================================================

if (length(all_plots) > 0) {
  
  zip(
    "All_RNAseq_Plots.zip",
    files = all_plots
  )
  
}


# ============================================================
# 27. COMBINE ALL PLOTS INTO ONE PDF
# ============================================================

if (length(all_plots) > 0) {
  
  pdf(
    "All_RNAseq_Plots.pdf",
    width = 10,
    height = 8
  )
  
  for (plot_file in all_plots) {
    
    img <- png::readPNG(plot_file)
    
    grid.newpage()
    
    grid.text(
      tools::file_path_sans_ext(
        basename(plot_file)
      ),
      y = unit(0.97, "npc"),
      gp = gpar(
        fontsize = 12,
        fontface = "bold"
      )
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

cat(
  "\nNumber of plots saved:",
  length(all_plots),
  "\n"
)

cat(
  "Created: All_RNAseq_Plots.zip\n"
)

cat(
  "Created: All_RNAseq_Plots.pdf\n"
)

