# ============================================================
# 1. LOAD PACKAGES
# ============================================================

library(tidyverse)
library(DESeq2)
library(pheatmap)


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


# Extract SRR accession numbers from sample names
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

# Remove genes with zero counts across ALL samples
counts <- counts[
  rowSums(counts[, sample_cols]) > 0,
]

# Check remaining genes
dim(counts)

# Check for remaining NA values
sum(is.na(counts[, sample_cols]))

# Check for zero-count genes
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

missing_samples <- setdiff(all_samples, colnames(counts))

if (length(missing_samples) > 0) {
  stop(
    paste(
      "Missing sample columns:",
      paste(missing_samples, collapse = ", ")
    )
  )
}

# Check for duplicate sample names
stopifnot(!anyDuplicated(colnames(counts)[3:ncol(counts)]))


# ============================================================
# 7. CREATE FOUR SEPARATE DATASETS
# ============================================================

microglia_LPS_counts <- counts[
  , c("gene", "length", microglia_LPS)
]

microglia_Saline_counts <- counts[
  , c("gene", "length", microglia_Saline)
]

astrocytes_LPS_counts <- counts[
  , c("gene", "length", astrocytes_LPS)
]

astrocytes_Saline_counts <- counts[
  , c("gene", "length", astrocytes_Saline)
]


# ============================================================
# 8. EXPORT FOUR CSV FILES
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
# 9. LOAD FOUR COUNT FILES
# ============================================================

microglia_LPS <- read.csv(
  "Microglia_LPS.csv",
  check.names = FALSE
)

microglia_Saline <- read.csv(
  "Microglia_Saline.csv",
  check.names = FALSE
)

astrocytes_LPS <- read.csv(
  "Astrocytes_LPS.csv",
  check.names = FALSE
)

astrocytes_Saline <- read.csv(
  "Astrocytes_Saline.csv",
  check.names = FALSE
)


# ============================================================
# 10. CREATE FUNCTION TO RUN DESEQ2
# ============================================================

run_deseq <- function(LPS, Saline, cell_type) {
  
  # Extract gene IDs
  genes <- LPS$gene
  
  # Check that gene IDs match between conditions
  stopifnot(identical(LPS$gene, Saline$gene))
  
  # Extract raw count data
  LPS_counts <- LPS[, -c(1, 2)]
  Saline_counts <- Saline[, -c(1, 2)]
  
  # Combine LPS and Saline samples
  count_data <- cbind(
    LPS_counts,
    Saline_counts
  )
  
  # Convert to matrix
  count_matrix <- as.matrix(count_data)
  
  # Check that counts are valid integers
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
  
  # Assign sample names
  rownames(sample_info) <- colnames(count_matrix)
  
  # Create DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = sample_info,
    design = ~ condition
  )
  
  # Filter genes with low counts
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
  
  # Convert results into a data frame
  res_df <- as.data.frame(res) |>
    rownames_to_column("gene") |>
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
  ) |>
    as.data.frame() |>
    rownames_to_column("gene")
  
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
  LPS = microglia_LPS,
  Saline = microglia_Saline,
  cell_type = "Microglia"
)


# ============================================================
# 12. RUN DESEQ2 FOR ASTROCYTES
# ============================================================

astrocytes_results <- run_deseq(
  LPS = astrocytes_LPS,
  Saline = astrocytes_Saline,
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
# 14. EXAMINE RESULTS
# ============================================================

head(microglia_dat)

head(astrocytes_dat)


# ============================================================
# 15. IDENTIFY SIGNIFICANT DEGS
# ============================================================

microglia_DEGs <- microglia_dat |>
  filter(
    !is.na(padj),
    padj < 0.05
  )

astrocytes_DEGs <- astrocytes_dat |>
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
# 16. COUNT UPREGULATED AND DOWNREGULATED GENES
# ============================================================

# Positive log2FoldChange = upregulated in LPS
# Negative log2FoldChange = downregulated in LPS

microglia_DEGs |>
  summarise(
    Upregulated = sum(log2FoldChange > 0),
    Downregulated = sum(log2FoldChange < 0)
  )

astrocytes_DEGs |>
  summarise(
    Upregulated = sum(log2FoldChange > 0),
    Downregulated = sum(log2FoldChange < 0)
  )

# ============================================================
# 17. LOAD VISUALIZATION PACKAGES
# ============================================================

library(ggplot2)
library(pheatmap)


# Extract DESeq2 objects and results

microglia_dds <- microglia_results$dds
astrocytes_dds <- astrocytes_results$dds

microglia_dat <- microglia_results$results
astrocytes_dat <- astrocytes_results$results


# Create output folders

dir.create(
  "Microglia_plots",
  showWarnings = FALSE
)

dir.create(
  "Astrocytes_plots",
  showWarnings = FALSE
)

plot_PCA <- function(dds, cell_type, output_folder) {
  
  # Transform counts
  rld <- rlog(dds, blind = FALSE)
  
  # Extract PCA data
  pcaData <- plotPCA(
    rld,
    intgroup = "condition",
    returnData = TRUE
  )
  
  # Calculate percentage of variance
  percentVar <- round(
    100 * attr(pcaData, "percentVar")
  )
  
  # Generate PCA plot
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
      ),
      plot.subtitle = element_text(
        size = 10
      )
    ) +
    labs(
      title = paste(
        "Global Gene Expression Profiles of",
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
  
  # Remove NA values
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
      
      # Prevent infinite values if adjusted p-value is zero
      neg_log10_padj = -log10(
        pmax(padj, .Machine$double.xmin)
      )
    )
  
  # Generate volcano plot
  p <- ggplot(
    labelled_dat,
    aes(
      x = log2FoldChange,
      y = neg_log10_padj,
      colour = up_down
    )
  ) +
    geom_point(
      alpha = 0.7,
      size = 1.5
    ) +
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed"
    ) +
    geom_vline(
      xintercept = c(-1, 1),
      linetype = "dashed"
    ) +
    theme_classic() +
    labs(
      title = paste(
        "Differential Gene Expression in",
        cell_type,
        "Following LPS Treatment"
      ),
      subtitle = "LPS vs Saline; dashed lines indicate adjusted p-value = 0.05 and log2FC = ±1",
      x = "Log2 Fold Change (LPS vs Saline)",
      y = "-Log10 Adjusted P-value",
      colour = "Expression Direction"
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
      "_Volcano.png"
    ),
    plot = p,
    width = 7,
    height = 6,
    dpi = 300
  )
}

plot_DEG_bar <- function(dat, cell_type, output_folder) {
  
  # Filter significant genes
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
  
  # Check whether significant genes exist
  if (nrow(labelled_dat) == 0) {
    message(
      paste("No significant DEGs found for", cell_type)
    )
    return(invisible(NULL))
  }
  
  # Set factor order
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
  
  # Generate bar plot
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

plot_top10 <- function(dat, cell_type, output_folder) {
  
  # Filter significant genes
  sig_dat <- dat %>%
    filter(
      !is.na(padj),
      !is.na(log2FoldChange),
      !is.na(lfcSE),
      padj < 0.05
    )
  
  # Select top 10 upregulated genes
  top10Genes <- sig_dat %>%
    filter(log2FoldChange > 0) %>%
    arrange(desc(log2FoldChange)) %>%
    slice_head(n = 10)
  
  # Select top 10 downregulated genes
  bot10Genes <- sig_dat %>%
    filter(log2FoldChange < 0) %>%
    arrange(log2FoldChange) %>%
    slice_head(n = 10)
  
  # Combine datasets
  joined_top10 <- bind_rows(
    top10Genes,
    bot10Genes
  ) %>%
    mutate(
      up_down = if_else(
        log2FoldChange > 0,
        "UP",
        "DOWN"
      )
    ) %>%
    arrange(log2FoldChange)
  
  # Check whether significant genes exist
  if (nrow(joined_top10) == 0) {
    message(
      paste("No significant DEGs found for", cell_type)
    )
    return(invisible(NULL))
  }
  
  # Set gene order
  joined_top10$gene <- factor(
    joined_top10$gene,
    levels = unique(joined_top10$gene)
  )
  
  # Generate plot
  p <- ggplot(
    joined_top10,
    aes(
      x = gene,
      y = log2FoldChange,
      fill = up_down
    )
  ) +
    geom_errorbar(
      aes(
        ymin = log2FoldChange - lfcSE,
        ymax = log2FoldChange + lfcSE
      ),
      width = 0.4
    ) +
    geom_col(
      colour = "black",
      width = 0.8
    ) +
    coord_flip() +
    geom_hline(
      yintercept = 0
    ) +
    theme_classic() +
    labs(
      title = paste(
        "Top 10 Upregulated and Downregulated Genes in",
        cell_type
      ),
      subtitle = "Significant DEGs following LPS treatment (adjusted p-value < 0.05)",
      x = "Gene ID",
      y = "Log2 Fold Change (LPS vs Saline)",
      fill = "Expression Direction"
    ) +
    theme(
      axis.text = element_text(
        colour = "black"
      ),
      axis.title = element_text(
        face = "bold"
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
      "_Top10_DEGs.png"
    ),
    plot = p,
    width = 8,
    height = 7,
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
  
  # Add condition information
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
  
  # Check whether significant genes exist
  if (length(gene_of_interest) == 0) {
    message(
      paste("No significant DEGs found for", cell_type)
    )
    return(invisible(NULL))
  }
  
  # Extract selected gene
  gene_data <- normalized_counts %>%
    filter(
      gene == gene_of_interest
    )
  
  # Generate normalized count plot
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
        gene_of_interest,
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
  
  # Apply variance stabilizing transformation
  vsd <- vst(
    dds,
    blind = FALSE
  )
  
  # Extract transformed counts
  vsd_matrix <- assay(vsd)
  
  # Calculate Z-scores for each gene
  Z <- t(
    scale(
      t(vsd_matrix)
    )
  )
  
  # Remove genes with NA or non-finite Z-scores
  Z_no_NA <- Z[
    apply(
      Z,
      1,
      function(x) all(is.finite(x))
    ),
    ,
    drop = FALSE
  ]
  
  # Check whether variable genes exist
  if (nrow(Z_no_NA) == 0) {
    message(
      paste("No suitable genes for heatmap in", cell_type)
    )
    return(invisible(NULL))
  }
  
  # Calculate variance of each gene
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
  
  # Create sample annotation
  annotation_col <- as.data.frame(
    colData(dds)[, "condition", drop = FALSE]
  )
  
  # Generate heatmap
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
# 24. MICROGLIA VISUALIZATIONS
# ============================================================

# PCA
plot_PCA(
  microglia_dds,
  "Microglia",
  "Microglia_plots"
)

# Volcano plot
plot_volcano(
  microglia_dat,
  "Microglia",
  "Microglia_plots"
)

# Significant DEG distribution
plot_DEG_bar(
  microglia_dat,
  "Microglia",
  "Microglia_plots"
)

# Top 10 upregulated and downregulated genes
plot_top10(
  microglia_dat,
  "Microglia",
  "Microglia_plots"
)

# Normalized counts
plot_normalized_counts(
  microglia_dds,
  microglia_dat,
  "Microglia",
  "Microglia_plots"
)

# Heatmap
plot_heatmap(
  microglia_dds,
  "Microglia",
  "Microglia_plots"
)


# ============================================================
# 25. ASTROCYTES VISUALIZATIONS
# ============================================================

# PCA
plot_PCA(
  astrocytes_dds,
  "Astrocytes",
  "Astrocytes_plots"
)

# Volcano plot
plot_volcano(
  astrocytes_dat,
  "Astrocytes",
  "Astrocytes_plots"
)

# Significant DEG distribution
plot_DEG_bar(
  astrocytes_dat,
  "Astrocytes",
  "Astrocytes_plots"
)

# Top 10 upregulated and downregulated genes
plot_top10(
  astrocytes_dat,
  "Astrocytes",
  "Astrocytes_plots"
)

# Normalized counts
plot_normalized_counts(
  astrocytes_dds,
  astrocytes_dat,
  "Astrocytes",
  "Astrocytes_plots"
)

# Heatmap
plot_heatmap(
  astrocytes_dds,
  "Astrocytes",
  "Astrocytes_plots"
)

# ============================================================
# 26. SAVE ALL PLOTS
# ============================================================

# Install packages if needed
if (!requireNamespace("png", quietly = TRUE)) {
  install.packages("png")
}

if (!requireNamespace("gridExtra", quietly = TRUE)) {
  install.packages("gridExtra")
}

if (!requireNamespace("grid", quietly = TRUE)) {
  install.packages("grid")
}

library(png)
library(gridExtra)
library(grid)


# ============================================================
# 27. COLLECT ALL PNG FILES
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

# Check number of saved plots
length(all_plots)

# Display plot filenames
basename(all_plots)


# ============================================================
# 28. SAVE ALL PLOTS INTO ONE ZIP FILE
# ============================================================

zip(
  "All_RNAseq_Plots.zip",
  files = all_plots
)


# ============================================================
# 29. COMBINE ALL PLOTS INTO ONE PDF
# ============================================================

pdf(
  "All_RNAseq_Plots.pdf",
  width = 10,
  height = 8
)

for (plot_file in all_plots) {
  
  # Read PNG image
  img <- png::readPNG(plot_file)
  
  # Display image on PDF page
  grid.newpage()
  
  grid.raster(
    img,
    width = unit(0.95, "npc"),
    height = unit(0.90, "npc"),
    interpolate = TRUE
  )
  
  # Add filename as page title
  grid.text(
    basename(plot_file),
    y = unit(0.97, "npc"),
    gp = gpar(
      fontsize = 12,
      fontface = "bold"
    )
  )
}

dev.off()


# ============================================================
# 30. CONFIRM OUTPUT
# ============================================================

cat(
  "Number of plots saved:",
  length(all_plots),
  "\n"
)

cat("Created: All_RNAseq_Plots.zip\n")
cat("Created: All_RNAseq_Plots.pdf\n")

