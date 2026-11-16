# ==================================================================================================
# 0. Load Required R Packages
# ==================================================================================================

# Suppress startup messages so the console stays clean while loading libraries
suppressPackageStartupMessages({
  library(tidyverse)  # Data manipulation (dplyr, tibble) and plotting (ggplot2)
  library(DESeq2)     # RNA-seq count normalization and variance-stabilizing transformation (vst)
  library(pheatmap)   # Generates clustered heatmaps with sample/gene annotation bars
  library(ggrepel)    # Prevents SRR sample text labels from overlapping on the PCA plot
  library(grid)       # Low-level graphics unit sizing
  library(gridExtra)  # Wraps heatmap gtables with uniform outer padding so edges never clip
})

# ==================================================================================================
# 1. Load & Clean Raw Counts
# ==================================================================================================

# Read the raw featureCounts table, ignoring metadata comment lines starting with '#'
counts <- read.delim("gene_counts.tsv", comment.char = "#", header = TRUE, check.names = FALSE)

# Drop columns 2 through 5 (Chr, Start, End, Strand), keeping Geneid, Length, and sample columns
counts <- counts[, -c(2:5)]

# Rename the first two columns to clean, lowercase names
colnames(counts)[1:2] <- c("gene", "length")

# Extract just the SRR accession ID (e.g., "SRR2938958") from long BAM file path column names
colnames(counts)[3:ncol(counts)] <- sub(".*(SRR[0-9]+).*", "\\1", colnames(counts)[3:ncol(counts)])

# Define the SRR sample IDs for each experimental group (from Srinivasan et al., 2016)
microglia_LPS     <- c("SRR2938958", "SRR2938964", "SRR2938974", "SRR2938976", "SRR2938981") # n = 5
microglia_Saline  <- c("SRR2938960", "SRR2938961", "SRR2938967", "SRR2938971", "SRR2938980") # n = 5
astrocytes_LPS    <- c("SRR2938963", "SRR2938973", "SRR2938977", "SRR2938985")               # n = 4 (1 contaminated excluded in paper)
astrocytes_Saline <- c("SRR2938962", "SRR2938968", "SRR2938969", "SRR2938972")               # n = 4 (1 contaminated excluded in paper)

# Combine all 18 sample IDs into a single ordered vector (Microglia LPS -> Microglia Saline -> Astrocytes LPS -> Astrocytes Saline)
all_samples       <- c(microglia_LPS, microglia_Saline, astrocytes_LPS, astrocytes_Saline)

# Ensure all count columns are stored as numeric values rather than text/characters
counts[all_samples] <- lapply(counts[all_samples], as.numeric)

# Keep only genes with no missing (NA) values AND at least 1 total read across the 18 samples
counts <- counts[complete.cases(counts[, all_samples]) & rowSums(counts[, all_samples]) > 0, ]

# ==================================================================================================
# 2. Map Ensembl Gene IDs to Gene Symbols Using GENCODE M39 GTF
# ==================================================================================================

# Read the GENCODE M39 mouse reference annotation file (tab-separated, skipping '#' header lines)
gtf <- read.delim("gencode.vM39.primary_assembly.annotation.gtf", header = FALSE, sep = "\t",
                  comment.char = "#", quote = "", stringsAsFactors = FALSE)

# Assign standard GTF column names to the 9 columns
colnames(gtf) <- c("chromosome", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")

# Build a clean 2-column lookup table (Ensembl gene_id -> readable gene_name like "Itgam" or "Gfap")
gene_annotation <- gtf %>%
  filter(feature == "gene") %>%                                                                # Keep only whole-gene rows (ignore exon/transcript rows)
  mutate(
    gene_id = sub("\\.[0-9]+$", "", stringr::str_match(attributes, 'gene_id "([^"]+)"')[, 2]), # Extract Ensembl ID and strip version suffix (e.g., ".15")
    gene_name = stringr::str_match(attributes, 'gene_name "([^"]+)"')[, 2]                     # Extract official mouse gene symbol (e.g., "Itgam")
  ) %>%
  filter(!is.na(gene_id)) %>%                                                                  # Remove any rows missing a valid Ensembl ID
  dplyr::select(gene_id, gene_name) %>%                                                        # Keep only the ID and symbol columns
  distinct(gene_id, .keep_all = TRUE)                                                          # Remove duplicate entries so each gene_id appears once

# ==================================================================================================
# 3. Build Combined 18-Sample DESeqDataSet & Run Variance-Stabilizing Transformation (VST)
# ==================================================================================================

# Extract the raw count columns for all 18 samples into a numeric matrix
all_mat <- as.matrix(counts[, all_samples])

# Convert matrix values to whole integers (required by DESeq2)
storage.mode(all_mat) <- "integer"

# Set the matrix row names to the Ensembl gene IDs
rownames(all_mat) <- counts$gene

# Create sample metadata table matching the exact column order of all_mat
all_meta <- data.frame(
  cell_type = factor(c(rep("Microglia", 10), rep("Astrocytes", 8))),                            # Cell type label for each sample
  condition = factor(c(rep("LPS", 5), rep("Saline", 5), rep("LPS", 4), rep("Saline", 4))),      # Treatment label for each sample
  group     = factor(c(rep("Microglia_LPS", 5), rep("Microglia_Saline", 5),
                       rep("Astrocytes_LPS", 4), rep("Astrocytes_Saline", 4))),                 # Combined cell_type + treatment group
  row.names = all_samples                                                                       # Match metadata row names to sample column names
)

# Create a unified DESeq2 object containing all 18 samples together
dds_all <- DESeqDataSetFromMatrix(countData = all_mat, colData = all_meta, design = ~ group)

# Pre-filter low-count noise: keep genes with >= 10 reads in at least 2 samples
dds_all <- dds_all[rowSums(counts(dds_all) >= 10) >= 2, ]

# Calculate library size factors to normalize for sequencing depth differences across samples
dds_all <- estimateSizeFactors(dds_all)

# Apply Variance-Stabilizing Transformation (VST) so high-count genes don't dominate PCA/distances
# blind = TRUE ensures the transformation is unbiased by sample group labels (best practice for QC)
vsd_all <- vst(dds_all, blind = TRUE)

# Create an output folder in your working directory to store the QC figures
dir.create("Sample_QC_plots", showWarnings = FALSE)

# ==================================================================================================
# 4. Define Unified Colorblind-Safe Palette & Clean Annotation Bars
# ==================================================================================================

# Define consistent colors for treatment groups, cell types, and marker categories across all plots
treat_colors  <- c("LPS" = "#D55E00", "Saline" = "#0072B2")        # Colorblind-safe Vermilion & Blue
cell_colors   <- c("Microglia" = "#E69F00", "Astrocytes" = "#009E73") # Amber & Bluish-Green
marker_colors <- c(
  "Microglia"          = "#E69F00",
  "Astrocytes"         = "#009E73",
  "Neurons (Neg Ctrl)" = "#7F7F7F"
)

# Combine palettes into a named list formatted for pheatmap's annotation_colors argument
annot_colors <- list(
  Treatment   = treat_colors,
  Cell_Type   = cell_colors,
  Marker_Type = marker_colors
)

# Create a clean metadata table with polished column names for heatmap annotation legends
clean_annot <- data.frame(
  Cell_Type = all_meta$cell_type,
  Treatment = all_meta$condition,
  row.names = rownames(all_meta)
)

# Define one shared blue-yellow-red color gradient for both heatmaps
shared_colors <- colorRampPalette(c("#313695", "#74ADD1", "#FFFFBF", "#F46D43", "#A50026"))(100)

# ==================================================================================================
# 5. QC Plot 1: Variance-Stabilized PCA
# ==================================================================================================

# Compute PCA coordinates from the top 500 most variable VST-transformed genes
pca_all <- plotPCA(vsd_all, intgroup = c("cell_type", "condition"), returnData = TRUE)

# Extract the percentage of total variance explained by PC1 and PC2 for axis labels
pct_var <- round(100 * attr(pca_all, "percentVar"))

# Flag SRR2938964 in its plot label as the weak-responder LPS sample identified in QC
pca_all <- pca_all %>%
  mutate(sample_label = if_else(name == "SRR2938964", paste0(name, " (Weak Responder)"), name))

# Build the polished PCA scatter plot using ggplot2
p_pca <- ggplot(pca_all, aes(x = PC1, y = PC2, colour = condition, shape = cell_type)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey85") +                        # Subtle horizontal reference line at y = 0
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey85") +                        # Subtle vertical reference line at x = 0
  geom_point(size = 4.5, alpha = 0.95) +                                                      # Plot samples colored by treatment and shaped by cell type
  geom_label_repel(
    aes(label = sample_label),
    size = 3.0,
    fontface = "bold",
    fill = alpha("white", 0.9),                                                               # Semi-transparent white pill background behind text
    box.padding = 0.75,                                                                       # Push label boxes further away from points and each other
    point.padding = 0.55,                                                                     # Prevent boxes from overlapping neighboring dots
    force = 4,                                                                                # Stronger repulsion force for tight replicate clusters
    force_pull = 0.5,
    min.segment.length = 0,                                                                   # Force EVERY label to have a connecting pointer line
    segment.color = "grey45",                                                                 # Thin pointer line connecting label pill to sample dot
    segment.size = 0.4,
    seed = 42,                                                                                # Reproducible label positions every time you run the script
    show.legend = FALSE,
    max.overlaps = Inf
  ) +
  scale_colour_manual(values = treat_colors) +                                                # Apply unified colorblind-safe treatment palette
  scale_shape_manual(values = c("Microglia" = 17, "Astrocytes" = 16)) +                       # Triangles for Microglia, circles for Astrocytes
  scale_x_continuous(expand = expansion(mult = 0.14)) +                                       # Extra left/right breathing room for labels
  scale_y_continuous(expand = expansion(mult = 0.14)) +                                       # Extra top/bottom breathing room for labels
  xlab(paste0("PC1: ", pct_var[1], "% variance (Cell Type Identity)")) +                      # Label x-axis with PC1 % variance and biological driver
  ylab(paste0("PC2: ", pct_var[2], "% variance (LPS Inflammatory Response)")) +               # Label y-axis with PC2 % variance and biological driver
  theme_classic(base_size = 12) +                                                             # Clean white background with crisp axis lines
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),                # Full black frame around the plot area
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(colour = "grey30", size = 10.5),
    legend.position = "right",
    legend.box.background = element_rect(colour = "grey80", linewidth = 0.5)
  ) +
  labs(
    title = "Sample QC: Variance-Stabilized PCA (All 18 Samples)",
    subtitle = "Primary separation by cell lineage (PC1) and LPS treatment (PC2); SRR2938964 shows attenuated LPS response",
    colour = "Treatment",
    shape = "Cell Type"
  )

# Display the PCA plot in RStudio's Plots pane
print(p_pca)

# Save a high-resolution (300 DPI) PNG image of the PCA plot
ggsave("Sample_QC_plots/01_VST_PCA.png", plot = p_pca, width = 9.5, height = 6.8, dpi = 300)

# ==================================================================================================
# 6. QC Plot 2: Sample-to-Sample Euclidean Distance Heatmap 
# ==================================================================================================

# Transpose the VST expression matrix (so samples are rows) and calculate pairwise Euclidean distances
sampleDists <- dist(t(assay(vsd_all)))

# Build the heatmap silently (silent = TRUE) so we can align the legends and viewport margins first
p_dist <- pheatmap(
  as.matrix(sampleDists),                                            # Convert distance object into an 18x18 matrix
  clustering_distance_rows = sampleDists,                            # Hierarchical tree clusters rows by whole-transcriptome similarity
  clustering_distance_cols = sampleDists,                            # Hierarchical tree clusters columns by whole-transcriptome similarity
  col = shared_colors,                                               # Unified blue-yellow-red palette matching Plot 3
  annotation_col = clean_annot,                                      # Add Cell_Type and Treatment color bars on top
  annotation_row = clean_annot,                                      # Add matching Cell_Type and Treatment color bars on left
  annotation_colors = annot_colors,                                  # Apply unified colorblind-safe annotation colors
  border_color = "white",                                            # Crisp white gridlines between cells
  fontsize = 10,
  fontsize_row = 9,
  fontsize_col = 9,
  angle_col = 45,                                                    # Angle bottom sample IDs at 45 degrees
  main = "Sample QC: Sample-to-Sample Euclidean Distance (VST Counts)",
  silent = TRUE
)

# Align the Treatment & Cell_Type legend vertically with the top of the numeric distance colorbar
gt_dist <- p_dist$gtable
leg_idx <- which(gt_dist$layout$name == "legend")                    # Locate the numeric distance colorbar
ann_idx <- which(gt_dist$layout$name == "annotation_legend")         # Locate the Treatment & Cell_Type legend

if (length(leg_idx) == 1 && length(ann_idx) == 1) {
  # Match both the top (t) and bottom (b) grid rows of the annotation legend to the colorbar
  gt_dist$layout[ann_idx, c("t", "b")] <- gt_dist$layout[leg_idx, c("t", "b")]
  # Nudge the internal legend grobs upward so "Treatment" sits flush with the top "300" tick mark
  ann_grob <- gt_dist$grobs[[ann_idx]]
  for (i in seq_along(ann_grob$children)) {
    if (!is.null(ann_grob$children[[i]]$y)) {
      ann_grob$children[[i]]$y <- ann_grob$children[[i]]$y + unit(0.18, "npc")
    }
  }
  gt_dist$grobs[[ann_idx]] <- ann_grob
}

# Save inside a scaled viewport (92% width/height) so outer labels and legends never touch the edge
png("Sample_QC_plots/02_Sample_Distance_Heatmap.png", width = 10.5, height = 8.5, units = "in", res = 300)
grid.newpage()
pushViewport(viewport(width = 0.92, height = 0.92))
grid.draw(gt_dist)
popViewport()
dev.off()

# ==================================================================================================
# 7. QC Plot 3: FACS Sort-Purity Marker Genes Heatmap
# ==================================================================================================

# Define the 12 canonical cell-type marker genes validated in Srinivasan et al. (2016) Figure 4b
marker_df <- data.frame(
  gene_name = c("Itgam", "Aif1", "Cx3cr1", "Ptprc",       # Rows 1-4: Microglial markers (CD11b+ sorted cells)
                "Gfap", "Aldh1l1", "Aqp4", "Slc1a2",      # Rows 5-8: Astrocytic markers (GFAP+ sorted cells)
                "Rbfox3", "Grin1", "Cx3cl1", "Nos1"),     # Rows 9-12: Neuronal markers (NeuN+ negative control for contamination)
  Marker_Type = factor(
    c(rep("Microglia", 4), rep("Astrocytes", 4), rep("Neurons (Neg Ctrl)", 4)),
    levels = c("Microglia", "Astrocytes", "Neurons (Neg Ctrl)")
  ),
  stringsAsFactors = FALSE
)

# Extract library-size-normalized counts from dds_all and attach readable mouse gene symbols
norm_all <- counts(dds_all, normalized = TRUE) %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  mutate(gene_id = sub("\\.[0-9]+$", "", gene)) %>%
  left_join(gene_annotation, by = "gene_id")

# Join from marker_df on the left so all 12 markers stay in strict biological order
marker_mat <- marker_df %>%
  left_join(
    norm_all %>% filter(gene_name %in% marker_df$gene_name) %>% distinct(gene_name, .keep_all = TRUE),
    by = "gene_name"
  ) %>%
  filter(!is.na(gene))

# Log2-transform the normalized counts (+1 pseudocount to avoid log2(0)) for clear visual contrast
expr_mat <- log2(as.matrix(marker_mat[, all_samples]) + 1)

# Label heatmap rows with official gene symbols
rownames(expr_mat) <- marker_mat$gene_name

# Create row annotation bar showing which cell type each marker belongs to
row_annot <- data.frame(Marker_Type = marker_mat$Marker_Type, row.names = marker_mat$gene_name)

# Build the sort-purity heatmap without clustering so LPS & Saline groups stay cleanly ordered
p_markers <- pheatmap(
  expr_mat,
  cluster_rows = FALSE,                                                        # Keep markers grouped: Microglia (1-4) -> Astrocytes (5-8) -> Neurons (9-12)
  cluster_cols = FALSE,                                                        # Keep samples ordered: Microglia LPS | Saline | Astrocytes LPS | Saline
  col = shared_colors,                                                         # Apply the same unified blue-yellow-red expression gradient
  annotation_row = row_annot,                                                  # Show marker category bar on the left
  annotation_col = clean_annot,                                                # Show sample Cell_Type and Treatment bars on top
  annotation_colors = annot_colors,                                            # Apply unified colorblind-safe annotation colors
  gaps_row = c(4, 8),                                                          # White horizontal dividers separating marker categories
  gaps_col = c(5, 10, 14),                                                     # White vertical dividers separating all 4 experimental groups
  border_color = "white",                                                      # Crisp white gridlines between cells
  fontsize = 10,
  fontsize_row = 10,
  fontsize_col = 9,
  angle_col = 45,                                                              # Angle sample IDs at 45 degrees for clean readability
  main = "Sample QC: FACS Sort-Purity Markers (log2 Normalized Counts + 1)",
  silent = TRUE
)

# Scale the gtable inside a centered 88% viewport (shifted slightly right) so both "Marker_Type" and "Neurons (Neg Ctrl)" fit completely
png("Sample_QC_plots/03_Sort_Purity_Markers.png", width = 11.5, height = 7.5, units = "in", res = 300)
grid.newpage()
pushViewport(viewport(x = 0.52, y = 0.5, width = 0.88, height = 0.90))
grid.draw(p_markers$gtable)
popViewport()
dev.off()

# ==================================================================================================
# 8. QC Plot 4: 3x4 Purity Marker Grid 
# ==================================================================================================

# Calculate gene lengths in kilobases (kb) matching the filtered genes in dds_all
gene_lengths_kb <- counts$length[match(rownames(dds_all), counts$gene)] / 1000

# Calculate total sequencing depth per sample in millions of mapped reads
lib_sizes_m     <- colSums(counts(dds_all)) / 1e6

# Compute FPKM (Fragments Per Kilobase Million) so linear y-axis scales match Srinivasan et al.'s nRPKM
fpkm_mat        <- sweep(sweep(counts(dds_all), 1, gene_lengths_kb, "/"), 2, lib_sizes_m, "/")

# Convert FPKM matrix to a data frame and attach readable mouse gene symbols from GENCODE M39
fpkm_df <- as.data.frame(fpkm_mat) %>%
  rownames_to_column("gene") %>%
  mutate(gene_id = sub("\\.[0-9]+$", "", gene)) %>%
  left_join(gene_annotation, by = "gene_id")

# Define the 12 canonical marker genes in the 3-row layout of Srinivasan et al. Fig. 4b:
# Row 1: Neuronal markers (Negative Control) | Row 2: Astrocytic markers | Row 3: Microglial markers
fig4b_markers <- data.frame(
  gene_name = c("Rbfox3", "Grin1", "Cx3cl1", "Nos1",      # Row 1: Neuronal markers (Neg Ctrl)
                "Gfap", "Aldh1l1", "Aqp4", "Slc1a2",      # Row 2: Astrocytic markers
                "Itgam", "Aif1", "Cx3cr1", "Ptprc"),      # Row 3: Microglial markers
  Marker_Category = c(
    rep("Neuronal\nmarkers\n(Neg Ctrl)", 4),
    rep("Astrocytic\nmarkers", 4),
    rep("Microglial\nmarkers", 4)
  ),
  stringsAsFactors = FALSE
)

# Reshape the 12 markers from wide to tidy long format and label Sorted Cell Types & Treatments explicitly
fig4b_long <- fig4b_markers %>%
  left_join(
    fpkm_df %>% filter(gene_name %in% fig4b_markers$gene_name) %>% distinct(gene_name, .keep_all = TRUE),
    by = "gene_name"
  ) %>%
  pivot_longer(cols = all_of(all_samples), names_to = "sample", values_to = "fpkm") %>%
  left_join(all_meta %>% rownames_to_column("sample"), by = "sample") %>%
  mutate(
    # Explicitly label the FACS-sorted cell populations (GFAP+ Astrocytes vs CD11b+ Microglia)
    Sorted_Cell_Type = factor(
      if_else(cell_type == "Astrocytes", "Astrocytes (GFAP+ cells)", "Microglia (CD11b+ cells)"),
      levels = c("Astrocytes (GFAP+ cells)", "Microglia (CD11b+ cells)")
    ),
    # Order Treatment factor with Saline baseline first, followed by LPS stimulation
    Treatment = factor(condition, levels = c("Saline", "LPS"))
  )

# Define a helper function that builds one clean Fig. 4b-style panel per marker gene
make_gene_panel <- function(gene_sym, show_ylab = FALSE) {
  # Subset expression values for the single target marker gene
  sub_df <- fig4b_long %>% filter(gene_name == gene_sym)
  
  # Calculate overall Mean and +/- 1 SD across all replicates per sorted cell population
  sum_df <- sub_df %>%
    group_by(Sorted_Cell_Type) %>%
    summarise(
      mean_val = mean(fpkm),
      sd_val   = sd(fpkm),
      .groups  = "drop"
    )
  
  # Set the upper y-axis limit 18% above the highest dot or error bar, starting strictly at 0
  y_max <- max(c(sub_df$fpkm, sum_df$mean_val + sum_df$sd_val), na.rm = TRUE) * 1.18
  if (y_max == 0) y_max <- 1
  
  ggplot() +
    # Draw vertical +/- 1 SD error bars per sorted cell population
    geom_errorbar(
      data = sum_df,
      aes(x = Sorted_Cell_Type, ymin = pmax(mean_val - sd_val, 0), ymax = mean_val + sd_val),
      width = 0.22,
      linewidth = 0.45,
      colour = "black"
    ) +
    # Draw horizontal crossbar indicating the population mean
    geom_crossbar(
      data = sum_df,
      aes(x = Sorted_Cell_Type, y = mean_val, ymin = mean_val, ymax = mean_val),
      width = 0.42,
      fatten = 1.4,
      colour = "black"
    ) +
    # Overlay individual mouse replicates: Color = Sorted Cell Type, Shape = Treatment (Saline vs LPS)
    geom_jitter(
      data = sub_df,
      aes(x = Sorted_Cell_Type, y = fpkm, colour = Sorted_Cell_Type, shape = Treatment),
      width = 0.14,
      size = 2.5,
      alpha = 0.95
    ) +
    # Apply unified colorblind-safe cell-type colors and distinct shapes for Saline (circle) vs LPS (triangle)
    scale_colour_manual(values = c("Astrocytes (GFAP+ cells)" = "#009E73", "Microglia (CD11b+ cells)" = "#E69F00")) +
    scale_shape_manual(values  = c("Saline" = 16, "LPS" = 17)) +
    scale_y_continuous(limits = c(0, y_max), expand = c(0, 0)) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold.italic", hjust = 0.5, size = 11.5),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),                                           # Hide repeated x-axis text (shown in legend)
      axis.ticks.x = element_blank(),
      axis.title.y = if (show_ylab) element_text(face = "bold", size = 10) else element_blank(),
      axis.text.y = element_text(colour = "black", size = 9),
      legend.position = "none",
      plot.margin = margin(6, 10, 6, 6)
    ) +
    labs(title = gene_sym, y = "FPKM")
}

# Generate all 12 individual gene panels ("FPKM" y-axis title only on column 1 of each row)
panel_list <- lapply(seq_along(fig4b_markers$gene_name), function(i) {
  make_gene_panel(fig4b_markers$gene_name[i], show_ylab = (i %% 4 == 1))
})

# Create bold row-header text grobs identifying the 3 Marker Gene Categories on the left
row_label_1 <- textGrob("Neuronal\nmarkers\n(Neg Ctrl)", gp = gpar(fontface = "bold", fontsize = 11))
row_label_2 <- textGrob("Astrocytic\nmarkers",           gp = gpar(fontface = "bold", fontsize = 11))
row_label_3 <- textGrob("Microglial\nmarkers",           gp = gpar(fontface = "bold", fontsize = 11))

# Build a dummy plot to extract a clean, boxed two-part legend (Color = Cell Population, Shape = Treatment)
legend_plot <- ggplot(fig4b_long, aes(x = Sorted_Cell_Type, y = fpkm, colour = Sorted_Cell_Type, shape = Treatment)) +
  geom_point(size = 3.2) +
  scale_colour_manual(
    values = c("Astrocytes (GFAP+ cells)" = "#009E73", "Microglia (CD11b+ cells)" = "#E69F00"),
    name = "FACS-Sorted Cell Population"
  ) +
  scale_shape_manual(
    values = c("Saline" = 16, "LPS" = 17),
    name = "Treatment Condition"
  ) +
  theme_void() +
  theme(
    legend.position = "left",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 10),
    legend.spacing.y = unit(0.35, "cm"),
    legend.box.background = element_rect(colour = "black", fill = "white", linewidth = 0.6),
    legend.box.margin = margin(8, 10, 8, 10)
  )

# Safely extract the combined legend grob
leg_grob_list <- ggplotGrob(legend_plot)$grobs
shared_legend <- leg_grob_list[[grep("guide-box", sapply(leg_grob_list, function(x) x$name))[1]]]

# Combine the 3 Marker Category row headers and 12 gene panels into a 3-row x 5-column layout
grid_body <- arrangeGrob(
  row_label_1, panel_list[[1]],  panel_list[[2]],  panel_list[[3]],  panel_list[[4]],
  row_label_2, panel_list[[5]],  panel_list[[6]],  panel_list[[7]],  panel_list[[8]],
  row_label_3, panel_list[[9]],  panel_list[[10]], panel_list[[11]], panel_list[[12]],
  ncol = 5,
  widths = c(1.15, 1, 1, 1, 1)
)

# Place the two-part legend on the left beside the 3x4 marker grid with a clear title
final_fig4b <- arrangeGrob(
  shared_legend,
  grid_body,
  ncol = 2,
  widths = c(1.75, 5.2),
  top = textGrob(
    "Sample QC: Cell-Type Marker Expression Across FACS-Sorted Populations (All 18 Samples)",
    gp = gpar(fontface = "bold", fontsize = 13)
  ),
  padding = unit(0.35, "in")
)

# Open a crisp 300 DPI PNG graphics device and draw the complete figure
png("Sample_QC_plots/04_Sort_Purity_Fig4b_Style.png", width = 12.2, height = 7.2, units = "in", res = 300)
grid.newpage()
grid.draw(final_fig4b)
dev.off()

# Print confirmation message to the RStudio Console when all 4 figures have finished saving
cat("\nDONE! Plots Saved.\n")