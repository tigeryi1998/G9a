# install
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

# Load DESeq2 library
library(DESeq2)

library(dplyr)

library(tidyr)

# 1. Read counts matrix (tab-delimited)
counts <- read.table("gene_counts_cleaned.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

# 2. Read metadata (CSV)
meta <- read.csv("metadata.csv", header = TRUE, stringsAsFactors = FALSE)

# 3. Reorder metadata rows to match count columns
# meta <- meta[match(colnames(counts), meta$sample), ]

# 3. Load gene map table, GeneID - GeneName
gene_map <- read.csv("ensembl_gene_map.csv", stringsAsFactors = FALSE)


# Check that sample names match exactly
if (!all(colnames(counts) == meta$sample)) {
  stop("Sample names in counts and metadata do not match!")
}


# 4. Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = meta,
                              design = ~ condition)

# set "WT" as the reference (baseline) level
dds$condition <- relevel(dds$condition, ref = "WT")


# 5. Pre-filter genes with low counts (keep genes with total counts >= 10)
dds <- dds[rowSums(counts(dds)) >= 10, ]


# 6. Run DESeq analysis
dds <- DESeq(dds)


# 7: Define full contrast vectors as a named list

## contrast = c("condition", "A", "B")
## log2 fold change of A relative to B, or log2(A / B)
## how much more A is than B
# KO vs WT shows the overall gene knockout effect.
# WT+M vs WT isolates the effect of morpholino treatment alone.
# KO isoforms vs WT show how isoforms rescue or alter gene expression compared to normal.
# KO isoforms vs KO show the differences between knockout and specific isoforms.
# KO isoforms vs each other reveal detailed differences between isoforms.

contrasts <- list(
  M_vs_WT    = c("condition", "WT+M", "WT"),   # Morpholino vs WT
  KO_vs_WT   = c("condition", "KO", "WT"),     # KO vs WT
  
  KO1_vs_WT  = c("condition", "KO-1", "WT"), # short FE, +E10
  KO2_vs_WT  = c("condition", "KO-2", "WT"), # short FE, -E10
  KO3_vs_WT  = c("condition", "KO-3", "WT"), # long FE, +E10
  KO4_vs_WT  = c("condition", "KO-4", "WT"), # long FE, -E10
  
  KO1_vs_KO  = c("condition", "KO-1", "KO"), # short FE, +E10
  KO2_vs_KO  = c("condition", "KO-2", "KO"), # short FE, -E10
  KO3_vs_KO  = c("condition", "KO-3", "KO"), # long FE, +E10
  KO4_vs_KO  = c("condition", "KO-4", "KO"), # long FE, -E10
  
  KO1_vs_KO2 = c("condition", "KO-1", "KO-2"), # +E10 vs -E10, short FE
  KO3_vs_KO4 = c("condition", "KO-3", "KO-4"), # +E10 vs -E10, long FE
  
  KO1_vs_KO3 = c("condition", "KO-1", "KO-3"), # short vs long FE, +E10
  KO2_vs_KO4 = c("condition", "KO-2", "KO-4") # short vs long FE, -E10
)




# 8: Extract log2 fold changes for each contrast vs WT
lfc_results <- data.frame(GeneID = rownames(dds))

for (name in names(contrasts)) {
  res <- results(dds, contrast = contrasts[[name]])
  lfc_results[[name]] <- res$log2FoldChange
}

# Merge with gene names
lfc_results <- merge(lfc_results, gene_map, by = "GeneID", all.x = TRUE)

# Reorder columns for readability
lfc_results <- lfc_results[, c("GeneID", "GeneName", names(contrasts))]

# 9. Replace NA values with 0 (optional)
lfc_viz <- lfc_results
lfc_viz[is.na(lfc_viz)] <- 1e-300

# Create output directory if it doesn't exist
output_dir <- "DE_CSV"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# 10. Save results as CSV
write.csv(lfc_viz, file = file.path(output_dir, "DE_log2FC.csv"), row.names = FALSE)




# 11. Extract both log2FC and padj
## Initialize a list to store results per contrast
full_results <- list()

for (name in names(contrasts)) {
  res <- results(dds, contrast = contrasts[[name]])
  df <- data.frame(
    GeneID = rownames(res),
    log2FC = res$log2FoldChange,
    padj = res$padj
  )
  # Merge in gene names
  df <- merge(df, gene_map, by = "GeneID", all.x = TRUE)
  
  # Cap padj=0 to avoid Inf from log10
  df$padj[df$padj == 0] <- 1e-300
  
  # Compute -log10 adjusted p-value
  df$neg_log_padj <- ifelse(is.na(df$padj), NA, -log10(df$padj))
  
  # Reorder columns
  df <- df[, c("GeneID", "GeneName", "log2FC", "padj", "neg_log_padj")]
  full_results[[name]] <- df
}


# 12. Filter for significant genes

for (name in names(full_results)) {
  # Filter: exclude NAs first!
  df <- full_results[[name]]
  df <- df[!is.na(df$padj) & !is.na(df$log2FC), ]
  
  # Apply significance filter: padj < 0.05 and abs(log2FC) > 0.3
  sig <- subset(df, padj < 0.05 & abs(log2FC) > 0.3)
  
  # Sort both full and sig by absolute log2FC descending
  df <- df[order(-abs(df$log2FC)), ]
  sig <- sig[order(-abs(sig$log2FC)), ]
  
  # Save to ./DE_CSV/
  write.csv(df, file = file.path(output_dir, paste0("DE_", name, ".csv")), row.names = FALSE)
  write.csv(sig, file = file.path(output_dir, paste0("DE_", name, "_sig.csv")), row.names = FALSE)
}







# 13. G9a gene
## EHMT2 
## ENSG00000204371.12


## KO1 > KO3 > KO4 > KO2 > KO > WT > WT+M

counts(dds)["ENSG00000204371.12", ]
# KO-0-1   KO-0-2   KO-1-1   KO-1-2   KO-2-1   KO-2-2   KO-3-1   KO-3-2   KO-4-1   KO-4-2     WT-1     WT-2 
# 7266     5066  1156775  1196388   768695   617530   978497   907427   646337   786207     7152     5915 
# WT-3 WTPosM-1 WTPosM-2 WTPosM-3 
# 4941     3474     3736     3057 


## KO1 > KO3 > KO2 > KO4 > KO > WT > WT+M

lfc_results[lfc_results$GeneName == "EHMT2", ]

#                   GeneID GeneName    M_vs_WT  KO_vs_WT KO1_vs_WT KO2_vs_WT KO3_vs_WT KO4_vs_WT KO1_vs_KO KO2_vs_KO
# 15852 ENSG00000204371.12    EHMT2 -0.6117749 0.3314718  7.947728  7.531955  7.733433  7.416948  7.616257  7.200484
#         KO3_vs_KO KO4_vs_KO KO1_vs_KO2 KO3_vs_KO4 KO1_vs_KO3 KO2_vs_KO4
# 15852  7.401961  7.085476  0.4157729  0.3164852  0.2142955  0.1150078




# Save plot to variable
p1 <- plotCounts(dds, gene = "ENSG00000204371.12", intgroup = "condition", returnData = TRUE)

# Plot manually using ggplot2 for full control
library(ggplot2)

gg <- ggplot(p1, aes(x = condition, y = count)) +
  geom_jitter(width = 0.2, size = 3, alpha = 0.8, color = "steelblue") +
  scale_y_log10() +
  labs(title = "Normalized counts for EHMT2 (G9a)",
       x = "Condition", y = "Normalized count (log10)") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Display plot
ggsave("EHMT2_counts.png", gg, width = 6, height = 4, dpi = 300)




# Filter the row for EHMT2
ehmt2_row <- lfc_results[lfc_results$GeneName == "EHMT2", ]

# Remove the GeneID and GeneName columns
ehmt2_long <- data.frame(
  Comparison = names(ehmt2_row)[-(1:2)],
  log2FC = as.numeric(ehmt2_row[ , -(1:2)])
)

# Add direction for coloring
ehmt2_long$Direction <- ifelse(ehmt2_long$log2FC > 0, "Up", "Down")

# Load ggplot2 and plot
library(ggplot2)
gg2 <- ggplot(ehmt2_long, aes(x = reorder(Comparison, log2FC), y = log2FC, fill = Direction)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  scale_fill_manual(values = c("Up" = "firebrick", "Down" = "steelblue")) +
  coord_flip() +
  labs(
    title = "Log2 Fold Changes of EHMT2 (G9a)",
    x = "Comparison",
    y = "log2 Fold Change"
  ) +
  theme_minimal(base_size = 14)


# Display plot
ggsave("EHMT2_log2FC.png", gg2, width = 6, height = 4, dpi = 300)




# Pearson Correlation 

library(ggplot2)
library(reshape2)
library(ggcorrplot)

# Remove GeneID and GeneName
lfc_matrix <- lfc_viz[ , -(1:2)]

# Compute correlation matrix
cor_mat <- cor(lfc_matrix, method = "pearson", use = "pairwise.complete.obs")

# Plot correlation matrix
gg3 <- ggcorrplot(cor_mat, 
           type = "lower", 
           lab = TRUE, 
           lab_size = 3, 
           colors = c("steelblue", "white", "firebrick"), 
           title = "Correlation of log2FC Across Conditions", 
           tl.cex = 10)





# signature genes

signature_genes <- read.csv("signatures.csv")

signature_genes_name <- unique(signature_genes$Gene)



# Heat Map 

library(pheatmap)

# Extract matrix of log2FC values
lfc_mat <- as.matrix(lfc_viz[, -(1:2)])
rownames(lfc_mat) <- lfc_viz$GeneName  # or use GeneID

# Compute variance across each gene
gene_vars <- apply(lfc_mat, 1, var)

# Select top 10 variable genes
top10_var_genes <- names(sort(gene_vars, decreasing = TRUE))[1:10]

# Add EHMT2 if not already in top10
genes_to_label <- unique(c(top10_var_genes, "EHMT2", signature_genes_name))

# Create vector of row labels: only label selected genes, blank for others
row_labels <- ifelse(rownames(lfc_mat) %in% genes_to_label, rownames(lfc_mat), "")

# Plot heatmap
gg4 <- pheatmap::pheatmap(lfc_mat,
                cluster_rows = TRUE,
                cluster_cols = TRUE,
                labels_row = row_labels,          # selectively label
                fontsize_row = 6,
                fontsize_col = 10,
                color = colorRampPalette(c("blue", "white", "red"))(100),
                breaks = seq(-5, 5, length.out = 101),  # symmetric color scale around 0
                main = "Heatmap of log2 Fold Changes (EHMT2 + Top 10 Variable Genes)",
                filename = "Heatmap_log2FC.png",     # ✅ Save directly here
                width = 10,
                height = 15
                )

# Subset lfc_results by GeneName
lfc_subset <- lfc_results[lfc_results$GeneName %in% genes_to_label, ]





library(ComplexHeatmap)
library(circlize)

# Prepare matrix: all genes
lfc_mat <- as.matrix(lfc_viz[, -(1:2)])
rownames(lfc_mat) <- lfc_viz$GeneName

# Compute gene variance
gene_vars <- apply(lfc_mat, 1, var)

# Get top 10 most variable genes
top10_var_genes <- names(sort(gene_vars, decreasing = TRUE))[1:10]

# Add EHMT2
genes_to_label <- unique(c(top10_var_genes, "EHMT2", signature_genes_name))

# Create annotation for row labels
row_anno <- rowAnnotation(
  labels = anno_mark(
    at = which(rownames(lfc_mat) %in% genes_to_label),
    labels = rownames(lfc_mat)[rownames(lfc_mat) %in% genes_to_label],
    labels_gp = gpar(fontsize = 8),
    side = "left",
    link_width = unit(1, "cm")
  )
)

# Plot heatmap
ht <- Heatmap(
  lfc_mat,
  name = "log2FC",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = FALSE,  # Hide all labels in row
  show_column_names = TRUE,
  col = colorRamp2(c(-5, 0, 5), c("blue", "white", "red")),
  column_names_gp = gpar(fontsize = 10),
  heatmap_legend_param = list(title = "log2FC"),
  left_annotation = row_anno,
  row_names_max_width = unit(10, "cm"),
  row_dend_width = unit(1, "cm")
)

# Save to file
png("Heatmap_log2FC_ComplexHeatmap.png", width = 1000, height = 1200, res = 150)
draw(ht)
dev.off()







# Volcano Plots

library(ggplot2)
library(ggrepel)

plot_volcano <- function(df, contrast_name, label_genes = "EHMT2", output_dir = "volcano_plots", N=30) {
  dir.create(output_dir, showWarnings = FALSE)
  
  # Filter valid rows
  df <- df[!is.na(df$padj) & !is.na(df$log2FC), ]
  
  # Compute sig points
  padj_cutoff <- 0.05
  log2FC_cutoff <- 0.5
  
  df$sig <- ifelse(df$padj < padj_cutoff & abs(df$log2FC) > log2FC_cutoff, "Significant", "Not Significant")
  
  # Significant genes
  sig_genes <- df[df$sig == "Significant", ]
  
  # Top N genes by absolute log2FC
  top_genes_df <- sig_genes[order(-abs(sig_genes$log2FC)), c("GeneName", "log2FC", "padj")]
  top_genes_df <- top_genes_df[!duplicated(top_genes_df$GeneName), ]  # remove any duplicates
  top_genes_df <- head(top_genes_df, N)
  
  # For labeling, combine top N genes with input label_genes
  label_genes_combined <- unique(c(top_genes_df$GeneName, label_genes))
  df$label <- ifelse(df$GeneName %in% label_genes_combined, df$GeneName, NA)
  
  # Compute -log10 adjusted p-value
  df$neg_log_padj <- ifelse(is.na(df$padj), NA, -log10(df$padj))
  
  # Label top genes with padj and log2fc
  # top_by_fc <- head(sig_genes[order(-abs(sig_genes$log2FC)), "GeneName"], N)
  # top_by_fc <- unique(top_by_fc)
  # label_genes <- unique(c(top_by_fc, label_genes)) 
  # df$label <- ifelse(df$GeneName %in% label_genes, df$GeneName, NA)
  
  # Plot
  p <- ggplot(df, aes(x = log2FC, y = neg_log_padj)) +
    geom_point(aes(color = sig), alpha = 0.3, size = 0.6) +
    scale_color_manual(values = c("Significant" = "firebrick", "Not Significant" = "lightblue")) +
    geom_vline(xintercept = c(-log2FC_cutoff, log2FC_cutoff), linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "gray50") +
    geom_text_repel(
      data = df[!is.na(df$label), ],
      aes(label = label),
      size = 2.5,
      box.padding = 0.25,
      max.overlaps = 200,
      segment.color = "gray30"
    ) +
    labs(
      title = paste("Volcano Plot:", contrast_name),
      x = "log2 Fold Change",
      y = expression(-log[10]~"(adjusted p-value)")
    ) +
    # ylim(0, y_limit) +
    theme_minimal(base_size = 13) +
    theme(panel.background = element_rect(fill = "white", color = NA))
  
  # Save
  ggsave(
    filename = file.path(output_dir, paste0("volcano_", contrast_name, ".png")),
    plot = p, width = 15, height = 10, dpi = 300
  )
  
  # Return top N genes with log2FC and padj
  return(top_genes_df)
}


# Loop through and plot all relevant contrasts

top_genes_by_contrast <- list()

for (name in names(full_results)) {
  df <- full_results[[name]]
  top_genes <- plot_volcano(df, contrast_name = name, label_genes = "EHMT2", N=50)
  top_genes_by_contrast[[name]] <- top_genes
}








# Look at how many timese EHMT2 G9a appears in the top genes in each contrasts
sum(sapply(top_genes_by_contrast, function(df) "EHMT2" %in% df$GeneName))








# Look at gene rank

get_gene_ranks <- function(rank_gene, full_results) {
  gene_ranks <- data.frame(
    Contrast = character(),
    Rank = numeric(),
    log2FC = numeric(),
    padj = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (name in names(full_results)) {
    df <- full_results[[name]]
    sig_df <- df[!is.na(df$padj) & df$padj < 0.05 & !is.na(df$log2FC) & !is.na(df$GeneName), ]
    
    if (!rank_gene %in% sig_df$GeneName) {
      rank <- NA
    } else {
      ordered_genes <- sig_df[order(-abs(sig_df$log2FC)), "GeneName"]
      rank <- match(rank_gene, ordered_genes)
    }
    
    gene_row <- df[df$GeneName == rank_gene, ]
    fc <- if (nrow(gene_row) > 0) gene_row$log2FC[1] else NA
    pval <- if (nrow(gene_row) > 0) gene_row$padj[1] else NA
    
    gene_ranks <- rbind(
      gene_ranks,
      data.frame(
        Contrast = name,
        Rank = rank,
        log2FC = fc,
        padj = pval
      )
    )
  }
  
  return(gene_ranks)
}

print(get_gene_ranks("EHMT2", full_results))








# KO vs WT down, then up in KO# vs KO

rescue <- function(
    full_results,
    query = "KO_vs_WT",
    keys = c("KO1_vs_KO", "KO2_vs_KO", "KO3_vs_KO", "KO4_vs_KO"),
    query_log2fc_threshold = 1.0,
    key_log2fc_threshold=1.0,
    padj_threshold = 0.05,
    is_query_down = TRUE,  # TRUE if query wants negative log2fc
    is_key_down = FALSE    # FALSE if key wants positive log2fc
) {
  
  # Vectorized helper to check if gene name is valid
  is_valid_gene <- function(name) {
    !is.na(name) & !grepl("^ENSG", name)
  }
  
  # Filter query contrast
  query_df <- full_results[[query]]
  query_df <- subset(
    query_df,
    is_valid_gene(GeneName) &
      !is.na(padj) & padj < padj_threshold &
      !is.na(log2FC) & (
        if (is_query_down == TRUE) log2FC < -query_log2fc_threshold else log2FC > query_log2fc_threshold
      )
  )
  
  # Start result with GeneName and log2FC from query
  result <- data.frame(
    GeneName = query_df$GeneName,
    check.names = FALSE
  )
  result[[query]] <- query_df$log2FC
  
  # Track matching key contrasts
  matched_any_key <- rep(FALSE, nrow(result))
  
  for (key in keys) {
    key_df <- full_results[[key]]
    key_df <- subset(
      key_df,
      is_valid_gene(GeneName) &
        !is.na(padj) & padj < padj_threshold &
        !is.na(log2FC) & (
          if (is_key_down == FALSE) log2FC > key_log2fc_threshold else log2FC < -key_log2fc_threshold
        )
    )
    
    # Match gene names
    log2fc_vals <- key_df$log2FC[match(result$GeneName, key_df$GeneName)]
    result[[key]] <- log2fc_vals
    
    # Mark if this key had a valid (non-NA) match
    matched_any_key <- matched_any_key | !is.na(log2fc_vals)
  }
  
  # Keep genes with at least one matched key contrast
  result <- result[matched_any_key, ]
  
  return(result)
}




rescue_genes <- rescue(full_results,
                       query = "KO_vs_WT",
                       keys = c("KO1_vs_KO", "KO2_vs_KO", "KO3_vs_KO", "KO4_vs_KO"),
                       query_log2fc_threshold = 5.0,
                       key_log2fc_threshold = 5.0,
                       padj_threshold = 0.05,
                       is_query_down = TRUE,  # TRUE if query wants negative log2fc
                       is_key_down = FALSE    # FALSE if key wants positive log2fc
                       )

print(rescue_genes, row.names = FALSE)

write.csv(rescue_genes, file = "rescue.csv", row.names = FALSE)





# signature genes

signature_genes <- read.csv("signatures.csv")

signature_genes_name <- unique(signature_genes$Gene)

lfc_signature <- lfc_results %>% filter(GeneName %in% signature_genes_name)

signature_list <- lapply(full_results, function(df) {
  subset(df, GeneName %in% signature_genes_name)
})



rescue_signature <- rescue(signature_list,
                       query = "KO_vs_WT",
                       keys = c("KO1_vs_KO", "KO2_vs_KO", "KO3_vs_KO", "KO4_vs_KO"),
                       query_log2fc_threshold = 0.5,
                       key_log2fc_threshold = 0.5,
                       padj_threshold = 0.05,
                       is_query_down = TRUE,  # TRUE if query wants negative log2fc
                       is_key_down = FALSE    # FALSE if key wants positive log2fc
)

print(rescue_signature, row.names = FALSE)


