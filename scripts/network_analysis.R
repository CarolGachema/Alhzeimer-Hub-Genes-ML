############################################
# Alzheimer's Disease Network Analysis
# Dataset: GSE528
############################################


# Load Required Libraries 

BiocManager::install("GEOquery")
library(GEOquery)
library(limma)          # For differential expression analysis
library(hgu133plus2.db) # For Affymetrix probe annotation
library(AnnotationDbi) 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("STRINGdb")# For gene ID mapping
BiocManager::install(c(
  "clusterProfiler",   # For KEGG pathway enrichment
  "org.Hs.eg.db",      # Human gene annotation database
  "enrichplot",        # Enrichment visualisation
  "pheatmap",          # Beautiful heatmaps
  "ComplexHeatmap"     # Advanced heatmap with annotations
))

install.packages(c("RColorBrewer", "viridis", "factoextra", "ggdendro"))

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(factoextra)
library(AnnotationDbi)

library(STRINGdb)       # For protein-protein interaction network
library(igraph)         # For network analysis and visualization
library(ggplot2)        # For enhanced visualizations
library(ggraph)         # For network visualization using ggplot2


# Set Working Directory
setwd("C:/Users/Gachema/Documents/MY PROJECTS/Alhzeimer-Hub-Genes-ML")


# Load GEO Dataset
# GSE5281 contains gene expression data from multiple brain regions of Alzheimer's Disease patients and healthy controls Using Affymetrix Human Genome U133 Plus 2.0 Array (GPL570)
gse <- getGEO(
  filename = "data/GSE5281_series_matrix.txt.gz",
  getGPL = FALSE
)


# Extract Expression Data and Metadata 
expression_data <- exprs(gse)    # Gene expression matrix
metadata <- pData(gse)           # Sample metadata


# Preview available characteristics to identify group labels
unique(metadata$characteristics_ch1)
head(metadata)


# Define Experimental Groups

# Classify samples as Alzheimer's Disease (AD) or Normal Control (NC)
# Use Disease State column for group classification
group <- ifelse(
  grepl("normal|control", metadata$`Disease State:ch1`, ignore.case = TRUE),
  "Control",
  "AD"
)


# Convert to factor
group <- factor(group, levels = c("Control", "AD"))


# Check counts
table(group)


# Convert to factor with Control as reference level
group <- factor(group, levels = c("Control", "AD"))


# Verify sample counts
table(group)
message(paste("Total samples:", length(group)))
message(paste("AD samples:", sum(group == "AD")))
message(paste("Control samples:", sum(group == "Control")))


# Build Design Matrix 
# No-intercept design matrix for two-group comparison
design <- model.matrix(~0 + group)
colnames(design) <- c("Control", "AD")
design


# Define Contrast 

# Compare AD vs Control to identify disease-associated gene expression changes
contrast <- makeContrasts(
  AD - Control,
  levels = design
)


# Fit Linear Model and Apply Empirical Bayes Statistics

fit  <- lmFit(expression_data, design)    # Fit linear model
fit2 <- contrasts.fit(fit, contrast)      # Apply contrast
fit2 <- eBayes(fit2)                      # Empirical Bayes moderation


# Extract Differential Expression Results 
results <- topTable(
  fit2,
  adjust.method = "fdr",    # Benjamini-Hochberg FDR correction
  number = Inf              # Return all genes
)

head(results)

# Filter Significant DEGs

# Apply significance thresholds:
# FDR-adjusted p-value < 0.05 AND |log2FC| > 1 (2-fold change)
sig_genes <- results[
  results$adj.P.Val < 0.05 & abs(results$logFC) > 1,
]

message(paste("Total significant DEGs:", nrow(sig_genes)))

# Annotate Probes with Gene Symbols

# Map Affymetrix probe IDs to HGNC gene symbols
gene_symbols <- mapIds(
  hgu133plus2.db,
  keys = rownames(sig_genes),
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)


# Add gene symbols to significant genes table
sig_genes$GeneSymbol <- gene_symbols


# Remove probes with no gene symbol annotation
sig_genes_annotated <- sig_genes[!is.na(sig_genes$GeneSymbol) & sig_genes$GeneSymbol != "", ]

message(paste("Annotated significant DEGs:", nrow(sig_genes_annotated)))

# View top 10 annotated DEGs
top10 <- sig_genes_annotated[order(sig_genes_annotated$adj.P.Val), ][1:10, ]
print(top10[, c("GeneSymbol", "logFC", "adj.P.Val")])

# Separate Upregulated and Downregulated Genes
up_genes   <- sig_genes_annotated[sig_genes_annotated$logFC > 1, ]
down_genes <- sig_genes_annotated[sig_genes_annotated$logFC < -1, ]

message(paste("Upregulated in AD:", nrow(up_genes)))
message(paste("Downregulated in AD:", nrow(down_genes)))

# Generate Volcano Plot

# Check the range of logFC values
summary(results$logFC)


# Filter out extreme outliers before plotting
results_filtered <- results[abs(results$logFC) < 10, ]

# Plot with filtered data
results_filtered$colour <- "Not Significant"
results_filtered$colour[results_filtered$adj.P.Val < 0.05 & results_filtered$logFC > 1]  <- "Upregulated"
results_filtered$colour[results_filtered$adj.P.Val < 0.05 & results_filtered$logFC < -1] <- "Downregulated"

volcano2 <- ggplot(results_filtered, aes(x = logFC, y = -log10(adj.P.Val), color = colour)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c(
    "Not Significant" = "grey70",
    "Upregulated"     = "#E74C3C",
    "Downregulated"   = "#3498DB"
  )) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  labs(
    title    = "Alzheimer's Disease vs Control: Differential Gene Expression",
    subtitle = paste("Dataset: GSE5281 | Filtered for |logFC| < 10"),
    x        = "Log2 Fold Change",
    y        = "-log10 Adjusted P-value",
    color    = "Expression Status"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

ggsave("figures/volcano_plot_AD_filtered.png", volcano2, width = 10, height = 7, dpi = 300)


# Build Protein-Protein Interaction Network

# Initialize STRING database for Homo sapiens (taxon ID: 9606)
# Version 11.5, minimum interaction score = 400 (medium confidence)
string_db <- STRINGdb$new(
  version        = "11.5",
  species        = 9606,
  score_threshold = 400,
  input_directory = "data/"
)


# Prepare gene list for STRING mapping
# Use top 200 most significant DEGs for network clarity
top_genes <- sig_genes_annotated[order(sig_genes_annotated$adj.P.Val), ][1:200, ]
gene_list <- data.frame(GeneSymbol = top_genes$GeneSymbol)


# Map gene symbols to STRING IDs
mapped_genes <- string_db$map(
  gene_list,
  "GeneSymbol",
  removeUnmappedRows = TRUE
)


message(paste("Genes mapped to STRING:", nrow(mapped_genes)))


# Retrieve Interaction Network

# Get all interactions between our mapped genes from STRING database
interactions <- string_db$get_interactions(mapped_genes$STRING_id)

message(paste("Total interactions retrieved:", nrow(interactions)))


# Build igraph Network Object

# Create network from interaction edges
network <- graph_from_data_frame(
  interactions[, c("from", "to", "combined_score")],
  directed = FALSE
)


# Map STRING IDs back to gene symbols for readable labels
id_to_symbol <- setNames(mapped_genes$GeneSymbol, mapped_genes$STRING_id)
V(network)$name <- id_to_symbol[V(network)$name]


# Remove any nodes with NA names
network <- delete_vertices(network, which(is.na(V(network)$name)))


message(paste("Network nodes (genes):", vcount(network)))
message(paste("Network edges (interactions):", ecount(network)))


# Calculate Network Statistics

# Degree centrality: number of connections each gene has
# Hub genes = highly connected genes = potential disease drivers!!
degree_centrality     <- degree(network)
betweenness_centrality <- betweenness(network, normalized = TRUE)


# Add centrality measures to network
V(network)$degree      <- degree_centrality
V(network)$betweenness <- betweenness_centrality


# Identify top 20 hub genes by degree centrality
hub_genes <- sort(degree_centrality, decreasing = TRUE)[1:20]
message("Top 20 Hub Genes by Degree Centrality:")
print(hub_genes)


# Visualize the Network

# Create a sub network of top 50 hub genes for cleaner visualization
top50_hubs  <- names(sort(degree_centrality, decreasing = TRUE)[1:50])
sub_network <- induced_subgraph(network, top50_hubs)


# Add fold change information for colouring nodes
fc_values <- setNames(sig_genes_annotated$logFC, sig_genes_annotated$GeneSymbol)
V(sub_network)$logFC <- fc_values[V(sub_network)$name]
V(sub_network)$logFC[is.na(V(sub_network)$logFC)] <- 0


# Generate network visualisation
network_plot <- ggraph(sub_network, layout = "fr") +
  geom_edge_link(aes(alpha = 0.3), color = "grey70") +
  geom_node_point(aes(
    size  = degree(sub_network),
    color = V(sub_network)$logFC
  )) +
  geom_node_text(
    aes(label = name),
    repel    = TRUE,
    size     = 3,
    max.overlaps = 20
  ) +
  scale_color_gradient2(
    low      = "#3498DB",
    mid      = "white",
    high     = "#E74C3C",
    midpoint = 0,
    name     = "Log2FC"
  ) +
  scale_size_continuous(name = "Degree", range = c(3, 12)) +
  labs(
    title    = "Protein-Protein Interaction Network: Top 50 Hub Genes in Alzheimer's Disease",
    subtitle = "Node size = connectivity | Node colour = expression direction"
  ) +
  theme_graph(base_family = "Arial") +
  theme(plot.title = element_text(face = "bold", size = 12))


# Save network plot
ggsave("figures/network_plot_AD.png", network_plot, width = 14, height = 10, dpi = 300)
message("Network plot saved!!")


# Save All Results

# Save DEG results
top_downregulated <- down_genes[order(down_genes$adj.P.Val), ][1:500, ]
write.csv(top_downregulated, "results/top_500_downregulated_genes.csv", row.names = FALSE)
top_significant <- results[order(results$adj.P.Val), ][1:300, ]
write.csv(top_significant, "results/top_300_significant_genes.csv", row.names = FALSE)
write.csv(down_genes,           "results/downregulated_genes_AD.csv")


# Save hub genes
hub_genes_df <- data.frame(
  GeneSymbol = names(hub_genes),
  Degree     = as.numeric(hub_genes)
)
write.csv(hub_genes_df, "results/hub_genes_AD.csv", row.names = FALSE)


# Save network statistics
network_stats <- data.frame(
  GeneSymbol    = V(network)$name,
  Degree        = degree_centrality,
  Betweenness   = betweenness_centrality
)
write.csv(network_stats, "results/network_statistics_AD.csv", row.names = FALSE)

message("All results saved!!")

# Save Workspace

# Save entire workspace for Part 2 Python ML analysis
save.image("results/alzheimers_workspace.RData")


#################################
# Prepare ML Dataset for Python
##############################



# Select top 20 hub genes
top20_gene_names <- names(hub_genes)[1:20]

# Match expression matrix rows using gene symbols
probe_to_gene <- mapIds(
  hgu133plus2.db,
  keys = rownames(expression_data),
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)

# Create expression dataframe
expression_df <- as.data.frame(expression_data)

# Add gene symbols
expression_df$GeneSymbol <- probe_to_gene

# Remove NA symbols
expression_df <- expression_df[!is.na(expression_df$GeneSymbol), ]

# Keep only hub genes
ml_data <- expression_df[
  expression_df$GeneSymbol %in% top20_gene_names,
]

# Remove duplicate genes
ml_data <- ml_data[!duplicated(ml_data$GeneSymbol), ]

# Set gene names as rownames
rownames(ml_data) <- ml_data$GeneSymbol
ml_data$GeneSymbol <- NULL

# Transpose so samples become rows
ml_data_t <- as.data.frame(t(ml_data))

# Add disease labels
ml_data_t$Group <- group

# Save ML-ready dataset
write.csv(
  ml_data_t,
  "results/alzheimers_ml_dataset.csv",
  row.names = TRUE
)

message("ML dataset exported successfully!")


#############################################
# HEATMAP OF TOP 20 HUB GENES
# Purpose: Visualise expression of hub genes across all 161 samples
#############################################


message("Generating hub gene heatmap...")

#Extract hub gene expression matrix 
# Get the probe IDs for our top 20 hub genes
# hub_genes_df was created in the main script
top20_gene_names <- names(hub_genes)[1:20]
top_hub_symbols <- hub_genes $GeneSymbol[1:20]

# Map gene symbols back to probe IDs using annotation database
hub_probe_ids <- mapIds(
  hgu133plus2.db,
  keys = top20_gene_names,
  column = "PROBEID",
  keytype = "SYMBOL",
  multiVals = "first"
)

# Remove any NAs
hub_probe_ids <- hub_probe_ids[!is.na(hub_probe_ids)]

# Get the probe IDs that exist in our expression data
valid_probes <- hub_probe_ids[hub_probe_ids %in% rownames(expression_data)]

message(paste("Hub gene probes found in expression data:", length(valid_probes)))

# Extract expression data for hub gene probes
hub_expression <- expression_data[valid_probes, ]

# Replace probe IDs with gene symbols as row names
rownames(hub_expression) <- names(valid_probes)

# ── A2. Prepare sample annotations ────────────────────────────────────────

# Create annotation dataframe for columns (samples)
# Shows which samples are AD vs Control
sample_annotation <- data.frame(
  Group = as.character(group),
  row.names = colnames(hub_expression)
)

# Define colours for the annotation
annotation_colors <- list(
  Group = c(
    "AD"      = "#E74C3C",   # Red for Alzheimer's
    "Control" = "#3498DB"    # Blue for healthy controls
  )
)

# ── A3. Generate heatmap ───────────────────────────────────────────────────

# Scale expression by row (gene) so patterns are visible across samples
# Each gene gets z-scored so high/low expression is relative to that gene's mean
hub_scaled <- t(scale(t(hub_expression)))

# Create the heatmap
png("figures/hub_genes_heatmap.png", width = 1200, height = 900, res = 150)

heatmap(
  hub_scaled,
  annotation_col   = sample_annotation,
  annotation_colors = annotation_colors,
  color            = colorRampPalette(c("#2980B9", "white", "#E74C3C"))(100),
  cluster_rows     = TRUE,     # Cluster genes by expression similarity
  cluster_cols     = TRUE,     # Cluster samples by expression similarity
  show_colnames    = FALSE,    # Too many samples to show names
  show_rownames    = TRUE,     # Show gene names
  fontsize_row     = 11,
  fontsize         = 10,
  border_color     = NA,
  main             = "Top 20 Hub Genes: Expression Across AD and Control Samples\nGSE5281 | Row-scaled expression (z-score)",
  # Add gaps between AD and Control samples
  gaps_col         = sum(group == "Control")
)

dev.off()
message("Heatmap saved to figures/hub_genes_heatmap.png!!")



#########################################################################
# KEGG PATHWAY ENRICHMENT ANALYSIS
# Purpose: Identify biological pathways enriched among significant DEGs
#########################################################################


# Convert gene symbols to Entrez IDs 
# KEGG requires Entrez gene IDs (numbers) not gene symbols
# We convert our significant DEG symbols to Entrez IDs
entrez_ids <- mapIds(
  org.Hs.eg.db,
  keys    = sig_genes_annotated$GeneSymbol,
  column  = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first"
)


# Remove NAs
entrez_ids <- entrez_ids[!is.na(entrez_ids)]
message(paste("Genes successfully converted to Entrez IDs:", length(entrez_ids)))

# Run KEGG enrichment

# enrichKEGG tests whether specific KEGG pathways are over-represented
# among our significant DEGs compared to the entire human genome
kegg_results <- enrichKEGG(
  gene          = entrez_ids,
  organism      = "hsa",         # hsa = Homo sapiens
  pAdjustMethod = "fdr",         # Benjamini-Hochberg FDR correction
  pvalueCutoff  = 0.05,          # Significance threshold
  qvalueCutoff  = 0.05           # FDR threshold
)


options(timeout = 300)


# Preview top 10 enriched KEGG pathways
print("Top 10 Enriched KEGG Pathways:")
print(head(as.data.frame(kegg_results), 10))


message(paste("Total significant KEGG pathways:", nrow(as.data.frame(kegg_results))))


# Visualise KEGG results

# Barplot: top 15 most enriched pathways
# Bar length = gene count | Colour = adjusted p-value
png("figures/kegg_barplot.png", width = 1100, height = 800, res = 150)

barplot(
  kegg_results,
  showCategory = 15,
  title        = "Top 15 KEGG Pathways Enriched in Alzheimer's Disease DEGs\nGSE5281 | FDR-corrected enrichment",
  font.size    = 10
)

dev.off()
message("KEGG barplot saved!!")


# Dotplot: richer visualisation showing gene ratio and significance
# Dot size = gene ratio | Dot colour = adjusted p-value
png("figures/kegg_dotplot.png", width = 1100, height = 800, res = 150)

dotplot(
  kegg_results,
  showCategory = 15,
  title        = "KEGG Pathway Enrichment Dotplot — Alzheimer's Disease\nGSE5281 | Dot size = gene ratio | Colour = significance"
)

dev.off()
message("KEGG dotplot saved!!")




####################################################
# MULTI-BRAIN-REGION COMPARISON
# Purpose: Compare DEGs and hub genes across different brain regions
# ##################################################


# Identify brain regions in metadata

# Check what brain regions are available in the dataset
print("Available brain regions:")
print(unique(metadata$`Organ Region:ch1`))


# Assign names to your group vector so R can align them by sample ID
names(group) <- rownames(metadata)

# Clean up trailing spaces in the brain region names (e.g., "hippocampus ")
metadata$`Organ Region:ch1` <- trimws(metadata$`Organ Region:ch1`)

# Define the regions variable cleanly
regions <- unique(metadata$`Organ Region:ch1`)
regions <- regions[!is.na(regions)]


# Create region-specific DEG analyses

# Define the regions to compare
# Adjust these based on what unique(metadata$`Organ Region:ch1`) shows you!!
regions <- unique(metadata$`Organ Region:ch1`)
regions <- regions[!is.na(regions)]


message(paste("Brain regions found:", paste(regions, collapse=", ")))


# Store results for each region
region_results <- list()
region_sig_counts <- data.frame(
  Region = character(),
  Total_DEGs = integer(),
  Upregulated = integer(),
  Downregulated = integer(),
  stringsAsFactors = FALSE
)


for(region in regions) {
  
  message(paste("Analysing region:", region))
  
  # Get samples from this brain region
  region_samples <- rownames(metadata)[
    !is.na(metadata$`Organ Region:ch1`) &
      metadata$`Organ Region:ch1` == region
  ]
  
  # Skip regions with too few samples
  if(length(region_samples) < 6) {
    message(paste("Skipping", region, "— too few samples:", length(region_samples)))
    next
  }
  
  # Subset expression data and group labels for this region
  region_expression <- expression_data[, region_samples]
  region_group <- group[region_samples]
  
  # Skip if only one group present
  if(length(unique(region_group)) < 2) {
    message(paste("Skipping", region, "— only one group present!!"))
    next
  }
  
  # Build design matrix for this region
  region_design <- model.matrix(~0 + region_group)
  colnames(region_design) <- c("Control", "AD")
  
  # Define contrast
  region_contrast <- makeContrasts(AD - Control, levels = region_design)
  
  # Fit model
  region_fit  <- lmFit(region_expression, region_design)
  region_fit2 <- contrasts.fit(region_fit, region_contrast)
  region_fit2 <- eBayes(region_fit2)
  
  # Extract results
  region_res <- topTable(
    region_fit2,
    adjust.method = "fdr",
    number = Inf
  )
  
  # Filter significant genes
  region_sig <- region_res[
    region_res$adj.P.Val < 0.05 & abs(region_res$logFC) > 1,
  ]
  
  # Store results
  region_results[[region]] <- region_sig
  
  # Count up and downregulated
  region_sig_counts <- rbind(region_sig_counts, data.frame(
    Region        = region,
    Total_DEGs    = nrow(region_sig),
    Upregulated   = sum(region_sig$logFC > 1),
    Downregulated = sum(region_sig$logFC < -1)
  ))
}


print("DEG counts by brain region:")
print(region_sig_counts)


# Visualise region comparison

if(nrow(region_sig_counts) > 0) {
  
  # Reshape for plotting
  region_long <- tidyr::pivot_longer(
    region_sig_counts,
    cols = c("Upregulated", "Downregulated"),
    names_to = "Direction",
    values_to = "Count"
  )
  
  # Plot
  region_plot <- ggplot(region_long, aes(x = Region, y = Count, fill = Direction)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.85) +
    scale_fill_manual(values = c("Upregulated" = "#E74C3C", "Downregulated" = "#3498DB")) +
    labs(
      title    = "Differentially Expressed Genes by Brain Region",
      subtitle = "Alzheimer's Disease vs Control | GSE5281",
      x        = "Brain Region",
      y        = "Number of Significant DEGs",
      fill     = "Expression Direction"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title    = element_text(face = "bold", size = 14),
      axis.text.x   = element_text(angle = 30, hjust = 1),
      legend.position = "top"
    )
  
  ggsave("figures/brain_region_comparison.png", region_plot,
         width = 12, height = 7, dpi = 300)
  
  message("Brain region comparison plot saved!!")
  
  # Save region results
  write.csv(region_sig_counts, "results/brain_region_DEG_counts.csv", row.names = FALSE)
}



###################################################################
# SUBTYPE CLUSTERING OF AD SAMPLES
# Purpose: Do AD samples naturally cluster into molecular subtypes??
###################################################################

message("Running AD subtype clustering analysis...")

# Extract AD samples only

# Subset expression data to AD samples only
ad_samples  <- names(group)[group == "AD"]
ad_expression <- expression_data[, ad_samples]

message(paste("AD samples for clustering:", ncol(ad_expression)))


# Select most variable genes for clustering

# Use top 500 most variable genes across AD samples
# High variance genes drive the biological differences between samples
gene_variance <- apply(ad_expression, 1, var)
top_var_genes <- names(sort(gene_variance, decreasing = TRUE))[1:500]
ad_var_expression <- ad_expression[top_var_genes, ]


# Dimensionality Reduction

# PCA reduces 500 genes to a few principal components
# Helps us visualise the structure of AD samples
ad_pca <- prcomp(t(ad_var_expression), scale. = TRUE)


# How much variance does each PC explain??
pca_variance <- summary(ad_pca)$importance[2, 1:10] * 100
message("Variance explained by first 10 PCs:")
print(round(pca_variance, 2))


# Determine optimal number of clusters

# Elbow plot: find where adding more clusters stops helping
png("figures/clustering_elbow_plot.png", width = 900, height = 600, res = 150)


fviz_nbclust(
  ad_pca$x[, 1:10],   # Use first 10 PCs
  kmeans,
  method = "wss",
  k.max  = 8
) +
  labs(
    title    = "Optimal Number of AD Subtypes",
    subtitle = "Within-cluster sum of squares — look for the 'elbow'"
  ) +
  theme_minimal()

dev.off()
message("Elbow plot saved — check it to decide number of clusters!!")


# K-means clustering

# Based on elbow plot — try k=3 first (adjust if needed!!)
K <- 3

set.seed(42)
kmeans_result <- kmeans(
  ad_pca$x[, 1:10],
  centers  = K,
  nstart   = 25,
  iter.max = 100
)

# Add cluster labels to PCA results
pca_df <- data.frame(
  PC1     = ad_pca$x[, 1],
  PC2     = ad_pca$x[, 2],
  Cluster = factor(kmeans_result$cluster)
)


message("AD Subtype Cluster Sizes:")
print(table(kmeans_result$cluster))


# Visualise clusters

cluster_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(size = 3.5, alpha = 0.8) +
  stat_ellipse(level = 0.75, linewidth = 1) +
  scale_color_manual(values = c("#E74C3C", "#3498DB", "#2ECC71",
                                "#F39C12", "#9B59B6")) +
  labs(
    title    = paste0("Alzheimer's Disease Molecular Subtypes (k=", K, ")"),
    subtitle = "PCA of top 500 variable genes | k-means clustering",
    x        = paste0("PC1 (", round(pca_variance[1], 1), "% variance)"),
    y        = paste0("PC2 (", round(pca_variance[2], 1), "% variance)"),
    color    = "Subtype"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title  = element_text(face = "bold", size = 14),
    legend.position = "right"
  )

ggsave("figures/ad_subtypes_pca.png", cluster_plot,
       width = 10, height = 7, dpi = 300)

message("AD subtype PCA plot saved!!")


# Save cluster assignments

cluster_assignments <- data.frame(
  Sample  = rownames(pca_df),
  Subtype = kmeans_result$cluster
)


write.csv(cluster_assignments, "results/ad_subtype_assignments.csv", row.names = FALSE)
message("Cluster assignments saved!!")


# Save workspace

save.image("results/alzheimers_workspace.RData")
message("Workspace updated and saved!!")
