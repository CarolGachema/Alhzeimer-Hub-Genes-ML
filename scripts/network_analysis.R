############################################
# Alzheimer's Disease Network Analysis
# Dataset: GSE5281
# Part 1: Differential Expression + Protein-Protein Interaction Network Analysis
# Author: Caroline Wambui Gachema
############################################


# Load Required Libraries 

BiocManager::install("GEOquery")
library(GEOquery)
library(limma)          # For differential expression analysis
library(hgu133plus2.db) # For Affymetrix probe annotation
library(AnnotationDbi)  if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("STRINGdb")# For gene ID mapping

library(STRINGdb)       # For protein-protein interaction network
library(igraph)         # For network analysis and visualization
library(ggplot2)        # For enhanced visualizations
library(ggraph)         # For network visualization using ggplot2


# Set Working Directory
setwd("C:/Users/Caroline/Documents/MY PROJECTS/Alhzeimer-Hub-Genes-ML")


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
