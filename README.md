#  Mapping the Molecular Landscape of Alzheimer's Disease
## A Multi-Omics Bioinformatics Analysis: Transcriptomics · Network Biology · Machine Learning


## Overview

This project employs a three stage multi-omics analytical framework to investigate the molecular landscape of Alzheimer's disease (AD), using publicly available gene expression data from the NCBI Gene Expression Omnibus. By combining differential gene expression analysis, protein-protein interaction network analysis, and machine learning classification, this project moves from identifying which genes are dysregulated in Alzheimer's disease, to mapping how those genes interact, to asking whether their expression patterns can *predict* disease status computationally.

Alzheimer's disease affects over 55 million people globally. It's been projected to reach 139 million by 2050. Despite decades of research, no disease-modifying treatment exists. Understanding the molecular mechanisms that drive neurodegeneration is one of the most urgent frontiers in modern biomedical science.


## Dataset

| Parameter | Details |
|-----------|---------|
| **GEO Accession** | [GSE5281](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5281) |
| **Platform** | Affymetrix Human Genome U133 Plus 2.0 Array (GPL570) |
| **Comparison** | Alzheimer's Disease vs Healthy Controls |
| **Brain Regions** | Multiple regions including Entorhinal Cortex, Hippocampus, and others |
| **Samples** | 161 total (100 AD, 61 Control) |
| **Source** | NCBI Gene Expression Omnibus (GEO) |


## Analytical Framework

```
Raw GEO Data → Expression Matrix → Group Classification (AD/Control)
       ↓
Stage 1: TRANSCRIPTOMICS
limma DEG Analysis → FDR Correction → Significant DEGs
       ↓
Stage 2: NETWORK ANALYSIS
STRINGdb PPI Network → igraph Centrality → Hub Gene Identification
       ↓
Stage 3: MACHINE LEARNING
Hub Gene Features → Python Classifier → Disease Prediction Model
```


## Methods

### Stage 1 — Transcriptomics: Differential Gene Expression Analysis
- **Package:** `limma` (Linear Models for Microarray Analysis)
- **Moderation:** Empirical Bayes (eBayes)
- **Multiple Testing Correction:** Benjamini-Hochberg FDR
- **Significance Thresholds:** adj.P.Val < 0.05 AND |log2FC| > 1
- **Annotation:** hgu133plus2.db (Affymetrix probe-to-gene mapping)
- **Visualisation:** ggplot2 volcano plot with ggrepel gene labels

### Stage 2 — Network Analysis: Protein-Protein Interaction Mapping
- **Database:** STRINGdb (version 11.5, score threshold ≥ 400)
- **Species:** Homo sapiens (taxon ID: 9606)
- **Network Analysis:** igraph (degree centrality, betweenness centrality)
- **Hub Gene Identification:** Top genes by degree centrality
- **Visualisation:** ggraph network plot (node size = connectivity, node colour = expression direction)

### Stage 3 — Machine Learning: Disease Classification
- **Language:** Python
- **Features:** Hub genes identified from network analysis
- **Task:** Binary classification (AD vs Control)
- **Framework:** scikit-learn
- *(In active development)*


## Key Results

| Metric | Value |
|--------|-------|
| **Total Samples** | 161 (100 AD, 61 Control) |
| **Significant DEGs** | See results/significant_genes_AD.csv |
| **Network Nodes** | See results/network_statistics_AD.csv |
| **Top Hub Genes** | See results/hub_genes_AD.csv |

### Visualisations

**Volcano Plot** — Colour-coded differential expression (red = upregulated in AD, blue = downregulated in AD)
 `figures/volcano_plot_AD.png` | `figures/volcano_plot_AD_filtered.png`

**PPI Network Plot** — Top 50 hub genes with interaction edges (node size = degree centrality, node colour = log2FC)
 `figures/network_plot_AD.png`


## Repository Structure

```
Alzheimers-Network-ML-Analysis/
│
├── data/
│   └── README.md                          # Instructions to download GSE5281
│
├── scripts/
│   ├── alzheimers_network.R               # Full annotated R pipeline (Stages 1 & 2)
│   └── ml_classifier.py                   # Python ML classifier (Stage 3 — in development)
│
├── results/
│   ├── all_genes_AD.csv                   # All genes with full DEG statistics
│   ├── significant_genes_AD.csv           # Filtered significant DEGs
│   ├── upregulated_genes_AD.csv           # Genes upregulated in AD
│   ├── downregulated_genes_AD.csv         # Genes downregulated in AD
│   ├── hub_genes_AD.csv                   # Top hub genes by degree centrality
│   └── network_statistics_AD.csv          # Full network centrality statistics
│
└── figures/
    ├── volcano_plot_AD.png                # Full volcano plot
    ├── volcano_plot_AD_filtered.png       # Filtered volcano plot (|logFC| < 10)
    └── network_plot_AD.png               # PPI network visualisation
```


## How to Reproduce This Analysis

### Prerequisites
Install required R packages:

```r
install.packages("BiocManager")
BiocManager::install(c(
  "GEOquery",
  "limma",
  "hgu133plus2.db",
  "AnnotationDbi",
  "STRINGdb"
))

install.packages(c("igraph", "ggplot2", "ggrepel", "ggraph"))
```

Install required Python packages:
```python
pip install scikit-learn pandas numpy matplotlib seaborn
```

### Steps
1. Clone this repository
2. Download GSE5281 from [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5281) and place in `data/`
3. Open `scripts/alzheimers_network.R` in RStudio
4. Run section by section — workspace saves automatically after each run!!
5. Results save to `results/` and figures to `figures/`


## Biological Context

**Why Multiple Brain Regions??**
Alzheimer's disease does not affect all brain regions equally. Neurodegeneration follows a characteristic spatiotemporal pattern, beginning in the entorhinal cortex and hippocampus before spreading to other regions. Analyzing gene expression across multiple brain regions provides a more complete picture of the molecular heterogeneity of the disease.

**Why Network Analysis??**
Differentially expressed genes do not act in isolation. Network analysis reveals the topology of molecular interactions, identifying hub genes that occupy critical positions and whose dysregulation may propagate dysfunction throughout entire biological pathways. Hub genes identified through network analysis are strong candidates for therapeutic targeting and biomarker development.

**Why Machine Learning??**
If the molecular signature of Alzheimer's disease is sufficiently distinct and consistent, it should be possible to train a computational model to recognise it, with implications for early detection, patient stratification, and precision medicine approaches to disease management.


## Limitations

- Exploratory analysis — findings require functional validation
- Microarray technology offers less precise quantification than RNA-sequencing
- Single dataset analysis — multi-dataset validation would strengthen findings
- Machine learning component in active development
- Future work should incorporate pathway enrichment (KEGG, Reactome) and survival analysis


## License

This project is licensed under the **MIT License** — free to use, modify, and build upon with attribution.

