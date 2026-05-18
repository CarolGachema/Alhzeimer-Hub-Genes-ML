# Identifying Hub Genes and Building a Machine Learning Classifier for Alzheimer's Disease Using Multi Omic Analysis



## Overview

This project employs a three-stage multi-omics analytical framework to investigate the molecular landscape of Alzheimer's disease (AD), using publicly available gene expression data from NCBI Gene Expression Omnibus. By combining differential gene expression analysis, protein-protein interaction network analysis, and machine learning classification, this project moves from identifying *which* genes are dysregulated in Alzheimer's disease, to mapping *how* those genes interact, to asking whether their expression patterns can *predict* disease status computationally.

Alzheimer's disease affects over 55 million people globally — a number projected to reach 139 million by 2050. Despite decades of research, no disease modifying treatment exists. Understanding the molecular mechanisms that drive neurodegeneration is one of the most urgent frontiers in modern biomedical science.


## Dataset

| Parameter | Details |
|-----------|---------|
| **GEO Accession** | [GSE5281](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5281) |
| **Platform** | Affymetrix Human Genome U133 Plus 2.0 Array (GPL570) |
| **Comparison** | Alzheimer's Disease vs Healthy Controls |
| **Total Samples** | 161 (100 AD, 61 Control) |
| **Source** | NCBI Gene Expression Omnibus (GEO) |


## Analytical Framework

```
Raw GEO Data → Expression Matrix → Group Classification (AD/Control)
       ↓
Stage 1: TRANSCRIPTOMICS (R)
limma DEG Analysis → FDR Correction → Significant DEGs → Annotation
       ↓
Stage 2: NETWORK ANALYSIS (R)
STRINGdb PPI Network → igraph Centrality Analysis → Hub Gene Identification
       ↓
Stage 3: MACHINE LEARNING (Python)
Hub Gene Features → 4 Classifiers → Random Forest Best Model → Disease Prediction
```


## Key Results

### Stage 2 — Network Analysis

| Metric | Value |
|--------|-------|
| **Hub Genes Identified** | 20 top hub genes by degree centrality |
| **Top Hub Gene** | ACTB (Beta-actin) |
| **Other Key Hubs** | NDUFS7, IDH3G, CDK7, SARS1, MDH2, CDK5 |

### Stage 3 — Machine Learning

| Metric | Value |
|--------|-------|
| **Best Model** | Random Forest |
| **Cross-Validation Accuracy** | 90.7% |
| **Test Set Accuracy** | 87.88% |
| **ROC AUC Score** | 0.9038 |
| **Training Samples** | 128 |
| **Testing Samples** | 33 |
| **Features Used** | 20 hub genes |

#### Model Comparison

| Model | CV Accuracy |
|-------|------------|
| **Random Forest**  | **0.907** |
| Support Vector Machine | 0.895 |
| Logistic Regression | 0.889 |
| Gradient Boosting | 0.863 |

#### Top Hub Genes by Feature Importance

| Rank | Gene | Biological Role |
|------|------|----------------|
| 1 | **ACTB** | Beta-actin — cytoskeletal integrity, neuronal structure |
| 2 | **NDUFS7** | Mitochondrial complex I — energy metabolism |
| 3 | **IDH3G** | Isocitrate dehydrogenase — TCA cycle |
| 4 | **CDK7** | Cyclin-dependent kinase — cell cycle regulation |
| 5 | **SARS1** | Seryl-tRNA synthetase — protein synthesis |
| 6 | **MDH2** | Malate dehydrogenase — metabolic regulation |
| 7 | **CDK5** | Neuronal kinase — tau phosphorylation in AD pathology |

> CDK5 is a well-established driver of aberrant tau phosphorylation — one of the hallmark pathological features of Alzheimer's disease. Its identification as a key hub gene directly validates the biological relevance of the network-driven feature selection approach.



## Visualisations

| Figure | Description |
|--------|------------|
| `figures/volcano_plot_AD_filtered.png` | Colour-coded DEG volcano plot |
| `figures/network_plot_AD.png` | PPI network — top 50 hub genes |
| `figures/hub_gene_distributions.png` | Hub gene expression: AD vs Control |
| `figures/model_comparison.png` | 4-model CV accuracy comparison |
| `figures/feature_importance.png` | Random Forest hub gene importance |
| `figures/confusion_matrix.png` | Test set classification results |
| `figures/roc_curve.png` | ROC curve (AUC = 0.9038) |


## Repository Structure

```
Alzheimers-Network-ML-Analysis/
│
├── data/
│   └── README.md                          # Download instructions for GSE5281
│
├── scripts/
│   ├── alzheimers_network.R               # Stage 1 & 2: DEG + PPI network (R)
│   └── alzheimers_ml.py                   # Stage 3: ML classification (Python)
│
├── results/
│   ├── significant_genes_AD.csv
│   ├── upregulated_genes_AD.csv
│   ├── downregulated_genes_AD.csv
│   ├── hub_genes_AD.csv
│   ├── network_statistics_AD.csv
│   ├── feature_importance.csv
│   └── model_comparison_summary.csv
│
└── figures/
    ├── volcano_plot_AD_filtered.png
    ├── network_plot_AD.png
    ├── hub_gene_distributions.png
    ├── model_comparison.png
    ├── feature_importance.png
    ├── confusion_matrix.png
    └── roc_curve.png
```


## How to Reproduce

### R Packages
```r
BiocManager::install(c("GEOquery","limma","hgu133plus2.db","AnnotationDbi","STRINGdb"))
install.packages(c("igraph","ggplot2","ggrepel","ggraph"))
```

### Python Packages
```python
pip install pandas numpy matplotlib seaborn scikit-learn
```

### Steps
1. Clone this repository
2. Download GSE5281 from [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5281) → place in `data/`
3. Run `alzheimers_network.R` in RStudio
4. Export hub gene expression as `results/alzheimers_ml_dataset.csv`
5. Run `alzheimers_ml.py` in Jupyter Notebook


## Biological Context and Interpretation

Alzheimer's disease does not affect all brain regions equally. Neurodegeneration follows a characteristic spatiotemporal pattern, beginning in the entorhinal cortex and hippocampus before spreading to other regions. Analyzing gene expression across multiple brain regions provides a more complete picture of the molecular heterogeneity of the disease.
Differentially expressed genes do not act in isolation. Network analysis reveals the topology of molecular interactions, identifying hub genes that occupy critical positions and whose dysregulation may propagate dysfunction throughout entire biological pathways. Hub genes identified through network analysis are strong candidates for therapeutic targeting and biomarker development.

**ACTB** — Beta-actin dysregulation reflects widespread neuronal structural breakdown characteristic of neurodegeneration!!

**CDK5** — A neuronal kinase with established roles in aberrant tau phosphorylation, one of Alzheimer's defining pathological hallmarks!!

**Mitochondrial genes (NDUFS7, IDH3G, MDH2)** - Consistent with growing evidence of bioenergetic failure as an early and targetable event in Alzheimer's disease pathogenesis!!


## Limitations

- Exploratory analysis requiring functional validation
- Microarray less precise than RNA-sequencing
- Small test set (n=33)
- Single dataset — multi-cohort validation would strengthen findings



