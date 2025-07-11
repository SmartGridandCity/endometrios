# 🧬 Multi-Method Feature Selection Pipeline for Endometrial Gene Expression

This repository contains a reproducible R pipeline for identifying discriminative gene signatures between **eutopic** and **ectopic endometrial tissues**. It is designed to tackle small-N, high-dimensional datasets and includes comprehensive model evaluation, visualization, and guidance for biological interpretation.

---

## Contents

-   `data/`: Raw expression matrix and phenotype labels.
-   `scripts/`: R scripts for running the full selection and analysis pipeline.
-   `results/`: All output tables, plots, and performance metrics.
-   `README.md`: This project overview and guide.

---

## Pipeline Overview

The pipeline executes the following workflow:

1.  **Data Preprocessing**: Loads the gene expression matrix and filters out low-variance genes.
2.  **Feature Selection**: Applies six distinct methods to identify robust gene signatures.
3.  **Model Training**: Trains a Random Forest classifier on each resulting feature set.
4.  **Evaluation & Visualization**: Compares model performance using ROC curves, OOB error, and accuracy metrics.
5.  **Interpretation**: Identifies consensus genes selected across multiple methods for downstream analysis.

---

## 🧠 Feature Selection Methods

This pipeline implements six diverse strategies to ensure a robust and comprehensive feature analysis:

| Method | Type | Description |
| :--- | :--- | :--- |
| **Limma Filter** | Univariate | Filters based on `logFC` from differential expression analysis. |
| **RFECV** | Wrapper | Performs recursive feature elimination with cross-validation. |
| **Gini Importance** | Embedded | Ranks features using Random Forest's native Gini impurity score. |
| **Lasso** | Embedded | Employs L1-penalized regression (`glmnet`) to force sparsity. |
| **Boruta** | Wrapper | Uses a stability selection algorithm with shadow features. |
| **DE Top Filter** | Ranked | Selects the top 20 genes based on `logFC` magnitude. |

> 📝 **Note:** Due to instability with a small sample size (`n = 18`), SHAP and permutation importance were excluded from the final pipeline.

---

## 📈 Results & Key Visualizations

### Model Performance Summary

The Boruta-selected feature set achieved the highest accuracy and lowest Out-of-Bag (OOB) error, indicating a robust and generalizable signature.

| Method | Accuracy | OOB Error |
| :--- | :---: | :---: |
| **Boruta** | **0.94** | **0.10** |
| RFECV | 0.89 | 0.12 |
| Lasso | 0.86 | 0.15 |
| Gini | 0.83 | 0.17 |
| Filtered | 0.78 | 0.22 |

*Full details are available in [`results/ModelComparison.csv`](results/ModelComparison.csv).*

### ROC Curve Comparison

ROC curves illustrate the superior classification performance of the Boruta and RFECV models.

![ROC Curve Comparison](results/ROC_Curve_Comparison.pdf)

### Feature Inclusion Heatmap

This heatmap shows the consensus of selected genes across different methods. Genes selected by multiple methods (darker cells) are stronger candidates for being true biomarkers.

![Feature Inclusion Heatmap](results/FeatureInclusion_Heatmap.pdf)

---

## 🔬 Shared Genes for Downstream Analysis

Genes selected by two or more methods are compiled into a single file to facilitate biological interpretation.

📁 **Find the list here: [`results/SharedGenes_MultiSelection.csv`](results/SharedGenes_MultiSelection.csv)**

These consensus genes are ideal candidates for:
-   Gene Ontology (GO) and pathway enrichment analysis.
-   Prioritization for experimental validation (e.g., qPCR, Western Blot).
-   Building a final, parsimonious biomarker panel.

---

## 🚀 Usage Instructions

### 1. Clone the Repository
```bash
git clone https://github.com/yourusername/endometrial-feature-selection.git
cd endometrial-feature-selection

### 2. Install Dependencies
Open R or RStudio and run the following command to install all required packages:

install.packages(c("caret", "randomForest", "glmnet", "Boruta",
                   "pheatmap", "pROC", "ggplot2", "ComplexUpset"))

### 3. Run the Pipeline
Execute the main script to run the entire analysis from start to finish. All results will be saved to the results/ directory.

source("scripts/run_pipeline.R")

## 📚 Acknowledgments
It draws on established best practices for machine learning in bioinformatics. Special thanks to the creators and maintainers of the open-source R packages used throughout this project.

## 🧵 Contact
For questions, suggestions, or collaboration opportunities, please open an issue on this repository or reach out via email.

## ⚠️ Disclaimer
This pipeline is intended for research and exploratory purposes only. Any biological conclusions drawn from its output should be validated with further experimental work.
