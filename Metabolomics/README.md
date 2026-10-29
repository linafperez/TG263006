# Metabolomics Analysis
This directory contains the workflow used to analyze untargeted metabolomic datasets associated with Parkinson’s disease cohorts. The pipeline reproduces the computational procedures used in the study to identify serum metabolites associated with Parkinson’s disease and prioritize metabolite signatures through statistical and machine learning approaches.

The workflow processes metabolomic feature tables derived from mass spectrometry experiments and performs statistical testing, multivariate modeling, and metabolite prioritization. The resulting metabolite signatures are subsequently integrated with microbial genomic and transcriptional signals to investigate microbiome–host metabolic interactions.

The metabolomic analysis focuses on identifying metabolites that differentiate Parkinson’s disease patients from healthy controls and evaluating their association with microbial pathways reconstructed from metagenomic and metatranscriptomic data.

---

## Pipeline Overview

The metabolomic analysis pipeline consists of the following stages:

1. Import metabolomic feature tables and metadata
2. Data preprocessing and normalization
3. Exploratory data analysis
4. Differential metabolite analysis
5. Multivariate statistical modeling
6. Machine learning classification
7. Metabolite prioritization
8. Generation of metabolite signatures for multi-omic integration

---

## Software Requirements

The following software versions were used for the analysis.

- R 4.5.2

### R Packages

- tidyverse
- data.table
- ggplot2
- mixOmics
- randomForest
- caret
- pheatmap
- ropls

---

## Directory Structure

```
Metabolomics
│
└── Metabolomics_analysis.Rmd
```

---

## Workflow

### Step 1 — Import Metabolomic Data

Metabolomic datasets were obtained from publicly available repositories associated with the Pereira et al. cohort.

The feature table contains metabolite intensities measured using untargeted mass spectrometry in serum samples.

Input files typically include:

```
metabolite_feature_table.tsv
sample_metadata.tsv
```

The analysis begins by loading the data into the R environment.

```r
Rscript Metabolomics_analysis.Rmd
```

During this step, metabolite intensity matrices and associated sample metadata are imported and formatted for downstream analyses.

---

### Step 2 — Data Preprocessing

Metabolomic data undergo preprocessing steps to improve comparability across samples.

These procedures include:

- removal of metabolites with excessive missing values
- normalization of metabolite intensities
- log transformation of feature abundances
- scaling of variables for multivariate analysis

The processed dataset is stored as a normalized metabolite matrix suitable for statistical modeling.

Output:

```
normalized_metabolite_matrix.tsv
```

---

### Step 3 — Exploratory Data Analysis

Exploratory analyses are performed to evaluate the global metabolic structure of the dataset.

These analyses include:

- Principal Component Analysis (PCA)
- clustering of samples based on metabolite profiles
- visualization of metabolite abundance distributions

Outputs:

```
PCA_plot.pdf
metabolite_heatmap.pdf
```

These analyses provide an overview of the metabolic differences between Parkinson’s disease and control samples.

---

### Step 4 — Differential Metabolite Analysis

Statistical tests are performed to identify metabolites significantly associated with Parkinson’s disease.

Differential analysis includes:

- group comparison tests between PD and healthy controls
- multiple testing correction
- effect size estimation

Output:

```
differential_metabolites.tsv
```

This table contains metabolites significantly enriched or depleted in Parkinson’s disease samples.

---

### Step 5 — Multivariate Modeling

Multivariate statistical approaches are applied to identify metabolic patterns that differentiate disease and control groups.

The primary multivariate method used in the analysis is:

- sparse Partial Least Squares Discriminant Analysis (sPLS-DA)

This method identifies metabolites contributing most strongly to group separation.

Output:

```
sPLSDA_model.pdf
```

---

### Step 6 — Machine Learning Classification

Random Forest models are trained to identify metabolites with strong discriminatory power between Parkinson’s disease and control samples.

The analysis includes:

- model training
- cross-validation
- variable importance ranking

Output:

```
random_forest_importance.tsv
```

This table ranks metabolites based on their contribution to classification performance.

---

### Step 7 — Metabolite Prioritization

Final metabolite signatures are defined by integrating results from:

- differential metabolite analysis
- sPLS-DA variable selection
- Random Forest feature importance

Metabolites consistently identified across methods are prioritized as robust disease-associated metabolic signatures.

Output:

```
final_metabolite_signature.tsv
```

---

## Final Outputs

The metabolomic pipeline produces the following outputs:

- normalized metabolite abundance matrix
- differential metabolite analysis results
- multivariate model outputs
- machine learning feature importance rankings
- prioritized metabolite signatures

These metabolite signatures are subsequently integrated with microbial genomic and transcriptomic features in the multi-omic integration analysis to characterize microbiome–host metabolic interactions associated with Parkinson’s disease.

