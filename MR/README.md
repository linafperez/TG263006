# Mendelian Randomization (MR) Analysis
This directoy contains the workflow used to perform Mendelian Randomization (MR) analyses investigating potential causal relationships between gut microbiome taxa and Parkinson’s disease (PD). The pipeline reproduces the analytical procedures used in the study to evaluate genetically instrumented microbial abundances and their association with PD risk.

The workflow integrates microbiome genome-wide association study (GWAS) summary statistics with Parkinson’s disease GWAS datasets from independent cohorts. Instrumental variables are derived from microbiome-associated genetic variants and tested against PD GWAS summary statistics using multiple MR estimators to ensure robust causal inference.

Two population-level analyses are implemented:

- **EUR cohort** — European ancestry PD GWAS datasets (FinnGen and GP2)
- **LATAM cohort** — Latin American PD GWAS dataset (LARGE-PD)

The analysis evaluates whether genetic variants associated with microbial taxa abundance exert a causal effect on Parkinson’s disease risk.

---

## Pipeline Overview

The MR analysis pipeline consists of the following stages:

1. Download microbiome GWAS summary statistics
2. Download Parkinson’s disease GWAS summary statistics
3. Select instrumental variables associated with microbiome taxa
4. Harmonize exposure and outcome datasets
5. Perform Mendelian Randomization analyses
6. Conduct sensitivity analyses
7. Generate summary statistics and plots

---

## Software Requirements

The analysis was performed using the following software versions.

- R 4.5.2

### R Packages

- TwoSampleMR
- MendelianRandomization
- MRPRESSO
- data.table
- tidyverse
- ggplot2

---

## Directory Structure

```
MR
│
├── Exposure
│   └── SNPs_exposure.R
│
├── Outcome
│   └── SNP_outcome.R
│
└── MR
    └── MR.Rmd
```

---

## Workflow

### Step 1 — Obtain Microbiome GWAS Summary Statistics

Exposure data consist of genetic variants associated with gut microbiome taxa abundance derived from the MiBioGen consortium meta-analysis.

The script below extracts SNPs significantly associated with microbial taxa and prepares the exposure dataset.

```r
Rscript Exposure/SNPs_exposure.R
```

Output:

```
exposure_snps.tsv
```

---

### Step 2 — Obtain Parkinson’s Disease GWAS Summary Statistics

Outcome data correspond to Parkinson’s disease GWAS summary statistics from independent cohorts.

European population datasets:

- FinnGen
- GP2

Latin American population dataset:

- LARGE-PD

The following script retrieves and formats the outcome dataset.

```r
Rscript Outcome/SNP_outcome.R
```

Output:

```
outcome_snps.tsv
```

---

### Step 3 — Harmonization of Exposure and Outcome Datasets

Exposure and outcome datasets are harmonized to ensure consistent allele orientation and effect direction.

This step removes:

- ambiguous SNPs (A/T or C/G)
- palindromic variants
- variants with inconsistent allele orientation

Harmonization is performed automatically during the MR analysis.

---

### Step 4 — Mendelian Randomization Analysis

The core MR analysis evaluates the causal effect of microbiome taxa abundance on Parkinson’s disease risk.

The following estimators are implemented:

- Inverse Variance Weighted (IVW)
- Weighted Median
- Weighted Mode
- Wald Ratio

The MR analysis is executed using the following script:

```r
Rscript MR/MR.Rmd
```

Outputs:

```
mr_results.tsv
forest_plot.pdf
```

---

### Step 5 — Sensitivity Analyses

To evaluate the robustness of causal estimates, several sensitivity analyses are performed.

These include:

- Cochran’s Q test for heterogeneity
- MR-Egger intercept to detect horizontal pleiotropy
- MR-PRESSO to identify outlier SNPs
- Leave-one-out analysis to assess the influence of individual variants

Outputs:

```
heterogeneity_results.tsv
pleiotropy_test.tsv
leave_one_out_results.tsv
```

---

## Final Outputs

The MR pipeline generates the following outputs:

- harmonized exposure–outcome datasets
- causal effect estimates for microbiome taxa
- sensitivity analysis results
- visualization of MR effect sizes

These results provide the statistical framework used to infer potential causal relationships between gut microbiome composition and Parkinson’s disease risk.

