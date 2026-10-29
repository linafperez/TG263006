# Multi-omic Integration
This directory contains the workflow used to integrate microbial genomic, transcriptional, and metabolomic datasets generated across the study. The integration framework connects strain-resolved microbial genomes, microbial gene functions, and host metabolite profiles in order to reconstruct functional microbiome–host interactions associated with Parkinson’s disease.

The integration analysis combines results from the metagenomic, metatranscriptomic, 16S rRNA gene, and metabolomic pipelines. Microbial genomes reconstructed from metagenomic assemblies are linked to amplicon sequence variants (ASVs) obtained from 16S datasets, while functional annotations and transcriptional activity are associated with serum metabolite signatures.

Two complementary integration analyses are implemented:

- **MAG–ASV Integration** – linking metagenome-assembled genomes with amplicon sequence variants
- **MAG–Gene–Metabolite Integration** – connecting microbial genes and metabolic pathways with host metabolite signatures

Together, these analyses reconstruct functional relationships connecting microbial taxa, encoded genes, metabolic pathways, and host metabolic responses.

---

## Pipeline Overview

The integration analysis pipeline consists of the following stages:

1. Import microbial genome and ASV datasets
2. Harmonize taxonomic assignments across sequencing platforms
3. Link ASVs to reconstructed microbial genomes
4. Import metabolomic signatures
5. Integrate microbial genes and metabolic pathways with host metabolites
6. Perform correlation and association analyses
7. Generate network representations of microbiome–host interactions

---

## Software Requirements

The following software versions were used for the analysis.

- R 4.5.2

### R Packages

- tidyverse
- data.table
- igraph
- ggplot2
- pheatmap
- reshape2
- mixOmics

---

## Directory Structure

```
Integration
│
├── MAG-ASV
│   └── MAG_ASV-integration.Rmd
│
└── MAG-genes-metabolite
    └── Microbe-gene-metabolite.Rmd
```

---

## Workflow

### Step 1 — Import Microbial Genome and ASV Data

Genome-resolved microbial datasets obtained from the metagenomic pipeline and ASV tables generated from the 16S analysis are imported into the R environment.

Inputs typically include:

```
MAG_taxonomy.tsv
MAG_abundance.tsv
ASV_taxonomy.tsv
ASV_table.tsv
```

The integration analysis begins by loading these datasets and preparing them for taxonomic harmonization.

```r
Rscript MAG-ASV/MAG_ASV-integration.Rmd
```

---

### Step 2 — Taxonomic Harmonization

Because metagenomic and amplicon sequencing pipelines generate taxonomic assignments using different methods and reference databases, a harmonization step is required.

During this step:

- MAG taxonomy assignments derived from GTDB-Tk are standardized
- ASV taxonomy assignments derived from SILVA are mapped to comparable taxonomic ranks
- shared taxonomic identifiers are generated to allow cross-platform comparison

This harmonization enables linking microbial genomes reconstructed from metagenomics with taxa identified in 16S datasets.

---

### Step 3 — MAG–ASV Linking

Amplicon sequence variants are linked to metagenome-assembled genomes using shared taxonomic identifiers and sequence similarity relationships.

This analysis identifies microbial taxa that are consistently detected across sequencing platforms.

Outputs include:

```
MAG_ASV_links.tsv
MAG_ASV_overlap.tsv
```

These tables define correspondences between genome-resolved taxa and amplicon-derived taxa.

---

### Step 4 — Import Metabolomic Signatures

Metabolite signatures identified in the metabolomics pipeline are imported and prepared for integration with microbial genomic and functional datasets.

Input files include:

```
final_metabolite_signature.tsv
normalized_metabolite_matrix.tsv
```

These datasets contain metabolites significantly associated with Parkinson’s disease.

---

### Step 5 — Integration of Microbial Genes and Metabolites

Functional annotations derived from microbial genomes and metatranscriptomic expression profiles are integrated with metabolomic datasets.

The integration analysis connects:

- microbial genes predicted from MAGs
- transcriptional activity measured in metatranscriptomic data
- host metabolite abundance profiles

This step is executed using:

```r
Rscript MAG-genes-metabolite/Microbe-gene-metabolite.Rmd
```

---

### Step 6 — Correlation Analysis

Associations between microbial genes and host metabolites are evaluated using correlation-based analyses.

The analysis identifies microbial pathways whose activity is associated with changes in host metabolic profiles.

Outputs include:

```
gene_metabolite_correlations.tsv
significant_microbe_metabolite_links.tsv
```

---

### Step 7 — Network Reconstruction

Significant associations between microbial genes and metabolites are visualized as interaction networks.

These networks illustrate potential microbiome–host functional interactions linked to Parkinson’s disease.

Outputs:

```
microbiome_metabolite_network.pdf
integration_heatmap.pdf
```

---

## Final Outputs

The integration pipeline produces the following outputs:

- harmonized taxonomic mappings between MAGs and ASVs
- microbial genome–ASV correspondence tables
- microbial gene–metabolite association datasets
- correlation statistics linking microbial functions to host metabolites
- network representations of microbiome–host interactions
