# Genome-resolved multi-omic workflow for Parkinson’s disease
Parkinson’s disease (PD) is increasingly hypothesized to originate in the gut, yet the mechanistic links between microbiome dysbiosis and systemic neurodegeneration remain poorly resolved. This study addresses that gap by performing a genome-resolved multi-omic meta-analysis integrating Mendelian randomization (MR), shotgun metagenomics, metatranscriptomics, 16S rRNA gene profiling, and untargeted metabolomics across European and Latin American cohorts. By reconstructing 356 high-quality metagenome-assembled genomes (MAGs) and linking strain-level microbial activity to host lipid remodeling, the work defines a functional immunometabolic interface connecting gut microbial programs to systemic lipid dysregulation in PD. This repository contains the fully reproducible computational workflow underlying that analysis, organized by omic layer and methodological order, enabling stepwise reprocessing, genome reconstruction, differential profiling, and cross-layer functional integration from raw data to mechanistic interpretation.

## Table of Contents
1. [Study Design](#study-design)  
2. [Repository Structure](#repository-structure)  
3. [Workflow Description](#workflow-description)  
   - [3.1. Mendelian Randomization (MR)](#31-mendelian-randomization-mr)  
   - [3.2. Shotgun Metagenomics](#32-shotgun-metagenomics)  
   - [3.3. Metatranscriptomics](#33-metatranscriptomics)  
   - [3.4. 16S rRNA Gene Analysis](#34-16s-rrna-gene-analysis)  
   - [3.5. Metabolomics](#35-metabolomics)  
   - [3.6. Multi-omic Integration](#36-multi-omic-integration)  
4. [Software and Dependencies](#software-and-dependencies)  
5. [Data Availability](#data-availability)  
6. [Reproducibility Notes](#reproducibility-notes)  
7. [Citation](#citation)

## Study Design

This study performs a genome-resolved multi-omic meta-analysis to investigate the role of the gut microbiome in Parkinson’s disease (PD). The workflow integrates multiple omic layers to reconstruct microbial genomic structure, functional activity, and host metabolic interactions.

The analysis combines datasets from independent cohorts in Europe and Latin America and integrates five complementary analytical layers:

1. Mendelian Randomization (MR)
2. Shotgun metagenomics
3. Metatranscriptomics
4. 16S rRNA gene sequencing
5. Untargeted metabolomics

All analyses were performed independently for each omic layer and subsequently integrated to identify strain-resolved microbial signatures associated with Parkinson’s disease.

### Included cohorts

The study includes publicly available datasets from three studies:

| Study | Region | Omics | Samples |
|------|------|------|------|
| Villette et al. | Europe | Shotgun metagenomics + metatranscriptomics | 46 PD / 49 HC |
| Pereira et al. | Europe | 16S rRNA sequencing + metabolomics | 63 PD / 61 HC |
| Forero-Rodríguez et al. | Latin America | 16S rRNA sequencing | 25 PD / 25 HC |

Microbiome data were derived from fecal samples, while metabolomic profiles were obtained from serum samples.

### Study workflow

The analytical workflow was designed to systematically reconstruct the functional landscape of the Parkinson’s disease gut microbiome through a multi-omic integration framework. The repository is organized so that each omic layer can be executed independently, while also allowing sequential execution of the complete analysis pipeline. Each directory contains a dedicated `README.md` describing the exact scripts required to reproduce the analysis for that specific component.

The recommended execution order of the repository is the following:

1. **Shotgun metagenomics (`Metagenomics/`)**  
   The workflow begins with the reconstruction of microbial genomes from shotgun metagenomic sequencing data. This pipeline performs quality control, host read removal, metagenome assembly, genome binning, and genome refinement to obtain a catalog of metagenome-assembled genomes (MAGs). Functional annotation and taxonomic classification of these genomes are also generated in this step.

2. **Metatranscriptomics (`Metatranscriptomics/`)**  
   After microbial genomes have been reconstructed, metatranscriptomic sequencing data are processed to quantify microbial gene expression. This pipeline performs quality filtering, host removal, rRNA depletion, transcriptome assembly, coding sequence prediction, and transcript quantification to generate microbial gene expression matrices.

3. **16S rRNA gene analysis (`16S_Metagenomics/`)**  
   Amplicon sequencing datasets from European and Latin American cohorts are processed to generate Amplicon Sequence Variants (ASVs) and taxonomic profiles. These results provide community-level microbial composition information that complements the genome-resolved metagenomic reconstruction.

4. **Metabolomics (`Metabolomics/`)**  
   Untargeted metabolomic datasets are analyzed to identify metabolites associated with Parkinson’s disease. Statistical and multivariate analyses are performed to define metabolite signatures that differentiate Parkinson’s disease patients from healthy controls.

5. **Mendelian Randomization (`MR/`)**  
   Mendelian Randomization analyses are conducted to evaluate potential causal relationships between microbiome-associated genetic variants and Parkinson’s disease risk. This step integrates microbiome GWAS summary statistics with Parkinson’s disease GWAS datasets from independent cohorts.

6. **Multi-omic integration (`Integration/`)**  
   The final stage integrates outputs generated by all previous pipelines. Microbial genomes reconstructed from metagenomics are linked to ASVs derived from 16S analysis, while microbial gene functions and transcriptional activity are integrated with metabolomic signatures. This step reconstructs functional relationships connecting microbial genomes, encoded genes, metabolic pathways, and host metabolites.

Together, these steps enable a systems-level reconstruction of microbiome–host interactions associated with Parkinson’s disease.

## Repository Structure

The repository is organized according to the different omic layers analyzed in this study and the final multi-omic integration steps. Each top-level directory corresponds to a specific analytical component of the workflow, and contains the scripts required to reproduce the analyses described in the manuscript.

The full directory structure of the repository is shown below.

```text
.
├── 16S_Metagenomics
│   ├── Classifier-training
│   │   ├── Curate
│   │   │   └── qiime_curate.sh
│   │   └── Download
│   │       └── qiime_download.sh
│   ├── EUR
│   │   ├── Denoising
│   │   │   ├── Cutadapt
│   │   │   │   └── qiime_cutadapt.sh
│   │   │   ├── DADA
│   │   │   │   └── qiime_denoise.sh
│   │   │   └── Summarize
│   │   │       └── qiime_summarize.sh
│   │   ├── Import
│   │   │   ├── Demux
│   │   │   │   └── qiime_demux.sh
│   │   │   └── qiime_import.sh
│   │   ├── Raw_data
│   │   │   └── sra_tools.sh
│   │   └── Taxonomy
│   │       ├── Classification
│   │       │   └── qiime_classification.sh
│   │       ├── Extract
│   │       │   └── qiime_extract.sh
│   │       ├── Filtering
│   │       │   └── qiime_filtering.sh
│   │       └── Train
│   │           └── qiime_train.sh
│   └── LATAM
│       ├── Denoising
│       │   ├── Cutadapt
│       │   │   └── qiime_cutadapt.sh
│       │   ├── DADA
│       │   │   └── qiime_denoise.sh
│       │   └── Summarize
│       │       └── qiime_summarize.sh
│       ├── Import
│       │   ├── Demux
│       │   │   └── qiime_demux.sh
│       │   └── qiime_import.sh
│       ├── Raw_data
│       │   └── sra-tools.sh
│       └── Taxonomy
│           ├── Classification
│           │   └── qiime_classification.sh
│           ├── Extract
│           │   └── qiime_extract.sh
│           ├── Filtering
│           │   └── qiime_filtering.sh
│           └── Train
│               └── qiime_train.sh
├── Integration
│   ├── MAG-ASV
│   │   └── MAG_ASV-integration.Rmd
│   └── MAG-genes-metabolite
│       └── Microbe-gene-metabolite.Rmd
├── MR
│   ├── Exposure
│   │   └── SNPs_exposure.R
│   ├── MR
│   │   └── MR.Rmd
│   └── Outcome
│       └── SNP_outcome.R
├── Metabolomics
│   └── Metabolomics_analysis.Rmd
├── Metagenomics
│   ├── Assembly
│   │   ├── MEGAHIT
│   │   │   ├── MEGAHIT.sh
│   │   │   └── MetaQuast
│   │   │       └── metaquast.sh
│   │   └── metaSPAdes
│   │       ├── MetaQuast
│   │       │   └── metaquast.sh
│   │       └── metaspades.sh
│   ├── Binning
│   │   ├── MEGAHIT
│   │   │   ├── Binner
│   │   │   │   ├── COMEBin
│   │   │   │   │   └── comebin.sh
│   │   │   │   ├── MetaBAT2
│   │   │   │   │   └── metabat2.sh
│   │   │   │   ├── SemiBin2
│   │   │   │   │   └── semibin2.sh
│   │   │   │   └── Vamb
│   │   │   │       └── vamb.sh
│   │   │   ├── CoverM
│   │   │   │   ├── coverm.sh
│   │   │   │   └── fix-vamb.py
│   │   │   └── Refinement
│   │   │       ├── DAS_Tool
│   │   │       │   └── das_tool.sh
│   │   │       └── Quality_control
│   │   │           ├── Clean_bin
│   │   │           │   ├── CheckM2
│   │   │           │   │   └── checkm2.sh
│   │   │           │   └── GUNC
│   │   │           │       └── gunc.sh
│   │   │           ├── Dereplication_filtering
│   │   │           │   └── dRep
│   │   │           │       └── drep.sh
│   │   │           └── Raw_bin
│   │   │               ├── CheckM2
│   │   │               │   └── checkm2.sh
│   │   │               └── GUNC
│   │   │                   └── gunc.sh
│   │   └── metaSPAdes
│   │       ├── Binner
│   │       │   ├── COMEBin
│   │       │   │   └── comebin.sh
│   │       │   ├── MetaBAT2
│   │       │   │   └── metabat2.sh
│   │       │   ├── SemiBin2
│   │       │   │   └── semibin2.sh
│   │       │   └── Vamb
│   │       │       └── vamb.sh
│   │       ├── CoverM
│   │       │   ├── coverm.sh
│   │       │   └── fix-vamb.py
│   │       └── Refinement
│   │           ├── DAS_Tool
│   │           │   └── das_tool.sh
│   │           └── Quality_control
│   │               ├── Clean_bin
│   │               │   ├── CheckM2
│   │               │   │   └── checkm2.sh
│   │               │   └── GUNC
│   │               │       └── gunc.sh
│   │               ├── Dereplication_filtering
│   │               │   └── dRep
│   │               │       └── drep.sh
│   │               └── Raw_bin
│   │                   ├── CheckM2
│   │                   │   └── checkm2.sh
│   │                   └── GUNC
│   │                       └── gunc.sh
│   ├── Features_MAGs
│   │   └── CoverM
│   │       └── coverm.sh
│   ├── Functional_annotation
│   │   ├── MAG_eggnog
│   │   │   └── merged_annotations.sh
│   │   └── eggmapper.sh
│   ├── Gene_prediction
│   │   └── GeneMarkS2
│   │       └── genemarks2.sh
│   ├── Host_removal
│   │   └── bowtie2.sh
│   ├── QC_MAGs
│   │   └── Final_MAGs
│   │       ├── CheckM2
│   │       │   └── checkm2.sh
│   │       ├── GUNC
│   │       │   └── gunc.sh
│   │       └── dRep
│   │           └── drep.sh
│   ├── Quality_control
│   │   ├── FastP
│   │   │   └── fastp.sh
│   │   ├── FastQC_clean
│   │   │   └── fastqc_clean.sh
│   │   └── FastQC_raw
│   │       └── fastqc_raw.sh
│   ├── Raw_data
│   │   └── sra-tools.sh
│   └── Taxonomy
│       ├── GTDB
│       │   ├── DB
│       │   └── taxonomy.sh
│       └── Phylogenomics
│           └── phylo.sh
└── Metatranscriptomics
    ├── Assembly
    │   └── rnaSPAdes
    │       ├── rnaQuast
    │       │   └── rnaquast.sh
    │       └── rnaspades.sh
    ├── Expression
    │   └── Meta-PD_DESeq2.Rmd
    ├── Host_removal
    │   └── STAR
    │       ├── Data
    │       │   └── symlinks.sh
    │       └── STAR.sh
    ├── Quality_control
    │   ├── FastP
    │   │   └── fastp.sh
    │   ├── FastQC_clean
    │   │   └── fastqc_clean.sh
    │   └── FastQC_raw
    │       └── fastqc_raw.sh
    ├── Quantification
    │   └── Salmon
    │       ├── Counts
    │       │   ├── import.sh
    │       │   └── tximport_salmon.R
    │       ├── Data
    │       │   └── symlinks.sh
    │       ├── MetaT
    │       └── meta_salmon.sh
    ├── Raw_data
    │   └── sra-tools.sh
    ├── Transcript_prediction
    │   ├── Diamond
    │   │   ├── DB
    │   │   │   └── MAG-protein_copy.sh
    │   │   └── magxtrans.sh
    │   └── TransDecoder
    │       └── transdecoder.sh
    └── rRNA_removal
        └── sortmerna.sh
```

## Workflow Description

The analytical framework implemented in this repository integrates multiple omic layers to characterize the relationship between the gut microbiome and Parkinson’s disease. Each omic dataset was processed using dedicated pipelines designed to ensure reproducibility and methodological consistency across cohorts. The workflow combines genome-resolved microbiome reconstruction, transcriptional profiling, metabolomic analysis, and causal inference approaches to obtain a systems-level view of microbiome–host interactions.

### Mendelian Randomization (MR)

Mendelian randomization analysis was performed to evaluate potential causal relationships between gut microbiome taxa and Parkinson’s disease risk. Genetic instruments associated with microbial taxa were obtained from large-scale microbiome genome-wide association studies. These variants were then tested against Parkinson’s disease GWAS summary statistics from independent cohorts.

The analysis includes multiple estimators to ensure robust causal inference, including inverse-variance weighted (IVW), weighted median, weighted mode, and Wald ratio methods. Sensitivity analyses were also conducted to evaluate pleiotropy, heterogeneity among instruments, and the stability of causal estimates.

### Shotgun Metagenomics

Shotgun metagenomic sequencing datasets were processed to reconstruct the genomic composition of the gut microbiome at genome resolution. After downloading raw sequencing reads from public repositories, reads were subjected to quality filtering and host contamination removal. Clean reads were assembled using MEGAHIT and metaSPAdes to generate metagenomic contigs.

Genome binning was performed using an ensemble strategy integrating multiple binning algorithms, including MetaBAT2, COMEBin, SemiBin2, and VAMB. Resulting bins were refined using DAS Tool and evaluated using CheckM2 and GUNC to assess completeness and contamination. High-quality genomes were dereplicated using dRep to obtain a non-redundant set of metagenome-assembled genomes (MAGs). These genomes were then taxonomically classified using the Genome Taxonomy Database and functionally annotated using eggNOG-based pipelines.

### Metatranscriptomics

Metatranscriptomic sequencing data were analyzed to characterize the functional activity of the gut microbiome. Raw RNA-seq reads were first quality-filtered and subjected to host read removal using genome alignment approaches. Ribosomal RNA sequences were removed to enrich for messenger RNA.

The remaining reads were assembled using RNA-SPAdes to reconstruct transcripts. Open reading frames were predicted using TransDecoder, and functional annotation was performed using DIAMOND against protein databases derived from reconstructed microbial genomes. Gene expression levels were quantified using Salmon, and differential expression analysis was conducted using DESeq2 to identify genes with significant transcriptional changes between Parkinson’s disease and control samples.

### 16S rRNA Gene Analysis

Amplicon sequencing datasets from European and Latin American cohorts were processed using QIIME2. Raw reads were imported, demultiplexed, and denoised using the DADA2 algorithm to generate high-resolution amplicon sequence variants (ASVs).

Taxonomic classification was performed using a Naïve Bayes classifier trained on curated reference sequences derived from the SILVA database. Downstream analyses included microbial community composition profiling, diversity analyses, and identification of taxa differentially abundant between Parkinson’s disease and healthy control groups.

### Metabolomics

Untargeted metabolomic datasets were analyzed to identify metabolites associated with Parkinson’s disease. Raw mass spectrometry data were processed for peak detection, alignment, and metabolite annotation. Statistical analyses were then performed to identify metabolites showing significant differences between disease and control groups.

Machine learning approaches, including sparse partial least squares discriminant analysis (sPLS-DA) and random forest models, were used to prioritize metabolites with strong discriminatory power. The final metabolite signature was defined by intersecting significant metabolites identified across statistical methods.

### Multi-omic Integration

To connect microbial genomic, taxonomic, transcriptional, and metabolic signals, a multi-layer integration strategy was implemented. Amplicon sequence variants obtained from 16S rRNA analysis were linked to reconstructed metagenome-assembled genomes through sequence-based taxonomic harmonization.

Functional annotations derived from microbial genomes and metatranscriptomic expression profiles were then integrated with metabolomic signatures using correlation-based analyses. This approach enabled the reconstruction of functional relationships linking microbial genomes, encoded genes, metabolic pathways, and host metabolites, providing a systems-level view of microbiome–host interactions associated with Parkinson’s disease.

## Software and Dependencies

The analyses implemented in this repository rely on a combination of widely used bioinformatics tools and statistical software for microbiome, transcriptomic, genomic, and metabolomic data processing. Most command-line workflows were executed in a Linux environment using Bash scripts, while downstream statistical analyses and multi-omic integration were performed in R.

Major software tools used across the different pipelines include FastQC and fastp for sequencing quality control, Bowtie2 and STAR for host read removal, MEGAHIT and metaSPAdes for metagenomic assembly, and RNA-SPAdes for transcriptome assembly. Genome binning was performed using MetaBAT2, COMEBin, SemiBin2, and VAMB, followed by bin refinement with DAS Tool. Genome quality was assessed using CheckM2 and GUNC, and genome dereplication was performed using dRep. Taxonomic classification of metagenome-assembled genomes was performed using the Genome Taxonomy Database toolkit (GTDB-Tk).

Functional annotation of predicted genes was performed using eggNOG-based annotation pipelines. Metatranscriptomic reads were quantified using Salmon and differential expression analyses were conducted using DESeq2 in R. Amplicon sequencing data were processed using QIIME2 with denoising performed by the DADA2 algorithm. Statistical and machine learning analyses for metabolomics were performed in R using established multivariate analysis packages.

All scripts included in this repository correspond to the exact commands used to generate the results described in the manuscript.

## Data Availability

All sequencing datasets used in this study were obtained from publicly available repositories. Raw sequencing reads and associated metadata can be accessed through the following accession numbers:

- **Villette et al. dataset** (shotgun metagenomics and metatranscriptomics): NCBI BioProject PRJNA782492  
- **Pereira et al. dataset** (16S rRNA sequencing and metabolomics): European Nucleotide Archive PRJEB27564  
- **Forero-Rodríguez et al. dataset** (16S rRNA sequencing, Latin American cohort): NCBI BioProject PRJNA975118  

Genome-wide association summary statistics used for Mendelian randomization analyses were obtained from publicly available GWAS resources, including MiBioGen, FinnGen, GP2, and LARGE-PD.

Due to repository size limitations, raw sequencing datasets and intermediate large files are not included in this GitHub repository. The scripts provided here allow full reproduction of the analysis once the original datasets are downloaded from their respective repositories.

## Reproducibility Notes

All analyses in this repository are organized as modular pipelines corresponding to the different omic layers included in the study. Each directory contains scripts representing the commands used at each stage of the workflow, from raw data acquisition to downstream statistical analysis.

The scripts were designed to be executed sequentially within each pipeline. Users attempting to reproduce the analyses should ensure that the required software tools are installed and that reference databases (such as GTDB and SILVA) are available locally.

Paths to input data, reference databases, and output directories may need to be adapted to match the local computational environment. The repository is intended to document the exact analytical workflow used in the study and facilitate reproducibility of the results.

## Citation
