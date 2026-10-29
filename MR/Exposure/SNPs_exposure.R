# Clear workspace
rm(list = ls())

# Load libraries
library(tidyverse)
library(data.table)
library(stringr)
library(TwoSampleMR)
library(MRInstruments)
library(ieugwasr)

# ===============================================================
# STEP 1: Define paths
# ===============================================================
base_dir <- "/opt/data/HPC01A/alexis_rojasc2/Alexis/GWAS/Exposure"
zip_file <- file.path(base_dir, "MiBioGen_QmbQTL_summary_genus.zip")
output_dir <- file.path(base_dir, "results")
unzipped_dir <- file.path(base_dir, "unzipped_mibiogen")
ref_geno_dir <- file.path(base_dir, "Genotype")
tmp_dir <- file.path(output_dir, "plink_tmp")

# Create output directories
invisible(lapply(c(output_dir, unzipped_dir, tmp_dir), dir.create, showWarnings = FALSE, recursive = TRUE))

# ===============================================================
# STEP 2: List input .gz files
# ===============================================================
gz_files <- list.files(unzipped_dir, pattern = "\\.gz$", full.names = TRUE)

# ===============================================================
# STEP 3: Palindromic SNP checker function
# ===============================================================
is_palindromic <- function(a1, a2) {
  (a1 == "A" & a2 == "T") | (a1 == "T" & a2 == "A") |
    (a1 == "C" & a2 == "G") | (a1 == "G" & a2 == "C")
}

# ===============================================================
# STEP 4: Function to process each genus file
# ===============================================================
process_mibiogen_genus <- function(file_path) {
  message("Read file: ", file_path)
  
  df <- fread(cmd = paste("zcat", shQuote(file_path)))
  
  setnames(df,
           old = c("rsID", "ref.allele", "eff.allele", "P.weightedSumZ"),
           new = c("SNP", "other_allele", "effect_allele", "pval"))
  
  df <- df %>% filter(pval < 1e-5)
  message("SNPs with p < 1e-5: ", nrow(df))
  if (nrow(df) == 0) return(NULL)
  
  df <- df %>% distinct(SNP, .keep_all = TRUE)
  
  exposure_dat <- df %>%
    transmute(
      SNP,
      beta.exposure = beta,
      se.exposure = SE,
      effect_allele.exposure = effect_allele,
      other_allele.exposure = other_allele,
      eaf.exposure = NA,
      pval.exposure = pval,
      samplesize.exposure = N,
      chr.exposure = chr,
      pos.exposure = bp,
      id.exposure = bac,
      exposure = bac
    )
  
  # ===============================================================
  # STEP 5: Use PLINK2 for clumping per population
  # ===============================================================
  clumped_all <- list()
  populations <- c("EUR", "SAS", "EAS", "AFR", "AMR")
  
  for (pop in populations) {
    ref_path <- file.path(ref_geno_dir, pop, paste0("1KG_", pop))
    clump_out <- file.path(tmp_dir, paste0("clumped_", pop))
    assoc_path <- file.path(tmp_dir, paste0("assoc_", pop, ".txt"))
    clumped_file <- paste0(clump_out, ".clumps")
    
    # Remove invalid SNPs and duplicates before writing assoc file
    if (!file.exists(clumped_file)) {
      fwrite(
        exposure_dat %>%
          separate_rows(SNP, sep = ",") %>% 
          filter(!is.na(SNP) & SNP != ".") %>%
          group_by(SNP) %>%
          slice_min(order_by = pval.exposure, n = 1, with_ties = FALSE) %>%
          ungroup() %>%
          select(SNP, pval.exposure) %>%
          rename(ID = SNP, P = pval.exposure),
          file = assoc_path, sep = "\t", quote = FALSE
      )
      
      cmd <- paste(
        "plink2",
        "--pfile", shQuote(ref_path),
        "--clump", shQuote(assoc_path),
        "--clump-snp-field ID",
        "--clump-field P",
        "--clump-p1 1 --clump-p2 1 --clump-r2 0.001 --clump-kb 10000",
        "--out", shQuote(clump_out)
      )
      
      message("Running PLINK2 for population: ", pop)
      system(cmd)
    } else {
      message("Using existing clumped file for population: ", pop)
    }
    
    if (file.exists(clumped_file)) {
      clumped_data <- fread(clumped_file)
      
      id_col <- if ("ID" %in% colnames(clumped_data)) "ID" else if ("SNP" %in% colnames(clumped_data)) "SNP" else NA
      if (is.na(id_col)) {
        message("Missing 'ID' or 'SNP' column in: ", clumped_file)
        next
      }
      
      snps <- clumped_data[[id_col]]
      clumped_all[[pop]] <- exposure_dat %>% filter(SNP %in% snps)
      message("SNPs after clumping [", pop, "]: ", length(snps))
    } else {
      message("Clumped file not found: ", clumped_file)
    }
  }
  
  # Merge all population results
  if (length(clumped_all) == 0) return(NULL)
  exp_clumped <- bind_rows(clumped_all, .id = "population") %>% distinct(SNP, .keep_all = TRUE)
  
  # ===============================================================
  # STEP 7: Remove palindromic SNPs
  # ===============================================================
  exp_clumped <- exp_clumped %>%
    rowwise() %>%
    filter(!is_palindromic(effect_allele.exposure, other_allele.exposure)) %>%
    ungroup()
  
  message("SNPs after palindromic filter: ", nrow(exp_clumped))
  if (nrow(exp_clumped) == 0) return(NULL)
  
  # ===============================================================
  # STEP 8: Filter by F-statistic
  # ===============================================================
  exp_clumped <- exp_clumped %>%
    mutate(F_stat = (beta.exposure / se.exposure)^2) %>%
    filter(F_stat >= 10)
  
  message("SNPs after F filter: ", nrow(exp_clumped))
  if (nrow(exp_clumped) == 0) return(NULL)
  
  # ===============================================================
  # STEP 9: Save results
  # ===============================================================
  taxon_file <- gsub("\\.summary\\.txt\\.gz$", "", basename(file_path))
  out_path <- file.path(output_dir, paste0("instruments_", taxon_file, ".txt"))
  fwrite(exp_clumped, out_path, sep = "\t")
  
  # ===============================================================
  # STEP 10: Clean up PLINK tmp files (once per .gz file)
  # ===============================================================
  tmp_files <- list.files(tmp_dir, full.names = TRUE, pattern = "\\.(clumps|log|txt)$")
  file.remove(tmp_files)
  
  return(exp_clumped)
}

# ===============================================================
# STEP 11: Process all .gz files
# ===============================================================
results_list <- lapply(gz_files, process_mibiogen_genus)
results_list <- Filter(Negate(is.null), results_list)

if (length(results_list) == 0) stop("No valid instruments found.")

# ===============================================================
# STEP 12: Combine and save all instruments
# ===============================================================
all_instruments <- bind_rows(results_list) %>%
  mutate(taxon = str_extract(id.exposure, "(?<=genus\\.)[^\\.]+"))

combined_output <- file.path(output_dir, "all_instruments_combined.txt")
fwrite(all_instruments, combined_output, sep = "\t")

cat("Combined instrument SNPs saved to:", combined_output, "\n")

