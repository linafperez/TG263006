# ===============================================================
# STEP 0: Clear workspace
# ===============================================================
rm(list = ls())

# ===============================================================
# STEP 1: Load libraries
# ===============================================================
library(tidyverse)
library(data.table)
library(stringr)
library(TwoSampleMR)
library(MRInstruments)
library(ieugwasr)

# ===============================================================
# STEP 2: Define paths
# ===============================================================
base_dir <- "/opt/data/HPC01A/alexis_rojasc2/Alexis/GWAS/Outcome"
output_dir <- file.path(base_dir, "results")
ref_geno_dir <- file.path(base_dir, "Genotype")
tmp_dir <- file.path(output_dir, "plink_tmp")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)

# ===============================================================
# STEP 3: List input .gz files
# ===============================================================
gz_files <- c(
  file.path(base_dir, "GWAS_EUR-GP2_rsid.txt.gz"),
  file.path(base_dir, "GWAS_EUR-FINN.gz"),
  file.path(base_dir, "GWAS_Latin_rsid.tsv.gz")
)

# ===============================================================
# STEP 4: Palindromic SNP checker
# ===============================================================
is_palindromic <- function(a1, a2) {
  (a1 == "A" & a2 == "T") | (a1 == "T" & a2 == "A") |
    (a1 == "C" & a2 == "G") | (a1 == "G" & a2 == "C")
}

# ===============================================================
# STEP 5: Processing function
# ===============================================================
process_outcome <- function(file_path) {
  message("Read file: ", file_path)
  
  region_tag <- if (grepl("Latin", file_path)) "LATAM" else "EUR"
  pop_use <- if (region_tag == "LATAM") "AMR" else "EUR"
  
  df <- fread(cmd = paste("zcat", shQuote(file_path)))
  
  # Standardize columns by file
  if (grepl("GP2", file_path)) {
    df <- df %>%
      rename(
        SNP = rsID,
        effect_allele = effect_allele,
        other_allele = other_allele,
        EAF = effect_allele_frequency,
        pval = p_value,
        beta = beta,
        SE = standard_error,
        chr = chromosome,
        bp = base_pair_position
      )
  } else if (grepl("FINN", file_path)) {
    df <- df %>%
      rename(
        SNP = rsids,
        effect_allele = alt,
        other_allele = ref,
        EAF = af_alt,
        pval = pval,
        beta = beta,
        SE = sebeta,
        chr = "#chrom",
        bp = pos
      )
  } else if (grepl("Latin", file_path)) {
    df <- df %>%
      rename(
        SNP = rsID,
        effect_allele = alt.x,
        other_allele = ref.x,
        EAF = freq,
        pval = Score.pval,
        beta = beta,
        SE = Score.SE,
        chr = chr.x,
        bp = pos.x
      )
  } else {
    stop("Unknown file structure: ", file_path)
  }
  
  df <- df %>% filter(pval < 1e-5)
  message("SNPs with p < 1e-5: ", nrow(df))
  if (nrow(df) == 0) return(NULL)
  
  df <- df %>% distinct(SNP, .keep_all = TRUE)
  
  outcome_dat <- df %>%
    transmute(
      SNP,
      beta.outcome = beta,
      se.outcome = SE,
      effect_allele.outcome = effect_allele,
      other_allele.outcome = other_allele,
      eaf.outcome = EAF,
      pval.outcome = pval,
      chr.outcome = chr,
      pos.outcome = bp,
      id.outcome = "PD",
      outcome = "PD"
    )
  
  # PLINK clumping
  ref_path <- file.path(ref_geno_dir, pop_use, paste0("1KG_", pop_use))
  clump_out <- file.path(tmp_dir, paste0("clumped_", pop_use))
  assoc_path <- file.path(tmp_dir, paste0("assoc_", pop_use, ".txt"))
  clumped_file <- paste0(clump_out, ".clumps")
  
  if (!file.exists(clumped_file)) {
    fwrite(
      outcome_dat %>%
        separate_rows(SNP, sep = ",") %>%
        filter(!is.na(SNP) & SNP != ".") %>%
        group_by(SNP) %>%
        slice_min(order_by = pval.outcome, n = 1, with_ties = FALSE) %>%
        ungroup() %>%
        select(SNP, pval.outcome) %>%
        rename(ID = SNP, P = pval.outcome),
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
    
    message("Running PLINK2 for population: ", pop_use)
    system(cmd)
  }
  
  if (!file.exists(clumped_file)) {
    message("Clumped file not found: ", clumped_file)
    return(NULL)
  }
  
  clumped_data <- fread(clumped_file)
  id_col <- if ("ID" %in% colnames(clumped_data)) "ID" else if ("SNP" %in% colnames(clumped_data)) "SNP" else NA
  if (is.na(id_col)) return(NULL)
  
  snps <- clumped_data[[id_col]]
  clumped_df <- outcome_dat %>% filter(SNP %in% snps)
  message("SNPs after clumping [", pop_use, "]: ", length(snps))
  
  # ===============================================================
  # STEP 7: Remove palindromic SNPs
  # ===============================================================
  out_clumped <- clumped_df %>%
    rowwise() %>%
    filter(!is_palindromic(effect_allele.outcome, other_allele.outcome)) %>%
    ungroup()
  message("SNPs after palindromic filter: ", nrow(out_clumped))
  if (nrow(out_clumped) == 0) return(NULL)
  
  # ===============================================================
  # STEP 8: Filter by F-statistic
  # ===============================================================
  out_clumped <- out_clumped %>%
    mutate(F_stat = (beta.outcome / se.outcome)^2) %>%
    filter(F_stat >= 10) %>%
    filter(!is.na(eaf.outcome) & eaf.outcome > 0.01)
  message("SNPs after F and MAF filter: ", nrow(out_clumped))
  if (nrow(out_clumped) == 0) return(NULL)
  
  # ===============================================================
  # STEP 9: Save results
  # ===============================================================
  out_clumped <- out_clumped %>% mutate(Region = region_tag)
  out_file <- file.path(output_dir, paste0("instruments_", basename(file_path)))
  fwrite(out_clumped, out_file, sep = "\t")
  
  # ===============================================================
  # STEP 10: Clean up PLINK tmp files (once per .gz file)
  # ===============================================================
  tmp_files <- list.files(tmp_dir, full.names = TRUE, pattern = "\\.(clumps|log|txt)$")
  file.remove(tmp_files)
  
  return(out_clumped)
}

# ===============================================================
# STEP 11: Process all input files
# ===============================================================
results_list <- lapply(gz_files, process_outcome)
results_list <- Filter(Negate(is.null), results_list)

if (length(results_list) == 0) stop("No valid instruments found.")

# ===============================================================
# STEP 12: Combine and save all instruments
# ===============================================================
all_instruments <- bind_rows(results_list)
combined_output <- file.path(output_dir, "all_instruments_combined.txt")
fwrite(all_instruments, combined_output, sep = "\t")

cat("Combined instrument SNPs saved to:", combined_output, "\n")

