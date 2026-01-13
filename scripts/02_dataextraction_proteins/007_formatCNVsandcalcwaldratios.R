# Date: 12-01-2026
# Author: A.L.Hanson
# Purpose: Harmonise Milind et al. pLOF and gene duplication associations for available proteins and complex traits 
# and calculate Wald Ratios

# Output: 
# Harmonised summary statistics with Wald ratios for each protein-->complex trait and cnv mask, saved in /data/harmonised/cnvs:
# E.g. milind_duplications_mask_[protein]_[outcome].rda

args <- commandArgs(trailingOnly = TRUE)
# Test args:
# args <- c("cardiometabolic/acan.deletions.WES.regenie.gz","p102.deletions.WES.regenie.gz","ACAN","Pulse_rate__automated_reading","milind_deletions","milind_deletions","mask")
# args <- c("cardiometabolic/pcsk9.plofs.WES.regenie.gz","p30760.plofs.WES.regenie.gz","PCSK9","HDL_cholesterol","milind_plofs","milind_plofs","mask")

exposure_study <- args[1]
outcome_study <- args[2]
exposure_name <- args[3]
outcome_name <- args[4]
exposure_source <- args[5]
outcome_source <- args[6]
class <- args[7]

library(dotenv)
library(TwoSampleMR)
library(dplyr)

# Load environment variables
dotenv::load_dot_env(file = "config.env")
data_dir <- Sys.getenv("data_dir")
res_dir <- Sys.getenv("results_dir")
prot_dir <- file.path(data_dir,"sumstats","Milind_CNV_associations")
project_dir <- Sys.getenv("project_dir")

out_dir <- file.path(data_dir,"harmonised","cnvs")

message("Formatting studies:")
message("Exposure: ", exposure_name, " (", exposure_study, ")\n",
  "Outcome: ", outcome_name, " (", outcome_study, ")\n",
  "Source: ", exposure_source, "\nClass: ", class)

main <- function(exposure_study, outcome_study){
  # Exposure pQTLs
  exposure_dat <- data.table::fread(file.path(data_dir, "sumstats/Milind_CNV_associations/CNV_Burden_Tests_Proteins", exposure_study), data.table = FALSE)
  # Filter to 1kb padded test (for deletions and duplications)
  if(exposure_source %in% c("milind_duplications","milind_deletions")){
    exposure_dat <- exposure_dat |> dplyr::filter(ALLELE1 == "Pad1KB.0.01")
  }
  message("n genes: ", nrow(exposure_dat))

  # Gene mappings
  genes <- data.table::fread(file.path(data_dir, "ensembl_geneIDmapping.txt"), data.table = F)
  
  if(nrow(exposure_dat) == 0){
    message("No features")
    harmonised_dat <- data.frame()
  }else{
    # Format exposure data
    exposure_dat <- exposure_dat |> dplyr::mutate(gene_ensembl = sub("\\..*", "", ID))
    exposure_dat <- exposure_dat |> dplyr::mutate(
      gene_symbol = genes[match(exposure_dat$gene_ensembl, genes$ensembl_gene_id), "hgnc_symbol"],
      protein_target = exposure_name,
      pval = 10^-LOG10P,
      effect_allele = "A", # Dummy column needed for TwoSampleMR
      other_allele = "C") # Dummy column needed for TwoSampleMR
    exposure_dat <- exposure_dat |> dplyr::filter(!(gene_symbol == ""))

    # Convert to TwoSampleMR format
    message("Formatting exposure data...")
    exposure_formatted <- TwoSampleMR::format_data(
      exposure_dat,
      type = "exposure",
      phenotype_col = "protein_target",
      chr_col= "CHROM",
      pos_col = "GENPOS",
      snp_col = "gene_symbol",
      beta_col = "BETA",
      se_col = "SE",
      eaf_col = "A1FREQ",
      effect_allele_col = "effect_allele",
      other_allele_col = "other_allele",
      pval_col = "pval")

    # Add cis-trans information
    exposure_formatted$cis_trans = ifelse(toupper(exposure_formatted$SNP) == exposure_formatted$exposure, "cis", "trans")
  
    # Extract instruments from outcome
    message("Formatting outcome data...")
    outcome_dat <- data.table::fread(file.path(data_dir, "sumstats/Milind_CNV_associations/CNV_Burden_Tests", outcome_study), data.table = FALSE)
  # Filter to 1kb padded test (for deletions and duplications)
    if(outcome_source %in% c("milind_duplications","milind_deletions")){
      outcome_dat <- outcome_dat |> dplyr::filter(ALLELE1 == "Pad1KB.0.01")
    }

    # Format outcome data
    outcome_dat <- outcome_dat |> dplyr::mutate(gene_ensembl = sub("\\..*", "", ID))
    outcome_dat <- outcome_dat |> dplyr::mutate(
      gene_symbol = genes[match(outcome_dat$gene_ensembl, genes$ensembl_gene_id), "hgnc_symbol"],
      outcome = outcome_name,
      pval = 10^-LOG10P,
      effect_allele = "A", # Dummy column needed for TwoSampleMR
      other_allele = "C") # Dummy column needed for TwoSampleMR
    outcome_dat <- outcome_dat |> dplyr::filter(!(gene_symbol == ""))

    outcome_formatted <- TwoSampleMR::format_data(
      outcome_dat,
      type = "outcome",
      phenotype_col = "outcome",
      chr_col= "CHROM",
      pos_col = "GENPOS",
      snp_col = "gene_symbol",
      beta_col = "BETA",
      se_col = "SE",
      eaf_col = "A1FREQ",
      effect_allele_col = "effect_allele",
      other_allele_col = "other_allele",
      pval_col = "pval")

    # Harmonise CNV data
    dat_harmonised <- list()
    dat_harmonised[[1]] <- harmonise_data(exposure_dat = exposure_formatted, outcome_dat = outcome_formatted)
    dat_harmonised[[1]]$test <- exposure_source
    names(dat_harmonised) <- exposure_source

    # Calculate Wald Ratios
    wald_ratios <- lapply(dat_harmonised, function(x){
      wrs <- apply(x, MARGIN = 1, FUN = function(x){
        TwoSampleMR::mr_wald_ratio(
        b_exp = as.numeric(x["beta.exposure"]),
        b_out = as.numeric(x["beta.outcome"]),
        se_exp = as.numeric(x["se.exposure"]),
        se_out = as.numeric(x["se.outcome"]))
      })

      wrs_df <- do.call("rbind",lapply(wrs, unlist))
      colnames(wrs_df) <- paste("wald_ratio", colnames(wrs_df), sep = "_")

      # Bind to harmonised data
      out <- cbind(x,wrs_df)
      return(out)
    })

  }

  # Write out harmonised studies
  message("Saving...")
  outname <- paste0(paste(exposure_source, class, exposure_name, outcome_name, sep = "_"),".rda")
  save(wald_ratios, file = file.path(out_dir, outname))
}

main(exposure_study = exposure_study, outcome_study = outcome_study)