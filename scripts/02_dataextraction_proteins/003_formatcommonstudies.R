# Date: 05-06-2025
# Author: A.L.Hanson
# Purpose: Harmonise OpenGWAS outcome summary statistics with Sun et al. protein exposure pQTLs
# Exposure variants are finemapped and do not require clumping

args <- commandArgs(trailingOnly = TRUE)
#args <- c("sun_pQTL_AAMDC_OID30236.tsv","ebi-a-GCST90000615_pQTLs_wproxies.vcf.bgz","AAMDC","VitD","sun_gwas","opengwas","variant")

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

out_dir <- file.path(data_dir,"harmonised","proteomics")

message("Formatting studies:")
message("Exposure: ", exposure_name, " (", exposure_study, ")\n",
  "Outcome: ", outcome_name, " (", outcome_study, ")\n",
  "Source: ", exposure_source, "\nClass: ", class)

main <- function(exposure_study, outcome_study){
  # Exposure pQTLs
  exposure_dat <- data.table::fread(file.path(data_dir, "sumstats/proteomics/common", exposure_study), data.table = FALSE)
  message("n instruments: ", nrow(exposure_dat))

  if(nrow(exposure_dat) == 0){
    message("No instruments for protein")
    harmonised_dat <- data.frame()
  }else{
    # Convert to TwoSampleMR format
    message("Formatting exposure data...")
    exposure_formatted <- TwoSampleMR::format_data(
      exposure_dat,
      type = "exposure",
      phenotype_col = "protein_target",
      snp_col = "VEP_RSID",
      beta_col = "BETA",
      se_col = "SE",
      eaf_col = "EAF",
      effect_allele_col = "EA",
      other_allele_col = "OA",
      pval_col = "P")

    # Add cis-trans information  
    exposure_formatted$cis_trans = exposure_dat[match(exposure_formatted$SNP, exposure_dat$VEP_RSID),]$cis_trans

    # Extract instruments from outcome
    message("Formatting outcome data...")
    outcome_dat <- VariantAnnotation::readVcf(file.path(data_dir, "sumstats/opengwas", outcome_study)) |>
      gwasvcf::vcf_to_granges() |> dplyr::as_tibble() |>
      dplyr::mutate(P = 10^-LP)

    outcome_formatted <- TwoSampleMR::format_data(
      outcome_dat,
      type = "outcome",
      snp_col = "ID",
      beta_col = "ES",
      se_col = "SE",
      eaf_col = "AF",
      effect_allele_col = "ALT",
      other_allele_col = "REF",
      pval_col = "P")
    outcome_formatted$outcome <- outcome_name

    # Harmonise SNP data
    dat_harmonised <- harmonise_data(exposure_dat = exposure_formatted, outcome_dat = outcome_formatted) 
    
    message("Final instrument count:", nrow(dat_harmonised))
  }

  # Write out harmonised studies
  message("Saving...")
  outname <- paste0(paste(exposure_source, class, exposure_name, outcome_name, sep = "_"),".rda")
  save(dat_harmonised, file = file.path(out_dir, outname))
}

main(exposure_study = exposure_study, outcome_study = outcome_study)