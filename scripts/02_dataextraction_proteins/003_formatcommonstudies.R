# Date: 05-06-2025
# Author: A.L.Hanson
# Purpose:
## Extract finemapped exposure variants (pQTLs) - no clumping required
## Harmonise OpenGWAS outcome summary statistics with Sun et al. protein exposure pQTLs

args <- commandArgs(trailingOnly = TRUE)
# Test args:
#args <- c("sun_pQTL_ABO_OID30675.tsv","GCST90013977_pQTLs_wproxies.vcf.bgz","ABO","RBCcount","sun_gwas","opengwas","variant")
#args <- c("sun_pQTL_APOE_OID30727.tsv","ieu-b-110_pQTLs_wproxies.vcf.bgz","APOE","LDL","sun_gwas","opengwas","variant")
#args <- c("sun_pQTL_AHSP_OID21078.tsv","GCST90014006_pQTLs_wproxies.vcf.bgz","AHSP","HbA1c","sun_gwas","opengwas","variant")

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
ld_dir <- Sys.getenv("ld_dir")

out_dir <- file.path(data_dir,"harmonised","proteomics")

message("Formatting studies:")
message("Exposure: ", exposure_name, " (", exposure_study, ")\n",
  "Outcome: ", outcome_name, " (", outcome_study, ")\n",
  "Source: ", exposure_source, "\nClass: ", class)

main <- function(exposure_study, outcome_study){
  # Exposure pQTLs
  exposure_dat <- data.table::fread(file.path(data_dir, "sumstats/proteomics/common", exposure_study), data.table = FALSE)
  message("n snps: ", nrow(exposure_dat), " n independent hits: ", sum(exposure_dat$finemapped_hit))

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

    # Add cis-trans information and chr:bp_A1_A2 ID
    exposure_formatted$cis_trans = exposure_dat[match(exposure_formatted$SNP, exposure_dat$VEP_RSID),]$cis_trans
    exposure_formatted$variant = exposure_dat[match(exposure_formatted$SNP, exposure_dat$VEP_RSID),]$RSID
    exposure_formatted$finemapped_hit = exposure_dat[match(exposure_formatted$SNP, exposure_dat$VEP_RSID),]$finemapped_hit
  
    # Generate two exposure SNP sets, one with finemapped variants, one clumped and pruned
    exposure_finalset <- list()
    # Extract finemapped variants
    exposure_finalset[[1]] <- exposure_formatted |> dplyr::filter(finemapped_hit == TRUE)
    # Perform clumping
    exposure_finalset[[2]] <- perform_clumping(exposure_df = exposure_formatted)
    names(exposure_finalset) <- c("exposure_common_finemapped", "exposure_common_clumped")

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
    dat_harmonised <- list()
    dat_harmonised[[1]] <- harmonise_data(exposure_dat = exposure_finalset$exposure_common_finemapped, outcome_dat = outcome_formatted) 
    dat_harmonised[[2]] <- harmonise_data(exposure_dat = exposure_finalset$exposure_common_clumped, outcome_dat = outcome_formatted)
    names(dat_harmonised) <- c("exposure_common_finemapped", "exposure_common_clumped")

    message("Final instrument count: ")
    print(unlist(lapply(dat_harmonised, nrow)))
  }

  # Write out harmonised studies
  message("Saving...")
  outname <- paste0(paste(exposure_source, class, exposure_name, outcome_name, sep = "_"),".rda")
  save(dat_harmonised, file = file.path(out_dir, outname))
}

# Functions
# Fixing bugs ieugwasr::ld_clump_local()
ld_clump_fix <- function (dat, clump_kb, clump_r2, clump_p, bfile, plink_bin) {
  shell <- ifelse(Sys.info()["sysname"] == "Windows", "cmd", "sh")
  fn <- tempfile()
  write.table(data.frame(SNP = dat[["rsid"]], P = dat[["pval"]]), 
              file = fn, row.names = FALSE, col.names = TRUE, quote = FALSE)
  fun2 <- paste0(shQuote(plink_bin, type = shell), " --bfile ", 
                 shQuote(bfile, type = shell), " --clump ", 
                 shQuote(fn, type = shell), " --clump-p1 ", clump_p, " --clump-r2 ", 
                 clump_r2, " --clump-kb ", clump_kb, " --out ", 
                 shQuote(fn, type = shell))
  system(fun2)
  res <- read.table(paste(fn, ".clumps", sep = ""), header = TRUE, comment.char = "")
  unlink(paste(fn, "*", sep = ""))
  y <- subset(dat, !dat[["rsid"]] %in% res[["ID"]])
  if (nrow(y) > 0) {
    message("Removing ", length(y[["rsid"]]), " of ", nrow(dat), 
            " variants due to LD with other variants or absence from LD reference panel")
  }
  return(subset(dat, dat[["rsid"]] %in% res[["ID"]]))
}

perform_clumping <- function(exposure_df){
  
  # Order SNP alleles alphabetically and capitalise to match reference
  df_in <- exposure_df 

  clumped_df <- ld_clump_fix(tibble(rsid=df_in$SNP, pval=df_in$pval.exposure, id=df_in$exposure),
                             clump_p = 1,
                             clump_r2 = 0.001,
                             clump_kb = 10000,
                             bfile = file.path(ld_dir, "full_rsid"),
                             plink_bin = "/usr/local/bin/plink2")
  
  positions_keep <- clumped_df$rsid
  
  return(subset(exposure_df, exposure_df$SNP %in% positions_keep))
}

main(exposure_study = exposure_study, outcome_study = outcome_study)