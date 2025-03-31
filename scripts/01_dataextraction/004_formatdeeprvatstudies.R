# Date: 27-03-2025
# Author: A.L.Hanson
# Purpose: Format Deeprvat gene score associations with exposure and outcome variables for use in MR

args <- commandArgs(trailingOnly = TRUE)

exposure_study <- args[1]
outcome_study <- args[2]
exposure_name <- args[3]
outcome_name <- args[4]
source <- args[5]
class <- args[6]

library(data.table)
library(dotenv)
library(TwoSampleMR)
library(dplyr)
library(biomaRt)

# Load environment variables
dotenv::load_dot_env(file = "config.env")
data_dir <- Sys.getenv("data_dir")
ld_dir <- Sys.getenv("ld_dir")

out_dir <- file.path(data_dir,"harmonised")

message("Formatting studies:")
message("Exposure: ", exposure_name, " (", exposure_study, ")\n",
  "Outcome: ", outcome_name, " (", outcome_study, ")\n",
  "Source: ", source, "\nClass: ", class)

main <- function(exposure_study, outcome_study, class) {
  infile <- c(
    exposure = file.path(data_dir, "sumstats", source, exposure_study),
    outcome = file.path(data_dir, "sumstats", source, outcome_study)
  )

  infile_formatted <- list()

  # Update column names and add gene chromosome and position
  for(i in 1:length(infile)){
    message("Formatting: ", basename(infile[i]), " class: ", class)
    infile_formatted[[i]] <- format_sumstats(path = infile[i])
  }
  names(infile_formatted) <- c("exposure","outcome")

  # Genes (instruments) to keep
  genes_keep <- infile_formatted$exposure |> 
    filter(pval_corrected < 0.05) |> # Bonferroni adjusted pval threshold
    pull(Gene)

  # Format exposure
  message("Formatting exposure...")
  exposure_formatted <- TwoSampleMR::format_data(
    infile_formatted$exposure |> filter(Gene %in% genes_keep),
    type = "exposure",
    phenotype_col = "description",
    snp_col = "Gene",
    beta_col = "beta",
    se_col = "SE",
    eaf_col = "AF",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    pval_col = "pval")

  # Format outcome
  message("Formatting outcome...")
  outcome_formatted <- TwoSampleMR::format_data(
    infile_formatted$outcome,
    type = "outcome",
    phenotype_col = "description",
    snp_col = "Gene",
    beta_col = "beta",
    se_col = "SE",
    eaf_col = "AF",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    pval_col = "pval")

  # Harmonise
  dat_harmonised <- TwoSampleMR::harmonise_data(exposure_dat = exposure_formatted, outcome_dat = outcome_formatted)

  message("Final instrument count:")
  print(nrow(dat_harmonised))

  # Write out list of harmonised studies split by variant frequency
  message("Saving...")
  outname <- paste0(paste(source, class, exposure_name, outcome_name, sep = "_"),".rda")
  save(dat_harmonised, file = file.path(out_dir, outname))
}

format_sumstats <- function(path) {
  sumstats <- fread(path,data.table = F) |> filter(cohort == "470k_eur")
  
  # Retrieve gene positions from Ensembl
  geneIDs <- sumstats$`Ensembl ID`
  ensembl <- biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

  genePos <- biomaRt::getBM(
    attributes = c("ensembl_gene_id", "chromosome_name", "start_position"),
    filters = "ensembl_gene_id",
    values = geneIDs,
    mart = ensembl)

  sumstats <- sumstats |> left_join(
    genePos, by = c(`Ensembl ID` = "ensembl_gene_id")) |> 
    mutate(
      effect_allele = "A", # Dummy column needed for TwoSampleMR
      other_allele = "C", # Dummy column needed for TwoSampleMR
      AF = 0.1,  # Dummy allele frequency
      SE = abs(beta/(qnorm(pval/2, lower.tail = FALSE)))) |> 
    rename(
      "description" = "Trait",
      "Gene" = "Gene symbol",
      "chr" = "chromosome_name",
      "pos" = "start_position")
    
  return(sumstats)
}

main(exposure_study = exposure_study, outcome_study = outcome_study, class = class)