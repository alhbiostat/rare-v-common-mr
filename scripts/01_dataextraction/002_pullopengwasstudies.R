# Date: 17-03-2025
# Author: A.L.Hanson
# Purpose: Extract common variant instruments for exposure and outcomes interest from OpenGWAS
# Output: Harmonised summary statistics for MR

args <- commandArgs(trailingOnly = TRUE)

exposure_study <- args[1]
outcome_study <- args[2]
exposure_name <- args[3]
outcome_name <- args[4]
source <- args[5]
class <- args[6]

library(dotenv)
library(TwoSampleMR)
library(dplyr)

# Load environment variables
dotenv::load_dot_env(file = "config.env")
data_dir <- Sys.getenv("data_dir")

out_dir <- file.path(data_dir,"harmonised")

message("Formatting studies:")
message("Exposure: ", exposure_name, " (", exposure_study, ")\n",
  "Outcome: ", outcome_name, " (", outcome_study, ")\n",
  "Source: ", source, "\nClass: ", class)

main <- function(exposure_study, outcome_study){
  # Get instruments for exposure from OpenGWAS
  exposure_dat <- TwoSampleMR::extract_instruments(exposure_study)
  message("n instruments: ", nrow(exposure_dat))

  # Extract instruments from outcome
  outcome_dat <- TwoSampleMR::extract_outcome_data(snps = exposure_dat$SNP, outcomes = outcome_study)

  # Harmonise SNP data (fill missing eafs from paired study if possible)
  harmonised_dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat) |>
    dplyr::mutate(eaf.exposure = ifelse(is.na(eaf.exposure) & !is.na(eaf.outcome), eaf.outcome, eaf.exposure),
      eaf.outcome = ifelse(is.na(eaf.outcome) & !is.na(eaf.exposure), eaf.exposure, eaf.outcome))

  message("Final instrument count:", nrow(harmonised_dat))

  # Write out harmonised studies
  message("Saving...")
  outname <- paste0(paste(source, class, exposure_name, outcome_name, sep = "_"),".rda")
  save(harmonised_dat, file = file.path(out_dir, outname))
}

main(exposure_study = exposure_study, outcome_study = outcome_study)