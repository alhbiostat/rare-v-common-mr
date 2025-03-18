# Date: 17-03-2025
# Author: A.L.Hanson
# Purpose: Extract common variant instruments for exposure and outcomes interest from OpenGWAS
# Output: Harmonised summary statistics for MR

# Studies:
# LDL direct: ieu-b-110 
# BMI: ukb-b-19953
# Vitamin D: ebi-a-GCST90000615
# Triglycerides: ieu-b-111
# Glycated Haemoglobin (HbA1c): ebi-a-GCST90014006
# Mean platelet volume: ebi-a-GCST90013981
# Mean corpuscular volume: ebi-a-GCST90013979
# IGF-1: ebi-a-GCST90014008
# WHRadjBMI: ieu-a-79
# RBC erythrocyte count: ebi-a-GCST90013977

# CAD: ebi-a-GCST90013864
# T2D: ebi-a-GCST90013892
# MS: ukb-b-17670
# Stroke: ukb-b-6358
# AF: ukb-b-964
# VTE: ukb-b-12040
# Prostate cancer: ukb-b-2160
# Hypertension: ukb-a-531

library(dotenv)
library(TwoSampleMR)

# Load environment variables
dotenv::load_dot_env(file = "config.env")
data_dir <- Sys.getenv("data_dir")

# Get instruments (nb some exposures are used twice for different outcomes)
exposure_studies <- c("ieu-b-110",
                      "ukb-b-19953",
                      "ebi-a-GCST90000615",
                      "ieu-b-111",
                      "ebi-a-GCST90014006",
                      "ebi-a-GCST90013981",
                      "ebi-a-GCST90014008",
                      "ieu-a-79",
                      "ebi-a-GCST90013977",
                      "ebi-a-GCST90013979",
                      "ebi-a-GCST90013981")

names(exposure_studies) <- c("LDL",
                             "BMI",
                             "Vitamin D",
                             "Triglycerides",
                             "HbA1c",
                             "Mean platelet volume",
                             "IGF-1",
                             "WHRadjBMI",
                             "RBC count",
                             "Mean corpuscular volume",
                             "Mean platelet volume")

# Extract instruments for exposures
exposure_dat <- lapply(exposure_studies, function(x){
  extract <- try(extract_instruments(x))
  if(class(extract) == "try-error"){
    return(NA)
  } else{
    return(extract)
  }
})

names(exposure_dat) <- names(exposure_studies)

# Extract instruments from corresponding paired outcome studies
outcome_studies <- c("ebi-a-GCST90013864",
                     "ebi-a-GCST90013892",
                     "ukb-b-17670",
                     "ukb-b-6358",
                     "ukb-b-964",
                     "ukb-b-12040",
                     "ukb-b-2160",
                     "ukb-a-531",
                     "ukb-b-6358",
                     "ukb-b-6358",
                     "ukb-b-6358")

names(outcome_studies) <- c("CAD",
                             "T2D",
                             "MS",
                             "Stroke",
                             "AF",
                             "VTE",
                             "Prostate cancer",
                             "Hypertension",
                             "Stroke",
                             "Stroke",
                             "Stroke")

outcome_dat <- list()
for (i in 1:length(outcome_studies)){
  outcome_dat[[i]] <- extract_outcome_data(snps = exposure_dat[[i]]$SNP, outcomes = outcome_studies[i])
}

names(outcome_dat) <- names(outcome_studies)

# Harmonise SNP data (fill missing eafs from paired study if possible)
harmonised_dat <- list()
for (i in 1:length(outcome_studies)){
  harmonised_dat[[i]] <- harmonise_data(exposure_dat = exposure_dat[[i]], outcome_dat = outcome_dat[[i]]) |>
    dplyr::mutate(eaf.exposure = ifelse(is.na(eaf.exposure) & !is.na(eaf.outcome), eaf.outcome, eaf.exposure),
           eaf.outcome = ifelse(is.na(eaf.outcome) & !is.na(eaf.exposure), eaf.exposure, eaf.outcome))
}

names(harmonised_dat) <- make.names(paste(names(exposure_studies), names(outcome_studies), sep = "_"))

# Save out paired and harmonised studies
save(harmonised_dat, file = file.path(data_dir,"harmonised","commonGWAS_harmonised.rda"))