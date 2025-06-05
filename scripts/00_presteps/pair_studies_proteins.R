# Date: 04-06-2025
# Author: A.L.Hanson
# Purpose: Generate list of protein (pQTL) - commplex trait study pairings for MR analysis 

library(here)

# Load environment variables
dotenv::load_dot_env(file = "config.env")
data_dir <- Sys.getenv("data_dir")
prot_dir <- file.path(data_dir,"sumstats","proteomics")

# pQTL summary statistics study files
prot_common <- list.files(file.path(prot_dir, "common"), pattern = "sun_pQTL_.*")
prot_exome <- list.files(file.path(prot_dir, "exome"), pattern = "dhindsa_pQTL_.*")
prot_aggregate <- list.files(file.path(prot_dir, "exome-aggregate"), pattern = "dhindsa_pQTL_aggregatetest_.*")

# Outcome studies (exposures from previous complex trait analysis)
out_studies <- read.csv(here("exposure_studies.csv"))

paired_studies_common <- expand.grid(
  exposure_study = prot_common,
  outcome_study = out_studies[out_studies$source == "opengwas", "exposure_study"]) |>
  dplyr::mutate(
    exposure_name = sub("sun_pQTL_(.*)_.*", "\\1", exposure_study),
    outcome_name = out_studies[match(outcome_study, out_studies$exposure_study), "exposure_name"],
    exposure_source = "sun_gwas",
    outcome_source = "opengwas",
    class = "variant") # 39,546

paired_studies_exome <- expand.grid(
  exposure_study = prot_exome,
  outcome_study = out_studies[(out_studies$source == "genebass" & out_studies$class == "variant"), "exposure_study"]) |>
  dplyr::mutate(
    exposure_name = sub("dhindsa_pQTL_(.*)_.*", "\\1", exposure_study),
    outcome_name = out_studies[match(outcome_study, out_studies$exposure_study), "exposure_name"],
    exposure_source = "dhindsa_exwas",
    outcome_source = "genebass",
    class = "variant") # 32,351

paired_studies_aggregate <- expand.grid(
  exposure_study = prot_aggregate,
  outcome_study = out_studies[(out_studies$source == "genebass" & out_studies$class == "mask"), "exposure_study"]) |>
  dplyr::mutate(
    exposure_name = sub("dhindsa_pQTL_aggregatetest_(.*)_.*", "\\1", exposure_study),
    outcome_name = out_studies[match(outcome_study, out_studies$exposure_study), "exposure_name"],
    exposure_source = "dhindsa_exwas",
    outcome_source = "genebass",
    class = "mask") # 32,351

paired_join <- rbind(paired_studies_common, paired_studies_exome, paired_studies_aggregate) # 94,248

write.csv(paired_join, here("allstudypairings_proteins.csv"), row.names = F, quote = F)