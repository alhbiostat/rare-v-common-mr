# Date: 13-01-2026
# Author: A.L.Hanson
# Purpose: Generate list of protein CNVs - commplex trait study pairings for MR analysis 
# Duplication and deletion associations taken from: https://zenodo.org/records/16800547

library(here)

# Load environment variables
dotenv::load_dot_env(file = "config.env")
data_dir <- Sys.getenv("data_dir")
prot_dir <- file.path(data_dir,"sumstats","Milind_CNV_associations")
dirs <- list.dirs(file.path(prot_dir, "CNV_Burden_Tests_Proteins"), recursive = F)

# CNV summary statistics study files
cnv_deletion <- unlist(lapply(dirs, list.files, pattern = "*.deletions.WES.regenie.gz", full = TRUE))
cnv_duplication <- unlist(lapply(dirs, list.files, pattern = "*.duplications.WES.regenie.gz", full = TRUE))
cnv_plof <- unlist(lapply(dirs, list.files, pattern = "*.plofs.WES.regenie.gz", full = TRUE))

# Outcome studies
out_deletion <- list.files(file.path(prot_dir, "CNV_Burden_Tests"), pattern = "*.deletions.WES.regenie.gz")
out_duplication <- list.files(file.path(prot_dir, "CNV_Burden_Tests"), pattern = "*.duplications.WES.regenie.gz")
out_plof <- list.files(file.path(prot_dir, "CNV_Burden_Tests"), pattern = "*.plofs.WES.regenie.gz")

# Outcome phenotype IDs
out_phenos <- read.csv(file.path(prot_dir, "selected_phenotypes.csv"))
## Exclude body composition by impedance, hand grip and spirometry measures to reduce tests w less clear biological interpretation
out_phenos <- out_phenos |> 
  dplyr::filter(!(category_title %in% c("Body composition by impedance", "Hand grip strength", "Spirometry")))

## Filter from file names
out_deletion <- out_deletion[grep(paste(paste0("^",out_phenos$field_id, "\\."), collapse = "|"), out_deletion)]
out_duplication <- out_duplication[grep(paste(paste0("^",out_phenos$field_id, "\\."), collapse = "|"), out_duplication)]
out_plof <- out_plof[grep(paste(paste0("^",out_phenos$field_id, "\\."), collapse = "|"), out_plof)]

paired_studies_deletion <- expand.grid(
  exposure_study = sub(".*_Proteins\\/", "", cnv_deletion),
  outcome_study = out_deletion) |>
  dplyr::mutate(
    exposure_name = toupper(sub("^.*\\/(.*)\\.deletions.*", "\\1", exposure_study)),
    outcome_name = out_phenos[match(sub("\\..*", "", outcome_study), out_phenos$field_id), "field_title"],
    exposure_source = "milind_deletions",
    outcome_source = "milind_deletions",
    class = "mask") # 166,611
# Format outcome name
paired_studies_deletion$outcome_name <- gsub("\\.","-",make.names(paired_studies_deletion$outcome_name))

paired_studies_duplication <- expand.grid(
  exposure_study = sub(".*_Proteins\\/", "", cnv_duplication),
  outcome_study = out_duplication) |>
  dplyr::mutate(
    exposure_name = toupper(sub("^.*\\/(.*)\\.duplications.*", "\\1", exposure_study)),
    outcome_name = out_phenos[match(sub("\\..*", "", outcome_study), out_phenos$field_id), "field_title"],
    exposure_source = "milind_duplications",
    outcome_source = "milind_duplications",
    class = "mask") # 166,611 
# Format outcome name
paired_studies_duplication$outcome_name <- gsub("\\.","-",make.names(paired_studies_duplication$outcome_name))

paired_studies_plof <- expand.grid(
  exposure_study = sub(".*_Proteins\\/", "", cnv_plof),
  outcome_study = out_plof) |>
  dplyr::mutate(
    exposure_name = toupper(sub("^.*\\/(.*)\\.plofs.*", "\\1", exposure_study)),
    outcome_name = out_phenos[match(sub("\\..*", "", outcome_study), out_phenos$field_id), "field_title"],
    exposure_source = "milind_plofs",
    outcome_source = "milind_plofs",
    class = "mask") # 166,611 
# Format outcome name
paired_studies_plof$outcome_name <- gsub("\\.","-",make.names(paired_studies_plof$outcome_name))

paired_join <- rbind(paired_studies_deletion, paired_studies_duplication, paired_studies_plof) # 499,833

write.csv(paired_join, here("allstudypairings_cnvs.csv"), row.names = F, quote = F)