# Date: 05-06-2025
# Author: A.L.Hanson
# Purpose: Add rsIDs for finemapped pQTLs from Sun et al. (as taken from GPMAP finemapped studies)

# Load environment variables
dotenv::load_dot_env(file = "config.env")
data_dir <- Sys.getenv("data_dir")
prot_dir <- file.path(data_dir,"sumstats","proteomics")

# Variant annotations
vep <- data.table::fread(file.path(data_dir, "variant_annotations", "vep_variantannotations_hg38_altered.txt"))

# Load summary statisitics
files <- list.files(file.path(prot_dir, "common"), pattern = "sun_pQTL_.*", full.names = T)

system(paste("mkdir -p", file.path(prot_dir, "common/archive")))

add_rsids <- function(file){

  dat <- data.table::fread(file)
  
  # Move original file to archive
  file_bn <- basename(file)
  file.rename(file, file.path(prot_dir, "common/archive", file_bn))

  vep_ids <- vep |> dplyr::filter(Uploaded_variation %in% dat$SNP) |>
    dplyr::select(Uploaded_variation, Existing_variation) |>
    dplyr::rename(SNP = Uploaded_variation) |>
    dplyr::mutate(VEP_RSID = gsub(",.*","",Existing_variation))
  
  dat <- dat |>
    dplyr::left_join(vep_ids, by = "SNP") |>
    dplyr::select(-Existing_variation)

  data.table::fwrite(dat, file, sep = "\t", quote = FALSE)
}

lapply(files, add_rsids)
