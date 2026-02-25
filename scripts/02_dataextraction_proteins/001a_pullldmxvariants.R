# Date: 20-02-2026
# Author: A.L.Hanson
# Purpose: Write out list of rare variant instruments from which to derive LD matrices using UKB WES data

# Load environment variables
dotenv::load_dot_env(file = "config.env")
data_dir <- Sys.getenv("data_dir")

# Read in extracted instruments for O-link proteins (threshold at p<1e-4)
files <- list.files(file.path(data_dir,"sumstats/proteomics/exome"), pattern = "dhindsa_pQTL",  full.names = T)

dat <- lapply(files, function(x){data.table::fread(x)})
names(dat) <- sub(".*pQTL_(.*)_.*.tsv","\\1",files)

# Filter to rare variants
dat_rare <- lapply(dat, function(x){
  x |> dplyr::filter(EAF <= 0.01 | 1-EAF <= 0.01)
})

# Filter to ExWAS P<=5e-8
dat_rare_top <- lapply(dat, function(x){
  x |> 
  dplyr::filter(EAF <= 0.01 | 1-EAF <= 0.01) |>
  dplyr::filter(P <= 5e-8)
})

# Merge variants, filter to unique and order
merged <- do.call("rbind", dat_rare) |>
  dplyr::select("Variant", "CHR", "BP", "EAF", "Situated_gene" = "GENE", "Instrumented_protein" = "protein_target", "cis_trans") |>
  dplyr::arrange(CHR, BP)

merged_top <- do.call("rbind", dat_rare_top) |>
  dplyr::select("Variant", "CHR", "BP", "EAF", "Situated_gene" = "GENE", "Instrumented_protein" = "protein_target", "cis_trans") |>
  dplyr::arrange(CHR, BP)

data.table::fwrite(merged, file = file.path(data_dir,"sumstats/proteomics/exome/rareinstruments_p1e-4.tsv"), sep = "\t")
data.table::fwrite(merged_top, file = file.path(data_dir,"sumstats/proteomics/exome/rareinstruments_p5e-8.tsv"), sep = "\t")