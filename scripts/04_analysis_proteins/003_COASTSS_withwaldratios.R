# Date: 24-11-2025
# Author: A.L.Hanson
# Purpose: Run rare variant allelic series analysis using COAST-SS with Wald Ratios derived for O'link protein cis coding variants
# and 11 outcomes of interest
# This essentially tests whether increading the deleterousness of a rare variant (e.g from benign missense -> deleterious missense -> protein truncating/loss of function)
# increases the strength of the causal association with the outcome
# Variants need to first be categoriesed into these groups using VEP annotations

# See:
# https://www.cell.com/ajhg/fulltext/S0002-9297(25)00368-4
# https://github.com/insitro/AllelicSeries

# Load environment variables
dotenv::load_dot_env(file = "config.env")
data_dir <- Sys.getenv("data_dir")
res_dir <- Sys.getenv("results_dir")
gpmap_data_dir <- Sys.getenv("gpmap_data_dir")

# Read in Wald Ratios (see 02_dataextraction_proteins/006_calcwaldratios.R)
load(file.path(res_dir, "results_molecular_waldratios_ciscodingvariants.rda")) #wald_ratios

# Load VEP variant annotations
vep <- data.table::fread(file.path(gpmap_data_dir,"variant_annotation","genebass_vep_rare_variantannotations_hg38_altered.txt"))

# Extract SIFT and PolyPhen calls
annot_col <- "Extra"
#vep[, REF_allele := stringr::str_extract(get(annot_col), "REF_ALLELE=[^;]+")] FIX
vep[, SIFT_token := stringr::str_extract(get(annot_col), "SIFT=[^;]+")]
vep[, PolyPhen_token := stringr::str_extract(get(annot_col), "PolyPhen=[^;]+")]
vep[, SIFT_label := ifelse(is.na(SIFT_token), NA_character_,
                           stringr::str_replace(SIFT_token, "^SIFT=([^\\(]+)\\(.*$", "\\1"))]
vep[, PolyPhen_label := ifelse(is.na(PolyPhen_token), NA_character_,
                               stringr::str_replace(PolyPhen_token, "^PolyPhen=([^\\(]+)\\(.*$", "\\1"))]

# Assign consequences to classes based on McCaw et al. classifications
# PTV = protein truncating variant
# DMV = damaging missense variant
# BTV = benign missense variant
concequences <- unique(unlist(sapply(vep$Consequence, function(x){unlist(strsplit(x,","))})))

ptv <- c("splice_acceptor_variant","splice_donor_variant","stop_gained","frameshift_variant")
missense <- c("missense_variant","inframe_deletion","inframe_insertion","stop_lost","start_lost","protein_altering_variant")
dmv <- c("possibly_damaging","probably_damaging","deleterious","deleterious_low_confidence")
bmv <- c("benign","tolerated","tolerated_low_confidence")

# Bind variant annotations
vep_annots <- vep |>
  dplyr::select(Uploaded_variation, Allele, Consequence, SIFT_label, PolyPhen_label)
### Continue from here









