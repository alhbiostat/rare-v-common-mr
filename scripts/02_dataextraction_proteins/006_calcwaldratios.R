# Date: 14-11-2025
# Author: A.L.Hanson
# Purpose: Extract protein expression and outcome associations for all rare coding variants in O'link protein coding genes (irrespective of pQTL significance)
# to derive Wald ratios for alleleic series analysis (COAST-SS)

# Output: 
# List of 11 outcomes with Wald ratios calculated for each cis-coding variant (for O'link proteins)
# See results/results_molecular_waldratios_ciscodingvariants.rda

# Load environment variables
dotenv::load_dot_env(file = "config.env")
data_dir <- Sys.getenv("data_dir")
res_dir <- Sys.getenv("results_dir")
prot_dir <- file.path(data_dir,"sumstats","proteomics")
project_dir <- Sys.getenv("project_dir")

# List rare variant proteomics studies
studies <- data.table::fread(file.path(prot_dir, "gpmap_ukbproteomics_studies.tsv"), data.table = F) |> 
  dplyr::filter(variant_type == "rare_exome") |>
  dplyr::mutate(protein_target = unlist(lapply(strsplit(trait_name, " "),function(x){x[1]}))) |>
  dplyr::arrange(protein_target)

# Loop through each set of O'link protein summary statistics and extract rare variants within target protein
# Note only associations of p<0.01 are published on the AzPheWAS portal

rare_coding <- mapply(FUN = function(x,y){
  dat <- data.table::fread(x) |> 
    dplyr::select(Variant, CHR, BP, OA, EA, ANNOTATION, GENE, `cDNA change`, `Amino acid change`, EAF, P, BETA, SE) |>
    dplyr::mutate(
      protein_target = y,
      GENE=gsub("'","",GENE)
    )
  target <- y
  dat <- dat |> dplyr::filter((GENE == protein_target & EAF <= 0.01) | (GENE == protein_target & EAF >= 1-0.01)) 
  return(dat)
}, x = studies$study_location, y = studies$protein_target, SIMPLIFY = F)

names(rare_coding) <- studies$protein_target

# Retain proteins with at least 2 cis-variants
rare_coding <- rare_coding[lapply(rare_coding, nrow) >= 2]

# Variants to extract from outcome studies
snps <- unique(unlist(lapply(rare_coding, function(x){x$Variant})))
snps <- gtools::mixedsort(snps)
#snps_pos <- sub("^([0-9].*-[0-9].*)-.*-.*$","\\1",snps)

# Outcome studies (Genebass)
# N.B these were considered exposures for complex trait MR
studies_outcomes <- data.table::fread(file.path(project_dir, "exposure_studies.csv")) |>
  dplyr::filter(source == "genebass", class == "variant")

# Loop through outcome studies,format SNP IDs to match those in `snps` and extract required variant
rare_coding_outcome <- lapply(file.path(data_dir, "sumstats/genebass", studies_outcomes$exposure_study), function(x){
  message("Reading file:", x)
  dat <- data.table::fread(x) |>
    dplyr::mutate(snp = sub("chr","",gsub(":|_|/","-",markerID))) |>
    dplyr::filter(snp %in% snps)
})
names(rare_coding_outcome) <- studies_outcomes$exposure_name

# Format protein and outcome associations, harmonise and and calculate Wald ratios using the TwoSampleMR package:
# Format exposure
rare_coding_proteins <- do.call("rbind", rare_coding)

pQTLs_formatted <- TwoSampleMR::format_data(
  data.frame(rare_coding_proteins),
  type = "exposure",
  phenotype_col = "protein_target",
  snp_col = "Variant",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "EAF",
  effect_allele_col = "EA",
  other_allele_col = "OA",
  pval_col = "P")

# Format outcome and harmonise
harmonised_dat <- lapply(rare_coding_outcome, function(x){
  dat <- x |>
    dplyr::mutate(
      OA = sapply(alleles, function(x){jsonlite::fromJSON(x)[1]}),
      EA = sapply(alleles, function(x){jsonlite::fromJSON(x)[2]}))

  outcome_formatted <- TwoSampleMR::format_data(
    data.frame(dat),
    type = "outcome",
    phenotype_col = "description",
    snp_col = "snp",
    beta_col = "BETA",
    se_col = "SE",
    eaf_col = "AF",
    effect_allele_col = "EA",
    other_allele_col = "OA",
    pval_col = "Pvalue")

  harmonised_dat <- TwoSampleMR::harmonise_data(exposure_dat = pQTLs_formatted, outcome_dat = outcome_formatted)
  return(harmonised_dat)
})

names(harmonised_dat) <- names(rare_coding_outcome)

# Calculate Wald Ratios
wald_ratios <- lapply(harmonised_dat, function(x){
  wrs <- apply(x, MARGIN = 1, FUN = function(x){
    TwoSampleMR::mr_wald_ratio(
    b_exp = as.numeric(x["beta.exposure"]),
    b_out = as.numeric(x["beta.outcome"]),
    se_exp = as.numeric(x["se.exposure"]),
    se_out = as.numeric(x["se.outcome"]))
  })

  wrs_df <- do.call("rbind",lapply(wrs, unlist))
  colnames(wrs_df) <- paste("wald_ratio", colnames(wrs_df), sep = "_")

  # Bind to harmonised data
  out <- cbind(x,wrs_df)
  return(out)
})

names(wald_ratios) <- names(harmonised_dat)

save(wald_ratios, file = file.path(res_dir, "results_molecular_waldratios_ciscodingvariants.rda"))
