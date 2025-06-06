# Date: 14-05-2025
# Author: A.L.Hanson
# Purpose: Extract common and rare variant instruments for UKB O-link protein exposures from GPMAP processed studies 
# Output: 
# Common (finemapped) variant summary statistics (top variant per finemapped region, p-value thesholded, for each protein)
# Rare variant summary statistics (p-value thesholded, for each protein)

library(EnsDb.Hsapiens.v86)

# Load environment variables
dotenv::load_dot_env(file = "config.env")
data_dir <- Sys.getenv("data_dir")
rawdata_dir <- Sys.getenv("rawdata_dir")
gpmap_data_dir <- Sys.getenv("gpmap_data_dir")
prot_dir <- file.path(data_dir,"sumstats","proteomics")

# Load list of UKB O-link single variant studies (and correct gene name)
studies <- data.table::fread(file.path(prot_dir, "gpmap_ukbproteomics_studies.tsv")) |>
  dplyr::mutate(
    gene = ifelse(variant_type == "common", gsub(":.*","",probe), gsub(" .*", "", trait_name)),
    probe = ifelse(variant_type == "rare_exome", gsub(" ",":",trait_name), probe),
    olink_id = sub(".*:(OID[0-9]+):.*", "\\1", probe),
  )

olink_proteins <- unique(studies$gene)

# Load list of gene aggregate studies
studies_aggregate <- data.table::fread(file.path(prot_dir, "azphewas_ukbproteomics_aggregatetests.tsv"))

# Get cis regions for each protein coding gene (from ST3 Sun et al.)
gene_ranges <- data.table::fread(file.path(prot_dir, "olink_generanges_Sunetal.csv"), check.names = TRUE)
# Take first symbol/position for proteins with multiple symbols
gene_ranges <- gene_ranges |> dplyr::mutate(
  Gene.CHROM = as.numeric(Gene.CHROM),
  Gene.symbol = ifelse(grepl(";", Gene.symbol), sub(";.*", "", Gene.symbol), Gene.symbol),
  Gene.start = ifelse(grepl(";", Gene.start), as.numeric(sub(";.*", "", Gene.start)), as.numeric(Gene.start)),
  Gene.end = ifelse(grepl(";", Gene.end), as.numeric(sub(";.*", "", Gene.end)), as.numeric(Gene.end))
)

gene_ranges$cis.start <- gene_ranges$Gene.start - 1000000
gene_ranges$cis.end <- gene_ranges$Gene.start + 1000000

# For each single variant study, read in finemapped (common) or extracted (rare) variants, filter on p-value threhold and combine across loci

message("Single variant studies:")

for(i in 1:nrow(studies)){
  study <- studies[i,]
  sample_size <- study$sample_size
  protein_target <- study$gene
  probe <- study$probe
  oid <- study$olink_id

  cis_region <- gene_ranges |> dplyr::filter(Olink.ID == oid)

  message("Processing study: ", study$study_name, " for protein: ", protein_target, " (", oid, ")")

  if(study$variant_type == "common"){

    output_file <- file.path(prot_dir, "common", paste0("sun_pQTL_",protein_target, "_", oid, ".tsv"))
    if(file.exists(output_file)){
      message("Output file already exists: ", output_file)
      next
    }

    # Finemapped regions for study
    study_dir <- study$extracted_location
    finemapped_dir <- paste0(study_dir, "finemapped")
    finemapped_files <- list.files(finemapped_dir, full.names = TRUE)
    
    if(length(finemapped_files) == 0){
      message("No finemapped files found for study: ", study$study_name)
      next
    }
    
    finemapped_regions <- sub(".*\\/(EUR_.*)\\.tsv.gz", "\\1", finemapped_files)

    # Read in finemapped variants across ld-blocks
    finemapped <- lapply(finemapped_files, function(x) {
      dt <- data.table::fread(x, colClasses = c("RSID" = "character"))

      if(!("IMPUTED" %in% colnames(dt))){
        dt$IMPUTED <- NA
      }

      dt |> dplyr::filter(IMPUTED == FALSE) |>
        dplyr::slice_max(n=1, abs(Z), with_ties = FALSE) # Keep top variant per region based on abs Z-score
    })

    finemapped <- finemapped |> 
      dplyr::bind_rows() |>
      dplyr::mutate(
        region = finemapped_regions[lapply(finemapped, nrow) == 1],
        protein_target = protein_target,
        probe = probe,
        olink_id = oid,
        sample_size = sample_size,
        variant_type = "common",
        cis_trans = ifelse(CHR == cis_region$Gene.CHROM & BP >= cis_region$cis.start & BP <= cis_region$cis.end, "cis", "trans")) |>
      dplyr::filter(P <= 5e-8 & EAF > 0.01 & 1-EAF > 0.01) |>
      dplyr::arrange(CHR, BP)

    # Save finemapped variants to file
    message("Num GWAS variants: ", nrow(finemapped))
    data.table::fwrite(finemapped, output_file, sep = "\t", quote = FALSE)
  }

  if(study$variant_type == "rare_exome"){

    output_file <- file.path(prot_dir, "exome", paste0("dhindsa_pQTL_",protein_target, "_", oid, ".tsv"))
    if(file.exists(output_file)){
      message("Output file already exists: ", output_file)
      next
    }

    # Summary statistics
    stats <- study$study_location

    # Read in summary statistics for exome variants (all frequencies)
    standardised <- data.table::fread(stats) |> 
      dplyr::mutate(
        protein_target = protein_target,
        probe = probe,
        olink_id = oid,
        sample_size = sample_size,
        variant_type = "exome",
        cis_trans = ifelse(CHR == cis_region$Gene.CHROM & BP >= cis_region$cis.start & BP <= cis_region$cis.end, "cis", "trans")) |>
      dplyr::filter(P <= 1e-4) |> # Dhindsa et al. suggestive p-value threshold
      dplyr::arrange(CHR, BP)

    # Save exome variants to file
    message("Num ExWAS variants: ", nrow(standardised))
    data.table::fwrite(standardised, output_file, sep = "\t", quote = FALSE)
  }
}

message("Gene aggregate studies:")

for(i in 1:nrow(studies_aggregate)){
  study <- studies_aggregate[i,]
  sample_size <- study$sample_size
  protein_target <- study$gene
  probe <- study$probe
  oid <- sub(".*_","",probe)

  message("Processing study: ", study$study_name, " for protein: ", protein_target, " (", oid, ")")

  # Summary statistics
  stats <- study$study_location

  # Read in summary statistics for exome variants (all frequencies)
  aggregates <- data.table::fread(stats) |> 
    dplyr::mutate(
      protein_target = protein_target,
      probe = probe,
      olink_id = oid,
      sample_size = sample_size,
      variant_type = "exome_aggregate",
      cis_trans = ifelse(Gene == protein_target, "cis", "trans")) |>
    dplyr::filter(pValue <= 1e-8 & CollapsingModel %in% c("syn","ptv","raredmg","ptvraredmg","UR")) # Dhindsa et al. aggregate test p-value threshold

    # Save exome variants to file
    message("Num ExWAS aggregates: ", nrow(aggregates))
    output_file <- file.path(prot_dir, "exome-aggregate", paste0("dhindsa_pQTL_aggregatetest_",protein_target, "_", oid, ".tsv"))
    data.table::fwrite(aggregates, output_file, sep = "\t", quote = FALSE)
}