# Date: 24-10-2025
# Author: A.L.Hanson
# Purpose: Lookup instruments for proteins (pQTLs) within each IV set in the summary statistics of all other proteins
# (for an assessment of pleiotropy)
# Output: 
# Compiled proteome-wide summary statitics for pQTLs significant for at least one protein

dotenv::load_dot_env(file = "config.env")
data_dir <- Sys.getenv("data_dir")
prot_dir <- file.path(data_dir,"sumstats","proteomics")

load(file = file.path(prot_dir, "protein_pQTLstolookup_5e-08.rda")) #snp_lookup

IV_sets <- names(snp_lookup)
IV_source <- c(rep("exome-aggregate",4), rep("exome",5), rep("common",2))
names(IV_source) <- IV_sets

# File locations:
sumstats <- c(
  "exome-aggregate" = file.path(prot_dir,"exome-aggregate"),
  "exome" = file.path(prot_dir,"exome"),
  "common" = file.path(prot_dir, "common"))

for (set in IV_sets){
  snps <- snp_lookup[[set]]
  sumstats_dir <- sumstats[IV_source[set]]

  # Write SNPs to temporary file for lookup
  temp_file <- file.path(prot_dir, paste0("tmp_snps_to_lookup_",set, ".txt"))
  write(snps, temp_file)

  out_file <- file.path(prot_dir, paste0("pQTL_lookup_", set, ".tsv"))

  if(IV_source[set] == "exome-aggregate"){
    header <- "Gene\tPhenotype\tCollapsingModel\tType\tpValue\tCategory\tnSamples\tYesQV\tNoQV\tConCaseMedian\tConCtrMedian\tbeta\tConBetaSe\tLCI\tUCI\tprotein_target\tprobe\tolink_id\tsample_size\tvariant_type\tcis_trans"
    writeLines(header, con = out_file)
    masktype <- sub(".*_","",set)
    # Run grep on appropriate summary statistics directory
    cmd <- paste("grep -w -F -f", temp_file, file.path(sumstats_dir, "*.tsv"), "| grep -w", masktype, ">>", out_file)
    system(cmd)
  }
  if(IV_source[set] == "exome"){
    header <- "Variant\tCHR\tBP\tOA\tEA\tVariantType\tPhenotype\tCategory\tModel\tANNOTATION\tGENE\tTranscript\tcDNAchange\tAminoacidchange\tExonrank\tNo.samples\tNo.AAgenotypes\tNo.ABgenotypes\tNo.BBgenotypes\tEAF\t%ABorBBgenotypes\t%BBgenotypes\tP\tBETA\tSE\tEffectsizeLCI\tEffectsizeUCI\tMedianvalue(cases)\tMedianvalue(controls)\tMAF\tprotein_target\tprobe\tolink_id\tsample_size\tvariant_type\tcis_trans"
    writeLines(header, con = out_file)
    # Run grep on appropriate summary statistics directory
    cmd <- paste("grep -w -F -f", temp_file, file.path(sumstats_dir, "*.tsv"), ">>", out_file)
    system(cmd)
  }
  if(IV_source[set] == "common"){
    header <- "RSID\tCHR\tBP\tEA\tOA\tEAF\tBETA\tSE\tP\tSNP\tZ\tregion\tprotein_target\tprobe\tolink_id\tsample_size\tvariant_type\tcis_trans\tfinemapped_hit\tVEP_RSID"
    writeLines(header, con = out_file)
    # Run grep on appropriate summary statistics directory
    cmd <- paste("grep -w -F -f", temp_file, file.path(sumstats_dir, "*.tsv"), ">>", out_file)
    system(cmd)
  }

  # Cleanup
  file.remove(temp_file)
}
