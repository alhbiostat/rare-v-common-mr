# Date: 18-03-2025
# Author: A.L.Hanson
# Purpose: Format Genebass extracted exposure and outcome data for use in MR
# Extract instruments of varying frequencies and clump

args <- commandArgs(trailingOnly = TRUE)

exposure_study <- args[1]
outcome_study <- args[2]
exposure_name <- args[3]
outcome_name <- args[4]
source <- args[5]
class <- args[6]

library(data.table)
library(dotenv)
library(TwoSampleMR)
library(dplyr)

# Load environment variables
dotenv::load_dot_env(file = "config.env")
data_dir <- Sys.getenv("data_dir")
ld_dir <- Sys.getenv("ld_dir")

out_dir <- file.path(data_dir,"harmonised")

message("Formatting studies:")
message("Exposure: ", exposure_name, " (", exposure_study, ")\n",
  "Outcome: ", outcome_name, " (", outcome_study, ")\n",
  "Source: ", source, "\nClass: ", class)

main <- function(exposure_study, outcome_study, class) {
  infile <- c(
    exposure = file.path(data_dir, "sumstats", source, exposure_study),
    outcome = file.path(data_dir, "sumstats", source, outcome_study)
  )

  infile_formatted <- list()

  # Update column names and format
  for(i in 1:length(infile)){
    message("Formatting: ", basename(infile[i]), " class: ", class)
    infile_formatted[[i]] <- format_sumstats(path = infile[i], class = class)
  }
  names(infile_formatted) <- c("exposure","outcome")

  # Split instruments by MAF or mask type and extract from exposure and outcome data for MR
  if(class == "variant"){
    maf_ranges <- list(c(0.05,0.95),c(0.01,0.05),c(0,0.01),c(0,0.001))
    exposure_split <- list()

    message("Formatting exposure...")
    for (maf_range in 1:length(maf_ranges)) {
      message("MAF: ", maf_ranges[[maf_range]][1], "-", maf_ranges[[maf_range]][2])
      
      exposure_split[[maf_range]] <- exp_format(
        exposure_df = infile_formatted$exposure, 
        maf_range = maf_ranges[[maf_range]],
        pval = 8e-9) # Genebass variant threshold
      
      message("n instruments: ", nrow(exposure_split[[maf_range]]))
    }
  
    names(exposure_split) <- c("exposure_common", "exposure_lowfreq", "exposure_rare", "exposure_ultrarare")

    # Perform LD clumping for common and low frequency variants
    # For rare variats lacking ld infomation, keep top hit per gene
    exposure_finalset <- list()
  
    message("LD clumping for common and low frequency variants...")
    exposure_finalset[[1]] <- perform_clumping(exposure_df = exposure_split$exposure_common)
    exposure_finalset[[2]] <- perform_clumping(exposure_df = exposure_split$exposure_lowfreq)
    
    exposure_finalset[[3]] <- exposure_split$exposure_rare
    exposure_finalset[[4]] <- exposure_split$exposure_ultrarare

    message("Keep top variant per gene for rare and ultra-rare variants...")
    exposure_finalset[[5]] <- keeptop_rare(exposure_df = exposure_split$exposure_rare)
    exposure_finalset[[6]] <- keeptop_rare(exposure_df = exposure_split$exposure_ultrarare)

    names(exposure_finalset) <- c("exposure_common", "exposure_lowfreq", "exposure_rare", "exposure_ultrarare", "exposure_rare_filt", "exposure_ultrarare_filt")

    # Format outcome
    message("Formatting outcome...")
    outcome_formatted <- TwoSampleMR::format_data(
      infile_formatted$outcome,
      type = "outcome",
      phenotype_col = "description",
      snp_col = "SNP",
      beta_col = "BETA",
      se_col = "SE",
      eaf_col = "AF",
      effect_allele_col = "effect_allele",
      other_allele_col = "other_allele",
      pval_col = "Pvalue")

    # Extract variants from outcome study and harmonise
    dat_harmonised <- lapply(exposure_finalset,
      function(x){
        TwoSampleMR::harmonise_data(exposure_dat = x, outcome_dat = outcome_formatted)
      })
    names(dat_harmonised) <- c("common", "lowfreq", "rare", "ultrarare", "rare_filt", "ultrarare_filt")

    message("Final instrument count:")
    print(unlist(lapply(dat_harmonised, nrow)))
  }

  if(class == "mask"){
    masks <- list("missense|LC","pLoF","pLoF|missense|LC","synonymous")
    exposure_split <- list()

    message("Formatting exposure...")
    for (mask in 1:length(masks)) {
      message("Mask: ", masks[[mask]])
      
      exposure_split[[mask]] <- exp_format_mask(
        exposure_df = infile_formatted$exposure, 
        mask = masks[[mask]],
        pval = 6.7e-7) # Genebass burden mask threshold
      
      message("n instruments: ", nrow(exposure_split[[mask]]))
    }

    names(exposure_split) <- paste("exposure", masks, sep = "_")

    # No LD clumping for masks

    # Format outcome
    message("Formatting outcome...")
    outcome_formatted <- TwoSampleMR::format_data(
      infile_formatted$outcome,
      type = "outcome",
      phenotype_col = "description",
      snp_col = "SNP",
      beta_col = "BETA_Burden",
      se_col = "SE_Burden",
      eaf_col = "AF",
      effect_allele_col = "effect_allele",
      other_allele_col = "other_allele",
      pval_col = "Pvalue_Burden")

    # Extract variants from outcome study and harmonise
    dat_harmonised <- lapply(exposure_split,
      function(x){
        TwoSampleMR::harmonise_data(exposure_dat = x, outcome_dat = outcome_formatted)
      })
    names(dat_harmonised) <- unlist(masks)

    message("Final instrument count:")
    print(unlist(lapply(dat_harmonised, nrow)))
  }

  # Write out list of harmonised studies split by variant frequency
  message("Saving...")
  outname <- paste0(paste(source, class, exposure_name, outcome_name, sep = "_"),".rda")
  save(dat_harmonised, file = file.path(out_dir, outname))

}

format_sumstats <- function(path, class) {
  
  if(class == "variant"){
    sumstats <- fread(
      path,
      drop = c(
        "call_stats","saige_version","description_more","category",
        "coding", "modifier","n_cases_defined","n_cases_both_sexes",
        "n_cases_females","n_cases_males","coding_description"),
      data.table = F)

    sumstats <- sumstats |> mutate(
      effect_allele = gsub('.*,"([A-Z].*)".*', "\\1", alleles),
      other_allele = gsub('.*"([A-Z].*)",.*', "\\1", alleles),
      chr = sub("chr(.*):.*", "\\1", locus),
      pos = sub("chr.*:(.*)", "\\1", locus),
      SNP = paste0(chr,":",pos,"_",other_allele,"_",effect_allele),
      description = ifelse(is.na(description), phenocode, description))
      
      return(sumstats)
  }
  
  if(class == "mask"){
    sumstats <- fread(
      path,
      drop = c(
        "markerIDs", "markerAFs", "Nmarker_MACCate_1", "Nmarker_MACCate_2", "Nmarker_MACCate_3",
        "Nmarker_MACCate_4", "Nmarker_MACCate_5", "Nmarker_MACCate_6", "Nmarker_MACCate_7","Nmarker_MACCate_8",
        "heritability", "saige_version","coding", "modifier","n_cases_defined","n_cases_both_sexes",
        "n_cases_females","n_cases_males","coding_description","description_more", "category"),
      data.table = F)
  
  sumstats <- sumstats |> mutate(
      effect_allele = "A", # Dummy column needed for TwoSampleMR
      other_allele = "C", # Dummy column needed for TwoSampleMR
      AF = 0.1, # Dummy allele frequency
      chr = sub(".chr(.*):.*", "\\1", interval),
      pos = sub(".chr.*:(.*)\\)", "\\1", interval),
      SNP = paste(gene_symbol, annotation, interval, sep = "_"),
      description = ifelse(is.na(description), phenocode, description))
      
      return(sumstats)
  }
}

exp_format <- function(exposure_df, maf_range, pval){
  
  snps_keep <- exposure_df |> 
    filter((AF > maf_range[1] & AF <= maf_range[2]) | (1-AF > maf_range[1] & 1-AF <= maf_range[2])) |>
    filter(Pvalue < pval) |> # Theshold from Genebass paper
    pull(markerID)
  
  exposure_formatted <- TwoSampleMR::format_data(exposure_df |> filter(markerID %in% snps_keep),
                                   type = "exposure",
                                   phenotype_col = "description",
                                   snp_col = "SNP",
                                   beta_col = "BETA",
                                   se_col = "SE",
                                   eaf_col = "AF",
                                   effect_allele_col = "effect_allele",
                                   other_allele_col = "other_allele",
                                   pval_col = "Pvalue")
  return(exposure_formatted)
}

exp_format_mask <- function(exposure_df, mask, pval){
  
  masks_keep <- exposure_df |> 
    filter(annotation == mask)|>
    filter(Pvalue_Burden < pval) |> # Theshold from Genebass paper
    pull(SNP)
  
  exposure_formatted <- TwoSampleMR::format_data(exposure_df |> filter(SNP %in% masks_keep),
                                   type = "exposure",
                                   phenotype_col = "description",
                                   snp_col = "SNP",
                                   beta_col = "BETA_Burden",
                                   se_col = "SE_Burden",
                                   eaf_col = "AF",
                                   effect_allele_col = "effect_allele",
                                   other_allele_col = "other_allele",
                                   pval_col = "Pvalue_Burden")
  return(exposure_formatted)
}

# Fixing bugs ieugwasr::ld_clump_local()
ld_clump_fix <- function (dat, clump_kb, clump_r2, clump_p, bfile, plink_bin) {
  shell <- ifelse(Sys.info()["sysname"] == "Windows", "cmd", "sh")
  fn <- tempfile()
  write.table(data.frame(SNP = dat[["rsid"]], P = dat[["pval"]]), 
              file = fn, row.names = FALSE, col.names = TRUE, quote = FALSE)
  fun2 <- paste0(shQuote(plink_bin, type = shell), " --bfile ", 
                 shQuote(bfile, type = shell), " --clump ", 
                 shQuote(fn, type = shell), " --clump-p1 ", clump_p, " --clump-r2 ", 
                 clump_r2, " --clump-kb ", clump_kb, " --out ", 
                 shQuote(fn, type = shell))
  system(fun2)
  res <- read.table(paste(fn, ".clumps", sep = ""), header = TRUE, comment.char = "")
  unlink(paste(fn, "*", sep = ""))
  y <- subset(dat, !dat[["rsid"]] %in% res[["ID"]])
  if (nrow(y) > 0) {
    message("Removing ", length(y[["rsid"]]), " of ", nrow(dat), 
            " variants due to LD with other variants or absence from LD reference panel")
  }
  return(subset(dat, dat[["rsid"]] %in% res[["ID"]]))
}

perform_clumping <- function(exposure_df){
  
  # Order SNP alleles alphabetically and capitalise to match reference
  df_in <- exposure_df 
  df_in$SNP <- apply(df_in, 1, function(x){
    allele_arrange <- order(as.character(x[c("effect_allele.exposure","other_allele.exposure")]))
    a1 <- x[c("effect_allele.exposure","other_allele.exposure")][allele_arrange][[1]]
    a2 <- x[c("effect_allele.exposure","other_allele.exposure")][allele_arrange][[2]]
    toupper(paste0(x['chr.exposure'], ":", x['pos.exposure'], "_", a1, "_", a2))
  })
  
  clumped_df <- ld_clump_fix(tibble(rsid=df_in$SNP, pval=df_in$pval.exposure, id=df_in$exposure),
                             clump_p = 1,
                             clump_r2 = 0.001,
                             clump_kb = 10000,
                             bfile = file.path(ld_dir, "full"),
                             plink_bin = "/usr/local/bin/plink2")
  
  positions_keep <- sub("_.*","",clumped_df$rsid)
  
  return(subset(exposure_df, sub("_.*","",exposure_df$SNP) %in% positions_keep))
}

keeptop_rare <- function(exposure_df){
  
  exposure_filt <- exposure_df |> 
    arrange(gene.exposure, pval.exposure) |> 
    filter(duplicated(gene.exposure) == FALSE) |> 
    arrange(chr.exposure, as.numeric(pos.exposure))
  
  return(exposure_filt)
}

main(exposure_study = exposure_study, outcome_study = outcome_study, class = class)