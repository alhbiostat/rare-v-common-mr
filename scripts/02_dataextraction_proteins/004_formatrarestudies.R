# Date: 05-06-2025
# Author: A.L.Hanson
# Purpose: Harmonise ExWAS pQTLs and gene-based aggregate tests with complex trait outcome studies

args <- commandArgs(trailingOnly = TRUE)
# Test args
#args <- c("dhindsa_pQTL_MAPK9_OID20557.tsv","genebass_ukbwes_p4080_Systolic_blood_pressure_automated_reading.tsv","MAPK9","SBP","dhindsa_exwas","genebass","variant")
#args <- c("dhindsa_pQTL_aggregatetest_WFIKKN2_OID20785.tsv","genebass_ukbwes_p4080_Systolic_blood_pressure_automated_reading_genetest.tsv","WFIKKN2","SBP","dhindsa_exwas","genebass","mask")
#args <- c("dhindsa_pQTL_ZP3_OID30265.tsv","genebass_ukbwes_p30870_Triglycerides.tsv","ZP3","Trig","dhindsa_exwas","genebass","variant")

exposure_study <- args[1]
outcome_study <- args[2]
exposure_name <- args[3]
outcome_name <- args[4]
exposure_source <- args[5]
outcome_source <- args[6]
class <- args[7]

library(data.table)
library(dotenv)
library(TwoSampleMR)
library(dplyr)

# Load environment variables
dotenv::load_dot_env(file = "config.env")
data_dir <- Sys.getenv("data_dir")
ld_dir <- Sys.getenv("ld_dir")

out_dir <- file.path(data_dir,"harmonised","proteomics")

message("Formatting studies:")
message("Exposure: ", exposure_name, " (", exposure_study, ")\n",
  "Outcome: ", outcome_name, " (", outcome_study, ")\n",
  "Source: ", exposure_source, "\nClass: ", class)

main <- function(exposure_study, outcome_study, class) {

  if(class == "variant"){
    infile <- c(
      exposure = file.path(data_dir, "sumstats/proteomics/exome", exposure_study),
      outcome = file.path(data_dir, "sumstats/genebass", outcome_study)
    )
    infile_formatted <- list()
    
    # Format exposure and outcome
    message("Formatting: ", basename(infile[1]), " class: ", class)
    infile_formatted[[1]] <- format_exposure(path = infile[1], class = class)
    snps <- infile_formatted[[1]]$SNP
    message("Formatting: ", basename(infile[2]), " class: ", class)
    infile_formatted[[2]] <- format_outcome(path = infile[2], class = class, snps = snps)

    names(infile_formatted) <- c("exposure","outcome")
  
    # Split instruments by MAF or mask type and extract from exposure and outcome data for MR
    maf_ranges <- list(c(0.01,0.99),c(0,0.01),c(0,0.001))
    exposure_split <- list()

    message("Splitting instruments by MAF...")
    for (maf_range in 1:length(maf_ranges)) {
      message("MAF: ", maf_ranges[[maf_range]][1], "-", maf_ranges[[maf_range]][2])
      
      exposure_split[[maf_range]] <- split_exposure(
        exposure_df = infile_formatted$exposure, 
        maf_range = maf_ranges[[maf_range]]) # Genebass variant threshold
      
      message("n instruments: ", nrow(exposure_split[[maf_range]]))
    }
  
    names(exposure_split) <- c("exposure_common", "exposure_rare", "exposure_ultrarare")

    # Perform LD clumping for common and variants
    # For rare variats lacking ld information, keep all, and top hit per gene (filt) post harmonisation
    exposure_finalset <- list()
  
    message("LD clumping for common variants...")
    exposure_finalset[[1]] <- perform_clumping(exposure_df = exposure_split$exposure_common)
    
    exposure_finalset[[2]] <- exposure_split$exposure_rare
    exposure_finalset[[3]] <- exposure_split$exposure_ultrarare

    names(exposure_finalset) <- c("exposure_common", "exposure_rare", "exposure_ultrarare")

    # Format outcome
    message("Formatting outcome for harmonisation...")
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
        if(nrow(x)==0){
          return(data.frame())
        }else{
          # Reappend cis/trans and variant annotation information
          x <- x|> left_join(infile_formatted$exposure[,c("SNP","cis_trans","ANNOTATION")] |> 
            mutate(SNP = tolower(SNP)))
          TwoSampleMR::harmonise_data(exposure_dat = x, outcome_dat = outcome_formatted)
        }})
    names(dat_harmonised) <- c("common", "rare", "ultrarare")

    message("Keep top variant per gene for rare and ultra-rare variants...")
    dat_harmonised[[4]] <- keeptop_rare(harmonised_df = dat_harmonised$rare)
    dat_harmonised[[5]] <- keeptop_rare(harmonised_df = dat_harmonised$ultrarare)

    names(dat_harmonised) <- c("common", "rare", "ultrarare", "rare_filt", "ultrarare_filt")

    message("Final instrument count:")
    print(unlist(lapply(dat_harmonised, nrow)))
  }

  if(class == "mask"){
    infile <- c(
      exposure = file.path(data_dir, "sumstats/proteomics/exome-aggregate", exposure_study),
      outcome = file.path(data_dir, "sumstats/genebass", outcome_study)
    )
    infile_formatted <- list()

    # Format exposure and outcome
    message("Formatting: ", basename(infile[1]), " class: ", class)
    infile_formatted[[1]] <- format_exposure(path = infile[1], class = class)
    message("Formatting: ", basename(infile[2]), " class: ", class)
    infile_formatted[[2]] <- format_outcome(path = infile[2], class = class)

    names(infile_formatted) <- c("exposure","outcome")
        
    # Split exposures and outcome studies by aggregate mask type 
    # (which differs slightly between Dhindsa and Genebass studies)
    masks_exposure <- list("raredmg","ptv","ptvraredmg","syn")
    masks_outcome <- list("missense|LC","pLoF","pLoF|missense|LC","synonymous")

    exposure_split <- list()
    outcome_split <- list()

    message("Splitting exposure masks...")
    for (mask in 1:length(masks_exposure)) {
      message("Mask: ", masks_exposure[[mask]])
      
      exposure_split[[mask]] <- split_exposure_mask(
        exposure_df = infile_formatted$exposure, 
        mask = masks_exposure[[mask]])
      
      message("n instruments: ", nrow(exposure_split[[mask]]))
    }

    names(exposure_split) <- paste("exposure", masks_exposure, sep = "_")

    message("Splitting outcome masks...")
    for (mask in 1:length(masks_outcome)) {
      message("Mask: ", masks_outcome[[mask]])
      
      outcome_split[[mask]] <- split_outcome_mask(
        outcome_df = infile_formatted$outcome, 
        mask = masks_outcome[[mask]])
    }

    names(outcome_split) <- paste("outcome", masks_outcome, sep = "_")

    # No LD clumping for masks

    # Extract genes from outcome study and harmonise between mask pairs
    # raredmg=missense|LC ptv=pLoF ptvraredmg=pLoF|missense|LC syn=synonymous

    dat_harmonised <- mapply(x = exposure_split, y = outcome_split,
      function(x,y){
        if(nrow(x)==0){
          return(data.frame())
        }else{
          TwoSampleMR::harmonise_data(exposure_dat = x, outcome_dat = y) |>
          mutate(
            test.exposure = unique(x$test),
            test.outcome = unique(y$test),
            cis_trans = ifelse(SNP == tolower(exposure), "cis", "trans")
          )
        }
      }, SIMPLIFY = FALSE)
    
    names(dat_harmonised) <- unlist(masks_exposure)

    message("Final instrument count:")
    print(unlist(lapply(dat_harmonised, nrow)))
  }

  # Write out list of harmonised studies split by variant frequency
  message("Saving...")
  outname <- paste0(paste(exposure_source, class, exposure_name, outcome_name, sep = "_"),".rda")
  save(dat_harmonised, file = file.path(out_dir, outname))
}

# FUNTIONS #

format_exposure <- function(path, class) {
  
  if(class == "variant"){
    sumstats <- fread(
      path,
      drop = c("Variant type", "Phenotype", "Category", "Model",
      "No. AA genotypes", "No. AB genotypes", "No. BB genotypes",
      "% AB or BB genotypes", "% BB genotypes",
      "Median value (cases)", "Median value (controls)"),
      sep = "\t",
      data.table = F)

    sumstats <- sumstats |> mutate(
      SNP = paste0(CHR,":",BP,"_",OA,"_",EA)
    )  
    return(sumstats)
  }

  if(class == "mask"){
    sumstats <- fread(
      path,
      drop = c(
        "Phenotype","Type","Category","YesQV","NoQV","ConCaseMedian","ConCtrMedian"),
      data.table = F)
  
    sumstats <- sumstats |> mutate(
        effect_allele = "A", # Dummy column needed for TwoSampleMR
        other_allele = "C", # Dummy column needed for TwoSampleMR
        AF = 0.1, # Dummy allele frequency
        chr = 1, # Dummy chromosome
        pos = 1, # Dummy position
        SNP = paste(Gene, CollapsingModel, sep = "_"),
        description = protein_target)
      
      return(sumstats)
  }
}

format_outcome <- function(path, class, snps) {
  
  if(class == "variant"){
    snps_search <- paste(c("locus",gsub("_.*","",snps)), collapse = "|")

    sumstats <- fread(
      cmd = paste0("rg ","'",snps_search,"' ",path),
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

split_exposure <- function(exposure_df, maf_range){
  
  snps_keep <- exposure_df |> 
    filter((EAF > maf_range[1] & EAF <= maf_range[2]) | (1-EAF > maf_range[1] & 1-EAF <= maf_range[2])) |>
    pull(SNP)

  if(length(snps_keep) == 0){
    exposure_formatted <- data.frame()
  }else{
    exposure_formatted <- TwoSampleMR::format_data(exposure_df |> filter(SNP %in% snps_keep),
                                    type = "exposure",
                                    phenotype_col = "protein_target",
                                    snp_col = "SNP",
                                    chr_col = "CHR",
                                    pos_col = "BP",
                                    beta_col = "BETA",
                                    se_col = "SE",
                                    eaf_col = "EAF",
                                    effect_allele_col = "EA",
                                    other_allele_col = "OA",
                                    pval_col = "P")
  }
  return(exposure_formatted)
}

split_exposure_mask <- function(exposure_df, mask){
  
  masks_keep <- exposure_df |> 
    filter(CollapsingModel == mask)|>
    pull(SNP)

  if(length(masks_keep) == 0){
    return(data.frame())
  }else{
    exposure_formatted <- TwoSampleMR::format_data(exposure_df |> filter(SNP %in% masks_keep),
                                    type = "exposure",
                                    phenotype_col = "description",
                                    snp_col = "Gene", #Pair exposure and outcome studies by gene for equivalent masks
                                    beta_col = "beta",
                                    se_col = "ConBetaSe",
                                    eaf_col = "AF",
                                    effect_allele_col = "effect_allele",
                                    other_allele_col = "other_allele",
                                    pval_col = "pValue") |>
                                    mutate(test = mask)
    return(exposure_formatted)
  }
}

split_outcome_mask <- function(outcome_df, mask){
  
  masks_keep <- outcome_df |> 
    filter(annotation == mask)|>
    pull(SNP)

  if(length(masks_keep) == 0){
    return(data.frame())
  }else{
    outcome_formatted <- TwoSampleMR::format_data(outcome_df |> filter(SNP %in% masks_keep),
                                    type = "outcome",
                                    phenotype_col = "description",
                                    snp_col = "gene_symbol", #Pair exposure and outcome studies by gene for equivalent masks
                                    beta_col = "BETA_Burden",
                                    se_col = "SE_Burden",
                                    eaf_col = "AF",
                                    effect_allele_col = "effect_allele",
                                    other_allele_col = "other_allele",
                                    pval_col = "Pvalue_Burden") |>
                                    mutate(test = mask)
    return(outcome_formatted)
  }
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
    chr <- as.character(as.numeric(x['chr.exposure']))
    pos <- as.character(as.numeric(x['pos.exposure']))
    toupper(paste0(chr, ":", pos, "_", a1, "_", a2))
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

keeptop_rare <- function(harmonised_df){

  if(nrow(harmonised_df) == 0){
    harmonised_filt <- data.frame()
  }else{
    harmonised_filt <- harmonised_df |> 
    # Nb. only genebass outcomes have gene annotated
      arrange(gene.outcome, pval.exposure) |> 
      filter(duplicated(gene.outcome) == FALSE) |> 
      arrange(chr.exposure, as.numeric(pos.exposure))
  }
  return(harmonised_filt)
}

main(exposure_study = exposure_study, outcome_study = outcome_study, class = class)