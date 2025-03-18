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

main <- function(exposure_study, outcome_study, class) {
  infile <- c(
    exposure = file.path(data_dir, "sumstats", source, exposure_study),
    outcome = file.path(data_dir, "sumstats", source, outcome_study)
  )

  infile_formatted <- list()

  for(i in 1:length(infiles)){
    message("Formatting: ", basename(infile[i]), " class: ", class)
    infile_formatted[[i]] <- format_sumstats(path = infile[i], class = class)
  }
  names(infile_formatted) <- c("exposure","outcome")

  # Extract sets of exposure variants of required frequency and format data for MR
  maf_ranges <- list(c(0.05,0.95),c(0.01,0.05),c(0,0.01),c(0,0.001))
  exposure_split <- list()
  
  message("Formatting exposure...")
  for (maf_range in 1:length(maf_ranges)) {
    message("MAF: ", maf_ranges[[maf_range]][1], "-", maf_ranges[[maf_range]][2])
    exposure_split[[maf_range]] <- exp_format(exposure_df = infile_formatted$exposure,  maf_range = maf_ranges[[maf_range]])
    message("n instruments: ", nrow(exposure_split[[maf_range]]))
  }
  
  names(exposure_split) <- c("exposure_common", "exposure_lowfreq", "exposure_rare", "exposure_ultrarare")

  # Perform LD clumping for common and low frequency variants
  message("LD clumping for common and low frequency variants...")
  exposure_common <- perform_clumping(exposure_df = exposure_split$exposure_common)
  exposure_lowfreq <- perform_clumping(exposure_df = exposure_split$exposure_lowfreq)
  


### UP TO HERE:
# Down filter rare variants to one per gene, but keep both sets of resutls
# Extract instruments from outcome GWAS
# Harmonise and save out

# Do the same with the masks

# Modify OpenGWAS common variant harmonisation to take studies listed in studypairings.csv

  # For rare variats not in reference panel keep top hit per gene
  exposure_rare <- exposure_split$exposure_rare |> 
    arrange(gene.exposure, pval.exposure) |> 
    filter(duplicated(gene.exposure) == FALSE) |> 
    arrange(chr.exposure, pos.exposure)

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
      SNP = paste0(chr,":",pos,"_",other_allele,"_",effect_allele))
      
      return(sumstats)
  }
  
  if(class == "mask"){}
}

exp_format <- function(exposure_df, maf_range){
  
  snps_keep <- exposure_df |> 
    filter((AF > maf_range[1] & AF <= maf_range[2]) | (1-AF > maf_range[1] & 1-AF <= maf_range[2])) |>
    filter(Pvalue < 5e-8) |> 
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
