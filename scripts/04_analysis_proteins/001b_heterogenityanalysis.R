# Date: 10-03-2026
# Author: A.L.Hanson
# Purpose: Run heterogenity analysis on IVs used for MR

# To run: Rscript 001b_heterogenityanalysis.R [HARMONISED OBJECT] [OUTCOME]
# e.g. Rscript 001b_heterogenityanalysis.R harmonised_studies_GW5e-08_EW1e-04_rarepruned_COMBINED.rda "BMI"

library(dotenv)
library(here)
library(dplyr)
library(TwoSampleMR)

dotenv::load_dot_env(file = here("config.env"))
data_dir <- file.path(Sys.getenv("data_dir"),"harmonised","proteomics")
res_dir <- Sys.getenv("results_dir")

args = commandArgs(trailingOnly=TRUE)

# Check for supplied arguments
if (length(args)!=2) {
  stop("A harmonised SNP object and an outcome must be specified", call.=FALSE)
} else{
  message("Input file: ", args[1])
  message("Outcome: ", args[2])
}

infile <- args[1]
outcome <- args[2]

# Load results
env <- new.env()
load(file.path(data_dir, infile), envir = env)
snps <- as.list(env)[[1]]

# Filter to the specified outcome
if (!outcome %in% names(snps)) {
  stop("Outcome not found in the data", call.=FALSE)
}
snps <- snps[outcome]

# Filter to variants only (remove gene mask IVs)
snps <- lapply(snps, function(x){
  lapply(x, function(y){
    y <- y[grep("variant",names(y))]
    return(y)
  })})

# For the specified outcome, perform instrument heterogenity analysis
all_het_results <- list()

ivs <- snps[[1]]
message("Performing heterogenity analysis for outcome: ", outcome)
for (j in 1:length(ivs)){
  message("Analysing: ", names(ivs)[j])
  set <- ivs[j][[1]]
  
  for (k in 1:length(set)) {
    df <- set[[k]]
    df <- df[!is.na(df$cis_trans) & df$cis_trans != "", ]
    set_name <- names(set)[k]
    
    # Check if dataframe is empty or has less than 2 rows
    if (is.null(df) || nrow(df) < 2) {
      message("Skipping ", set_name, ": insufficient data (", if(is.null(df)) 0 else nrow(df), " rows)")
      next
    }
    
    # Split dataframe by cis_trans
    if (!"cis_trans" %in% colnames(df)) {
      message("Skipping ", set_name, ": no cis_trans column")
      next
    }
    
    splits <- split(df, df$cis_trans)
    
    for (variant_type in names(splits)) {
      sub_df <- splits[[variant_type]]
      sub_df <- sub_df[sub_df$mr_keep == TRUE, ]
      
      # Check again for each split
      if (nrow(sub_df) < 2) {
        message("Skipping ", set_name, " - ", variant_type, ": insufficient data (", nrow(sub_df), " rows)")
        next
      }
      
      # Perform heterogeneity analysis
      het <- mr_heterogeneity(sub_df)
      
      # Add columns
      het$variant_type <- variant_type
      het$variant_set <- set_name
      het$n_snps <- nrow(sub_df)
      het$outcome <- outcome
      het$iv_set <- names(ivs)[j]
      
      # Collect results
      all_het_results <- c(all_het_results, list(het))
    }
  }
}

# Combine all results into a single dataframe
if (length(all_het_results) > 0) {
  combined_het <- do.call(rbind, all_het_results)
  message("Heterogeneity analysis completed for ", outcome)
  # Save to file
  snpset <- sub(".*(GW.*)_COMBINED.rda","\\1",infile)
  outfile <- file.path(res_dir, paste0("results_heterogenity_", snpset, "_", outcome))
  write.csv(combined_het, file = paste0(outfile, ".csv"), row.names = FALSE)
} else {
  message("No heterogeneity results for ", outcome, ".")
}
