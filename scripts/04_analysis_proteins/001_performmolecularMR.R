# Date: 07-08-2025
# Author: A.L.Hanson
# Purpose: Read in harmonised summary statistics for a given outcome (combined across exposures using 005_joinharmonisedstudies.R)
# perform MR and output results
# To run: Rscript 001_performmolecularMR.R [HARMONISED OBJECT]

# Common GWAS IVs (split cis, trans, combined): 
#1. MR with finemapped hits (accounting for LD)
#2. MR with clumped and pruned hits

# Exome IVs (split cis, trans, combined):
#1. MR with clumped and pruned common variants
#2. MR with rare and ultrarare unpruned variants
#3. MR with rare and ultrarare pruned variants (top hit per gene)

# Exome masks (split cis, trans, combined):
# MR with aggregate IVs:

library(dotenv)
library(here)
library(dplyr)
library(ggplot2)

dotenv::load_dot_env(file = here("config.env"))
data_dir <- file.path(Sys.getenv("data_dir"),"harmonised","proteomics")
res_dir <- Sys.getenv("results_dir")

args = commandArgs(trailingOnly=TRUE)

# Check for supplied outcome
if (length(args)!=1) {
  stop("A harmonised study data object must be specified", call.=FALSE)
} else{
  message("Input file: ", args[1])
}

harmonised_file <- args[1]
load(file.path(data_dir, harmonised_file)) # harmonised_studies

# Load harmonised study data

# Perform MR on each exposure-outcome pair, across each of 10 instrument sets (split by cis, trans and combined IVs)

# Check if MR results already exist
if(file.exists(file.path(res_dir, "results_complextrait_MR.rda"))){
  ## FIX
  load(file.path(res_dir, "results_complextrait_MR.rda"))
  results_MR_old <- results_MR

  # If new study pairs exist that are not in results:
  study_pairs <- names(results_MR_old)
  # New study pairs
  study_pairs <- names(harmonised_studies)[!(names(harmonised_studies) %in% study_pairs)]

  if(length(study_pairs) == 0){
    message("No new study pairs to analyse")
  }else{
    message("New studies to analyse: ", paste(study_pairs, collapse = ", "))
    harmonised_studies <- harmonised_studies[study_pairs]
    results_MR <- list()
  }
}else{
  results_MR <- list()
}

# Analyse each (new) study pair
for(i in 1:length(harmonised_studies)){
  
  res_list <- list()
  
  for(j in 1:length(harmonised_studies[[i]])){

    dat <- harmonised_studies[[i]][j]
    
    message("Performing MR for ", names(harmonised_studies)[i],": ", names(harmonised_studies[[i]])[j])
    
    skip_pair <- (nrow(dat[[1]]) == 1 & is.na(dat[[1]]$mr_keep[1])) | (nrow(dat[[1]]) == 1 & dat[[1]]$mr_keep[1] == FALSE)

    if(!is.na(dat[[1]]$SNP[1]) & skip_pair == FALSE){
      # Perform LD aware MR for common GWAS finemapped IV sets
      if(names(dat) == "sun_gwas_variant"){

        ## need to pull out clumped as well as finemapped IVs
        ### CODE HERE 

      } else {
        # Split cis, trans and combined instruments
        dat <- dat[[1]]
        dat_split <- split(dat, dat$cis_trans)
        dat_split$combined <- dat

        res_split <- lapply(dat_split, function(x){
          res <- TwoSampleMR::mr(x)
          res$pair <- names(harmonised_studies)[i]
          res$study <- names(harmonised_studies[[i]])[j]
          IV_set <- unique(x$cis_trans)
          res$IV_set <- paste(IV_set[order(IV_set)],collapse = "_")
          return(res)
        })
        res <- unique(do.call("rbind", res_split))
      }
    } else {
      res <- data.frame(
        id.exposure = NA,
        id.outcome = NA,
        outcome = NA,
        exposure = NA,
        method = NA,
        nsnp = NA,
        b = NA,
        se = NA,
        pval = NA,
        pair = NA,
        study = NA
      )
    }
    
    res_list <- c(res_list, list(res))
  }
  res_join <- do.call("rbind", res_list)
  results_MR <- c(results_MR, list(res_join))
}

names(results_MR) <- names(harmonised_studies)

## Append new results to existing results_MR object
if(exists("results_MR_old")){
  results_MR <- c(results_MR_old, results_MR)
  # Reorder studies
  order_studies <- names(results_MR)[order(names(results_MR))]
  results_MR <- results_MR[order_studies]
  
  ## Write out results_MR object
  save(results_MR, file = file.path(res_dir, "results_complextrait_MR.rda"))
}else{
  ## Write out harmonised_study object
  save(results_MR, file = file.path(res_dir, "results_complextrait_MR.rda"))
}
