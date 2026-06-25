# Date: 15-06-2026
# Author: A.L.Hanson
# Purpose: Modification of 001_performmolecularMR.R to repeat rare/ultrarare-variant MR with LD pruned true-cis variants in the protein coding gene
# Utilising saved object "harmonised_studies_GW5e-08_EW1e-04_rarepruned_COMBINED.rda" containing ld pruned harmonised SNPs at the
# relaxed p-value theshold for rare/ultrarare IVs

library(dotenv)
library(here)
library(dplyr)
library(ggplot2)

dotenv::load_dot_env(file = here("config.env"))
data_dir <- file.path(Sys.getenv("data_dir"),"harmonised")
res_dir <- Sys.getenv("results_dir")

harmonised_file <- "harmonised_studies_GW5e-08_EW1e-04_rarepruned_COMBINED.rda"
results_file <- "results_molecular_MR_GW5e-08_EW5e-08_rarepruned_truecis_COMBINED.rda"

# Load harmonised study data
load(file.path(data_dir, harmonised_file)) # harmonised_snps_relaxed

# Keep only the ultrarare IV sets and restrict rows to gene.outcome == exposure (true-cis IVs)
harmonised_snps_relaxed <- lapply(harmonised_snps_relaxed, function(outcome_list){
  lapply(outcome_list, function(exposure_list){
    
    keep_names <- intersect(names(exposure_list), c(
      "dhindsa_exwas_variant_ultrarare_pruned",
      "dhindsa_exwas_variant_rare_pruned"
    ))
    exposure_list <- exposure_list[keep_names]
    lapply(exposure_list, function(df){
      if(is.data.frame(df) && all(c("gene.outcome", "exposure") %in% names(df))){
        df[df$gene.outcome == df$exposure, , drop = FALSE]
      } else {
        df
      }
    })
  })
})

# Perform MR on each exposure-outcome pair, across each rare/ultrarare instrument set (true cis)
results_MR_all <- list()
for(exp in 1:length(harmonised_snps_relaxed)){
  results_MR <- list()
  
  for(i in 1:length(harmonised_snps_relaxed[[exp]])){
    res_list <- list()
    for(j in 1:length(harmonised_snps_relaxed[[exp]][[i]])){
      dat <- harmonised_snps_relaxed[[exp]][[i]][j]
      message("Performing MR for ", names(harmonised_snps_relaxed)[exp],": ", names(harmonised_snps_relaxed[[exp]])[i])
      skip_pair <- (all(is.na(dat[[1]]$mr_keep)) | all(dat[[1]]$mr_keep == FALSE))
      
      if(!is.na(dat[[1]]$SNP[1]) & skip_pair == FALSE){
        # Perform standard TwoSampleMR::mr()
        res <- TwoSampleMR::mr(dat[[1]])
        res$pair <- names(harmonised_snps_relaxed[[exp]][i])
        res$IV_set <- names(dat)
      } else {
        res <- data.frame(
          id.exposure = NA, id.outcome = NA, outcome = NA, exposure = NA, method = NA, nsnp = NA,
          b = NA, se = NA, pval = NA, pair = names(harmonised_snps_relaxed[[exp]][i]), IV_set = names(dat))
      }
      res_list <- c(res_list, list(res))
    }
    
    res_join <- do.call("rbind", res_list)
    rownames(res_join) <- NULL
    res_join <- list(res_join)
    names(res_join) <- names(harmonised_snps_relaxed[[exp]][i])
    results_MR <- c(results_MR, res_join)
  }
  results_MR <- list(results_MR)
  names(results_MR) <- names(harmonised_snps_relaxed)[exp]
  results_MR_all <- c(results_MR_all,results_MR)
}

## Write out results_MR object
message("Writing to results file: ", file.path(res_dir, results_file))
save(results_MR_all, file = file.path(res_dir, results_file))
