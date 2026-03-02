# Date: 25-02-2026
# Author: A.L.Hanson
# Purpose: Modification of 001_performmolecularMR.R to repeat rare-variant MR with LD pruned variants
# See 02_dataextraction_proteins/001b_prunerareIVs.sh

# To run: Rscript 001a_performmolecularMR_prunerareIVs.R [HARMONISED OBJECT]

# Exome IVs (split cis, trans, combined):
#1. MR with rare and ultrarare pruned variants (r2>0.1)

library(dotenv)
library(here)
library(dplyr)
library(ggplot2)

dotenv::load_dot_env(file = here("config.env"))
data_dir <- file.path(Sys.getenv("data_dir"),"harmonised","proteomics")
res_dir <- Sys.getenv("results_dir")
pruned_dir <- file.path(Sys.getenv("data_dir"),"sumstats","proteomics","exome","protein_pQTL_pruning")

args = commandArgs(trailingOnly=TRUE)

# Check for supplied outcome
if (length(args)!=1) {
  stop("A harmonised study data object must be specified", call.=FALSE)
} else{
  message("Input file: ", args[1])
}

harmonised_file <- args[1]
# Load harmonised study data
load(file.path(data_dir, harmonised_file)) # harmonised_studies
message("Number of exposure:outcome pairs: ", length(harmonised_studies))

# Load list of variants to keep (r2 < 0.1)
# prunein_files <- grep("prune.in", list.files(pruned_dir, pattern = "p5e-8_.*", full = TRUE), value = T)
prunein_files <- grep("pQTLs_chr[0-9].*\\.prune\\.in", list.files(pruned_dir, full = TRUE), value = T)
prunein <- do.call("rbind", lapply(prunein_files, function(x){data.table::fread(x, header = F, data.table = F)}))[,1]
# Correct variant name formatting
prunein <- unlist(lapply(strsplit(prunein, ":"), function(x){paste0(x[1],":",x[2],"_",tolower(x[3]),"_",tolower(x[4]))}))

# Subset to rare and ultrarare variant instrument sets
harmonised_studies <- lapply(harmonised_studies, function(x){
  x[c("dhindsa_exwas_variant_rare","dhindsa_exwas_variant_ultrarare")]
})

# Prune variants from instrument sets
harmonised_studies_pruned <- lapply(harmonised_studies, function(x){
  lapply(x, function(y){
    if(!(is.na(y$SNP[1]))){
      y |> dplyr::filter(SNP %in% prunein)
    } else {y}
  })
})

outcome <- sub(".*_","",sub("harmonised_studies_(.*)\\.rda","\\1",harmonised_file))
#extension <- sub("harmonised_studies_(.*)_.*\\.rda","\\1",harmonised_file)
results_file <- paste0("results_molecular_MR", "_rarepruned_", outcome, ".rda")
harmonised_file_pruned <- paste0("harmonised_studies", "_rarepruned_", outcome, ".rda")

# Write out pruned IVs
if(!(file.exists(file.path(data_dir, harmonised_file_pruned)))){
  save(harmonised_studies_pruned, file = file.path(data_dir, harmonised_file_pruned))
}

# Perform MR on each exposure-outcome pair, across each rare/ultrare instrument set (split by cis, trans and combined IVs)

# Analyse each (new) study pair
for(i in 1:length(harmonised_studies_pruned)){

  # Check if MR results already exist
  if(file.exists(file.path(res_dir, results_file))){
    load(file.path(res_dir, results_file))
    results_MR_old <- results_MR

    # If new study pairs exist that are not in results:
    study_pairs <- names(results_MR_old)
    # New study pairs
    study_pairs <- names(harmonised_studies_pruned)[!(names(harmonised_studies_pruned) %in% study_pairs)]

    if(length(study_pairs) == 0){
      message("No new study pairs to analyse")
      break()
    }else{
      message("Pairs remaining: ", length(study_pairs))
      harmonised_studies_pruned <- harmonised_studies_pruned[study_pairs]
      i <- 1
      results_MR <- list()
    }
  }else{
    results_MR <- list()
  }

  res_list <- list()
  
  for(j in 1:length(harmonised_studies_pruned[[i]])){
    
    dat <- harmonised_studies_pruned[[i]][j]
    
    message("Performing MR for ", names(harmonised_studies)[i],": ", names(harmonised_studies[[i]])[j])
    
    skip_pair <- (all(is.na(dat[[1]]$mr_keep)) | all(dat[[1]]$mr_keep == FALSE))

    if(!is.na(dat[[1]]$SNP[1]) & skip_pair == FALSE){
      # Split cis and trans instruments
      dat1 <- dat[[1]]
      if(all(c("cis","trans") %in% unique(dat1$cis_trans))){
        dat_split <- split(dat1, dat1$cis_trans)
        dat_split$combined <- dat1
      }else{
        dat_split <- split(dat1, dat1$cis_trans)
      }
  
      # Perform standard TwoSampleMR::mr()
      res_split <- lapply(dat_split, function(x){

        skip_split <- (all(is.na(x$mr_keep)) | all(x$mr_keep == FALSE))
        if(skip_split == TRUE){
          res <- data.frame(
            id.exposure = NA, id.outcome = NA, outcome = NA, exposure = NA, method = NA, nsnp = NA,
            b = NA, se = NA, pval = NA, pair = names(harmonised_studies_pruned)[i], cis_trans = NA, IV_set = names(dat))
          return(res)
        }
        cis_trans <- unique(x$cis_trans)
        res <- TwoSampleMR::mr(x)
        res$pair <- names(harmonised_studies_pruned)[i]
        res$cis_trans <- paste(cis_trans[order(cis_trans)],collapse = "_")
        res$IV_set <- names(dat)
        return(res)
      })
      res <- unique(do.call("rbind", res_split))
    
    } else {
      res <- data.frame(
              id.exposure = NA, id.outcome = NA, outcome = NA, exposure = NA, method = NA, nsnp = NA,
              b = NA, se = NA, pval = NA, pair = names(harmonised_studies)[i], cis_trans = NA, IV_set = names(dat))
    }
    res_list <- c(res_list, list(res))
  }

  res_join <- do.call("rbind", res_list)
  rownames(res_join) <- NULL
  res_join <- list(res_join)
  names(res_join) <- names(harmonised_studies_pruned[i])

  results_MR <- c(results_MR, res_join)

  ## Append new results to existing results_MR object
  if(exists("results_MR_old")){
    results_MR <- c(results_MR_old, results_MR)
    # Reorder studies
    order_studies <- names(results_MR)[order(names(results_MR))]
    results_MR <- results_MR[order_studies]
    
    ## Write out results_MR object
    message("Writing to results file: ", results_file)
    save(results_MR, file = file.path(res_dir, results_file))
  }else{
    ## Write out results_MR object
    message("Writing to results file: ", results_file)
    save(results_MR, file = file.path(res_dir, results_file))
  }
}
