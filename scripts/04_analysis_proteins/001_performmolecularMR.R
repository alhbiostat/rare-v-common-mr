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
#1. MR with aggregate IVs: pLOF, missense, pLOF/missense, synonymous

library(dotenv)
library(here)
library(dplyr)
library(ggplot2)

dotenv::load_dot_env(file = here("config.env"))
data_dir <- file.path(Sys.getenv("data_dir"),"harmonised","proteomics")
res_dir <- Sys.getenv("results_dir")
ld_dir <- Sys.getenv("ld_dir")

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

outcome <- sub("harmonised_studies_(.*)\\.rda","\\1",harmonised_file)
results_file <- paste0("results_molecular_MR_", outcome, ".rda")

# Perform MR on each exposure-outcome pair, across each of 10 instrument sets (split by cis, trans and combined IVs)

# Check if MR results already exist
if(file.exists(file.path(res_dir, results_file))){

  load(file.path(res_dir, results_file))
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
      # Split cis and trans instruments
      dat1 <- dat[[1]]
      if(all(c("cis","trans") %in% unique(dat1$cis_trans))){
        dat_split <- split(dat1, dat1$cis_trans)
        dat_split$combined <- dat1
      }else{
        dat_split <- split(dat1, dat1$cis_trans)
      }
  
      # Perform LD aware MR for common GWAS finemapped IV sets
      if(names(dat) == "sun_gwas_variant_common_finemapped"){
        # Perform MR with correlated instruments using MendelianRandomisation::mr_ivw()
        res_split <- lapply(dat_split, function(x){
          
          skip_split <- (nrow(x) == 1 & is.na(x$mr_keep[1])) | (nrow(x) == 1 & x$mr_keep[1] == FALSE)
          if(skip_split == TRUE){
            res <- data.frame(
              id.exposure = NA, id.outcome = NA, outcome = NA, exposure = NA, method = NA, nsnp = NA,
              b = NA, se = NA, pval = NA, pair = names(harmonised_studies)[i], cis_trans = NA, IV_set = names(dat))
            return(res)
          }

          snplist <- x$SNP
          cis_trans <- unique(x$cis_trans)
          
          # Generate LD matrix for SNPs (alleles should be correctly oriented as using alphatised GPMAP data)
          ld_mat <- ieugwasr::ld_matrix(
            snplist,
            plink_bin = genetics.binaRies::get_plink_binary(),
            bfile = file.path(ld_dir, "full_rsid"),
            with_alleles = FALSE
          )

          if(length(snplist) == 1){
            # Remove any SNPs not in ld matrix
            snplist <- snplist[snplist %in% colnames(ld_mat)]
            x <- x |> dplyr::filter(SNP %in% snplist)
            rownames(ld_mat) <- paste("snp",1:nrow(ld_mat),sep = "_")
            colnames(ld_mat) <- rownames(ld_mat)
          }else{
            # Remove any SNPs not in ld matrix
            snplist <- snplist[snplist %in% colnames(ld_mat)]
            x <- x |> dplyr::filter(SNP %in% snplist)
            # Reorder ld matrix to match exposure data
            ld_mat <- ld_mat[snplist, snplist]
            rownames(ld_mat) <- paste("snp",1:nrow(ld_mat),sep = "_")
            colnames(ld_mat) <- rownames(ld_mat)
          }

          # Convert to MendelianRandomization package input with correlations
          dat2 <- MendelianRandomization::mr_input(
            bx = x$beta.exposure,
            bxse = x$se.exposure,
            by = x$beta.outcome,
            byse = x$se.outcome,
            corr=ld_mat)

          res <- MendelianRandomization::mr_ivw(object = dat2, correl = TRUE)
          # Reformat results to match TwoSampleMR::mr()
          res <- data.frame(
            id.exposure = unique(x$id.exposure),
            id.outcome = unique(x$id.outcome),
            outcome = unique(x$outcome),
            exposure = unique(x$exposure),
            method = ifelse(length(snplist) == 1,"Wald ratio", "Inverse variance weighted"),
            nsnp = res@SNPs,
            b = res@Estimate,
            se = res@StdError,
            pval = res@Pvalue,
            pair = names(harmonised_studies)[i],
            cis_trans = paste(cis_trans[order(cis_trans)],collapse = "_"),
            IV_set = names(dat))
          
          return(res)
        })
      } else {
        # For all other instrument sets use standard TwoSampleMR::mr()
        res_split <- lapply(dat_split, function(x){

          skip_split <- (nrow(x) == 1 & is.na(x$mr_keep[1])) | (nrow(x) == 1 & x$mr_keep[1] == FALSE)
          if(skip_split == TRUE){
            res <- data.frame(
              id.exposure = NA, id.outcome = NA, outcome = NA, exposure = NA, method = NA, nsnp = NA,
              b = NA, se = NA, pval = NA, pair = names(harmonised_studies)[i], cis_trans = NA, IV_set = names(dat))
            return(res)
          }

          cis_trans <- unique(x$cis_trans)
          res <- TwoSampleMR::mr(x)
          res$pair <- names(harmonised_studies)[i]
          cis_tans <- unique(x$cis_trans)
          res$cis_trans <- paste(cis_trans[order(cis_trans)],collapse = "_")
          res$IV_set <- names(dat)
          return(res)
        })
      }
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
  save(results_MR, file = file.path(res_dir, results_file))
}else{
  ## Write out harmonised_study object
  save(results_MR, file = file.path(res_dir, results_file))
}