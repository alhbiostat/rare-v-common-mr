# Date: 07-08-2025
# Author: A.L.Hanson
# Purpose: Read in harmonised CNV mask summary statistics and wald ratios for a given set of molecular exposures and outcomes, and combine 
# into a single object across instrument sets

# N.b Wald ratios have been derived for all available CVN masks, not just those with exposure associations

# To run: Rscript 008_joinharmoinisedstudiesCNV.R [OUTCOME] 

library(dotenv)
library(here)
library(dplyr)
library(ggplot2)

dotenv::load_dot_env(file = here("config.env"))
data_dir <- file.path(Sys.getenv("data_dir"),"harmonised","cnvs")
res_dir <- Sys.getenv("results_dir")

args = commandArgs(trailingOnly=TRUE)

# Check for supplied outcome
if (length(args) == 1){
  message("Outcome: ", args[1])
} else{
  stop("A single outcome must be specified", call.=FALSE)
}

outcome <- args[1]
harmonised_file <- paste0("harmonised_studies_", outcome, ".rda") # Amalgamated output file
results_file <- paste0("results_molecular_waldratios_CNVs_", outcome, ".rda") # Copy to results

filelist <- list.files(data_dir, pattern = paste0("*_", outcome, ".rda"))

# Load harmonised study data and reformat into single list of lists containing each protein exposure and CNV mask (deletion, duplication, plof)

if(harmonised_file %in% filelist){

  # Load previous harmonised_studies object to append to
  load(file.path(data_dir, harmonised_file))
  harmonised_studies_old <- harmonised_studies

  filelist <- filelist[!(grepl("harmonised_studies", filelist))]

  # If new study pairs exist that are not in harmonised_studies_old:
  # Exposure - outcome study pairs
  study_pairs <- sub("^[^_]+_[^_]+_[^_]+_(.*)\\.rda$","\\1",filelist) |> unique()
  # New study pairs
  study_pairs <- study_pairs[!(study_pairs %in% names(harmonised_studies_old))]

  if(length(study_pairs) == 0){
    message("No new study pairs to harmonise")
  }else{
    message("Num. new studies to harmonise: ", length(study_pairs))
    # Copy original harmonised_studies object to backup
    file.copy(file.path(data_dir, harmonised_file), 
            file.path(data_dir, "backup", paste0("harmonised_studies_", outcome, "_backup_", Sys.Date(), ".rda")))
    # New harmonised_studies object
    harmonised_studies <- as.list(rep(NA, length(study_pairs)))
    names(harmonised_studies) <- study_pairs
  }
}else{
  # Exposure - outcome study pairs
  filelist <- filelist[!(grepl("harmonised_studies", filelist))]

  study_pairs <- sub("^[^_]+_[^_]+_[^_]+_(.*)\\.rda$","\\1",filelist) |> unique()
  study_pairs <- study_pairs[order(study_pairs)]

  harmonised_studies <- as.list(rep(NA, length(study_pairs)))
  names(harmonised_studies) <- study_pairs
}

message("Joining harmonised study data across instrument sets")
#study_pairs <- study_pairs[1:5] #test

for (pair in study_pairs){
  message("exposure:outcome: ", pair)
  files <- file.path(data_dir, 
                  grep(paste0("_",pair,".rda"), list.files(path = data_dir), value = T))
  file_source <- sub(paste0("_mask_",pair,".rda"),"",basename(files))
  names(files) <- file_source
  
  for (file in file_source){
    env <- new.env()
    
    if(file == "milind_deletions"){
      load(files[file], envir = env)
      obj <- as.list(env)
      obj1 <- do.call(c, obj)
      names(obj1) <- names(obj$wald_ratios)
    }
    
    else if(file == "milind_duplications"){
      load(files[file], envir = env)
      obj <- as.list(env)
      obj2 <- do.call(c, obj)
      names(obj2) <- names(obj$wald_ratios)
    }
    
    else if(file == "milind_plofs"){
      load(files[file], envir = env)
      obj <- as.list(env)
      obj3 <- do.call(c, obj)
      names(obj3) <- names(obj$wald_ratios)
    }
  }

  # If no results for given instrument set, set to empty list
  if(!exists("obj1")){
    obj1 <- replicate(1, data.frame(), simplify = F)
    names(obj1) <- "milind_deletions"
  }
  if(!exists("obj2")){
    obj2 <- replicate(1, data.frame(), simplify = F)
    names(obj2) <- "milind_duplications"
  }
  if(!exists("obj3")){
    obj3 <- replicate(1, data.frame(), simplify = F)
    names(obj3) <- "milind_plofs"
  }

  out_studies <- c(obj1,obj2,obj3)
  
  # Add column names to empty dataframes
  col_names <-  c("SNP", "effect_allele.exposure", "other_allele.exposure",
  "effect_allele.outcome", "other_allele.outcome", "beta.exposure", "beta.outcome",
  "eaf.exposure", "eaf.outcome", "remove", "palindromic", "ambiguous",
  "id.outcome", "chr.outcome", "pos.outcome", "se.outcome", "outcome", "pval.outcome",
  "mr_keep.outcome", "pval_origin.outcome", "chr.exposure",
  "pos.exposure", "se.exposure", "pval.exposure", "exposure", "mr_keep.exposure",
  "pval_origin.exposure", "id.exposure", "cis_trans", "action", "SNP_index", "mr_keep", "samplesize.outcome",
  "test", "wald_ratio_b", "wald_ratio_se", "wald_ratio_pval", "wald_ratio_nsnp")
  
  if(any(lapply(out_studies,nrow) == 0)){
    empty_df <- which(lapply(out_studies,nrow) == 0)
    for(df in empty_df){
      out_studies[[df]] <- data.frame(matrix(data = NA, nrow = 1, ncol = length(col_names), dimnames = list(NULL,col_names)))}
  }
    
  harmonised_studies[[pair]] <- out_studies

  rm(obj,obj1,obj2,obj3,out_studies) 
}

## Append new harmonised studies to existing object
if(exists("harmonised_studies_old")){

  harmonised_studies <- c(harmonised_studies_old, harmonised_studies)
  # Reorder studies
  order_studies <- names(harmonised_studies)[order(names(harmonised_studies))]
  harmonised_studies <- harmonised_studies[order_studies]

  ## Write out harmonised_study object
  message("Writing output to: ", harmonised_file)
  save(harmonised_studies, file = file.path(data_dir, harmonised_file))
  wald_ratios <- harmonised_studies
  save(wald_ratios, file = file.path(res_dir, results_file))
}else{
  ## Write out harmonised_study object
  message("Writing output to: ", harmonised_file)
  save(harmonised_studies, file = file.path(data_dir, harmonised_file))
  wald_ratios <- harmonised_studies
  save(wald_ratios, file = file.path(res_dir, results_file))
}