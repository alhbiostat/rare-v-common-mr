# Date: 07-08-2025
# Author: A.L.Hanson
# Purpose: Read in harmonised summary statistics for a given set of molecular exposures and outcome, and combine 
# into a single object across instrument sets for MR
# To run: Rscript 005_joinharmoinisedstudies.R [OUTCOME]

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
  stop("A single outcome must be specified", call.=FALSE)
} else{
  message("Outcome: ", args[1])
}

outcome <- args[1]
harmonised_file <- paste0("harmonised_studies_", outcome, ".rda")

# Load harmonised study data and reformat into single list of lists (each study pairing split into instrument classes)

filelist <- list.files(path = data_dir, pattern = paste0("_",outcome,".rda"))

if(harmonised_file %in% filelist){

  # Load previous harmonised_studies object to append to
  load(file.path(data_dir, harmonised_file))
  harmonised_studies_old <- harmonised_studies

  filelist <- filelist[!(grepl(harmonised_file, filelist))]

  # If new study pairs exist that are not in harmonised_studies_old:
  # Exposure - outcome study pairs
  study_pairs <- sub("^.*_.*_.*_(.*_.*).rda","\\1",filelist) |> unique()
  # New study pairs
  study_pairs <- study_pairs[!(study_pairs %in% names(harmonised_studies_old))]

  if(length(study_pairs) == 0){
    message("No new study pairs to harmonise")
  }else{
    message("Num. new studies to harmonise: ", length(study_pairs))
    # Copy original harmonised_studies object to backup
    file.copy(file.path(data_dir, harmonised_file), 
            file.path(data_dir, paste0("harmonised_studies_", outcome, "_backup_", Sys.Date(), ".rda")))
    # New harmonised_studies object
    harmonised_studies <- as.list(rep(NA, length(study_pairs)))
    names(harmonised_studies) <- study_pairs
  }
}else{
  # Exposure - outcome study pairs
  study_pairs <- sub("^.*_.*_.*_(.*_.*).rda","\\1",filelist) |> unique()
  study_pairs <- study_pairs[order(study_pairs)]

  harmonised_studies <- as.list(rep(NA, length(study_pairs)))
  names(harmonised_studies) <- study_pairs
}

message("Joining harmonised study data across instrument sets")

for (pair in study_pairs){
  message("exposure:outcome: ", pair)
  files <- file.path(data_dir, 
                  grep(paste0("_",pair,".rda"), list.files(path = data_dir), value = T))
  file_source <- sub(paste0("_",pair,".rda"),"",basename(files))
  names(files) <- file_source
  
  for (file in file_source){
    env <- new.env()
    
    if(file == "dhindsa_exwas_mask"){
      load(files[file], envir = env)
      obj <- as.list(env)
      obj1 <- do.call(c, obj)
      names(obj1) <- paste0("dhindsa_exwas_mask_",names(obj$dat_harmonised))
    }
    
    else if(file == "dhindsa_exwas_variant"){
      load(files[file], envir = env)
      obj <- as.list(env)
      obj2 <- do.call(c, obj)
      names(obj2) <- paste0("dhindsa_exwas_variant_",names(obj$dat_harmonised))
    }
    
    else if(file == "sun_gwas_variant"){
      load(files[file], envir = env)
      obj <- as.list(env)
      obj3 <- do.call(c, obj)
      names(obj3) <- paste0("sun_gwas_variant",sub("exposure","",names(obj$dat_harmonised)))
    }
  }

  # If no results for given instrument set, set to empty list
  if(!exists("obj1")){
    obj1 <- replicate(4, data.frame(), simplify = F)
    names(obj1) <- paste0("dhindsa_exwas_mask_",c("raredmg","ptv","ptvraredmg","syn"))
  }
  if(!exists("obj2")){
    obj2 <- replicate(5, data.frame(), simplify = F)
    names(obj2) <- paste0("dhindsa_exwas_variant_",c("common","rare","ultrarare","rare_filt","ultrarare_filt"))
  }
  if(!exists("obj3")){
    obj3 <- replicate(2, data.frame(), simplify = F)
    names(obj3) <- paste0("sun_gwas_variant_", c("common_finemapped", "common_clumped"))
  }

  out_studies <- c(obj1,obj2,obj3)
  
  # Add column names to empty dataframes
  col_names <- lapply(out_studies, colnames) |> unlist() |> unique()
  
  if(any(lapply(out_studies,nrow) == 0)){
    empty_df <- which(lapply(out_studies,nrow) == 0)
    for(df in empty_df){
      out_studies[[df]] <- data.frame(matrix(data = NA, nrow = 1, ncol = length(col_names), dimnames = list(NULL,col_names)))}
  }
  
  # Remove instruments with MAF<1x10^-5 (some present in genebass studies)
  out_studies <- lapply(out_studies, function(x){
    if(!is.na(x$eaf.exposure[1])){
      x |> dplyr::filter(eaf.exposure > 1*10^-5 & 1-eaf.exposure > 1*10^-5)
    }else{
      x
    }
  })
  
  harmonised_studies[[pair]] <- out_studies 
}

# Annotations for study sets
instrument_type <- c(rep("mask",4), rep("variant",7))
names(instrument_type) <- names(harmonised_studies[[1]])

instrument_class <- c("missense|LC","pLOF","pLOF|missense|LC","synonymous","exome:common",
                      "exome:rare","exome:ultrarare","exome:rare","exome:ultrarare","genome:common","genome:common")
names(instrument_class) <- names(harmonised_studies[[1]])

rm(obj,obj1,obj2,obj3,out_studies)

## Append new harmonised studies to existing object
if(exists("harmonised_studies_old")){

  harmonised_studies <- c(harmonised_studies_old, harmonised_studies)
  # Reorder studies
  order_studies <- names(harmonised_studies)[order(names(harmonised_studies))]
  harmonised_studies <- harmonised_studies[order_studies]

  ## Write out harmonised_study object
  message("Writing output to: ", harmonised_file)
  save(harmonised_studies, file = file.path(data_dir, harmonised_file))
}else{
  ## Write out harmonised_study object
  message("Writing output to: ", harmonised_file)
  save(harmonised_studies, file = file.path(data_dir, harmonised_file))
}