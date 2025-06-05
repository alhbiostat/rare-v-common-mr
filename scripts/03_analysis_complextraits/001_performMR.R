# Date: 04-04-2025
# Author: A.L.Hanson
# Purpose: Read in harmonised summary statistics across exposure-outcome pairs, perform MR and output results

## CONTINUE:
# Check formatting of harmonised_studies object in the absence of deepRVAT results
# Add new harmonised results to harmonised_studies object and save
# Run MR or new harmonised studies and append to results_MR object

library(dotenv)
library(here)
library(dplyr)
library(ggplot2)

dotenv::load_dot_env(file = here("config.env"))
data_dir <- Sys.getenv("data_dir")
res_dir <- Sys.getenv("results_dir")

# Load harmonised study data and reformat into single list of lists (each study pairing split into instrument classes)

filelist <- list.files(path = file.path(data_dir, "harmonised"), pattern = ".rda")

if("harmonised_studies.rda" %in% filelist){

  # Load previous harmonised_studies object to append to
  load(file.path(data_dir, "harmonised", "harmonised_studies.rda"))
  harmonised_studies_old <- harmonised_studies

  filelist <- filelist[!(grepl("harmonised_studies", filelist))]

  # If new study pairs exist that are not in harmonised_studies_old:
  # Exposure - outcome study pairs
  study_pairs <- sub("^.*_.*_(.*_.*).rda","\\1",filelist) |> unique()
  # New study pairs
  study_pairs <- study_pairs[!(study_pairs %in% names(harmonised_studies_old))]

  if(length(study_pairs) == 0){
    message("No new study pairs to harmonise")
  }else{
    message("New studies to harmonise: ", paste(study_pairs, collapse = ", "))
    # Copy original harmonised_studies object to backup
    file.copy(file.path(data_dir, "harmonised", "harmonised_studies.rda"), 
            file.path(data_dir, "harmonised", paste0("harmonised_studies_backup_", Sys.Date(), ".rda")))
    # New harmonised_studies object
    harmonised_studies <- as.list(rep(NA, length(study_pairs)))
    names(harmonised_studies) <- study_pairs
  }
}else{
  # Exposure - outcome study pairs
  study_pairs <- sub("^.*_.*_(.*_.*).rda","\\1",filelist) |> unique()

  harmonised_studies <- as.list(rep(NA, length(study_pairs)))
  names(harmonised_studies) <- study_pairs
}

for (pair in study_pairs){
  files <- file.path(data_dir,"harmonised", 
                  grep(paste0("_",pair,".rda"), list.files(path = file.path(data_dir, "harmonised")), value = T))
  file_source <- sub(paste0("_",pair,".rda"),"",basename(files))
  names(files) <- file_source
  
  for (file in file_source){
    env <- new.env()
    
    if(file == "deeprvat_genescore"){
      load(files[file], envir = env)
      obj1 <- as.list(env)
      names(obj1) <- "deeprvat_genescore"
    }
    
    else if(file == "genebass_mask"){
      load(files[file], envir = env)
      obj <- as.list(env)
      obj2 <- do.call(c, obj)
      names(obj2) <- paste0("genebass_",names(obj$dat_harmonised))
    }
    
    else if(file == "genebass_variant"){
      load(files[file], envir = env)
      obj <- as.list(env)
      obj3 <- do.call(c, obj)
      names(obj3) <- paste0("genebass_",names(obj$dat_harmonised))
    }
    
    else if(file == "opengwas_variant"){
      load(files[file], envir = env)
      obj4 <- as.list(env)
      names(obj4) <- paste0("opengwas_common")
    }
  }

  # If no deepRVAT results for study pair, set to empty list
  if(!exists("obj1")){
    obj1 <- list(data.frame())
    names(obj1) <- "deeprvat_genescore"
  }
  
  out_studies <- c(obj1,obj2,obj3,obj4)
  
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
instrument_type <- c("genescore",rep("mask",4), rep("variant",7))
names(instrument_type) <- names(harmonised_studies[[1]])

instrument_class <- c("genescore","missense|LC","pLOF","pLOF|missense|LC","synonymous","exome:common","exome:lowfreq",
                      "exome:rare","exome:ultrarare","exome:rare","exome:ultrarare","genome:common")
names(instrument_class) <- names(harmonised_studies[[1]])

rm(obj,obj1,obj2,obj3,obj4,out_studies)

## Append new harmonised studies to existing object
if(exists("harmonised_studies_old")){
  harmonised_studies <- c(harmonised_studies_old, harmonised_studies)
  # Reorder studies
  order_studies <- names(harmonised_studies)[order(names(harmonised_studies))]
  harmonised_studies <- harmonised_studies[order_studies]

  ## Write out harmonised_study object
  save(harmonised_studies, file = file.path(data_dir, "harmonised", "harmonised_studies.rda"))
}else{
  ## Write out harmonised_study object
  save(harmonised_studies, file = file.path(data_dir, "harmonised", "harmonised_studies.rda"))
}

# Perform MR on each exposure-outcome pair, across each of 12 instrument sets

# Check if MR results already exist
if(file.exists(file.path(res_dir, "results_complextrait_MR.rda"))){
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
    
    skip_pair <- (nrow(dat[[1]]) == 1 & dat[[1]]$mr_keep[1] == FALSE)

    if(!is.na(dat[[1]]$SNP[1]) & skip_pair == FALSE){
      res <- TwoSampleMR::mr(dat[[1]])
      res$pair <- names(harmonised_studies)[i]
      res$study <- names(harmonised_studies[[i]])[j]
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
