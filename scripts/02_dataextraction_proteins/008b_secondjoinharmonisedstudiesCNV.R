# Date: 07-08-2025
# Author: A.L.Hanson
# Purpose: Combine per-trait results from 008_joinharmonisedstudiesCNV.R into a single list and filter on p-value

library(dotenv)
library(here)
library(dplyr)

dotenv::load_dot_env(file = here("config.env"))
res_dir <- Sys.getenv("results_dir")

files_cnvs <- list.files(res_dir, pattern = "results_molecular_waldratios_CNVs_[A-Z].*", full.names = T)

# Import results
res_CNVs <- list()
for(i in 1:length(files_cnvs)){
  message("Loading", files_cnvs[i])
  load(files_cnvs[i])
  res_CNVs[[i]] <- wald_ratios
}

names(res_CNVs) <- sub(".*_CNVs_(.*)\\.rda","\\1", files_cnvs)

# Determine number of tests run per protein for multiple testing adjustment
n_deletions <- unlist(lapply(res_CNVs, function(x){lapply(x, function(y){nrow(y$milind_deletions)})})) # max = 89
n_duplications <- unlist(lapply(res_CNVs, function(x){lapply(x, function(y){nrow(y$milind_duplications)})})) # max = 308
n_plofs <- unlist(lapply(res_CNVs, function(x){lapply(x, function(y){nrow(y$milind_plofs)})})) #max = 340
n_lofs <- unlist(lapply(res_CNVs, function(x){lapply(x, function(y){nrow(y$milind_lofs)})})) #max = 2104

# Adjust for the maximum number of genes with deletions, duplicatons and plofs available across traits
n_tests <- sum(max(n_deletions),max(n_duplications),max(n_plofs),max(n_lofs)) #737 gene-level tests/trait max

# Alpha threshold
alpha <- signif(0.05/n_tests,2) #6.8e-05
alpha_relaxed <- 0.05

# Filter pval.exposure and save objects
CNV_strict <- lapply(res_CNVs, lapply, lapply,
                     dplyr::filter,
                     pval.exposure <= alpha)

CNV_relaxed <- lapply(res_CNVs, lapply, lapply,
                     dplyr::filter,
                     pval.exposure <= alpha_relaxed)

save(res_CNVs, file = file.path(res_dir, "results_molecular_waldratios_CNVs_alloutcomes.rda"))
save(CNV_strict, file = file.path(res_dir, "results_molecular_waldratios_CNVs_alloutcomes_p7e-05.rda"))
save(CNV_relaxed, file = file.path(res_dir, "results_molecular_waldratios_CNVs_alloutcomes_p0.05.rda"))