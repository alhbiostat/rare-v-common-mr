# Date: 12-01-2026
# Author: A.L.Hanson
# Purpose: Harmonise Milind et al. pLOF and gene duplication associations for available proteins and complex traits 

# Output: 

args <- commandArgs(trailingOnly = TRUE)
# Test args:
args <- c("cardiometabolic/acan.deletions.WES.regenie.gz","p102.deletions.WES.regenie.gz","ACAN","Pulse rate, automated reading","milind_deletions","milind_deletions","mask")

exposure_study <- args[1]
outcome_study <- args[2]
exposure_name <- args[3]
outcome_name <- args[4]
exposure_source <- args[5]
outcome_source <- args[6]
class <- args[7]

library(dotenv)
library(TwoSampleMR)
library(dplyr)

# Load environment variables
dotenv::load_dot_env(file = "config.env")
data_dir <- Sys.getenv("data_dir")
res_dir <- Sys.getenv("results_dir")
prot_dir <- file.path(data_dir,"sumstats","Milind_CNV_associations")
project_dir <- Sys.getenv("project_dir")