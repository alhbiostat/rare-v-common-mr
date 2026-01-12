# Date: 12-01-2026
# Author: A.L.Hanson
# Purpose: Harmonise Milind et al. pLOF and gene duplication associations for available proteins and complex traits 

# Output: 

# Load environment variables
dotenv::load_dot_env(file = "config.env")
data_dir <- Sys.getenv("data_dir")
res_dir <- Sys.getenv("results_dir")
prot_dir <- file.path(data_dir,"sumstats","proteomics")
project_dir <- Sys.getenv("project_dir")