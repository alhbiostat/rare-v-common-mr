#!/bin/bash

# Run heterogeneity analysis in parallel for each outcome in outcomes.txt
# Usage: bash run_heterogeneity_parallel.sh [HARMONISED_FILE]
# e.g. bash run_heterogeneity_parallel.sh harmonised_studies_GW5e-08_EW1e-04_rarepruned_COMBINED.rda

harmonised_file=$1
scripts_dir="scripts/04_analysis_proteins"

echo "Running heterogeneity analysis in parallel for each outcome..."

while IFS= read -r outcome; do
  echo "Starting job for outcome: $outcome"
  Rscript ${scripts_dir}/001b_heterogenityanalysis.R "$harmonised_file" "$outcome" &
done < outcomes.txt

wait

echo "All heterogeneity analysis jobs completed."