#!/bin/bash

# Run 008_joinharmonisedstudiesCNV.R in parallel for each outcome in phenos_cnvs.txt

source "config.env"

scripts_dir=${project_dir}/scripts/02_dataextraction_proteins
phenos_file=${project_dir}/phenos_missing.txt

# Loop over each line in phenos_cnvs.txt
while IFS= read -r outcome; do
    echo "Starting job for outcome: $outcome"
    Rscript "${scripts_dir}/008_joinharmonisedstudiesCNV.R" "$outcome" &
done < "$phenos_file"

# Wait for all background jobs to finish
wait

echo "All jobs completed."