#!/bin/bash

# Split rare variant CNV mask exposure-outcome pairs to harmonise into chunks and run in parallel
# Usage: source runformatting_proteins_cnvs_parallel.sh [NUM_CHUNKS] [STUDY_PAIR_LIST]
# Example: source runformatting_proteins_cnvs_parallel.sh 20 allstudypairings_cnvs.csv

source "config.env"

chunks=${1:-4}  # Default to 4 chunks if not provided
studies="${project_dir}/${2}"
#studies=${project_dir}/allstudypairings_cnvs.csv
#studies=${project_dir}/testrun.txt
scripts_dir=${project_dir}/scripts/02_dataextraction_proteins
chunk_dir=${data_dir}/tmp_study_chunks

mkdir -p $chunk_dir

# --- Split file ---

echo "Splitting $studies into $chunks chunks..."

# Remove header and count lines
lines=$(cat $studies | wc -l)
lines_per_chunk=$(( (lines + chunks - 1) / chunks ))

# Split and add header to each
cat $studies | split -l $lines_per_chunk - ${chunk_dir}/chunk_
header=$(head -n 1 $studies)
for file in ${chunk_dir}/chunk_*; do
    sed -i "1s/^/$header\n/" "$file"
done

# --- Run chunks in parallel ---
echo "Running in parallel..."
for file in ${chunk_dir}/chunk_*; do
    source ${scripts_dir}/runformatting_proteins_cnvs.sh $file &
done

wait
echo "All chunks completed"

# --- Clean up ---
rm -r $chunk_dir

echo "DONE"