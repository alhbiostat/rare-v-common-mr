#!/bin/bash

# Split common variant GWAS exposure-outcome pairs to harmonise into chunks and run in parallel
# Usage: source run_parallel_formatting.sh [NUM_CHUNKS]
# Example: source run_parallel_formatting.sh 20

source "config.env"

chunks=${1:-4}  # Default to 4 chunks if not provided
studies=${project_dir}/allstudypairings_proteins.csv
scripts_dir=${project_dir}/scripts/02_dataextraction_proteins
chunk_dir=${data_dir}/tmp_study_chunks

mkdir -p $chunk_dir

grep "sun_gwas" $studies > ${chunk_dir}/allstudypairings_proteins_common.csv
studies_common=${chunk_dir}/allstudypairings_proteins_common.csv

# --- Split file ---

echo "Splitting $studies_common into $chunks chunks..."

# Remove header and count lines
lines=$(cat $studies_common | wc -l)
lines_per_chunk=$(( (lines + chunks - 1) / chunks ))

# Split and add header to each
cat $studies_common | split -l $lines_per_chunk - ${chunk_dir}/chunk_
header=$(head -n 1 $studies)
for file in ${chunk_dir}/chunk_*; do
    sed -i "1s/^/$header\n/" "$file"
done

# --- Run chunks in parallel ---
echo "Running in parallel..."
for file in ${chunk_dir}/chunk_*; do
    source ${scripts_dir}/runformatting_proteins.sh $file &
done

wait
echo "All chunks completed"

# --- Clean up ---
rm -r $chunk_dir

echo "DONE"