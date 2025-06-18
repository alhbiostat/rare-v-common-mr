#!/bin/bash

source config.env
stats_dir="${data_dir}/sumstats/opengwas"

cd $stats_dir

# Download available summary statistics from OpenGWAS in vcf format
for id in `cat studyids.txt`; do
  echo "Downloading:" $id
  wget https://gwas.mrcieu.ac.uk/files/${id}/${id}.vcf.gz
  wget https://gwas.mrcieu.ac.uk/files/${id}/${id}.vcf.gz.tbi
done

# Remaining studies downloaded directly from GWAS Catalogue and converted to GWAS .vcf format using
# convert_to_gwasvcf.R