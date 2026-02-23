#!/bin/bash

# Date: 20-02-2026
# Author: A.L.Hanson
# Purpose: Prune rare variants used as O'link protein instruments (pQTL) based on LD in UKB EUR reference

# This script interacts with the UK Biobank Research Analysis Platform (RAP) via the command line
# using dx-toolkit. To access UKB RAP data and software you must first be logged into the platform
# with <dx-login> and have selected the appropriate project with <dx-select>

# Run:
# run using: ./001b_prunerareIVs.sh

# Inputs:
# rareinstruments.tsv (list of rare protein pQTLs)
# ids_EUR.txt (European sample IDs - KING ancestry defined)

# Steps
# 1. Extract required variants (pQTLs) and individuals (EUR unrelated; see repo ukbrap-ieu-gwas for QC) from WES PLINK format files
# 2. Perform variant pruning based on LD (r2 > 0.001)

# Output:
# wes_pQTLs_chr${i}.prune.in .prune.out (independent variants to retain)

# ======================================================== #

source ../../config.env

for i in {21..22}; do
    run_prune="

    # Variant IDs
    awk -F'\t' 'NR>1 && \$2==${i} {print \$1}' rareinstruments.tsv > varset.txt
    sed -i 's/-/:/g' varset.txt

    # Copy plink files
    cp /mnt/project/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants,\ PLINK\ format\ -\ final\ release/ukb23158_c${i}_b0_v1.* .

    # Filter variants and samples
    plink --bfile ukb23158_c${i}_b0_v1 \
    --keep ukb_unrel_whitebritish_QCed.txt \
    --extract varset.txt \
    --make-bed \
    --out tmpset

    # Num vars
    echo Num. pQTLs:
    wc -l varset.txt
    echo Variants available:
    wc -l tmpset.bim

    # Prune variants (exclude variants within 1MB with r2 > 0.3)
    plink --bfile tmpset \
    --indep-pairwise 500kb 50 0.3 \
    --out wes_pQTLs_chr${i}

    # LD matrix
    plink2 --bfile tmpset \
    --r2-unphased square \
    --out wes_pQTLs_chr${i}_ldmx

    rm ukb23158_c${i}_b0_v1.* 
    rm tmpset.*
    rm varset.txt
    "

    dx run swiss-army-knife \
    -iin="${rap_data_dir}/ukb_unrel_whitebritish_QCed.txt" \
    -iin="${rap_data_dir}/WES/pqtls/rareinstruments.tsv" \
    -icmd="${run_prune}" \
    --tag="prune rare pqtls chr${i}" \
    --instance-type="mem2_ssd2_v2_x8" \
	  --destination="${rap_proj}:${rap_data_dir}/WES/pqtls/" \
	  --brief --yes
done