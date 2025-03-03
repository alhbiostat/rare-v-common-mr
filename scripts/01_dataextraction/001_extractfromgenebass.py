# Date: 03-03-2025
# Author: A.L.Hanson
# Purpose: Extract gene and variant based summary stats from Genebass Hail Matrix Table format data (stored on P1 server)
# for required complex traits

# Run on P1 within docker container

import pandas as pd
import hail as hl
import os
import sys

os.environ['PYSPARK_SUBMIT_ARGS'] = '--driver-memory 10g --executor-memory 10g pyspark-shell'
DATA_DIR = "/local-scratch/data/hg38/genebass/"
OUT_DIR = ""

# Read in hail matix table
mt = hl.read_matrix_table(os.path.join(DATA_DIR, "variant_results.mt"))
mt_gene = hl.read_matrix_table(os.path.join(DATA_DIR, "results.mt"))

# Extract required studies basesd on description
# target_studies = {'LDL direct',
#                   'Body mass index (BMI)',
#                   'Vitamin D',
#                   'Triglycerides',
#                   'Glycated haemoglobin (HbA1c)',
#                   'Cystatin C',
#                   'Mean platelet (thrombocyte) volume',
#                   'IGF-1',
#                   'Waist-to-hip ratio adjusted for body mass index (WHRadjBMI)',
#                   'Red blood cell (erythrocyte) count',
#                   'Mean corpuscular volume',
#                   'Mean platelet (thrombocyte) volume',
#                   'Coronary artery disease (CAD)',
#                   'Date E11 first reported (non-insulin-dependent diabetes mellitus)',
#                   'Multiple Sclerosis (MS)',
#                   'Atrial fibrillation',
#                   'Venous thromboembolism (VTE)',
#                   'Prostate cancer',
#                   'Date I10 first reported (essential (primary) hypertension)'}

# phenocode
target_studies = {'30780',
                  '21001',
                  '30890',
                  '30870',
                  '30750',
                  '30720',
                  '30100',
                  '30770',
                  'WHRadjBMI_custom',
                  '50',
                  '30010',
                  '30040',
                  '30100',
                  'CAD_custom',
                  '130708',
                  'MS_custom',
                  'Afib_custom',
                  'VTE_custom',
                  'ProstateCancer_custom',
                  '131286',
                  'Ischemic_stroke_custom'}

mt_set = mt.filter_cols(hl.literal(target_studies).contains(mt.phenocode))
mt_gene_set = mt_gene.filter_cols(hl.literal(target_studies).contains(mt_gene.phenocode))

for key in range(1,len(target_studies)+1):
    col_key = mt_set.cols().take(key)[key-1]

    mt_study = mt_set.filter_cols(
        (mt_set.trait_type == col_key['trait_type']) &
        (mt_set.phenocode == col_key['phenocode']) &
        (mt_set.pheno_sex == col_key['pheno_sex']) &
        (mt_set.coding == col_key['coding']) &
        (mt_set.modifier == col_key['modifier']))
    
    mt_study_gene = mt_gene_set.filter_cols(
        (mt_gene_set.description == col_key['description']) &
        (mt_gene_set.phenocode == col_key['phenocode']) &
        (mt_gene_set.modifier == col_key['modifier']))

    desc =  mt_study.cols().select('description').collect()[0].description
    print("pheno:", desc)

    outname = "".join(["genebass_ukbwes_","p",col_key['phenocode'], "_", desc.replace(" ","_").replace("(", "").replace(")", ""), ".tsv"])
    outname_gene = "".join(["genebass_ukbwes_","p",col_key['phenocode'], "_", desc.replace(" ","_").replace("(", "").replace(")", ""), "_genetest.tsv"])

    # Extract summary statistics (variant level)
    study_stats = mt_study.entries()
    # Extract summary statistics (gene level)
    study_stats_gene = mt_study_gene.entries()

    # Write out rare summary statistics
    print("Writing to:", os.path.join(OUT_DIR, outname))
    study_stats.export(os.path.join(OUT_DIR, outname))
    
    ## FIX PATHS AND MOVE TO P1