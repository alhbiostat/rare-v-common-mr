# Date: 16-03-2026
# Author: A.L.Hanson
# Purpose: LD prune all common-variant GWAS derived IVs using the UKB LD referene panel for use in pleiotropy analyses
# N.B These varaints are taken from those written to pQTL_lookup_sun_gwas_variant_common_clumped.tsv by 04_abalysis_proteins/002_lookuppQTLs.R
## which pulls all common GWAS pQTLS from /data/sumstats/proteomics/common.
## Although all instruments are pruned within each MR analysis, when all instruments are aggregated for assessment of pleiotropic
## effects across pQTLs, some LD likely remains between variants in differing instrument sets

source config.env

# Variant list:
snplist=${data_dir}/sumstats/proteomics/common/gwascommoninstruments_p5e-8.tsv
# Reference data set (UKB)
ref_set=${ld_dir}/full

# Pull required variants
plink2 --bfile $ref_set \
--extract $snplist \
--make-bed \
--out tmpset

# Pairwise pruning of variants (slightly more relaxed than IV selection pruning to retain sufficient coverage of proteins for analysis)
plink2 --bfile tmpset \
--indep-pairwise 500kb 1 0.1 \
--out ${data_dir}/sumstats/proteomics/common/protein_pQTLpruning/gwas_pQTLs_p5e-8

rm tmpset.*
