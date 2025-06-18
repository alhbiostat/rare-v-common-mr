# Date: 10-06-2025
# Author: A.L.Hanson
# Purpose: 
# 1. Convert EBI GWAS Catalog downloaded summary statistics to GWAS VCF format
# 2. Add proxy variants for pQTLs if necessary
# 3. Write out .vcf files to be used in data harmonisation
# Guide: https://mrcieu.github.io/gwasvcf/articles/guide.html

library(gwasvcf)
library(ieugwasr)
library(gwasglue2)
library(magrittr)
library(VariantAnnotation)
library(AnnotationHub)
library(rtracklayer)

# Load environment variables
dotenv::load_dot_env(file = "config.env")
data_dir <- Sys.getenv("data_dir")
ld_dir <- Sys.getenv("ld_dir")

set_plink()
set_bcftools("/home/kf23130/Software/bcftools")

# LD file for proxies
ld_file <- file.path(ld_dir, "full_rsid")

# Read in variant annotations
annots <- data.table::fread(file.path(data_dir,"/variant_annotations/vep_variantannotations_hg38_altered.txt"))

# Read in list of proteomics pQTLs
## THIS MAY NEED TO BE REPEATED WHEN FINEMAPPING IS FIXED (AS WILL SUN GWAS VARIANT EXTRACTION)
pQTLs <- data.table::fread(file.path(data_dir,"sumstats/proteomics/common/all_pQTLs.tsv"), fill = T, header = F)
colnames(pQTLs) <- c("SNP","RSID")

# Read in chain file for liftOver
ah <- AnnotationHub()
chain <- ah[["AH14150"]] #hg19ToHg38

## ---- STEP 1 ----- ##

# Read in .tsv.gz files from GWAS Catalogue and convert to GWAS vcf format
infiles <- list.files(file.path(data_dir,"sumstats/opengwas"), pattern = ".tsv.gz", full.names = T)
sumstats <- lapply(infiles, function(x){data.table::fread(x)})

study_names <- sub(".*(GCST.*)_.*","\\1",infiles)

# Set alphabetised SNP IDs and match rsIDs
# Note, ~1,000,000 SNPs without rsID matches here

sumstats <- lapply(sumstats, function(x){
  a1 <- x$effect_allele
  a2 <- x$other_allele
  
  x <- x |> dplyr::mutate(
    SNP = ifelse(a1 < a2, paste0(chromosome,":",base_pair_location,"_",a1,"_",a2),
      paste0(chromosome,":",base_pair_location,"_",a2,"_",a1)))

  rsid <- annots[match(x$SNP, annots$Uploaded_variation),]$Existing_variation
  rsid <- unlist(lapply(strsplit(rsid, ","), function(x){x[1]}))

  x$RSID <- rsid
  x$RSID <- ifelse(is.na(x$RSID), x$SNP, x$RSID)

  return(x)
})

# Convert to vcf
vcfs <- mapply(x = sumstats, y = study_names, FUN = function(x,y){
  x %$% gwasvcf::create_vcf(
  chrom = chromosome,
  pos = base_pair_location,
  nea = other_allele,
  ea = effect_allele,
  snp = RSID,
  effect = beta,
  se = standard_error,
  pval = p_value,
  name = y)
}, SIMPLIFY = F)

# Write out vcf formatted files
lapply(vcfs, function(x){
  file_name <- file.path(data_dir, "sumstats/opengwas", paste0(samples(header(x)), "_buildGRCh38.vcf"))
  writeVcf(x, file_name, index = T)})

## ---- STEP 2 ----- ##
# Read in all studies .vcf files
infiles_vcf <- grep("vcf.*z$", list.files(file.path(data_dir,"sumstats/opengwas"), full.names = T), value = T)

# Find proxies and add to data, extract pQTLs only
find_proxies <- function(vcf){
  
  message("Processing:", vcf)
  dat <- VariantAnnotation::readVcf(vcf)

  # If genomic build is GRCh37 (Hg19), update to 38
  if(grepl("19|37", unique(genome(dat)))){
    message("Lifting over genome build...")

    # Remove any variants on patches or sex chromosomes
    seq_out <- which(!(seqnames(dat) %in% 1:22))
    
    if(length(seq_out) != 0){
      dat_lift <- dat[-seq_out,]
    }else{dat_lift <- dat}
    
    # Convert to chrX format for liftover
    if(seqlevels(dat_lift)[1] != "chr1"){
      seqlevels(dat_lift) <- paste0("chr", seqlevels(dat_lift))
    }

    # Liftover variants
    gr <- rowRanges(dat_lift)
    lifted <- rtracklayer::liftOver(gr, chain)
    # Remove variants with multiple mappings
    unique_lift <- elementNROWS(lifted) == 1
    lifted_gr <- unlist(lifted[unique_lift])

    # Filter vcf and update coordiantes
    dat_lift <- dat_lift[unique_lift,]
    rowRanges(dat_lift) <- lifted_gr

    # Convert to chrX format back to numeric for proxy lookup
    if(seqlevels(dat_lift)[1] == "chr1"){
      seqlevels(dat_lift) <- sub("chr","",seqlevels(dat_lift))
    }

    dat <- dat_lift
  }

  snps <- names(rowRanges(dat))

  # Find proxies for missing variants
  message("Finding proxies...")
  missing <- pQTLs[which(!(pQTLs$RSID %in% snps)),]$RSID
  proxies <- gwasvcf::query_gwas(dat, rsid = missing, proxies = "yes", bfile = ld_file, tag_r2 = 0.6)
  geno(proxies)$ID <- matrix(rownames(geno(proxies)$ID), ncol = 1) # Set ID as missing SNP, PR as proxy ID

  # Add proxies to original vcf
  geno(dat)$PR <- matrix(NA, nrow = nrow(dat), ncol = ncol(dat))
  PR_field <- DataFrame(Number = 1, Type = "String", Description = "Proxy rsid")
  geno(header(dat))["PR",] <- PR_field

  merged <- rbind(dat, proxies)

  # Keep pQTLS only
  merged_pqtls <- merged[names(rowRanges(merged)) %in% pQTLs$RSID]
  message("Final variant count: ", length(rowRanges(merged_pqtls)))

  # Write out
  file_name <- file.path(data_dir, "sumstats/opengwas", paste0(samples(header(merged_pqtls)), "_pQTLs_wproxies.vcf"))
  writeVcf(merged_pqtls, file_name, index = T)
}

lapply(infiles_vcf, find_proxies)