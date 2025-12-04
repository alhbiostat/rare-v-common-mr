# Date: 24-11-2025
# Author: A.L.Hanson
# Purpose: Run rare variant allelic series analysis using COAST-SS with Wald Ratios derived for O'link protein cis coding variants
# and 11 outcomes of interest
# This essentially tests whether increading the deleterousness of a rare variant (e.g from benign missense -> deleterious missense -> protein truncating/loss of function)
# increases the strength of the causal association with the outcome
# Variants need to first be categoriesed into these groups using VEP annotations

# See:
# https://www.cell.com/ajhg/fulltext/S0002-9297(25)00368-4
# https://github.com/insitro/AllelicSeries

# Load environment variables
dotenv::load_dot_env(file = "config.env")
data_dir <- Sys.getenv("data_dir")
res_dir <- Sys.getenv("results_dir")
gpmap_data_dir <- Sys.getenv("gpmap_data_dir")

library(AllelicSeries)

# Read in Wald Ratios (see 02_dataextraction_proteins/006_calcwaldratios.R)
load(file.path(res_dir, "results_molecular_waldratios_ciscodingvariants.rda")) #wald_ratios

# Load VEP variant annotations
vep <- data.table::fread(file.path(gpmap_data_dir,"variant_annotation","genebass_vep_rare_variantannotations_hg38_altered.txt"))

# Extract SIFT and PolyPhen calls
annot_col <- "Extra"
#vep[, REF_allele := stringr::str_extract(get(annot_col), "REF_ALLELE=[^;]+")] FIX
vep[, SIFT_token := stringr::str_extract(get(annot_col), "SIFT=[^;]+")]
vep[, PolyPhen_token := stringr::str_extract(get(annot_col), "PolyPhen=[^;]+")]
vep[, SIFT_label := ifelse(is.na(SIFT_token), NA_character_,
                           stringr::str_replace(SIFT_token, "^SIFT=([^\\(]+)\\(.*$", "\\1"))]
vep[, PolyPhen_label := ifelse(is.na(PolyPhen_token), NA_character_,
                               stringr::str_replace(PolyPhen_token, "^PolyPhen=([^\\(]+)\\(.*$", "\\1"))]

# Assign consequences to classes based on McCaw et al. classifications
# PTV = protein truncating variant = 4
# DMV = damaging missense variant = 3
# BTV = benign missense variant = 2
# SV = synonymous varaint = 1
consequences <- unique(unlist(sapply(vep$Consequence, function(x){unlist(strsplit(x,","))})))

ptv <- c("splice_acceptor_variant|splice_donor_variant|stop_gained|frameshift_variant")
missense <- c("missense_variant|inframe_deletion|inframe_insertion|stop_lost|start_lost|protein_altering_variant")
synonymous <- c("synonymous_variant")
noncoding <- c("3_prime_UTR_variant|5_prime_UTR_variant|downstream_gene_variant|upstream_gene_variant|intron_variant")

# SIFT/PolyPhen (missense designation)
dmv <- c("possibly_damaging","probably_damaging","deleterious","deleterious_low_confidence")
bmv <- c("benign","tolerated","tolerated_low_confidence")

# Bind variant annotations for exposure variants pQTLs
vep_annots <- vep |>
  dplyr::select(Uploaded_variation, Allele, Consequence, SIFT_label, PolyPhen_label) |>
  dplyr::mutate(snp_join = tolower(gsub(":|_","-",Uploaded_variation)), .before = "Uploaded_variation")

# Assign numerical annotations (most deleterious last)
vep_annots$annots <- NA
# SV
vep_annots$annots[grepl(synonymous,vep_annots$Consequence)] <- 1
# BMV
vep_annots$annots[
  grepl(missense,vep_annots$Consequence) &
  (vep_annots$PolyPhen_label %in% bmv | is.na(vep_annots$PolyPhen_label) | vep_annots$PolyPhen_label == "unknown") &
  vep_annots$SIFT_label %in% bmv] <- 2
vep_annots$annots[
  grepl(missense,vep_annots$Consequence) &
  vep_annots$PolyPhen_label %in% c("benign","possibly_damaging") &
  (vep_annots$SIFT_label %in% c("deleterious_low_confidence","tolerated","tolerated_low_confidence") | is.na(vep_annots$SIFT_label))] <- 2
vep_annots$annots[
  grepl(missense,vep_annots$Consequence) &
  is.na(vep_annots$PolyPhen_label) &
  is.na(vep_annots$SIFT_label)] <- 2 #If SIFT and PolyPhen calls are missing but variant is missense, assume BMV
# DMV
vep_annots$annots[
  grepl(missense,vep_annots$Consequence) &
  (vep_annots$PolyPhen_label %in% dmv | is.na(vep_annots$PolyPhen_label) | vep_annots$PolyPhen_label == "unknown") &
  vep_annots$SIFT_label %in% dmv] <- 3
vep_annots$annots[
  grepl(missense,vep_annots$Consequence) &
  vep_annots$PolyPhen_label %in% dmv &
  vep_annots$SIFT_label == "tolerated_low_confidence"] <- 3
vep_annots$annots[
  grepl(missense,vep_annots$Consequence) &
  vep_annots$PolyPhen_label %in% dmv &
  (vep_annots$SIFT_label %in% dmv | is.na(vep_annots$SIFT_label))] <- 3
# PTV
vep_annots$annots[grepl(ptv,vep_annots$Consequence)] <- 4

# Append to data
wald_ratios <- lapply(wald_ratios, function(x){
  # Alphabetise SNP names for joining
  pos <- sub("^(.*-.*)-.*-.*","\\1",x$SNP)
  alleles <- x[,c("effect_allele.exposure", "other_allele.exposure")]
  a1 <- tolower(apply(alleles, 1, function(x){x[order(as.character(x))][1]}))
  a2 <- tolower(apply(alleles, 1, function(x){x[order(as.character(x))][2]}))
  x$snp_join = paste(pos,a1,a2,sep = "-")
  
  # Conservative consequence calling
  x <- dplyr::left_join(x,vep_annots,by="snp_join")
  
  # Relaxed consequence calling for any remaining unannotated SNPs
  # non coding variant = 1
  x$annots_relaxed <- x$annots
  x$annots_relaxed[grepl(noncoding, x$Consequence) & is.na(x$annots)] <- 1
  # dissagreement between SIFT and PolyPhen = 2
  x$annots_relaxed[grepl(missense, x$Consequence) & is.na(x$annots)] <- 2

  # Merge all missense variants (non coding/synonymous = 1, missense = 2, pLOF = 3)
  x$annots_missensemerged <- x$annots_relaxed
  x$annots_missensemerged[which(x$annots_relaxed == 3)] <- 2
  x$annots_missensemerged[which(x$annots_relaxed == 4)] <- 3

  # Check the effect alelle matches the VEP annotated allele
  x$allele_match = ifelse(
    nchar(x$effect_allele.exposure) == 1 & nchar(x$other_allele.exposure) == 1, x$effect_allele.exposure == x$Allele, ifelse(
      nchar(x$other_allele.exposure) > 1 & nchar(x$effect_allele.exposure) == 1 & x$Allele == "-", TRUE, ifelse(
        nchar(x$effect_allele.exposure) > 1 & nchar(x$other_allele.exposure) == 1, substr(x$effect_allele.exposure, 2, nchar(x$effect_allele.exposure)) == x$Allele, FALSE)))

  return(x)
  }
)

# How many variants missing annotations?
#lapply(wald_ratios, function(x){table(is.na(x$annots))})
#lapply(wald_ratios, function(x){table(is.na(x$annots_relaxed))})
#lapply(wald_ratios, function(x){table(x$allele_match)})

# Save annotated variant Wald ratios
#save(wald_ratios, file = file.path(res_dir, "results_molecular_waldratios_ciscodingvariants_VEPannotated.rda"))

# Run COAST-SS for using Wald ratios and variant annotations
run_COASTSS <- function(wrs, annot_col, weights){
  # Exclude mismatched allele annotations and split genes
  dat <- wrs |> dplyr::filter(allele_match == TRUE)
  dat <- split(dat, f = dat$exposure)

  res_list <- lapply(dat, function(x){
    message("Exposure: ", unique(x$exposure), " Outcome: ", unique(x$outcome))
    x <- x[!is.na(x[,annot_col]),]
    
    if(nrow(x) == 0){
      res <- data.frame(
        "exposure" = NA,
        "outcome" = NA,
        "annot_count" = NA,
        "annot_meanbeta" = NA,
        "annot_meanabsbeta" = NA,
        "alleleic_skat_p" = NA)
    }else if(length(unique(x[,annot_col])) == 1){
      count_vars <- table(x[,annot_col])
      count_vars_string <- paste0(names(count_vars), "=", count_vars)
      
      res <- data.frame(
        "exposure" = unique(x$exposure),
        "outcome" = unique(x$outcome),
        "annot_count" = count_vars_string, # number of variants in each annotation class
        "annot_meanbeta" = paste(unique(x[,annot_col]), "=", round(mean(x$wald_ratio_b),5)), # mean wald ratio of variants in class
        "annot_meanabsbeta" = paste(unique(x[,annot_col]), "=", round(mean(abs(x$wald_ratio_b)),5)),
        "alleleic_skat_p" = NA)
    }else{
      count_vars <- table(x[,annot_col])
      count_vars_string <- paste(paste0(names(count_vars), "=", count_vars),collapse=";")
      mean_betas <- x |> dplyr::group_by_at(annot_col) |> dplyr::summarise(mean(wald_ratio_b))
      mean_betas_string = paste(paste0(mean_betas[[1]], "=", round(mean_betas[[2]],5)), collapse=";")
      mean_absbetas <- x |> dplyr::group_by_at(annot_col) |> dplyr::summarise(mean(abs(wald_ratio_b)))
      mean_absbetas_string = paste(paste0(mean_absbetas[[1]], "=", round(mean_absbetas[[2]],5)), collapse=";")

      coastss <- AllelicSeries::COASTSS(
        anno = x[,annot_col],
        beta = x$wald_ratio_b,
        se = x$wald_ratio_se,
        maf = x$eaf.exposure,
        weights = weights # Weighting for annotation category
      )
      res <- data.frame(
        "exposure" = unique(x$exposure),
        "outcome" = unique(x$outcome),
        "annot_count" = count_vars_string, # number of variants in each annotation class
        "annot_meanbeta" = mean_betas_string, # mean wald ratio of variants in class
        "annot_meanabsbeta" = mean_absbetas_string,
        "alleleic_skat_p" = coastss@Pvals[coastss@Pvals$test == "allelic_skat","pval"])
    }
    return(res)
    })
  out <- do.call("rbind", res_list)
  return(out)
}

coastss_annots <- lapply(wald_ratios, function(x){
  res = run_COASTSS(wrs = x, annot_col = "annots", weights = c(1,2,3,4)); res |> dplyr::arrange(alleleic_skat_p)})
save(coastss_annots, file = file.path(res_dir, "results_molecular_COASTSS.rda"))

coastss_annots_relaxed <- lapply(wald_ratios, function(x){
  res = run_COASTSS(wrs = x, annot_col = "annots_relaxed", weights = c(1,2,3,4)); res |> dplyr::arrange(alleleic_skat_p)})
save(coastss_annots_relaxed, file = file.path(res_dir, "results_molecular_COASTSS_relaxed.rda"))

coastss_annots_missensemerged <- lapply(wald_ratios, function(x){
  res = run_COASTSS(wrs = x, annot_col = "annots_missensemerged", weights = c(1,2,3)); res |> dplyr::arrange(alleleic_skat_p)})
save(coastss_annots_missensemerged, file = file.path(res_dir, "results_molecular_COASTSS_missensemerged.rda"))
