# Functions for handling molecular exposure --> outcome MR results across rare and common variant instrument sets

IV_sets <- c(
  "dhindsa_exwas_mask_raredmg",
  "dhindsa_exwas_mask_ptv",
  "dhindsa_exwas_mask_ptvraredmg",
  "dhindsa_exwas_mask_syn",
  "dhindsa_exwas_variant_common",
  "dhindsa_exwas_variant_rare",
  "dhindsa_exwas_variant_ultrarare",
  "dhindsa_exwas_variant_rare_filt",
  "dhindsa_exwas_variant_ultrarare_filt",
  "sun_gwas_variant_common_finemapped",
  "sun_gwas_variant_common_clumped")

# Generate summary table of results detailing:
# min_p: the minimum p-value for IVW estimates across instrument sets
# min_p_beta: the corresponding effect estimate
# thresh_p: whether the p-value passes a given significance threshold

# Input: 
# results_MR: list containing MR results for all protein exposures --> outcome, 
# p_thresh: p-value threshold to call an IVW estimate significant 
# n_snps: number of instruments in set in order to consider results

summarise_results <- function(results_MR, p_thresh, n_snps = 0){
  
  min_p <- do.call("rbind", lapply(results_MR, function(x){
    df <- x |> na.exclude() |>
      dplyr::filter(method %in% c("Inverse variance weighted"), nsnp >= n_snps) |> 
      dplyr::group_by(IV_set) |> 
      dplyr::filter(pval == min(pval))
    
    out <- df$pval
    names(out) <- df$IV_set
    
    out2 <- out[IV_sets]
    names(out2) <- IV_sets
    
    return(out2)
  })) 
  
  min_p_beta <- do.call("rbind", lapply(results_MR, function(x){
    df <- x |> na.exclude() |>
      dplyr::filter(method %in% c("Inverse variance weighted"), nsnp >= n_snps) |> 
      dplyr::group_by(IV_set) |> 
      dplyr::filter(pval == min(pval))
    
    out <- df$b
    names(out) <- df$IV_set
    
    out2 <- out[IV_sets]
    names(out2) <- IV_sets
    
    return(out2)
  })) 
  
  # Which IVW estimates significant at threshopld
  thresh_p <- min_p < p_thresh
  
  return(list(min_p = min_p, min_p_beta = min_p_beta, thresh_p = thresh_p))
}

# Filter results for molecular exposure --> complex trait MR to retain pairs where
# the IVW causal estimate for at least one instrument set is significant (at p_thresh)
# and the number of SNPs used for each variant based cis-trans instrument set is >=n_snps

# Input: list containing MR results for all protein exposures --> outcome
# Output: filtered results containing exposure --> outcome pairs with sufficient instruments per set (for variants) and at least one significant causal estimate

filter_results <- function(results_MR, p_thresh, n_snps = 0){
  
  keep <- lapply(results_MR, function(x){
    # Instruments per instrument set
    total_snps <- x |>
      dplyr::summarise(n_snps = unique(nsnp), .by = c(IV_set, cis_trans)) |> 
      dplyr::filter(grepl("variant", IV_set), !grepl("_filt", IV_set)) |> 
      dplyr::group_by(IV_set) |>
      dplyr::summarise(max_snps = max(n_snps)) # Max across cis, trans or cis_trans
    ## Sufficient instruments in all IV sets to keep exposure-outcome pair?:
    keep_snps <- all(total_snps$max_snps >= n_snps)
    
    # Significant IVW estimates
    total_sig <- x |> dplyr::filter(
      grepl("variant", IV_set), 
      !grepl("_filt", IV_set),
      method == "Inverse variance weighted", 
      nsnp >= n_snps)
    ## At least one significant IVW estimate (for instrument sets >= n_snps)?
    keep_sig <- any(total_sig$pval <= p_thresh)
    
    keep <- all(keep_sig,keep_snps)
    return(keep)
  })
  
  keep <- unlist(keep)
  results_filtered <- results_MR[keep]
  return(results_filtered)
}


