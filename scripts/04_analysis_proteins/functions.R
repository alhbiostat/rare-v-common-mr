# Filter results for molecular exposure --> complex trait MR to retain pairs where
# the IVW causal estimate for at least one instrument set is significant (at p_thresh)
# and the number of SNPs used for each variant based cis-trans instrument set is >=n_snps

# Input: list containing MR results for all protein exposures --> outcome
# Output: filtered results containing exposure --> outcome pairs with sufficient instruments per set and at least one significant causal estimate

filter_results <- function(results_MR, p_thresh, n_snps){
  
  keep <- lapply(results_MR, function(x){
    # Instruments per instrument set
    total_snps <- x |> dplyr::group_by(IV_set, cis_trans) |> 
      dplyr::summarise(n_snps = unique(nsnp)) |> 
      dplyr::filter(cis_trans == "cis_trans", grepl("variant", IV_set))
    ## Sufficient instruments to keep exposure-outcome pair?:
    keep_snps <- all(total_snps$n_snps >= n_snps)
    
    # Significant IVW estimates
    total_sig <- x |> dplyr::filter(method == "Inverse variance weighted", nsnp >= n_snps)
    ## At least one significant IVW estimate (for instrument sets >= n_snps)?
    keep_sig <- any(total_sig$pval <= p_thresh)
    
    keep <- all(keep_sig,keep_snps)
    return(keep)
  })
  
  keep <- unlist(keep)
  results_filtered <- results_MR[keep]
  return(results_filtered)
}

# Generate summary table of results
