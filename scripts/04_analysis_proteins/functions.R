# Functions for handling molecular exposure --> outcome MR results across rare and common variant instrument sets

# IV_sets <- c(
#   "dhindsa_exwas_mask_raredmg",
#   "dhindsa_exwas_mask_ptv",
#   "dhindsa_exwas_mask_ptvraredmg",
#   "dhindsa_exwas_mask_syn",
#   "dhindsa_exwas_variant_common",
#   "dhindsa_exwas_variant_rare",
#   "dhindsa_exwas_variant_ultrarare",
#   "dhindsa_exwas_variant_rare_filt",
#   "dhindsa_exwas_variant_ultrarare_filt",
#   "sun_gwas_variant_common_finemapped",
#   "sun_gwas_variant_common_clumped")

library(ggplot2)
library(dplyr)
library(patchwork)

## =============================================================================================================== ##

# Generate summary table of results detailing:
# min_p: the minimum p-value for IVW estimates across instrument sets
# min_p_beta: the corresponding effect estimate
# thresh_p: whether the p-value passes a given significance threshold

# Input: 
# results_MR: list containing MR results for all protein exposures --> outcome, 
# p_thresh: p-value threshold to call an IVW estimate significant 
# n_snps: number of instruments in set in order to consider results

summarise_results <- function(results_MR, p_thresh, n_snps = 0, cis_trans = "any"){
  
  min_p <- do.call("rbind", lapply(results_MR, function(x){
    if(cis_trans == "any"){
      df <- x |> na.exclude() |>
        dplyr::filter(method %in% c("Inverse variance weighted"), nsnp >= n_snps) |> 
        dplyr::group_by(IV_set) |> 
        dplyr::filter(pval == min(pval))
    }
    if(cis_trans == "cis"){
      df <- x |> na.exclude() |>
        dplyr::filter(method %in% c("Inverse variance weighted"), nsnp >= n_snps, cis_trans == "cis") |> 
        dplyr::group_by(IV_set) |> 
        dplyr::filter(pval == min(pval))
    }
    if(cis_trans == "trans"){
      df <- x |> na.exclude() |>
        dplyr::filter(method %in% c("Inverse variance weighted"), nsnp >= n_snps, cis_trans == "trans") |> 
        dplyr::group_by(IV_set) |> 
        dplyr::filter(pval == min(pval))
    }
    out <- df$pval
    names(out) <- df$IV_set
    
    out2 <- out[IV_sets]
    names(out2) <- IV_sets
    
    return(out2)
  })) 
  
  min_p_beta <- do.call("rbind", lapply(results_MR, function(x){
    if(cis_trans == "any"){
      df <- x |> na.exclude() |>
        dplyr::filter(method %in% c("Inverse variance weighted"), nsnp >= n_snps) |> 
        dplyr::group_by(IV_set) |> 
        dplyr::filter(pval == min(pval))
    }
    if(cis_trans == "cis"){
      df <- x |> na.exclude() |>
        dplyr::filter(method %in% c("Inverse variance weighted"), nsnp >= n_snps, cis_trans == "cis") |> 
        dplyr::group_by(IV_set) |> 
        dplyr::filter(pval == min(pval))
    }
    if(cis_trans == "trans"){
      df <- x |> na.exclude() |>
        dplyr::filter(method %in% c("Inverse variance weighted"), nsnp >= n_snps, cis_trans == "trans") |> 
        dplyr::group_by(IV_set) |> 
        dplyr::filter(pval == min(pval))
    }
    
    out <- df$b
    names(out) <- df$IV_set
    
    out2 <- out[IV_sets]
    names(out2) <- IV_sets
    
    return(out2)
  })) 
  
  # Which IVW estimates significant at threshold
  thresh_p <- min_p < p_thresh
  
  return(list(min_p = min_p, min_p_beta = min_p_beta, thresh_p = thresh_p))
}

## =============================================================================================================== ##

# Filter results for molecular exposure --> complex trait MR to retain pairs where
# the IVW causal estimate for at least one instrument set is significant (at p_thresh)
# and the number of SNPs used for each variant based cis-trans instrument set is >=n_snps

# Input: list containing MR results for all protein exposures --> outcome
# Output: filtered results containing exposure --> outcome pairs with sufficient instruments per set (for variants) and at least one significant causal estimate

filter_results <- function(results_MR, p_thresh, n_snps = 0, cistrans = "cis"){
  
  keep <- lapply(results_MR, function(x){
    # Instruments per instrument set
    total_snps <- x |>
      dplyr::summarise(n_snps = unique(nsnp), .by = c(IV_set, cis_trans)) |> 
      dplyr::filter(grepl("variant", IV_set) & 
                      !grepl("_filt", IV_set) & 
                      cis_trans %in% cistrans) |> 
      dplyr::group_by(IV_set) |>
      dplyr::summarise(max_snps = max(n_snps)) # Max across cis, trans or cis_trans
    ## Sufficient instruments in all IV sets to keep exposure-outcome pair?:
    keep_snps <- all(total_snps$max_snps >= n_snps)
    
    # Significant IVW estimates
    total_sig <- x |> 
      dplyr::filter(grepl("variant", IV_set) & 
                      !grepl("_filt", IV_set) & 
                      cis_trans %in% cistrans &
                      method == "Inverse variance weighted" &
                      nsnp >= n_snps)
    ## At least one significant IVW estimate (for instrument sets >= n_snps)?
    keep_sig <- any(total_sig$pval <= p_thresh)
    
    keep <- all(keep_sig,keep_snps)
    return(keep)
  })
  
  keep <- names(which(keep == TRUE))
  results_filtered <- results_MR[keep]
  return(results_filtered)
}

## =============================================================================================================== ##

# Test for significant heterogeneity between a pair of effect estimates

effect_heterogeneity <- function(beta_vec, se_vec){
  # Weights
  w <- 1/(se_vec^2)
  
  # Mean effect
  beta <- sum(beta_vec*w)/sum(w)
  se <- sqrt(1/sum(w))
  pval <- pnorm(abs(beta/se), lower.tail = FALSE)
  
  # Heterogeneity
  Qj <- w * (beta-beta_vec)^2
  Q <- sum(Qj)
  Qdf <- length(beta_vec)-1
  
  if(Qdf == 0) Q <- 0
  Qjpval <- pchisq(Qj, 1, lower.tail=FALSE)
  Qpval <- pchisq(Q, Qdf, lower.tail=FALSE)
  return(c(Q=Q, Qpval=Qpval))
}

## =============================================================================================================== ##

# Plot IV positions on gene
library(GenomicRanges)

genome_axis <- GenomeAxisTrack()
# MANE-select transcripts
mane_ensembl <- rtracklayer::import(file.path(data_dir, "MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz"))

plot_genomefeatures <- function(chr, regionstart, regionstop, marker_positions){
  
  region <- GRanges(seqnames = chr, ranges = IRanges(regionstart, regionstop))
  # Extract reference annotations in region
  mane_region <- mane_ensembl[mane_ensembl %over% region]
  mane_region <- mane_region[mcols(mane_region)$type == "exon" & mcols(mane_region)$tag == "MANE_Select"]
  
  # Genome ranges
  grtrack <- GeneRegionTrack(
    mane_region,
    chromosome = chr,
    name = "Gene Model",
    transcript = "transcript_id",
    exon = "exon_id",
    gene = "gene_id",
    symbol = "gene_name",
    stacking = "squish",
    collapseTranscripts = "meta",
    showId = FALSE
  )
  
  displayPars(grtrack) <- list(
    transcriptAnnotation = "gene_name",
    just.group = "above",
    min.height = 20)
  
  # Convert to snp positons to GRanges object
  snp_gr <- GRanges(
    seqnames = marker_positions$chr,
    ranges = IRanges(start = marker_positions$pos, width = 1))
  
  # SNP annotation track
  snptrack <- AnnotationTrack(
    snp_gr,
    name = "SNPs",
    shape = "box",
    stacking = "dense")
  
  plotTracks(list(genome_axis, grtrack, snptrack), from = regionstart, to = regionstop,
             sizes = c(1, 3, 1))
}

## ggplot2 version 

plot_genomefeatures_gg <- function(chr, regionstart, regionstop, marker_positions) {
  
  region <- GRanges(seqnames = chr, ranges = IRanges(regionstart, regionstop))
  # Extract reference annotations in region
  #mane_region <- mane_ensembl[mane_ensembl %over% region]
  genes_in <- unique(mane_ensembl[mane_ensembl %over% region]$gene_name)
  mane_region <- mane_ensembl[mane_ensembl$gene_name %in% genes_in]
  mane_region <- mane_region[mcols(mane_region)$type == "exon" & mcols(mane_region)$tag == "MANE_Select"]
  
  exons <- as.data.frame(mane_region) |>
    select(start, end, strand, gene_name)
  
  # Order genes by their leftmost exon start; assign numeric y positions
  gene_extents <- exons |>
    group_by(gene_name, strand) |>
    summarise(gene_start = min(start), gene_end = max(end), .groups = "drop") |>
    arrange(gene_start) |>
    mutate(y_pos = row_number())
  gene_extents$gene_start <- sapply(gene_extents$gene_start, function(x){max(x, regionstart)})
  gene_extents$gene_end <- sapply(gene_extents$gene_end, function(x){min(x, regionstop)})
  
  exons <- exons |>
    left_join(gene_extents |> select(gene_name, y_pos), by = "gene_name")
  
  # --- Gene model panel ---------------------------------------------------------
  p_genes <- ggplot() +
    # Thin backbone spanning full gene extent
    geom_segment(
      data = gene_extents,
      aes(x = gene_start, xend = gene_end, y = y_pos, yend = y_pos),
      linewidth = 0.4, colour = "grey30"
    ) +
    # Exon boxes
    geom_rect(
      data = exons,
      aes(xmin = start, xmax = end, ymin = y_pos - 0.3, ymax = y_pos + 0.3),
      fill = "steelblue", colour = NA
    ) +
    scale_x_continuous(
      limits  = c(regionstart-1000, regionstop+1000),
      expand  = c(0, 0),
      labels  = scales::label_number(scale = 1e-6, suffix = " Mb")
    ) +
    scale_y_continuous(
      breaks = gene_extents$y_pos,
      labels = gene_extents$gene_name,
      expand = expansion(add = 0.7)
    ) +
    labs(x = NULL, y = NULL) +
    theme_bw() +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text.x        = element_blank(),
      axis.ticks.x       = element_blank()
    )
  
  # --- SNP track ----------------------------------------------------------------
  snp_df <- marker_positions |>
    filter(chr == !!chr, pos >= regionstart, pos <= regionstop)
  
  p_snps <- ggplot(snp_df, aes(x = pos)) +
    geom_segment(aes(xend = pos, y = 0, yend = 1, colour = IV_set), linewidth = 0.8) +
    scale_x_continuous(
      limits = c(regionstart-1000, regionstop+1000),
      expand = c(0, 0),
      labels = scales::label_number(scale = 1e-6, suffix = " Mb")
    ) +
    scale_y_continuous(breaks = NULL, expand = c(0, 0)) +
    scale_colour_manual(values = cols_IVsets) +
    labs(x = paste0(gsub("chr", "Chr ", chr), " position"), y = "IVs") +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.grid   = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y  = element_text(size = 8, colour = "grey40")
    )
  
  # --- Combine ------------------------------------------------------------------
  n_genes <- nrow(gene_extents)
  p_genes / p_snps + plot_layout(heights = c(max(n_genes, 2), 1))
}

## =============================================================================================================== ##

# Bind two nested lists with the same structure
bind_nested <- function(x, y) {
  if (is.data.frame(x) && is.data.frame(y)) {
    bind_rows(x, y)
  } else if (is.list(x) && is.list(y)) {
    Map(bind_nested, x, y)
  } else {
    stop("Structures do not match.")
  }
}


