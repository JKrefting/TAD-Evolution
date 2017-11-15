# Read in chain BED file
readChainFile <- function(species){
  ch_file <- paste0("data/chains/hg19.", species, ".chain.bed")
  return(import(ch_file, seqinfo = hum_seqinfo))
}

# Read in breakpoint BED file
readBPFile <- function(species, threshold){
  bp_file <- paste0("data/breakpoints/hg19.", species, ".", 
                    as.character(format(threshold, scientific = FALSE)), ".bp.bed")
  return(import(bp_file, seqinfo = hum_seqinfo))
}

# ' Extracts domain boundaries, each enlarged by 'boundary_dist' in both directions.
# '
# ' @param domains A GRanges object containing TAD ranges.
# ' @paream boundary_dist An integer indicating the length of boundary enlargement.
# ' @return GRanges containing only the boundaries. 
getBoundaries <- function(domains, boundary_dist){
  boundaries <- c(GRanges(seqnames(domains), IRanges(start(domains) - boundary_dist, start(domains) + boundary_dist), seqinfo = hum_seqinfo),
                  GRanges(seqnames(domains), IRanges(end(domains) - boundary_dist, end(domains) + boundary_dist), seqinfo = hum_seqinfo))
  return(boundaries)
}

# subdivide the regions into bins of equal size and find number of breakpoints that fall into each bin
calcHitsPerBin <- function(breakpoints, regions, n_bins){
  bins <- tile(regions, n_bins)  # GRangesList
  # check which breakpoint falls into which bin
  allBins <- unlist(bins)
  hitsListAll <- countOverlaps(allBins, breakpoints)
  # add up vectors as rows in a matrix
  hitsMatrix <- matrix(hitsListAll, ncol = n_bins, byrow = TRUE)
  hitsPerBin <- colSums(hitsMatrix)
  return(hitsPerBin)
}

# sample n_control random breakpoints for each actual breakpoint and chromosome
sampleBreakpoints <- function(breakpoints_per_chr, hum_seqinfo, n_controls){
  
  starts <- mapply(function(seqnames, count, n_controls) 
    sample(1:seqlengths(hum_seqinfo[as.character(seqnames)]), 
           count * n_controls, replace = TRUE), 
    breakpoints_per_chr$seqnames, breakpoints_per_chr$count, n_controls, SIMPLIFY = TRUE)
  
  seqnames <- mapply(function(seqnames, count, n_controls) 
    rep(seqnames, count * n_controls),
    breakpoints_per_chr$seqnames, breakpoints_per_chr$count, n_controls)
  
  rdm_breakpoints <- GRanges(unlist(seqnames), IRanges(unlist(starts), width = 1), seqinfo = hum_seqinfo)
  
  return(rdm_breakpoints)
}

# ' Querry the Transcription Start Sites (TSS) from Ensembl (via biomaRt).
# ' 
# ' @param ensembl Gene ensembl
# ' @param geneid_list List of gene identifiers to querry for (optional)
# ' @return TSS of <longest transcript> for each querried gene 