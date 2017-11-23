# Read in fill BED file
readFillFile <- function(species){
 fi_file <- paste0("data/fills/hg19.", species, ".fill.bed")
  return(import(fi_file, seqinfo = hum_seqinfo))
}

# Read in breakpoint BED file
readBPFile <- function(species, threshold){
  bp_file <- paste0("data/breakpoints/hg19.", species, ".", 
                    as.character(format(threshold, scientific = FALSE)), ".bp.bed")
  return(import(bp_file, seqinfo = hum_seqinfo))
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

# transforms a vector a p-values into asterisks depending on their value
# * = pVal < 0.05, ** = pVal < 0.01, *** = pVal < 0.001 
asterisks <- function(pVals){
  blank <- pVals >= 0.05 
  sig <- pVals < 0.05
  hsig <- pVals < 0.01
  hhsig <- pVals < 0.001
  pVals[blank] <- ""
  pVals[sig] <- "*"
  pVals[hsig] <- "**"
  pVals[hhsig] <- "***"
  return(pVals)
}