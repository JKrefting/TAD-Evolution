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

# ' Calculates the number of breakpoints that fall into each bin for each region
# ' 
# ' @param breakpoints GRanges with breakpoint coordinates
# ' @param region_bins GRangesList where every region is subdivided into n_bins
# ' @param n_bins Integer indicating the number of bins each region is subdivided into
# ' @return The hits per bin summed up over all regions
calcHitsPerBin <- function(breakpoints, region_bins, n_bins){
  # check which breakpoint falls into which bin
  all_bins <- unlist(region_bins)
  hits_list_all <- countOverlaps(all_bins, breakpoints)
  # add up vectors as rows in a matrix
  hits_matrix <- matrix(hits_list_all, ncol = n_bins, byrow = TRUE)
  # sum up over all regions
  hits_per_bin <- colSums(hits_matrix)
  
  return(hits_per_bin)
}

# ' Sample random breakpoints for each actual breakpoint per chromosome
# ' 
# ' @param breakpoints_per_chr Dataframe providing the number of actual breakpoints per chromosome
# ' @param region_bins hum_seqinfo Seqinfo object for the construction of random breakpoint GRanges 
# ' @return GRanges containing the sampled random breakpoints
sampleBreakpoints <- function(breakpoints_per_chr, hum_seqinfo){
  
  breakpoints_per_chr <- breakpoints_per_chr %>% 
    mutate(
      chrom_len = seqlengths(hum_seqinfo[seqnames])
      )
  
  # sample the coordinates of random breakpoints
  starts <- map2(breakpoints_per_chr$chrom_len, breakpoints_per_chr$count,
                 sample.int)
  
  seqnames <- rep(breakpoints_per_chr$seqnames, breakpoints_per_chr$count)  
  
  rdm_breakpoints <- GRanges(seqnames, IRanges(unlist(starts), width = 1), seqinfo = hum_seqinfo)
  
  return(rdm_breakpoints)
}

# ' A vector of p values is transformed into a vector of asterisks representing significance levels
# ' 
# ' @param p_vals Integer (vector) containing p values
# ' @return String (vector) containing asterisks for the respective significance level
# none = p_val > 0.05; * = p_val <= 0.05; ** = p_val <= 0.01; *** = p_val <= 0.001 
asterisks <- function(p_vals){
  blank <- p_vals > 0.05 
  sig <- p_vals <= 0.05
  hsig <- p_vals <= 0.01
  hhsig <- p_vals <= 0.001
  p_vals[blank] <- ""
  p_vals[sig] <- "*"
  p_vals[hsig] <- "**"
  p_vals[hhsig] <- "***"
  
  return(p_vals)
}

# ' Calculate the standard error of the mean
# ' source: https://stackoverflow.com/questions/2676554/in-r-how-to-find-the-standard-error-of-the-mean
# ' @param x Numeric vector
# ' @return Numeric indicating standard error of x
std <- function(x) sd(x)/sqrt(length(x))

# ' Capitalize first letters of all words in a string
# ' source: https://stackoverflow.com/questions/6364783/capitalize-the-first-letter-of-both-words-in-a-two-word-string
# ' @param x String
# ' @return String
simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}