# =============================================================================================================================
# Analysis of rearrangement breakpoint distributions with regard to stuctural domains of the human genome. 
# =============================================================================================================================

require(devtools)
devtools::install_version("tidyverse", version = "1.1.1", repos = "http://cran.us.r-project.org")
require(tidyverse)
require(BSgenome.Hsapiens.UCSC.hg19)

# read metadata for analysis
METADATA <- read_tsv("metadata.csv")
SPECIES <- select(METADATA, genome_assembly, trivial_name)
DOMAINS <- METADATA %>% 
  select(genomic_domain_type, genomic_domain_path) %>%
  filter(!is.na(genomic_domain_type))
THRESHOLDS <- unlist(METADATA %>% 
  filter(!is.na(min_size_threshold)) %>%
  select(min_size_threshold)
)

# number of bins to subdivide each structural domain
NBINS <- 20

# number of random control breakpoints per actual breakpoint
NCONTROLS <- 100

# size around a boundary that is investigated (only relevant for 2nd analysis)
BOUNDARY_SIZE <- 4*10^5

# Load human seqinfo
genome <- BSgenome.Hsapiens.UCSC.hg19
hum_seqinfo <- seqinfo(genome)

# Load function to read in breakpoint BED file
readBPFile <- function(species, threshold){
  bp_file <- paste0("data/breakpoints/hg19.", species, ".", 
                   as.character(format(threshold, scientific = FALSE)), ".bp.bed")
  return(import(bp_file, seqinfo=hum_seqinfo))
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

# -------------------------------------------------------------------------------------------------------------------
# Analysis of breakpoints around structural domains
# -------------------------------------------------------------------------------------------------------------------

# tibble to gather all results for domain, species, threshold
all_results <- tibble()

for (D in DOMAINS$genomic_domain_path){
  
  domains <- import(unlist(D), seqinfo = hum_seqinfo)
  
  # get domain type to store with every result
  domain_type <- unlist(DOMAINS %>%
                          filter(genomic_domain_path == D) %>%
                          select(genomic_domain_type)
  )
  
  
  # enlarge each domain by 50% of its width to each side
  domains_plus <- resize(domains, fix = "center", width=2*width(domains))
  
  # results fot all species and thresholds concerning a domain set
  results <- tibble()
  
  for (S in SPECIES$genome_assembly){
    
    for (THR in THRESHOLDS){
      
      # tibble to store results for this loop 
      results_loop <- tibble(bin = seq(1, NBINS))
      
      # load breakpoint and TAD bed files
      breakpoints <- readBPFile(S, THR)
      
      # handle case of no breakpoints for threshold
      if (length(breakpoints) < 1){
        results_loop <- cbind.data.frame(results_loop, hits = rep(0, NBINS))
        results_loop <- cbind.data.frame(results_loop, rdm_hits = rep(0, NBINS))
        results_loop <- cbind.data.frame(results_loop, species = rep(factor(S), NBINS))
        results_loop <- cbind.data.frame(results_loop, threshold = rep(factor(THR), NBINS))
        results_loop <- cbind.data.frame(results_loop, domains = rep(factor(domain_type), NBINS))
        results_loop <- cbind.data.frame(results_loop, n_breakpoints = rep(0, NBINS))
        results <- rbind.data.frame(results, results_loop)
        next
      }
      # determine distribution of breakpoints
      # subdivide the genomic regions into bins of equal size and find number of breakpoints that fall into each bin
      hits <- calcHitsPerBin(breakpoints, domains_plus, NBINS)
      results_loop <- cbind.data.frame(results_loop, hits = hits)
      
      # generate random breakpoints
      # count number of breakpoints for each chr in 'breakpoints'
      breakpoints_per_chr <- as.tibble(breakpoints) %>%
        group_by(seqnames) %>%
        summarise(count = n())
      # generate random breakpoints for each chromosome
      rdm_breakpoints <- sampleBreakpoints(breakpoints_per_chr, hum_seqinfo, NCONTROLS)
      
      # determine distribution of random breakpoints
      rdm_hits <- calcHitsPerBin(rdm_breakpoints, domains_plus, NBINS)
      results_loop <- cbind.data.frame(results_loop, rdm_hits = rdm_hits)
      
      # complete result_loop info columns
      results_loop <- cbind.data.frame(results_loop, species = rep(factor(S), NBINS))
      results_loop <- cbind.data.frame(results_loop, threshold = rep(factor(THR), NBINS))
      results_loop <- cbind.data.frame(results_loop, domains = rep(factor(domain_type), NBINS))
      results_loop <- cbind.data.frame(results_loop, n_breakpoints = rep(length(breakpoints), NBINS))
      results <- rbind.data.frame(results, results_loop)
    } 
  } 
  
  all_results <- rbind.data.frame(all_results, results)
  
}

dir.create("results/", showWarnings = FALSE)
all_results <- as.tibble(all_results)
saveRDS(all_results, file = "results/breakpoints_at_domains.rds")

# -------------------------------------------------------------------------------------------------------------------
# Analysis of breakpoints at BOUNDARIES of structural domains
# -------------------------------------------------------------------------------------------------------------------

# tibble to gather all results for domain, species, threshold
all_results <- tibble()

for (D in DOMAINS$genomic_domain_path){
  
  domains <- import(unlist(D), seqinfo = hum_seqinfo)
  
  # get domain type to store with ever result
  domain_type <- unlist(DOMAINS %>%
                          filter(genomic_domain_path == D) %>%
                          select(genomic_domain_type)
  )
  
  # extract boundaries
  boundaries <- c(GRanges(seqnames(domains), 
                          IRanges(start(domains), start(domains)), 
                          seqinfo = hum_seqinfo),
                  GRanges(seqnames(domains), 
                          IRanges(end(domains), end(domains)), 
                          seqinfo = hum_seqinfo)
                  )  
  
  # enlarge
  boundaries_plus <- resize(boundaries, fix = "center", BOUNDARY_SIZE)
  
  # results fot all species and thresholds concerning a domain set
  results <- tibble()
  
  for (S in SPECIES$genome_assembly){
    
    for (THR in THRESHOLDS){
      
      # tibble to store results for this loop 
      results_loop <- tibble(bin = seq(1, NBINS))
      
      # load breakpoint and TAD bed files
      breakpoints <- readBPFile(S, THR)
      
      # handle case of no breakpoints for threshold
      if (length(breakpoints) < 1){
        results_loop <- cbind.data.frame(results_loop, hits = rep(0, NBINS))
        results_loop <- cbind.data.frame(results_loop, rdm_hits = rep(0, NBINS))
        results_loop <- cbind.data.frame(results_loop, species = rep(factor(S), NBINS))
        results_loop <- cbind.data.frame(results_loop, threshold = rep(factor(THR), NBINS))
        results_loop <- cbind.data.frame(results_loop, boundaries = rep(factor(domain_type), NBINS))
        results_loop <- cbind.data.frame(results_loop, n_breakpoints = rep(0, NBINS))
        results <- rbind.data.frame(results, results_loop)
        next
      }
      
      # determine distribution of breakpoints
      # subdivide the genomic regions into bins of equal size and find number of breakpoints that fall into each bin
      hits <- calcHitsPerBin(breakpoints, boundaries_plus, NBINS)
      results_loop <- cbind.data.frame(results_loop, hits = hits)
      
      # generate random breakpoints
      # count number of breakpoints for each chr in 'breakpoints'
      breakpoints_per_chr <- as.tibble(breakpoints) %>%
        group_by(seqnames) %>%
        summarise(count = n())
      # generate random breakpoints for each chromosome
      rdm_breakpoints <- sampleBreakpoints(breakpoints_per_chr, hum_seqinfo, NCONTROLS)
      
      # determine distribution of random breakpoints
      rdm_hits <- calcHitsPerBin(rdm_breakpoints, boundaries_plus, NBINS)
      results_loop <- cbind.data.frame(results_loop, rdm_hits = rdm_hits)
      
      # complete result_loop info columns
      results_loop <- cbind.data.frame(results_loop, species = rep(factor(S), NBINS))
      results_loop <- cbind.data.frame(results_loop, threshold = rep(factor(THR), NBINS))
      results_loop <- cbind.data.frame(results_loop, boundaries = rep(factor(domain_type), NBINS))
      results_loop <- cbind.data.frame(results_loop, n_breakpoints = rep(length(breakpoints), NBINS))
      results <- rbind.data.frame(results, results_loop)
    }
  }
  
  all_results <- rbind.data.frame(all_results, results)
  
}

dir.create("results/", showWarnings = FALSE)
all_results <- as.tibble(all_results)
saveRDS(all_results, file = "results/breakpoints_at_boundaries.rds")

# -------------------------------------------------------------------------------------------------------------------
# Analysis of breakpoint distances to domain boundaries
# -------------------------------------------------------------------------------------------------------------------

# Calculates the distances of each breakpoint set (each species and threshold) to their
# next boundaries in domains. The same is done for randomly generated breakpoints. 

# tibble to gather all results for domain, species, threshold
results <- tibble()

for (D in DOMAINS$genomic_domain_path){
  
  domains <- import(unlist(D), seqinfo = hum_seqinfo)
  
  # get domain type to store with ever result
  domain_type <- unlist(DOMAINS %>%
                          filter(genomic_domain_path == D) %>%
                          select(genomic_domain_type)
  )
  
  print(domain_type)
  
  # extract boundaries
  boundaries <- c(GRanges(seqnames(domains), 
                          IRanges(start(domains), start(domains)), 
                          seqinfo = hum_seqinfo),
                  GRanges(seqnames(domains), 
                          IRanges(end(domains), end(domains)), 
                          seqinfo = hum_seqinfo)
                  )  
  
  
  for (S in SPECIES$genome_assembly){
    
    print(S)
    
    for (THR in THRESHOLDS){
      
      # load breakpoints, if empty -> skip
      breakpoints <- readBPFile(S, THR)
      if (length(breakpoints) < 1){
        tmp_1 <- tibble(
          boundaries = factor(domain_type),
          species = factor(S), 
          threshold= factor(THR),
          distance = NA,
          type = factor("breakpoint")
        )
        tmp_2 <- tibble(
          boundaries = factor(domain_type),
          species = factor(S), 
          threshold= factor(THR),
          distance = NA,
          type = factor("random")
        )
        tmp <- rbind.data.frame(tmp_1, tmp_2)
        results <- rbind.data.frame(results, tmp)
        next
      }
      
      # Returns the distance for each range in x to its nearest neighbor in the subject.
      bp_dist <- distanceToNearest(breakpoints, boundaries)
      
      # do the same for random breakpoints
      # count number of breakpoints for each chr in 'breakpoints'
      breakpoints_per_chr <- as.tibble(breakpoints) %>%
        group_by(seqnames) %>%
        summarise(count = n())
      # generate random breakpoints for each chromosome
      rdm_breakpoints <- sampleBreakpoints(breakpoints_per_chr, hum_seqinfo, NCONTROLS)
      # calculate distances    
      rdm_dist <- distanceToNearest(rdm_breakpoints, boundaries)
      
      # store distances in result tibble
      tmp_1 <- tibble(
        boundaries = factor(domain_type),
        species = factor(S), 
        threshold= factor(THR),
        distance = mcols(bp_dist)$distance,
        type = factor("breakpoint")
      )
      tmp_2 <- tibble(
        boundaries = factor(domain_type),
        species = factor(S), 
        threshold= factor(THR),
        distance = mcols(rdm_dist)$distance,
        type = factor("random")
      )
      
      tmp <- rbind.data.frame(tmp_1, tmp_2)
      results <- rbind.data.frame(results, tmp)
      
    }
  }
}

dir.create("results/", showWarnings = FALSE)
results <- as.tibble(results)
saveRDS(results, file = "results/distances_to_boundaries.rds")
