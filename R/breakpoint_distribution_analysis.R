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
  
  
  boundaries <- c(GRanges(seqnames(domains), 
                          IRanges(start(domains) - boundary_dist, start(domains) + boundary_dist), 
                          seqinfo = hum_seqinfo),
                  GRanges(seqnames(domains), 
                          IRanges(end(domains) - boundary_dist, end(domains) + boundary_dist), 
                          seqinfo = hum_seqinfo)
                  )  
  
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
    # calculate hits
    hits_at_boundaries <- calcHitsAtBoundaries(breakpoints, boundaries, n_bins, 'list')
    # add as col to plot_data_loop
    plot_data_loop <- cbind.data.frame(plot_data_loop, hits = hits_at_boundaries)
    
    rdm_breakpoints_n_controls <- generateRandomBreakpointSets(breakpoints, hum_seqinfo, n_controls)
    # generate hits per bin (domain bins) for each random breakpoint set
    rdm_hits_at_boundaries <- calcHitsAtBoundaries(rdm_breakpoints_n_controls, boundaries, n_bins, 'list')
    # add each rdm_hits_per_bin as column to plot_data_loop
    plot_data_loop <- cbind.data.frame(plot_data_loop, rdmHits = rdm_hits_at_boundaries)
    
    # complete plot_data_loop and add to plot_data
    plot_data_loop <- cbind.data.frame(plot_data_loop, vsSpecies = rep(species, n_bins))
    plot_data_loop <- cbind.data.frame(plot_data_loop, threshold = rep(threshold, n_bins))
    plot_data_loop <- cbind.data.frame(plot_data_loop, nBreakpoints = rep(length(breakpoints), n_bins))
    plot_data <- rbind.data.frame(plot_data, plot_data_loop)
  } # end for threshold
} # end for species

saveRDS(plot_data, file = "results/dataframes/breakpoints_boundaries_rao_large.rds")
