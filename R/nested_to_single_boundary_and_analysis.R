require(tidyverse)
require(stringr)
require(BSgenome.Hsapiens.UCSC.hg19)

setwd("/home/jan/projects/TAD-Evolution/GM12878_single_test/TAD-Evolution")

source("R/functions.R")

hg19_seqinfo <- seqinfo(BSgenome.Hsapiens.UCSC.hg19)

# ==============================================================================
# extract GM12878 boundaries, seperate them into nested (inside at least another
# TAD) and single (not overlapping any TADs) boundaries
# ==============================================================================

tad_path <- "data/TADs/GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.bed"

tads <- import.bed(tad_path, seqinfo = hg19_seqinfo)

# extract boundaries
bdies <- c(
  GRanges(seqnames = seqnames(tads), 
                 ranges = IRanges(start = start(tads), width = 1), 
                 strand = strand(tads)),
  GRanges(seqnames = seqnames(tads), 
          ranges = IRanges(start = end(tads), width = 1), 
          strand = strand(tads))
  )

# 'assign' start and end boundaries to their corresponding TADs
start_hits <- hits <- findOverlaps(bdies, tads, type = "start")
end_hits <- hits <- findOverlaps(bdies, tads, type = "end")


# find boundaries which are within a TAD (includes the start and end boundaries 
# as hits)
all_hits <- findOverlaps(bdies, tads, type = "within")

# filter out start and end hits from all hits to identify boundaries which are 
# located in at least one other TAD
tmp_hits <- setdiff(all_hits, start_hits)
final_hits <- setdiff(tmp_hits, end_hits)

# get boundaries overlapping TADs (= nested)
bdy_idx <- unique(queryHits(final_hits))
nested_bdies <- bdies[bdy_idx]

# get boundaries NOT overlapping another TAD (= single)
bdy_idx <- setdiff(seq(1:length(bdies)), queryHits(final_hits))
single_bdies <- bdies[bdy_idx]

# save as bed files
export.bed(nested_bdies,
           "data/TADs/GSE63525_GM12878_primary+replicate_Arrowhead_domainlist_nested_boundaries.bed",
           format = "bed")

export.bed(single_bdies,
           "data/TADs/GSE63525_GM12878_primary+replicate_Arrowhead_domainlist_single_boundaries.bed",
           format = "bed")


# ==============================================================================
# Breakpoint distribution at GM12878 single and nested boundaries
# ==============================================================================

# Read metadata for analysis
SPECIES <- read_tsv("species_meta.tsv")
METADATA <- read_tsv("metadata.tsv")
DOMAINS <- read_tsv("domains_meta.tsv")

# save as bed files
nested_bdies <- import.bed("data/TADs/GSE63525_GM12878_primary+replicate_Arrowhead_domainlist_nested_boundaries.bed",
           seqinfo = hg19_seqinfo)

single_bdies <- import.bed("data/TADs/GSE63525_GM12878_primary+replicate_Arrowhead_domainlist_single_boundaries.bed",
           seqinfo = hg19_seqinfo)

# minimum size threasholds for chains (syntenic blocks) to be considered 
# rearrangement blocks
THRESHOLDS <- unlist(METADATA %>% 
                       filter(!is.na(min_size_threshold)) %>%
                       dplyr::select(min_size_threshold)
)

# number of bins where breakpoint distribution is investigated
NBINS <- unlist(METADATA %>% 
                  filter(!is.na(n_bins)) %>%
                  dplyr::select(n_bins)
)

# number of random control breakpoints per actual breakpoint
NCONTROLS <- unlist(METADATA %>% 
                      filter(!is.na(n_controls)) %>%
                      dplyr::select(n_controls)
)

# read domain to GRB overlap assignent
domain_to_GRB <- read_rds("results/domain_to_GRB.rds")

# area around a boundary that is investigated for breakpoint enrichment
BOUNDARY_AREA <-  unlist(METADATA %>% 
                           filter(!is.na(boundary_plus_adjacence)) %>%
                           dplyr::select(boundary_plus_adjacence)
)

domain_result_list <- map(DOMAINS$domain_path, function(D){
  
  domains <- import(unlist(D), seqinfo = hg19_seqinfo)
  
  # get domain type to store with every result
  domain_type <- unlist(DOMAINS %>%
                          filter(domain_path == D) %>%
                          dplyr::select(domain_type)
  )
  
  # if domain type is not already just boundaries, extract boundaries
  if (!grepl("bdies", domain_type)){
    boundaries <- c(GRanges(seqnames(domains), 
                          IRanges(start(domains), start(domains)), 
                          seqinfo = hg19_seqinfo),
                  GRanges(seqnames(domains), 
                          IRanges(end(domains), end(domains)), 
                          seqinfo = hg19_seqinfo)
                  )
  } else {
    boundaries <- domains
  }
  
  # enlarge
  boundaries_plus <- resize(boundaries, fix = "center", BOUNDARY_AREA)
  
  # returns GRangesList in which every domain is subdivided into NBINS
  boundaries_plus_bins <- tile(boundaries_plus, NBINS)
  
  species_result_list <- map(SPECIES$genome_assembly, function(S){
    
    thr_result_list <- map(THRESHOLDS, function(THR){
      
      # load breakpoint and TAD bed files
      breakpoints <- readBPFile(S, THR)
      
      # handle case of no breakpoints for threshold
      if (length(breakpoints) < 1){
        # tibble to store results for this loop 
        # add to results
        actual_results <- tibble(
          bin = seq(NBINS),
          hits = 0,
          sample = "real",
          replicate = 1,
          species = S,
          threshold = THR,
          domains = domain_type,
          n_breakpoints = length(breakpoints)
        )
        
        return(actual_results)
      }
      
      # -------------------------------------------------------------------------------------------------------------------
      # Distribution of actual breakpoints
      # -------------------------------------------------------------------------------------------------------------------
      
      # find number of breakpoints that fall into each bin
      hits <- calcHitsPerBin(breakpoints, boundaries_plus_bins, NBINS)
      
      # add to results
      actual_results <- tibble(
        bin = seq(NBINS),
        hits = hits,
        sample = "real",
        replicate = 1,
        species = S,
        threshold = THR,
        domains = domain_type,
        n_breakpoints = length(breakpoints)
      )
      
      # -------------------------------------------------------------------------------------------------------------------
      # Distribution of random breakpoints
      # -------------------------------------------------------------------------------------------------------------------
      
      # count number of breakpoints for each chr in 'breakpoints'
      breakpoints_per_chr <- tibble(seqnames = as.character(seqnames(breakpoints))) %>%
        mutate(seqnames = as.character(seqnames)) %>%
        group_by(seqnames) %>%
        summarise(count = n())
      
      control_results_list <- map(seq(NCONTROLS), function(i){
        # generate the corresponding number of random breakpoints for each chromosome
        rdm_breakpoints <- sampleBreakpoints(breakpoints_per_chr, hg19_seqinfo)
        
        # determine distribution of random breakpoints
        rdm_hits <- calcHitsPerBin(rdm_breakpoints, boundaries_plus_bins, NBINS)
        
        # add to results
        results_loop <- tibble(
          bin = seq(NBINS),
          hits = rdm_hits,
          sample = "random",
          replicate = i,
          species = S,
          threshold = THR,
          domains = domain_type,
          n_breakpoints = length(breakpoints)
        )
        return(results_loop)
      })
      # return combined tibble for single threshold
      return(bind_rows(actual_results, control_results_list))
    }) 
    # return combined tibble for all thresholds / single species
    return(bind_rows(thr_result_list))
  })
  # return combined tibble for all species / single domain
  return(bind_rows(species_result_list))
})

results <- bind_rows(domain_result_list)

dir.create("results/", showWarnings = FALSE)
write_tsv(results, "results/breakpoints_at_boundaries.tsv")
write_rds(results, "results/breakpoints_at_boundaries.rds")