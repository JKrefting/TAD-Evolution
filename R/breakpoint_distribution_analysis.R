
# =============================================================================================================================
# Analysis of rearrangement breakpoint distributions with regard to stuctural domains of the human genome. 
# =============================================================================================================================

# require(BSgenome.Hsapiens.UCSC.hg19)
require(TxDb.Hsapiens.UCSC.hg38.knownGene)
require(rtracklayer)
source("R/functions.R")
require(tidyverse)

# Read metadata for analysis
SPECIES <- read_tsv("species_meta.tsv", col_types = "ccc")
DOMAINS <- read_tsv("domains_meta.tsv", col_types = "ccc")
METADATA <- read_tsv("metadata.tsv", col_types = "iiiiii")

# minimum size threasholds for chains (syntenic blocks) to be considered rearrangement blocks
THRESHOLDS <- METADATA %>% pull(min_size_threshold)

# number of bins where breakpoint distribution is investigated
NBINS <- METADATA %>% 
  filter(!is.na(n_bins)) %>%
  pull(n_bins)

# number of random control breakpoints per actual breakpoint
NCONTROLS <- METADATA %>%
  filter(!is.na(n_controls)) %>%
  pull(n_controls)

# read domain to GRB overlap assignent
domain_to_GRB <- read_rds("results/domain_to_GRB.rds")

hum_seqinfo <- seqinfo(TxDb.Hsapiens.UCSC.hg38.knownGene)

# ------------------------------------------------------------------------------
# Analysis of breakpoints around structural domains
# ------------------------------------------------------------------------------

# read domains
all_domains <- DOMAINS %>% 
  mutate(
    gr = map(domain_path, import.bed, seqinfo = hum_seqinfo),
    domain_subtype = "all"
  ) %>% 
  select(-domain_path, -domain_res)

# get GRB domains
grb_domains <- all_domains %>% 
  filter(domain_type != "GRB") %>% 
  mutate(
    sub_idx = map(domain_type, ~pull(filter(domain_to_GRB, domain_type == .x, GRB_class == "GRB"), domain_id)),
    gr = map2(gr, sub_idx, ~ .x[.y]),
    domain_subtype = "GRB"
  ) %>%
  select(domain_type, domain_subtype, gr)

# get non-GRB domains
nongrb_domains <- all_domains %>% 
  filter(domain_type != "GRB") %>% 
  mutate(
    sub_idx = map(domain_type, ~pull(filter(domain_to_GRB, domain_type == .x, GRB_class == "nonGRB"), domain_id)),
    gr = map2(gr, sub_idx, ~ .x[.y]),
    domain_subtype = "nonGRB"
  ) %>%
  select(domain_type, domain_subtype, gr)

# get screend domains
screened_domains <- all_domains %>% 
  filter(domain_type != "GRB") %>% 
  mutate(
    sub_idx = map(domain_type, ~pull(filter(domain_to_GRB, domain_type == .x, GRB_class == "screened"), domain_id)),
    gr = map2(gr, sub_idx, ~ .x[.y]),
    domain_subtype = "screened"
  ) %>%
  select(domain_type, domain_subtype, gr)

domains_combined = bind_rows(
  all_domains,
  grb_domains,
  nongrb_domains,
  screened_domains
  )

# parse domains and tile in bins
domain_gr_df <- domains_combined %>% 
  mutate(
    domains_plus = map(gr, ~ resize(.x, fix = "center", width = 2 * width(.x))),
    domain_plus_bin = map(domains_plus, tile, NBINS),
  ) %>% 
  select(-gr, -domains_plus)

# for each combination of threshold and species read breakpoints
bp_df <- crossing(
  threshold = THRESHOLDS, 
  species = SPECIES$genome_assembly
  ) %>% 
  mutate(
    bp_file = str_c("data/breakpoints/hg38.", species, ".", threshold, ".bp.flt.flt_adj_fill.bed"),
    breakpoints = map(bp_file, import.bed, seqinfo = hum_seqinfo),
    n_breakpoints = map_int(breakpoints, length),
    breakpoints_per_chr = map(breakpoints, ~rename(count(tibble(seqnames = as.character(seqnames(.x))), seqnames), count = n))
  )  
  
# for each species and threshold simulate brekapoint NCONTROLS time
rand_bp_df <- bp_df %>% 
  crossing(replicate = seq(NCONTROLS)) %>% 
  mutate(
    breakpoints = map(breakpoints_per_chr, sampleBreakpoints, hum_seqinfo),
    sample = "random"
  ) %>% 
  select(-bp_file, -n_breakpoints, -breakpoints_per_chr)
  
# combine real and randomized breakpoints  
all_bp_df <- rand_bp_df %>% 
  bind_rows(
    bp_df %>% 
      select(threshold, species, breakpoints) %>% 
      mutate(replicate = 1, sample = "real") 
  )

write_rds(all_bp_df, "results/all_bp_df.rds")
# all_bp_df <- read_rds("results/all_bp_df.rds")

# add domains and calculate breakpoint hits per bin arround domains
all_df <- crossing(domain_gr_df, all_bp_df) %>% 
  mutate(
    bin = list(seq(NBINS)),
    hits = map2(breakpoints, domain_plus_bin, calcHitsPerBin, NBINS)
  )


reakpoints_at_domains <- all_df %>% 
  select(-domain_plus_bin, -breakpoints) %>% 
  unnest(bin, hits)

write_rds(reakpoints_at_domains, "results/breakpoints_at_domains.rds")

# -------------------------------------------------------------------------------------------------------------------
# Analysis of breakpoints at BOUNDARIES of structural domains
# -------------------------------------------------------------------------------------------------------------------

# area around a boundary that is investigated for breakpoint enrichment
BOUNDARY_AREA <- METADATA %>%
  filter(!is.na(boundary_plus_adjacence)) %>%
   pull(boundary_plus_adjacence)

#' function to get area arounc boundares for each TAD
get_boundary_bins <- function(domains, boundary_area, nbins){
  
  # extract boundaries
  boundaries <- c(GRanges(seqnames(domains), 
                          IRanges(start(domains), start(domains)), 
                          seqinfo = hum_seqinfo),
                  GRanges(seqnames(domains), 
                          IRanges(end(domains), end(domains)), 
                          seqinfo = hum_seqinfo)
  )  
  
  # enlarge
  boundaries_plus <- resize(boundaries, fix = "center", boundary_area)
  
  # returns GRangesList in which every domain is subdivided into NBINS
  boundaries_plus_bins <- tile(boundaries_plus, nbins)
  
  return(boundaries_plus_bins)  
}

# get all rearrangement hits around boundaries
boundary_hits <- all_domains %>% 
  mutate(boundary_bins = map(gr, get_boundary_bins, BOUNDARY_AREA, NBINS)) %>% 
  crossing(all_bp_df) %>% 
  mutate(
    bin = list(seq(NBINS)),
    hits = map2(breakpoints, boundary_bins, calcHitsPerBin, NBINS)
  )

# remove larege objects and unpack bins and hits tables into one large table
boundary_hits_tidy <- boundary_hits %>% 
  select(-gr, -boundary_bins, -breakpoints) %>% 
  unnest(bin, hits)

write_rds(boundary_hits_tidy, "results/boundary_hits_tidy.rds")

# domain_result_list <- map(DOMAINS$domain_path, function(D){
#   
#   domains <- import(unlist(D), seqinfo = hum_seqinfo)
#   
#   # get domain type to store with ever result
#   domain_type <- unlist(DOMAINS %>%
#                           filter(domain_path == D) %>%
#                           dplyr::select(domain_type)
#   )
#   
#   # extract boundaries
#   boundaries <- c(GRanges(seqnames(domains), 
#                           IRanges(start(domains), start(domains)), 
#                           seqinfo = hum_seqinfo),
#                   GRanges(seqnames(domains), 
#                           IRanges(end(domains), end(domains)), 
#                           seqinfo = hum_seqinfo)
#   )  
#   
#   # enlarge
#   boundaries_plus <- resize(boundaries, fix = "center", BOUNDARY_AREA)
#   
#   # returns GRangesList in which every domain is subdivided into NBINS
#   boundaries_plus_bins <- tile(boundaries_plus, NBINS)
#   
#   species_result_list <- map(SPECIES$genome_assembly, function(S){
#     
#     thr_result_list <- map(THRESHOLDS, function(THR){
#       
#       # load breakpoint and TAD bed files
#       breakpoints <- readBPFile(S, THR)
#       
#       # handle case of no breakpoints for threshold
#       if (length(breakpoints) < 1){
#         # tibble to store results for this loop 
#         # add to results
#         actual_results <- tibble(
#           bin = seq(NBINS),
#           hits = 0,
#           sample = "real",
#           replicate = 1,
#           species = S,
#           threshold = THR,
#           domains = domain_type,
#           n_breakpoints = length(breakpoints)
#         )
#         
#         return(actual_results)
#       }
#       
#       # -------------------------------------------------------------------------------------------------------------------
#       # Distribution of actual breakpoints
#       # -------------------------------------------------------------------------------------------------------------------
#       
#       # find number of breakpoints that fall into each bin
#       hits <- calcHitsPerBin(breakpoints, boundaries_plus_bins, NBINS)
#       
#       # add to results
#       actual_results <- tibble(
#         bin = seq(NBINS),
#         hits = hits,
#         sample = "real",
#         replicate = 1,
#         species = S,
#         threshold = THR,
#         domains = domain_type,
#         n_breakpoints = length(breakpoints)
#       )
#       
#       # -------------------------------------------------------------------------------------------------------------------
#       # Distribution of random breakpoints
#       # -------------------------------------------------------------------------------------------------------------------
#       
#       # count number of breakpoints for each chr in 'breakpoints'
#       breakpoints_per_chr <- tibble(seqnames = as.character(seqnames(breakpoints))) %>%
#         mutate(seqnames = as.character(seqnames)) %>%
#         group_by(seqnames) %>%
#         summarise(count = n())
#       
#       control_results_list <- map(seq(NCONTROLS), function(i){
#         # generate the corresponding number of random breakpoints for each chromosome
#         rdm_breakpoints <- sampleBreakpoints(breakpoints_per_chr, hum_seqinfo)
#         
#         # determine distribution of random breakpoints
#         rdm_hits <- calcHitsPerBin(rdm_breakpoints, boundaries_plus_bins, NBINS)
#         
#         # add to results
#         results_loop <- tibble(
#           bin = seq(NBINS),
#           hits = rdm_hits,
#           sample = "random",
#           replicate = i,
#           species = S,
#           threshold = THR,
#           domains = domain_type,
#           n_breakpoints = length(breakpoints)
#         )
#         return(results_loop)
#       })
#       # return combined tibble for single threshold
#       return(bind_rows(actual_results, control_results_list))
#     }) 
#     # return combined tibble for all thresholds / single species
#     return(bind_rows(thr_result_list))
#   })
#   # return combined tibble for all species / single domain
#   return(bind_rows(species_result_list))
# })
# 
# results <- bind_rows(domain_result_list)
# 
# dir.create("results/", showWarnings = FALSE)
# write_tsv(results, "results/breakpoints_at_boundaries.tsv")
# write_rds(results, "results/breakpoints_at_boundaries.rds")
# results <- read_rds("results/breakpoints_at_boundaries.rds")

# -------------------------------------------------------------------------------------------------------------------
# Analysis of breakpoint distances to domain boundaries
# -------------------------------------------------------------------------------------------------------------------

# Calculates the distances of each breakpoint set (each species and threshold) to their
# next boundaries in domains. The same is done for randomly generated breakpoints. 

# tibble to gather all results for domain, species, threshold
results <- tibble()

for (D in DOMAINS$domain_path){
  
  domains <- import(unlist(D), seqinfo = hum_seqinfo)
  
  # get domain type to store with ever result
  domain_type <- unlist(DOMAINS %>%
                          filter(domain_path == D) %>%
                          dplyr::select(domain_type)
  )
  
  # extract boundaries
  boundaries <- c(GRanges(seqnames(domains), 
                          IRanges(start(domains), start(domains)), 
                          seqinfo = hum_seqinfo),
                  GRanges(seqnames(domains), 
                          IRanges(end(domains), end(domains)), 
                          seqinfo = hum_seqinfo)
                  )  
  
  
  for (S in SPECIES$genome_assembly){
    
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
      rdm_breakpoints <- sampleBreakpoints(breakpoints_per_chr, hum_seqinfo)
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
write_tsv(results, "results/distances_to_boundaries.tsv")
write_rds(results, "results/distances_to_boundaries.rds")
