require(tidyverse)
require(biomaRt)
require(stringr)
require(BSgenome.Hsapiens.UCSC.hg19)

source("R/functions.R")

# Read metadata for analysis
SPECIES <- read_tsv("species_meta.tsv")
DOMAINS <- read_tsv("domains_meta.tsv")
METADATA <- read_tsv("metadata.tsv")
THRESHOLDS <- unlist(METADATA %>% 
                       filter(!is.na(min_size_threshold)) %>%
                       dplyr::select(min_size_threshold)
)

# Load human seqinfo
genome <- BSgenome.Hsapiens.UCSC.hg19
hum_seqinfo <- seqinfo(genome)

# =============================================================================================================================
# Prerequisite for domain categorisation in rearranged and conserved domains. 
# For every domain determine if #1 it is enclosed by a syntenic block (chain) #2 a rearrangement breakpoint 
# occurs inside the domain boundaries (for each threshold)
# =============================================================================================================================

# tibble to store all classification info
domain_classes <- tibble()

for (D in DOMAINS$domain_path) {
  
  domains <- import(unlist(D), seqinfo = hum_seqinfo)
  
  # get domain type to store with every result
  domain_type <- unlist(DOMAINS %>%
                          filter(domain_path == D) %>%
                          dplyr::select(domain_type)
  )
  
  print(domain_type)
  
  # get domain resolution for domain shrinking
  domain_res <- unlist(DOMAINS %>%
                          filter(domain_path == D) %>%
                          dplyr::select(domain_res)
  )
  
  # shrink domains to accommodate data resolution
  domains_minus <- resize(domains, fix = "center", 
                          ifelse(width(domains) - (4 * domain_res) <= 0, 
                                  1, 
                                 width(domains) - (4 * domain_res)
                                 )
                          )
  
  for (S in SPECIES$genome_assembly) {
    
    print(S)
    
    chains <- readFillFile(S) 
    
    # check if domains inside chains
    enclosed_by_chain <- overlapsAny(domains, chains, type = "within")
    
    # tibble to gather rearranged results for all thresholds
    rearr_all_thr <- tibble()
    
    for (THR in THRESHOLDS){
      # for each syntenic block size threshold, check if a rearrangement breakpoint 
      # occurrs inside a domain with a safety margin of at least 40 kb to each TAD boundary
      
      breakpoints <- readBPFile(S, THR)
      
      # handle case of no breakpoints for threshold
      if (length(breakpoints) < 1){
        rearr_this_thr <- tibble(rearranged_by_breakpoint = rep(NA, length(domains)),
                                 threshold = rep(THR, length(domains))
        )
        rearr_all_thr <- rbind.data.frame(rearr_all_thr, rearr_this_thr)
        next
      }
      
      # find rearranged domains
      rearr_by_bp <- overlapsAny(domains_minus, breakpoints)
      
      # tibble to gather results for each threshold
      rearr_this_thr <- tibble(rearranged_by_breakpoint = rearr_by_bp,
                               threshold = THR)
      rearr_all_thr <- rbind.data.frame(rearr_all_thr, rearr_this_thr)
    }
    
    # tibble for domain classification by species
    domain_classes_by_species  <- tibble(domain_id = rep(1:length(domains), length(THRESHOLDS)),
                                         domain_type = domain_type,
                                         enclosed_by_chain = rep(enclosed_by_chain, length(THRESHOLDS)),
                                         rearranged_by_breakpoint = rearr_all_thr$rearranged_by_breakpoint,
                                         threshold = rearr_all_thr$threshold,
                                         species = S
    )
    
    domain_classes <- rbind.data.frame(domain_classes, domain_classes_by_species)
    
  }
  
}

saveRDS(domain_classes, "results/domain_classification_more_savely.rds")
