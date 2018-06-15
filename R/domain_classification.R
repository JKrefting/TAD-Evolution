# require(BSgenome.Hsapiens.UCSC.hg19)
require(TxDb.Hsapiens.UCSC.hg38.knownGene)
require(biomaRt)
require(tidyverse)
require(stringr)

source("R/functions.R")

# Read metadata for analysis
SPECIES <- read_tsv("species_meta.tsv")
DOMAINS <- read_tsv("domains_meta.tsv")
METADATA <- read_tsv("metadata.tsv")
THRESHOLDS <- unlist(METADATA %>% 
                       filter(!is.na(min_size_threshold)) %>%
                       dplyr::select(min_size_threshold)
)

#Conserved domains do NOT have breakpoints from all thresholds (smallest size threshold)
CONSV_BP_THR <- THRESHOLDS[1]
# Rearranged domains have breakpoints from the largest threshold
REARR_BP_THR <- THRESHOLDS[3]


# Load human seqinfo
# genome <- BSgenome.Hsapiens.UCSC.hg19
hum_seqinfo <- seqinfo(TxDb.Hsapiens.UCSC.hg38.knownGene)
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
    
    chains <- readFillFile(S, seqinfo = hum_seqinfo) 
    
    # check if domains inside chains
    enclosed_by_chain <- overlapsAny(domains, chains, type = "within")
    
    # tibble to gather rearranged results for all thresholds
    rearr_all_thr <- tibble()
    
    for (THR in THRESHOLDS){
      # for each syntenic block size threshold, check if a rearrangement breakpoint 
      # occurrs inside a domain with a safety margin of at least 40 kb to each TAD boundary
      print(THR)
      
      breakpoints <- readBPFile(S, THR)

      # find rearranged domains
      rearr_by_bp <- overlapsAny(domains_minus, breakpoints)
      
      # tibble to gather results for each threshold
      rearr_this_thr <- tibble(rearranged_by_breakpoint = rearr_by_bp,
                               threshold = THR)
      rearr_all_thr <- bind_rows(rearr_all_thr, rearr_this_thr)
    }
    
    # tibble for domain classification by species
    domain_classes_by_species  <- tibble(domain_id = rep(1:length(domains), length(THRESHOLDS)),
                                         domain_type = domain_type,
                                         enclosed_by_chain = rep(enclosed_by_chain, length(THRESHOLDS)),
                                         rearranged_by_breakpoint = rearr_all_thr$rearranged_by_breakpoint,
                                         threshold = rearr_all_thr$threshold,
                                         species = S
    )
    
    domain_classes <- bind_rows(domain_classes, domain_classes_by_species)
    
  }
  
}

write_tsv(domain_classes, "results/domain_classification.tsv")
write_rds(domain_classes, "results/domain_classification.rds")

#===============================================================================
# calssify each TAD for each species in either "conserved" or "rearraged" 
#===============================================================================

# get subset with conserved threshold
domains_conserved_th <- domain_classes %>% 
  filter(threshold == CONSV_BP_THR) %>% 
  select(domain_id, domain_type, species, rearranged_by_breakpoint)

# get subset with rearranged threshold
domains_rearranged_th <- domain_classes %>% 
  filter(threshold == REARR_BP_THR) %>% 
  select(domain_id, domain_type, species, rearranged_by_breakpoint)


# remove threshold column and get unique domain per type and species
tidy_domain_groups <- domain_classes %>% 
  select(-threshold) %>% 
  distinct(domain_id, domain_type, species, .keep_all = TRUE)

# add rearraged column from specific size thresholds
tidy_domain_groups <- tidy_domain_groups %>% 
  left_join(domains_conserved_th, suffix = c("", "_conserved_th"),
            by = c("domain_id", "domain_type", "species")) %>% 
  left_join(domains_rearranged_th, suffix = c("", "_rearranged_th"),
            by = c("domain_id", "domain_type", "species")) 


# add conserved and rearranged categories as columns
tidy_domain_groups <- tidy_domain_groups %>% 
  mutate(
    conserved = enclosed_by_chain & !rearranged_by_breakpoint_conserved_th,
    rearranged = !enclosed_by_chain & rearranged_by_breakpoint_rearranged_th,
    category = ifelse(conserved, 
                      "Conserved", 
                      ifelse(rearranged, "Rearranged", NA)
                      )
  )

# save to output file
write_rds(tidy_domain_groups, "results/tidy_domain_groups.rds")


