# =============================================================================================================================
# Comparing expression levels of orthologous genes with regard to conserved TADs and TAD rearrangements
# =============================================================================================================================

require(tidyverse)
require(biomaRt)
require(stringr)
require(BSgenome.Hsapiens.UCSC.hg19)

source("R/functions.R")

# Read metadata for analysis
METADATA <- read_tsv("metadata.csv")
SPECIES <- dplyr::select(METADATA, genome_assembly, trivial_name)
DOMAINS <- METADATA %>% 
  dplyr::select(genomic_domain_type, genomic_domain_path) %>%
  filter(!is.na(genomic_domain_type))
THRESHOLDS <- unlist(METADATA %>% 
                       filter(!is.na(min_size_threshold)) %>%
                       dplyr::select(min_size_threshold)
)

# Minimal distance of a breakpoint to (both) domain boundaries 
# for the conserved and rearranged domain classifications
MIN_BP_BOUNDARY_DIST <- 4*10^4

# Conserved domains do NOT have breakpoints from this size threshold
CONSV_BP_THR <- 10^4
# Rearranged domains have breakpoints from this size threshold
REARR_BP_THR <- 10^6

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

for (D in DOMAINS$genomic_domain_path) {
  
  domains <- import(unlist(D), seqinfo = hum_seqinfo)
  
  # get domain type to store with every result
  domain_type <- unlist(DOMAINS %>%
                          filter(genomic_domain_path == D) %>%
                          select(genomic_domain_type)
  )
  
  print(domain_type)
  
  for (S in SPECIES$genome_assembly) {
    
    print(S)
    
    chains <- readFillFile(S) 
    
    # check if domains inside chains
    enclosed_by_chain <- overlapsAny(domains, chains, type = "within")
    
    # tibble to gather rearranged results for all thresholds
    rearr_all_thr <- tibble()
    
    for (THR in THRESHOLDS){
      
      print(THR)
      
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
      
      # Calculate the distance to closest breakpoint for each domain
      dist_domains_to_bp <- distanceToNearest(domains, breakpoints)
      
      # Preselect 'rearranged' domains
      dist_rearr_domains <- dist_domains_to_bp[mcols(dist_domains_to_bp)$distance == 0] # bp inside domain, margin unknown
      # Select conserved domains
      dist_conserved_domains <- dist_domains_to_bp[mcols(dist_domains_to_bp)$distance > MIN_BP_BOUNDARY_DIST] # bp outside, margin verified
      
      # handle case of no domain rearrangement
      if (length(dist_rearr_domains) < 1){
        rearr_this_thr <- tibble(rearranged_by_breakpoint = rep(FALSE, length(domains)),
                                 threshold = rep(THR, length(domains))
        )
        rearr_all_thr <- rbind.data.frame(rearr_all_thr, rearr_this_thr)
        next
      }
      
      # Verify rearranged TADs
      savely_rearranged <- mapply(function(dm_idx, bp_idx){
        # Select rearranged domain
        dm <- domains[dm_idx, ]
        # Extract domain boundaries
        boundaries <- getBoundaries(dm, 0)
        # Select breakpoint that overlaps domain
        bp <- breakpoints[bp_idx, ]
        # Calculate distances between both boundaries and the breakpoint
        dist <- abs(distance(boundaries, bp))
        # If at least one of the distances is smaller than 'min_distance', return FALSE 
        if (sum(dist < MIN_BP_BOUNDARY_DIST) >= 1) {
          # TRUE = Breakpoint far enough from boundaries
          # FALSE = Breakpoint too close to boundary
          return(FALSE)
        }
        return(TRUE)
      }, queryHits(dist_rearr_domains), subjectHits(dist_rearr_domains))
      
      # Filter out savely rearranged domains from <<Hits object>>
      dist_rearr_domains <- dist_rearr_domains[savely_rearranged]
      
      # Write labels into vector
      rearr_by_bp <- vector(length=length(domains))
      rearr_by_bp[queryHits(dist_rearr_domains)] <- TRUE
      rearr_by_bp[queryHits(dist_conserved_domains)] <- FALSE
      
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

saveRDS(domain_classes, "results/domain_classification_fills.rds")


# =============================================================================================================================
# Gene expression correlation of orthologs
# =============================================================================================================================

# preprocess downloaded expression data
source("R/ortholog_expression_data.R")
# read expression data tables
human_exp_match <- read_tsv("data/ExpressionAtlas/formated.human_exp.matching_subset.tsv")
mouse_exp_match <- read_tsv("data/ExpressionAtlas/formated.mouse_exp.matching_subset.tsv")

# -------------------------------------------------------------------------------------------
# Find human-mouse orthologs
# ------------------------------------------------------------------------------------------- 

# Load gene ensembl human gene ids and orthologs
ensembl <- useMart(host="grch37.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", verbose=FALSE)

# Attributes that are returned
orth_attr = c("ensembl_gene_id", 
              "mmusculus_homolog_ensembl_gene", "mmusculus_homolog_orthology_type", "mmusculus_homolog_orthology_confidence")

# Query orthologs by human gene ID
orthologs = getBM(attributes=orth_attr,
                  filters=c("ensembl_gene_id", "biotype", "chromosome_name"),
                  values=list(ensemble_gene_id=human_exp_match$'Gene ID', biotype="protein_coding", chromosome_name=c(1:22, "X", "Y")),
                  mart=ensembl)

# --------------------------------------------------------------------------------------------
# Prepare data for correlation
# --------------------------------------------------------------------------------------------  

# For genes that do not have orthologs set empty string to NA
orthologs[orthologs[ , 2] == "", 2:4] <- NA

# Filter out one2one orthologs
orthologs <- orthologs %>% filter(mmusculus_homolog_orthology_type == "ortholog_one2one")

# Filter out orthologs with high confidence (http://www.ensembl.org/info/genome/compara/Ortholog_qc_manual.html)
orthologs <- orthologs %>% filter(mmusculus_homolog_orthology_confidence == 1)

# genes without ortholog should be eliminated, for safety
orthologs <- orthologs %>% filter(!is.na(mmusculus_homolog_ensembl_gene))

# Order genes in human_exp_match according to ordering in orthologs
human_expr <- dplyr::rename(human_exp_match,  
                     "ensembl_gene_id" = "Gene ID",
                     "human_gene_name" = "Gene Name")
human_expr <- left_join(orthologs, human_expr, by = "ensembl_gene_id")

# Order genes in mouse_exp_match according to ordering in orthologs
mouse_expr <- dplyr::rename(mouse_exp_match, 
                     "mmusculus_homolog_ensembl_gene" = "Gene ID",
                     "mouse_gene_name" = "Gene Name")
mouse_expr <- left_join(orthologs, mouse_expr, by = "mmusculus_homolog_ensembl_gene")

# --------------------------------------------------------------------------------------------
# Calculate correlations
# --------------------------------------------------------------------------------------------  

correlations  <- vector()
for (i in 1:nrow(human_expr)){
  # Loop through each row of the dataframes (i.e. loop through the orthologs) 
  # and convert the expression values of the tissues into a numeric vector for each gene to calculate
  # the correlations
  cor <- cor(as.numeric(human_expr[i, 6:ncol(human_expr)]), 
             as.numeric(mouse_expr[i, 6:ncol(mouse_expr)]),
             use="pairwise.complete.obs", method="pearson")
  correlations <- c(correlations, cor)
}

# make correlations identifiable by human gene id
correlations <- tibble(ensembl_gene_id = human_expr$ensembl_gene_id,
                       correlation = correlations)

# =============================================================================================================================
# Assign the ortholog pairs (with correlations) to their associated domains and categorise
# ============================================================================================================================= 

# --------------------------------------------------------------------------------------------
# Receive the transcription start sites for each human ortholog
# -------------------------------------------------------------------------------------------- 

# Get human transcription start sites
hum_tss_df <- as.tibble(getBM(attributes=c(
  'ensembl_gene_id', 'external_gene_name', 'chromosome_name', 'strand',
  'transcription_start_site', 'transcript_start', 'transcript_end'), 
  filters=c('biotype'), 
  values=list(biotype="protein_coding"), mart = ensembl)
)

# Filter out duplicated entries (take only TSS from largest transcript)
hum_tss_df$transcript_length <- hum_tss_df$transcript_end - hum_tss_df$transcript_start
hum_tss_df <- hum_tss_df[with(hum_tss_df, order(ensembl_gene_id, -transcript_length)), ]
hum_tss_df <- hum_tss_df[!duplicated(hum_tss_df$ensembl_gene_id), ]

# Delete the spare columns
hum_tss_df <- hum_tss_df[ ,1:4]

hum_tss_df <- hum_tss_df %>%
  # Adapt UCSC name style
  mutate(chr = str_c("chr", chromosome_name)) %>%
  # Filter for valid names
  filter(chr %in% seqnames(hum_seqinfo))

# Generate GRanges
hum_tss_gr <- GRanges(seqnames = hum_tss_df$chr,
                   strand = ifelse(hum_tss_df$strand == 1, "+", "-"),
                   ranges = IRanges(start = hum_tss_df$transcription_start_site,
                                    width = 1),
                   geneID = hum_tss_df$ensembl_gene_id,
                   seqinfo = hum_seqinfo
                   )
# --------------------------------------------------------------------------------------------
# Start building final result tibble
# --------------------------------------------------------------------------------------------  

# add ortholog gene info (if present) to ALL tss retrieved
inter_results <- hum_tss_df %>%
  left_join(dplyr::select(human_expr, 1:4), by = "ensembl_gene_id")

# add the correlation results to otrholog pairs
inter_results <- inter_results %>%
  left_join(correlations, by = "ensembl_gene_id")

# --------------------------------------------------------------------------------------------
# Assign the domains to ortholog pair in 'results' if overlapping with human TSS,
# then denote conserved or rerranged status of each domain
# -------------------------------------------------------------------------------------------- 


domain_classes <- readRDS("results/domain_classification_fills.rds")

tmp <- tibble()

for (D in DOMAINS$genomic_domain_path) {
  
  domains <- import(unlist(D), seqinfo = hum_seqinfo)
  
  # get domain type to store with every result
  domain_type <- unlist(DOMAINS %>%
                          filter(genomic_domain_path == D) %>%
                          dplyr::select(genomic_domain_type)
  )

  # to add columns to results, then rbind loop results
  results_loop <- inter_results
  
  # determine gene association to domain
  tss_domain_hits <- findOverlaps(hum_tss_gr, domains)
  
  # add associated domain (if present) to orthologs
  # NOTE: order of genes in hum_tss_gr and results_loop (inter_results) equal
  results_loop$domain_id <- NA
  results_loop[queryHits(tss_domain_hits), ]$domain_id <- subjectHits(tss_domain_hits)
  results_loop$domain_type <- domain_type
  
  # categorise domain as conserved, identifiable by domain id
  conserved <- domain_classes %>%
    filter(domain_type == domain_type,
           species == "mm10",
           threshold == CONSV_BP_THR
    ) %>%
    transmute(domain_id = domain_id,
              conserved = enclosed_by_chain & !rearranged_by_breakpoint)
  
  # categorise domain as rearranged
    rearranged <- domain_classes %>%
    filter(domain_type == domain_type,
           species == "mm10",
           threshold == REARR_BP_THR
    ) %>%
    transmute(domain_id = domain_id,
              rearranged = !enclosed_by_chain & rearranged_by_breakpoint)
  
  # add conserved / rearranged info of domains by merging by domain id
  results_loop$conserved <- conserved[match(results_loop$domain_id, conserved$domain_id), ]$conserved 
  results_loop$rearranged <- rearranged[match(results_loop$domain_id, rearranged$domain_id), ]$rearranged 
  
  
  # determine genes outside any domain
  results_loop$outside <- !overlapsAny(hum_tss_gr, domains)
  
  tmp <- rbind.data.frame(tmp, results_loop)
}

results <- tmp

# --------------------------------------------------------------------------------------------
# Finish categorisation of orthpairs localised in conserved, rearranged or outside domains
# -------------------------------------------------------------------------------------------- 

results$category <- NA
results$category <- ifelse(results$conserved, "Conserved", results$category)
results$category <- ifelse(results$rearranged, "Rearranged", results$category)
results$category <- ifelse(results$outside, "Outside", results$category)
results$category <- factor(results$category, levels = c("Conserved", "Rearranged", "Outside"))

saveRDS(results, "results/ortholog_expression_correlation.rds")
