# =============================================================================================================================
# Comparing expression levels of orthologous genes with regard to conserved TADs and TAD rearrangements
# =============================================================================================================================

require(tidyverse)
require(biomaRt)
require(stringr)
require(BSgenome.Hsapiens.UCSC.hg19)

source("R/functions.R")

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

# minimal distance of a breakpoint to (both) domain boundaries 
# for the conserved and rearranged domain classifications
MIN_BP_BOUNDARY_DIST <- 4*10^4

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
    
    chains <- readChainFile(S) 
    
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

saveRDS(domain_classes, "results/domain_classification.rds")


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
human_expr <- rename(human_exp_match,  
                     "ensembl_gene_id" = "Gene ID",
                     "human_gene_name" = "Gene Name")
human_expr <- left_join(orthologs, human_expr, by = "ensembl_gene_id")

# Order genes in mouse_exp_match according to ordering in orthologs
mouse_expr <- rename(mouse_exp_match, 
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

# --------------------------------------------------------------------------------------------
# Start building final result tibble
# --------------------------------------------------------------------------------------------  

results <- tibble(human_gene_id = human_expr$ensembl_gene_id,
                  human_gene_name = human_expr$human_gene_name,
                  mouse_gene_id = mouse_expr$mmusculus_homolog_ensembl_gene,
                  mouse_gene_name = mouse_expr$mouse_gene_name,
                  correlation = correlations)

# =============================================================================================================================
# Assign the ortholog pairs (with correlations) to their associated domains
# ============================================================================================================================= 

# --------------------------------------------------------------------------------------------
# Receive the transcription start sites for each human ortholog
# -------------------------------------------------------------------------------------------- 

# Get human transcription start sites, use gene ids in human expr as filter
hum_tss_df <- as.tibble(getBM(attributes=c(
  'ensembl_gene_id', 'chromosome_name', 'strand',
  'transcription_start_site', 'transcript_start', 'transcript_end'), 
  filters=c('ensembl_gene_id', 'biotype'), 
  values=list(human_expr$ensembl_gene_id, biotype="protein_coding"), mart = ensembl)
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

# Order as in expression tables
hum_tss_df <- hum_tss_df[match(human_expr$ensembl_gene_id, hum_tss_df$ensembl_gene_id), ]

# Generate GRanges
hum_tss_gr <- GRanges(seqnames = hum_tss_df$chr,
                   strand = ifelse(hum_tss_df$strand == 1, "+", "-"),
                   ranges = IRanges(start = hum_tss_df$transcription_start_site,
                                    width = 1),
                   geneID = hum_tss_df$ensembl_gene_id,
                   seqinfo = hum_seqinfo
                   )

# --------------------------------------------------------------------------------------------
# Assign the domains to ortholog pair if overlapping with human TSS
# -------------------------------------------------------------------------------------------- 

tmp <- results
results <- tibble()

for (D in DOMAINS$genomic_domain_path) {
  
  # tp add columns to results, then rbind loop results
  results_loop <- tmp
  
  domains <- import(unlist(D), seqinfo = hum_seqinfo)
  
  # get domain type to store with every result
  domain_type <- unlist(DOMAINS %>%
                          filter(genomic_domain_path == D) %>%
                          select(genomic_domain_type)
                        )
  
  tss_domain_hits <- findOverlaps(hum_tss_gr, domains)
  
  # add domain information to results df
  results_loop$domain_index <- NA
  results_loop[queryHits(tss_domain_hits), ]$domain_index <- subjectHits(tss_domain_hits)
  results_loop$domain_type <- domain_type
  
  results <- rbind.data.frame(results, results_loop)
}
  
# NOTE: jetzt  nur noch mit den domain indizes in results in domain classes die TAD category bestimmen
# get rearranged TADs
domain_classes <- readRDS("results/domain_classification.rds")





# get the tad idx for genes in rearranged and intact TADs
domains$intact <- consv
domains$rearr <- rearr
intact_tad_idx <- findOverlaps(hum_tss, domains[domains$intact], select = "first")
rearr_tad_idx <- findOverlaps(hum_tss, domains[domains$rearr], select = "first")
intact_tad_idx[is.na(intact_tad_idx)] <- rearr_tad_idx[is.na(intact_tad_idx)]
hum_tss$tad <- intact_tad_idx

# # Calculate distances of TSS to nearest breakpoint
# distances <- distanceToNearest(hum_tss, breakpoints)
# hum_tss$distGeneToBreakpoint <- -1
# hum_tss[queryHits(distances)]$distGeneToBreakpoint <- mcols(distances)$distance
# if (sum(hum_tss$distGeneToBreakpoint == -1) > 0) hum_tss[hum_tss$distGeneToBreakpoint == -1]$distGeneToBreakpoint <- NA

# ------------------------------------------------------------- Gather everything in dataframe
# Merge the "correlations" and "hum_tss" dataframes.

# Extract columns in the same row order as in results (not matched = NA)
class <- hum_tss[match(results$geneID, hum_tss$geneID), ]$class
tad_idx <- hum_tss[match(results$geneID, hum_tss$geneID), ]$tad
in_grb <- hum_tss[match(results$geneID, hum_tss$geneID), ]$grb
# dist_gene_to_breakpoint <- hum_tss[match(correlations$geneID, hum_tss$geneID), ]$distGeneToBreakpoint

# Gather everything in correlations df
results <- cbind.data.frame(results, class=class, tad=tad_idx, grb = in_grb)

# filter out na class
results <- filter(results, !is.na(class))

# Add species label to data, factorise, order
results$class <- as.factor(results$class)
results$class <- factor(results$class, levels = c("Conserved", "Rearranged", "Outside"))

results <- as.tibble(results)

# Save results
saveRDS(results, str_c("results/dataframes/expressions_part1_fantom_data_dixon_domains_fills_and_bp_double_condition_10k_1000k"))
# }
