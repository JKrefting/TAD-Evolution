# =============================================================================================================================
# Comparing expression levels of orthologous genes with regard to conserved TADs and TAD rearrangements
# =============================================================================================================================

require(tidyverse)
require(biomaRt)
require(stringr)
require(BSgenome.Hsapiens.UCSC.hg19)

source("functions.R")

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

# Load human seqinfo
genome <- BSgenome.Hsapiens.UCSC.hg19
hum_seqinfo <- seqinfo(genome)

# preprocess downloaded expression data
source("R/ortholog_expression_data.R")
# read expression data tables
human_exp_match <- read_tsv("data/ExpressionAtlas/formated.human_exp.matching_subset.tsv")
mouse_exp_match <- read_tsv("data/ExpressionAtlas/formated.mouse_exp.matching_subset.tsv")

# -------------------------------------------------------------------------------------------------------------------
# Find human-mouse orthologs
# ------------------------------------------------------------------------------------------------------------------- 

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

# -------------------------------------------------------------------------------------------------------------------
# Prepare data for correlation
# -------------------------------------------------------------------------------------------------------------------  

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

# -------------------------------------------------------------------------------------------------------------------
# Calculate correlations
# -------------------------------------------------------------------------------------------------------------------  

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

# -------------------------------------------------------------------------------------------------------------------
# Categorise TADs as conserved or rearranged
# -------------------------------------------------------------------------------------------------------------------

TADs_with_breakpoint <- tibble(domain_id = 1:length(domains))
TADs_enclosed_by_chain <- tibble(domain_id = 1:length(domains))

for (D in DOMAINS$genomic_domain_path) {
  
  domains <- import(unlist(D), seqinfo = hum_seqinfo)
  
  for (S in SPECIES$genome_assembly) {
    
    chains <- readChainFile(S) 
    
    # check if domains inside chains
    dm_enclosed_by_chain <- overlapsAny(domains, chains, type = "within")
    
    for (THR in THRESHOLDS){
      
      # for each syntenic block size threshold, check if a rearrangement breakpoint 
      # occurrs inside a domain with a safety margin of at least 40 kb to each TAD boundary
      breakpoints <- readBPFile(S, THR)
    
      # Calculate the distance to closest breakpoint
      dist_domains_to_bp <- distanceToNearest(domains, breakpoints)
      
      # Filter out 'rearranged' and 'Conserved' domains from <<Hits object>>
      # Note: Breakpoints closely outside of domain boundaries are ALSO NOT considered
      dist_rearr_domains <- dist_domains_to_bp[mcols(dist_domains_to_bp)$distance == 0]
      dist_intact_domains <- dist_domains_to_bp[mcols(dist_domains_to_bp)$distance > min_dist_bp_to_boundary]
      
      # Calculate a boolean vector that indicates which domains have a breakpoint NOT too close to a boundary
      # TRUE = Breakpoint far enough from boundaries
      # FALSE = Breakpoint too close to boundary
      savely_rearranged <- mapply(function(dm_idx, bp_idx){
        # Select rearranged domain
        dm <- domains[dm_idx, ]
        # Extract domain boundaries
        boundaries <- extractBoundaries(dm, 0)
        # Select breakpoint that overlaps domain
        bp <- breakpoints[bp_idx, ]
        # Calculate distances between both boundaries and the breakpoint
        dist <- distance(boundaries, bp)
        # If at least one of the distances is smaller than 'min_distance', return FALSE 
        if (sum(dist < min_dist_bp_to_boundary) >= 1) {
          return(FALSE)
        }
        return(TRUE)
      }, queryHits(dist_rearr_domains), subjectHits(dist_rearr_domains))
      
      # Filter out savely rearranged domains from <<Hits object>>
      dist_rearr_domains <- dist_rearr_domains[savely_rearranged]
      
      # Write labels into vector
      domain_rearranged <- vector(length=length(domains))
      domain_rearranged[queryHits(dist_rearr_domains)] <- TRUE
      domain_rearranged[queryHits(dist_intact_domains)] <- FALSE
    
    # Coerce to tibble
    is_rearranged <- tibble(tmp = is_rearranged)
    names(is_rearranged) <- unname(unlist(COLOQUIAL_NAMES[species]))
    are_TADs_rearranged <- cbind(are_TADs_rearranged, is_rearranged)
  }
  
}

# Coerce to tibble
are_TADs_rearranged <- as.tibble(are_TADs_rearranged)

# get rearranged TADs
are_TADs_rearranged_bp <- readRDS("results/dataframes/are_TADs_rearranged_by_breakpoints_dixon")

are_TADs_rearranged_bp_10k <- are_TADs_rearranged_bp %>% 
  filter(threshold == 10000) %>%
  select(-threshold, -rearrByTreeIndex) %>%
  gather("species", "Rearranged", 2:9) %>%
  mutate(species = factor(species))

are_TADs_rearranged_bp_100k <- are_TADs_rearranged_bp %>% 
  filter(threshold == 100000) %>%
  select(-threshold, -rearrByTreeIndex) %>%
  gather("species", "Rearranged", 2:9) %>%
  mutate(species = factor(species))

# get conserved TADs 
are_TADs_rearranged_fl <- readRDS("results/dataframes/are_TADs_rearranged_by_fills_dixon")
# "flip" boolean table
are_TADs_conserved_fl <- are_TADs_rearranged_fl %>%
  mutate_at(2:9, funs(!.)) %>%
  gather("species", "Conserved", 2:9) %>%
  mutate(species = factor(species))


# find conserved and rearranged TADs with stricter rules
data_df <- tibble(DomainIndex = are_TADs_conserved_fl$DomainIndex,
                  species = are_TADs_conserved_fl$species,
                  Conserved = are_TADs_conserved_fl$Conserved & !are_TADs_rearranged_bp_10k$Rearranged,
                  Rearranged = are_TADs_rearranged_bp_10k$Rearranged & !are_TADs_conserved_fl$Conserved
)

are_TADs_consv_and_rearr <- readRDS("results/dataframes/are_TADs_consv_fills_and_rearr_bp_double_condition_10k_1000k")

consv <- pull(are_TADs_consv_and_rearr %>% 
                filter(species == "Mouse"), Conserved)

rearr <- pull(are_TADs_consv_and_rearr %>%
                filter(species == "Mouse"), Rearranged)

intact_domains <- domains[consv]
rearr_domains <- domains[rearr]

# -------------------------------------------------------------------------------------------------------------------
# Assign ortholog pairs to TAD categories
# -------------------------------------------------------------------------------------------------------------------   

# Get human transcription start sites, gene ids in human expr as filter
hum_tss_raw <- getBM(attributes=c(
    'ensembl_gene_id', 'chromosome_name', 'strand',
    'transcription_start_site', 'transcript_start', 'transcript_end'), 
    filters=c('ensembl_gene_id', 'biotype'), 
    values=list(human_expr$ensembl_gene_id, biotype="protein_coding"), mart = ensembl)

# Filter out duplicated entries (take only TSS from largest transcript)
hum_tss_raw$transcript_length <- hum_tss_raw$transcript_end - hum_tss_raw$transcript_start
hum_tss_raw <- hum_tss_raw[with(hum_tss_raw, order(ensembl_gene_id, -transcript_length)), ]
hum_tss_raw <- hum_tss_raw[!duplicated(hum_tss_raw$ensembl_gene_id), ]

# Delete the spare columns
hum_tss_raw <- hum_tss_raw[ ,1:4]

# Adapt UCSC name style and filter names
hum_tss_raw <- hum_tss_raw %>%
  mutate(chr = str_c("chr", chromosome_name)) %>%
  filter(chr %in% seqnames(hum_seqinfo))

# Generate GRanges
hum_tss <- GRanges(seqnames = hum_tss_raw$chr,
                   strand = ifelse(hum_tss_raw$strand == 1, "+", "-"),
                   ranges = IRanges(start = hum_tss_raw$transcription_start_site,
                                    width = 1),
                   geneID = hum_tss_raw$ensembl_gene_id,
                   seqinfo = hum_seqinfo
)

# Combine with human gene IDs in data frame
results <- tibble(gene_id = human_expr$ensembl_gene_id, cor = correlations)

# Check if genes were subject to a TAD rearrangement in second species
# if breakpoint of second species within TAD -> rearranged, otherwise not rearranged or outside
gene_location <- geneLocationRegardingDomains(rearr_domains, intact_domains, domains, hum_tss)
hum_tss$class <- gene_location
# check also for genes in grb tads
hum_tss$grb <- overlapsAny(hum_tss, grb_domains)

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
