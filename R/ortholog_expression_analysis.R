# =============================================================================================================================
# Comparing expression levels of orthologous genes with regard to conserved TADs and TAD rearrangements
# =============================================================================================================================

require(BSgenome.Hsapiens.UCSC.hg19)
require(biomaRt)
require(tidyverse)
require(stringr)

source("R/functions.R")

# # Read metadata for analysis
# METADATA <- read_tsv("metadata.csv")
# SPECIES <- dplyr::select(METADATA, genome_assembly, trivial_name)
# DOMAINS <- METADATA %>% 
#   dplyr::select(genomic_domain_type, genomic_domain_path) %>%
#   filter(!is.na(genomic_domain_type))
# THRESHOLDS <- unlist(METADATA %>% 
#                        filter(!is.na(min_size_threshold)) %>%
#                        dplyr::select(min_size_threshold)
# )

# Read metadata for analysis
SPECIES <- read_tsv("species_meta.tsv")
DOMAINS <- read_tsv("domains_meta.tsv")
METADATA <- read_tsv("metadata.tsv")

# minimum size threasholds for chains (syntenic blocks) to be considered rearrangement blocks
THRESHOLDS <- unlist(METADATA %>% 
                       filter(!is.na(min_size_threshold)) %>%
                       dplyr::select(min_size_threshold)
)

#Conserved domains do NOT have breakpoints from all thresholds (smallest size threshold)
CONSV_BP_THR <- THRESHOLDS[1]
# Rearranged domains have breakpoints from the largest threshold
REARR_BP_THR <- THRESHOLDS[3]

# Load human seqinfo
genome <- BSgenome.Hsapiens.UCSC.hg19
hum_seqinfo <- seqinfo(genome)
# =============================================================================================================================
# read domains with classification in conserved and rearraged
# =============================================================================================================================
domain_classes <- read_rds("results/domain_classification.rds")

# =============================================================================================================================
# Gene expression correlation of orthologs
# =============================================================================================================================

# define column types
col_types_exp_data = cols(
  .default = col_double(),
  `Gene ID` = col_character(),
  `Gene Name` = col_character()
)

# read expression data of mathcing tissues
human_exp_match <- read_tsv(
  "data/ExpressionAtlas/formated.human_exp.matching_subset.tsv", 
  col_types = col_types_exp_data)
mouse_exp_match <- read_tsv(
  "data/ExpressionAtlas/formated.mouse_exp.matching_subset.tsv",
   col_types = col_types_exp_data)

# -------------------------------------------------------------------------------------------
# Find human-mouse orthologs
# ------------------------------------------------------------------------------------------- 

# Load gene ensembl human gene ids and orthologs
ensembl <- useMart(host = "grch37.ensembl.org", 
                   biomart = "ENSEMBL_MART_ENSEMBL", 
                   dataset = "hsapiens_gene_ensembl")

# Attributes that are returned
orth_attr = c("ensembl_gene_id",  "mmusculus_homolog_ensembl_gene", 
              "mmusculus_homolog_orthology_type", 
              "mmusculus_homolog_orthology_confidence")

# filters = c("biotype", "chromosome_name"),
# values = list(
#   biotype = "protein_coding", 
#   chromosome_name = c(1:22, "X", "Y")),


# Query orthologs by human gene ID
orthologsDF = getBM(attributes = orth_attr, mart = ensembl) 

# filter all human genes to have an "one2one" ortholog in mouse with confidance  = 1
orthologs <- as.tibble(orthologsDF) %>% 
  filter(
    mmusculus_homolog_ensembl_gene != "",
    mmusculus_homolog_orthology_type == "ortholog_one2one",
    mmusculus_homolog_orthology_confidence == 1) 


# filter orthologs for genes contained in human expression data
orthologs_with_exp <- orthologs %>%
  filter(ensembl_gene_id %in% human_exp_match$`Gene ID`) %>% 
  filter(mmusculus_homolog_ensembl_gene %in% mouse_exp_match$`Gene ID`) %>% 
  select(-mmusculus_homolog_orthology_type, -mmusculus_homolog_orthology_confidence)

# --------------------------------------------------------------------------------------------
# Prepare data for correlation
# --------------------------------------------------------------------------------------------  
# Order genes in human_exp_match according to ordering in orthologs
human_expr <- rename(human_exp_match,  
                     "ensembl_gene_id" = "Gene ID",
                     "human_gene_name" = "Gene Name")
human_expr <- left_join(orthologs_with_exp, human_expr, by = "ensembl_gene_id")

# Order genes in mouse_exp_match according to ordering in orthologs
mouse_expr <- dplyr::rename(mouse_exp_match, 
                     "mmusculus_homolog_ensembl_gene" = "Gene ID",
                     "mouse_gene_name" = "Gene Name")
mouse_expr <- left_join(orthologs_with_exp, mouse_expr, by = "mmusculus_homolog_ensembl_gene")

# --------------------------------------------------------------------------------------------
# Calculate correlations
# --------------------------------------------------------------------------------------------  

cor_values <- map_dbl(1:nrow(human_expr), 
                      ~ cor(as.numeric(human_expr[.x, 4:ncol(human_expr)]), 
                            as.numeric(mouse_expr[.x, 4:ncol(mouse_expr)]),
                            method = "pearson")
                      ) 

# # make correlations identifiable by human gene id
# correlations <- tibble(
#   ensembl_gene_id = human_expr$ensembl_gene_id,
#   correlation = cor_values)

# Build table with ortholy assignment, and correlation
correlationDF <- orthologs_with_exp %>% 
  mutate(
    correlation = cor_values
  )

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
  mart = ensembl)
)

tssDF <- as.tibble(hum_tss_df) %>% 
  # Filter out duplicated entries (take only TSS from largest transcript)
  mutate(transcript_length = transcript_end - transcript_start) %>% 
  arrange(ensembl_gene_id, desc(transcript_length)) %>% 
  distinct(ensembl_gene_id, .keep_all = TRUE) %>% 
  # Delete the spare columns
  select(1:5) %>% 
  # Adapt UCSC name style
  mutate(chr = str_c("chr", chromosome_name)) %>%
  # Filter for valid names
  filter(chr %in% seqnames(hum_seqinfo))


# 
# hum_tss_df <- hum_tss_df[ ,1:4]
# 
# hum_tss_df <- hum_tss_df %>%
#   # Adapt UCSC name style
#   mutate(chr = str_c("chr", chromosome_name)) %>%
#   # Filter for valid names
#   filter(chr %in% seqnames(hum_seqinfo))

# Generate GRanges
hum_tss_gr <- GRanges(seqnames = tssDF$chr,
                   strand = ifelse(tssDF$strand == 1, "+", "-"),
                   ranges = IRanges(start = tssDF$transcription_start_site,
                                    width = 1),
                   geneID = tssDF$ensembl_gene_id,
                   gene_name = tssDF$external_gene_name,
                   seqinfo = hum_seqinfo
                   )
write_rds(hum_tss_gr, "results/ortholog_expression.hum_tss_gr.rds")
# --------------------------------------------------------------------------------------------
# Start building final result tibble
# --------------------------------------------------------------------------------------------  

# add ortholog gene info (if present) to ALL tss retrieved
inter_results <- tssDF %>%
  left_join(dplyr::select(human_expr, 1:4), by = "ensembl_gene_id")

# add the correlation results to otrholog pairs
inter_results <- inter_results %>%
  left_join(correlationDF, by = "ensembl_gene_id")

# --------------------------------------------------------------------------------------------
# Assign the domains to ortholog pair in 'results' if overlapping with human TSS,
# then denote conserved or rerranged status of each domain
# -------------------------------------------------------------------------------------------- 

domain_classes <- read_rds("results/domain_classification_fills.rds")

tmp <- tibble()

for (D in DOMAINS$domain_path) {
  
  domains <- import(unlist(D), seqinfo = hum_seqinfo)
  
  # get domain type to store with every result
  domain_type <- DOMAINS %>%
    filter(domain_path == D) %>%
    dplyr::select(domain_type) %>% 
    pull()

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

write_rds(results, "results/ortholog_expression_correlation.rds")

