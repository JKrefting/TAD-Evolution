# =============================================================================================================================
# Comparing expression levels of orthologous genes with regard to conserved TADs and TAD rearrangements
# =============================================================================================================================

require(tidyverse)
require(biomaRt)
require(stringr)
require(BSgenome.Hsapiens.UCSC.hg19)

source("R/functions.R")

# Read metadata for analysis
SPECIES <- read_tsv("species_meta.tsv")
DOMAINS <- read_tsv("domains_meta.tsv")
METADATA <- read_tsv("metadata.tsv")

# minimum size threasholds for chains (syntenic blocks) to be considered rearrangement blocks
THRESHOLDS <- unlist(METADATA %>% 
                       filter(!is.na(min_size_threshold)) %>%
                       dplyr::select(min_size_threshold)
)

# Conserved domains do NOT have breakpoints from all thresholds (smallest size threshold)
CONSV_BP_THR <- THRESHOLDS[1]
# Rearranged domains have breakpoints from the largest threshold
REARR_BP_THR <- THRESHOLDS[3]

# Load human seqinfo
genome <- BSgenome.Hsapiens.UCSC.hg19
hum_seqinfo <- seqinfo(genome)

# =============================================================================================================================
# Gene expression correlation of orthologs
# =============================================================================================================================

# read expression data tables
# human_exp_match <- read_tsv("data/ExpressionAtlas/formated.human_exp.matching_subset.tsv")
# mouse_exp_match <- read_tsv("data/ExpressionAtlas/formated.mouse_exp.matching_subset.tsv")
human_exp_match <- read_tsv("../../data/gxa/FANTOM5.human_exp.matching_subset.tsv")
mouse_exp_match <- read_tsv("../../data/gxa/FANTOM5.mouse_exp.matching_subset.tsv")

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

# genes without ortholog should be eliminated (for safety)
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
hum_tss_df <- hum_tss_df %>% dplyr::select(-transcript_start, -transcript_end, -transcript_length)

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


domain_classes <- read_rds("results/domain_classification.rds")

# categorise domains as conserved, identifiable by domain id
consv_dm <- domain_classes %>%
  filter(species == "mm10",
         threshold == CONSV_BP_THR
  ) %>%
  mutate(conserved = enclosed_by_chain & !rearranged_by_breakpoint) %>%
  dplyr::select(domain_id, domain_type, conserved)

# categorise domains as rearranged
rearr_dm <- domain_classes %>%
  filter(species == "mm10",
         threshold == REARR_BP_THR
  ) %>%
  mutate(rearranged = !enclosed_by_chain & rearranged_by_breakpoint) %>%
  dplyr::select(domain_id, domain_type, rearranged)

results <- tibble()

for (D in DOMAINS$domain_path) {
  
  domains <- import(unlist(D), seqinfo = hum_seqinfo)
  
  this_domain_type <- unlist(DOMAINS %>%
                               filter(domain_path == D) %>%
                               dplyr::select(domain_type)
  )
  
  # --------------------------------------------------------------------------------------------
  # assign domains to orthologs
  # -------------------------------------------------------------------------------------------- 
  
  # replicate correlation results (inter) with every domain type
  results_loop <- inter_results
  
  # determine gene association to domain
  tss_domain_hits <- findOverlaps(hum_tss_gr, domains)
  
  # add associated domain (if present) to orthologs
  # NOTE: order of genes in hum_tss_gr and results_loop (inter_results) equal
  results_loop$domain_id <- NA
  results_loop[queryHits(tss_domain_hits), ]$domain_id <- subjectHits(tss_domain_hits)
  results_loop$domain_type <- this_domain_type
  
  # --------------------------------------------------------------------------------------------
  #  assign domain categories (conserved, rearranged, outside) to orthologs
  # -------------------------------------------------------------------------------------------- 

  # filter for current domain type
  this_consv_dm <- consv_dm %>% filter(domain_type == this_domain_type)
  this_rearr_dm <- rearr_dm %>% filter(domain_type == this_domain_type)
  
  # add conserved / rearranged info of domains by merging by domain id
  results_loop$conserved <- this_consv_dm[match(results_loop$domain_id, this_consv_dm$domain_id), ]$conserved 
  results_loop$rearranged <- this_rearr_dm[match(results_loop$domain_id, this_rearr_dm$domain_id), ]$rearranged 
  
  # determine genes outside any domain
  results_loop$outside <- !overlapsAny(hum_tss_gr, domains)
  
  results <- rbind.data.frame(results, results_loop)
}

# --------------------------------------------------------------------------------------------
# Finish categorisation of orthpairs localised in conserved, rearranged or outside domains
# -------------------------------------------------------------------------------------------- 

results$category <- NA
results$category <- ifelse(results$conserved, "Conserved", results$category)
results$category <- ifelse(results$rearranged, "Rearranged", results$category)
results$category <- ifelse(results$outside, "Outside", results$category)
results$category <- factor(results$category, levels = c("Conserved", "Rearranged", "Outside"))

write_tsv(results, "results/ortholog_expression_correlation_old.tsv")
write_rds(results, "results/ortholog_expression_correlation_old.rds")

