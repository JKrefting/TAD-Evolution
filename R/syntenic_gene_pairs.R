#'##############################################################################
#' Analyse syntenic gene pairs betwen human and other species to assess 
#' breakpoint identification accuracy.
#' 
#' Goal: Get FDR = FP / (TP + FP) for rearrangement identification, where
#' FP : gene pairs with rearragement breakpoint but adjacent in mouse
#' TP : gene pair without rearragement breakpoint and adjacent in mouse
#' 
#' Plan:
#' - get list of all human adjacent genes and annotate with the follwoing
#'   - Intergenic distance
#'   - adjacent in mouse?
#'      - same chrom
#'      - same strand combination
#'   - breakpoint between human genes
#'   
#' Issues: 
#'  - update readBPFile() funciton to include seqinfo as argument
#'  
#'##############################################################################

# require(BSgenome.Hsapiens.UCSC.hg19)
# require(BSgenome.Hsapiens.UCSC.hg38)
# require(BSgenome.Mmusculus.UCSC.mm10)
require(TxDb.Hsapiens.UCSC.hg38.knownGene)
require(TxDb.Mmusculus.UCSC.mm10.knownGene)
require(biomaRt)
source("R/functions.R")
require(tidyverse)
require(stringr)


ENSEMBL_URL = "aug2017.archive.ensembl.org"

SPECIES <- read_tsv("species_meta.tsv")
THRESHOLDS <- read_tsv("metadata.tsv") %>% pull(min_size_threshold)

# Get human genes as GRanges --------------------------------------------------#

# Load human seqinfo
# genome <- BSgenome.Hsapiens.UCSC.hg19
hum_seqinfo <- seqinfo(TxDb.Hsapiens.UCSC.hg38.knownGene)

# Load gene ensembl human gene ids and orthologs
# ensembl <- useMart(host = "grch37.ensembl.org", 
ensembl <- useMart(host = ENSEMBL_URL,
                   biomart = "ENSEMBL_MART_ENSEMBL", 
                   dataset = "hsapiens_gene_ensembl")

# Get human transcription start sites
hum_gene_df <- as.tibble(getBM(attributes = c(
  'ensembl_gene_id', 'external_gene_name', 'chromosome_name', 'strand',
  'start_position', 'end_position'), 
  mart = ensembl)
  ) 

geneDF <- hum_gene_df %>% 
  # Adapt UCSC name style
  mutate(chr = str_c("chr", chromosome_name)) %>%
  # Filter for valid names
  filter(chr %in% seqnames(hum_seqinfo))

# Generate GRanges
geneGR <- GRanges(seqnames = geneDF$chr,
                      strand = ifelse(geneDF$strand == 1, "+", "-"),
                      ranges = IRanges(start = geneDF$start_position,
                                       end = geneDF$end_position),
                      geneID = geneDF$ensembl_gene_id,
                      gene_name = geneDF$external_gene_name,
                      seqinfo = hum_seqinfo
)

# Load data of other species --------------------------------------------------#

# Load gene ensembl species gene ids and orthologs
mart <- useMart(host = ENSEMBL_URL, biomart = "ENSEMBL_MART_ENSEMBL")

martData <- listDatasets(mart) %>% as_tibble() %>% 
  mutate_all(as.character) %>% 
  write_tsv("results/ensemble_datasets.tsv")

speciesDF <- SPECIES %>% 
    mutate(
    ensembl_name = str_c(str_to_lower(str_sub(species_name, 1,1)),
                         str_split_fixed(species_name, " ", n = 3)[,2]),
    dataset = str_c(ensembl_name, "_gene_ensembl")
  ) %>% 
  left_join(martData, by = "dataset") %>% 
  write_tsv("results/species_ensemble_datasets.tsv")

# compute syntey and rearranged pairs -----------------------------------------#

speciesDF <- speciesDF %>% 
  filter(!is.na(description)) %>% 
  mutate(
    df = map2(ensembl_name, genome_assembly, get_syntenic_pairs, 
              geneGR = geneGR, 
              ensembl_url = ENSEMBL_URL,
              size_thresholds = THRESHOLDS)
  )

write_rds(speciesDF, "results/speciesDF.rds")
#speciesDF <- read_rds("results/speciesDF.rds")

syntenicDF <- speciesDF %>% 
  select(genome_assembly, trivial_name, df) %>% 
  unnest(df) %>% 
  gather("size_threshold", "breakpoint", starts_with("breakpoints_")) %>% 
  separate(size_threshold, c("bp", "size_threshold"), sep = "_") %>% 
  mutate(size_threshold = parse_integer(size_threshold)) %>% 
  select(-bp)

write_rds(syntenicDF, "results/syntenicDF.rds")

gene_pair_stats <- syntenicDF %>% 
  filter(
    # filter out gene pairs with larger distance than size threshold
    dist <= size_threshold, 
    # if syntenic filter also for species distance smaller than threshold
    !syntenic | dist_species <= size_threshold
  ) %>% 
  group_by(genome_assembly, size_threshold, syntenic, breakpoint) %>% 
  summarize(n = n()) %>% 
  mutate(percent = n / sum(n) * 100) %>% 
  ungroup() %>% 
  mutate(
    syntenic = ifelse(syntenic, "syntenic", "non-syntenic"),
    rearranged = ifelse(breakpoint, "rearranged", "non-rearranged")
    ) %>% 
  right_join(SPECIES, by = "genome_assembly") %>% 
  filter(!is.na(size_threshold)) %>% 
  write_tsv("results/syntenic_genes.gene_pair_stats.tsv")


# Get statistics about species, gene pairs. 
# We need: (for species and threshold)
# - adjacent_genes_with_orthologs
# - with_breakpoint
# - 
gene_pair_counts <- syntenicDF  %>% 
  group_by(genome_assembly, size_threshold) %>% 
  summarize(adjacent_genes_with_orthologs = n()) %>% 
  select(-size_threshold) %>% 
  distinct() 

gene_pair_counts_out <- SPECIES %>% 
  left_join(gene_pair_counts, by = "genome_assembly") %>% 
  write_tsv("results/gene_pair_counts.tsv")


syntenic_performance_DF <- syntenicDF %>% 
  filter(
    
    # filter out gene pairs with larger distance than size threshold
    dist <= size_threshold, 
    
    # if syntenic filter also for species distance smaller than threshold
    !syntenic | dist_species <= size_threshold
    
    ) %>%
  group_by(trivial_name, size_threshold) %>%
  summarize(
    n = n(),
    n_adjacent = sum(adjacent_species, na.rm = TRUE),
    n_same_strand = sum(same_strand, na.rm = TRUE),
    n_syntenic = sum(syntenic, na.rm = TRUE),
    n_breakpoint = sum(breakpoint, na.rm = TRUE),
    TP = sum(!syntenic & breakpoint, na.rm = TRUE),
    FP = sum(syntenic & breakpoint, na.rm = TRUE),
    TN = sum(syntenic & !breakpoint, na.rm = TRUE),
    FN = sum(!syntenic & !breakpoint, na.rm = TRUE)
  ) %>% 
  mutate(
    FDR = FP / (FP + TP),
    FPR = (FP / (FP + TN)),
    PPV = TP / (TP + FP),
    sensitivity = TP / (TP + FN),
    specificity = TN / (TN + FP)
  )

syntenic_performance_DF <- SPECIES %>% 
  left_join(syntenic_performance_DF, by = "trivial_name")

write_rds(syntenic_performance_DF, "results/syntenic_performance_DF.rds")
write_tsv(syntenic_performance_DF, "results/syntenic_performance_DF.tsv")
