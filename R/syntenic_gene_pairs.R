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
#'   - TSS distance
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
require(tidyverse)
require(stringr)

source("R/functions.R")

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
hum_tss_df <- as.tibble(getBM(attributes = c(
  'ensembl_gene_id', 'external_gene_name', 'chromosome_name', 'strand',
  'transcription_start_site', 'transcript_start', 'transcript_end'), 
  mart = ensembl)
  ) 

tssDF <- hum_tss_df %>% 
  # Filter out duplicated entries (take only TSS from largest transcript)
  mutate(transcript_length = transcript_end - transcript_start) %>% 
  arrange(ensembl_gene_id, desc(transcript_length)) %>% 
  distinct(ensembl_gene_id, .keep_all = TRUE) %>% 
  # Adapt UCSC name style
  mutate(chr = str_c("chr", chromosome_name)) %>%
  # Filter for valid names
  filter(chr %in% seqnames(hum_seqinfo))

# Generate GRanges
tssGR <- GRanges(seqnames = tssDF$chr,
                      strand = ifelse(tssDF$strand == 1, "+", "-"),
                      ranges = IRanges(start = tssDF$transcription_start_site,
                                       width = 1),
                      geneID = tssDF$ensembl_gene_id,
                      gene_name = tssDF$external_gene_name,
                      seqinfo = hum_seqinfo
)

# Load data of other species --------------------------------------------------#

# Load gene ensembl mouse gene ids and orthologs
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
  left_join(martData, by = "dataset")

# compute syntey and rearranged pairs -----------------------------------------#

speciesDF <- speciesDF %>% 
  filter(!is.na(description)) %>% 
  mutate(
    df = map2(ensembl_name, genome_assembly, get_syntenic_pairs, 
              tssGR = tssGR, 
              ensembl_url = "aug2017.archive.ensembl.org",
              size_thresholds = THRESHOLDS)
  )

write_rds(speciesDF, "results/speciesDF.rds")


syntenicDF <- speciesDF %>% 
  select(genome_assembly, trivial_name, df) %>% 
  unnest(df) %>% 
  gather("size_threshold", "breakpoint", starts_with("breakpoints_")) %>% 
  separate(size_threshold, c("bp", "size_threshold"), sep = "_") %>% 
  mutate(size_threshold = parse_integer(size_threshold)) %>% 
  select(-bp)

write_rds(syntenicDF, "results/syntenicDF.rds")


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
    PPV = TP / (TP + FP),
    sensitivity = TP / (TP + FN),
    specificity = TN / (TN + FP)
  )

write_rds(syntenic_performance_DF, "results/syntenic_performance_DF.rds")
write_tsv(syntenic_performance_DF, "results/syntenic_performance_DF.tsv")


# plot FDR ---------------------------------------------------------------------

p <- ggplot(syntenic_performance_DF, aes(x = trivial_name, y = FDR, fill = as.factor(size_threshold))) + 
  geom_bar(stat = "identity", position = position_dodge(1)) +
  geom_text(aes(label = round(FDR, 2)), position = position_dodge(1), hjust = "inward") +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "bottom")
ggsave("results/syntenic_genes.performance.FDR_by_species.pdf", w = 6, h = 6)

# plot ROC curve ---------------------------------------------------------------

p <- ggplot(syntenic_performance_DF, aes(x = 1- specificity, y = sensitivity, color = trivial_name)) + 
  geom_line() + 
  theme_bw() + 
  lims(x = c(0, 1), y = c(0, 1))
p
# plot PRC curve ---------------------------------------------------------------

p <- ggplot(syntenic_performance_DF, aes(x = sensitivity, y = PPV, color = trivial_name)) + 
  geom_line() + 
  theme_bw() + 
  lims(x = c(0, 1), y = c(0, 1))
p


# example testing -------------------------------------------------------------#

# false positiv breakpoints:
false_positives <- syntenicDF %>% 
  filter(
    trivial_name == "chimpanzee", 
    size_threshold == 100000, 
    dist <= size_threshold, 
    syntenic, 
    breakpoint
    ) %>% 
  mutate(
    gene_name_1 = tssGR$gene_name[match(g1, tssGR$geneID)],
    gene_name_2 = tssGR$gene_name[match(g2, tssGR$geneID)]
  ) %>%
  select(gene_name_1, gene_name_2, everything()) %>% 
  write_tsv("results/syntenic_gene_pairs.hg38_chimpanzee.false_positive_pairs_100000.tsv")


