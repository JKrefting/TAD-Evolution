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
require(TxDb.Hsapiens.UCSC.hg38.knownGene)
# require(BSgenome.Mmusculus.UCSC.mm10)
require(TxDb.Mmusculus.UCSC.mm10.knownGene)
require(biomaRt)
require(tidyverse)
require(stringr)

source("R/functions.R")

ENSEMBL_URL = "aug2017.archive.ensembl.org"

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


# Get mouse genes as GRanges --------------------------------------------------#

# Load mouse seqinfo
# mouse_seqinfo <- seqinfo(BSgenome.Mmusculus.UCSC.mm10)

MOUSE_ENSEMBL_STR = "mmusculus"

# Load gene ensembl mouse gene ids and orthologs
mouse_ensembl <- useMart(host = ENSEMBL_URL, 
                   biomart = "ENSEMBL_MART_ENSEMBL", 
                   dataset = "mmusculus_gene_ensembl")

# Get mouse transcription start sites
mouse_tss_df <- as.tibble(getBM(attributes=c(
  'ensembl_gene_id', 'external_gene_name', 'chromosome_name', 'strand',
  'transcription_start_site', 'transcript_start', 'transcript_end'), 
  mart = mouse_ensembl)
) 

mouseTssDF <- mouse_tss_df %>% 
  # Filter out duplicated entries (take only TSS from largest transcript)
  mutate(transcript_length = transcript_end - transcript_start) %>% 
  arrange(ensembl_gene_id, desc(transcript_length)) %>% 
  distinct(ensembl_gene_id, .keep_all = TRUE) %>% 
  # Adapt UCSC name style
  mutate(chr = str_c("chr", chromosome_name)) 
# %>%
#   # Filter for valid names
#   filter(chr %in% seqnames(mouse_seqinfo))

# Generate GRanges
mouseTssGR <- GRanges(seqnames = mouseTssDF$chr,
                 strand = ifelse(mouseTssDF$strand == 1, "+", "-"),
                 ranges = IRanges(start = mouseTssDF$transcription_start_site,
                                  width = 1),
                 geneID = mouseTssDF$ensembl_gene_id,
                 gene_name = mouseTssDF$external_gene_name
                 # ,
                 # seqinfo = mouse_seqinfo
)


# Get human genes with orthologs in mouse -------------------------------------#

# Attributes that are returned
orth_attr = c("ensembl_gene_id",  "mmusculus_homolog_ensembl_gene", 
              "mmusculus_homolog_orthology_type", 
              "mmusculus_homolog_orthology_confidence")

# Query orthologs by human gene ID
orthologsDF = getBM(attributes = orth_attr, mart = ensembl) 

# filter all human genes to have an "one2one" ortholog in mouse
orthologs <- as.tibble(orthologsDF) %>% 
  filter(
    mmusculus_homolog_ensembl_gene != "",
    mmusculus_homolog_orthology_type == "ortholog_one2one"
  ) %>% 
  select(ensembl_gene_id, mmusculus_homolog_ensembl_gene)


# Get adjacent gene pairs  ----------------------------------------------------#

human_ortholog_GR <- tssGR[tssGR$geneID %in% orthologs$ensembl_gene_id]
mouse_ortholog_GR <- mouseTssGR[mouseTssGR$geneID %in% orthologs$mmusculus_homolog_ensembl_gene]

getAdjacentPairs <- function(tssGR){

  # get the next tss for each gene along the chromsome
  # using the precede() function from GenomicRanges package
  nextGene = precede(tssGR, ignore.strand = TRUE)
  
  firstGR <- tssGR[!is.na(nextGene)]
  nextGR <- tssGR[nextGene[!is.na(nextGene)]]

  gP = tibble(
      g1 = firstGR$geneID, 
      g2 = nextGR$geneID,
      dist = abs(start(firstGR) - start(nextGR)),
      strand = case_when(
        as.logical(strand(firstGR) == strand(nextGR)) ~ "same",
        as.logical(strand(firstGR) == "+" & strand(nextGR) == "-") ~ "convergent",
        as.logical(strand(firstGR) == "-" & strand(nextGR) == "+") ~ "divergent"
      )
    )

}

human_adjacent <- getAdjacentPairs(human_ortholog_GR)
mouse_adjacent <- getAdjacentPairs(mouse_ortholog_GR)

# build gene pair data set  ---------------------------------------------------#

# douplicate mouse adjacent pairs to have A-B and B-A pairs
mouse_adjacent_all <- bind_rows(
    mouse_adjacent,
    rename(mouse_adjacent, g1 = g2, g2 = g1)
  ) %>% 
  rename(dist_mouse = dist, strand_mouse = strand) %>% 
  mutate(adjacent_mouse = TRUE)

df <- human_adjacent %>% 
  # add ortholog to first gene
  left_join(orthologs, by = c("g1" = "ensembl_gene_id")) %>% 
  rename(g1_mouse = mmusculus_homolog_ensembl_gene) %>% 
  # add ortholog to second gene
  left_join(orthologs, by = c("g2" = "ensembl_gene_id")) %>% 
  rename(g2_mouse = mmusculus_homolog_ensembl_gene) %>% 
  # test if mouse orthologs are adjacent
  left_join(mouse_adjacent_all, by = c("g1_mouse" = "g1", "g2_mouse" = "g2")) %>% 
  mutate(
    adjacent_mouse = !is.na(adjacent_mouse),
    same_strand = strand == strand_mouse,
    syntenic = adjacent_mouse & same_strand
    )

# same_chrom = 
#   strand(mouse_ortholog_GR[match(g1_mouse, mouse_ortholog_GR$geneID)]) == strand(mouse_ortholog_GR)[match(g2_mouse, mouse_ortholog_GR$geneID)],

# dDF <- df %>% 
#   count(adjacent_mouse)
# 
# p <- ggplot(df, aes(x = dist, y = dist_mouse)) +
#   geom_point(alpha = .25) +
#   scale_x_log10() + scale_y_log10()

# Add breakpoints to gene pairs  ----------------------------------------------#

# build GRanges for span of adjacent paris

gr1 = human_ortholog_GR[match(df$g1, human_ortholog_GR$geneID)]
gr2 = human_ortholog_GR[match(df$g2, human_ortholog_GR$geneID)]

adjacentGR <- GRanges(
  seqnames = seqnames(gr1),
  IRanges(start(gr1), start(gr2)),
  seqinfo = seqinfo(gr1)
)

# load breakpoint data 
# Read metadata for analysis
SPECIES <- read_tsv("species_meta.tsv")
# DOMAINS <- read_tsv("domains_meta.tsv")
THRESHOLDS <- read_tsv("metadata.tsv") %>% 
  pull(min_size_threshold)

# DOMAINS <- read_tsv("domains_meta.tsv")

# iterate over all size thresholds for rearragement breakpoints
# for (THR in THRESHOLDS){

# read_breakpoints <- function(species, threshold){
#   bp_file <- paste0("data/breakpoints/hg19.", species, ".", 
#                     as.character(format(threshold, scientific = FALSE)), ".bp.bed")
#   return(import(bp_file, seqinfo = hum_seqinfo))
# }

rearraged = THRESHOLDS %>%
  set_names(paste0("breakpoints_", THRESHOLDS)) %>%

  # build path to breakpoint file
  map(~ paste0("data/breakpoints/hg38.", "mm10", ".",
               as.character(format(.x, scientific = FALSE)),  ".bp.bed")) %>%
  # read breakopoints
  map(import.bed, seqinfo = hum_seqinfo) %>% 
  
  # get rearraged pairs by overlap with adjacent region
  map(~ overlapsAny(adjacentGR, .x)) %>% 
  as.tibble()

# add rearraged columns to df
df <- df %>% 
  bind_cols(rearraged)

# Write adjacent pairs data to ouptut files -----------------------------------#

write_tsv(df, "results/syntenic_gene_pairs.mouse.tsv")  
write_rds(df, "results/syntenic_gene_pairs.mouse.rds")  

# df <- read_rds("results/syntenic_gene_pairs.mouse.rds")

# Analyse syntenic gene pairs and rearragements ===============================#

tidyDF <- df %>% 
  gather("size_threshold", "breakpoint", starts_with("breakpoints_")) %>% 
  separate(size_threshold, c("bp", "size_threshold"), sep = "_") %>% 
  mutate(
    size_threshold = parse_integer(size_threshold)
  ) %>% 
  select(-bp)

countDF <- tidyDF %>% 
  group_by(size_threshold) %>% 
  count(breakpoint)

# test dist and dist in mouse
p <- tidyDF %>% 
  filter(dist <= size_threshold) %>% 
  ggplot(aes(x = dist, y = dist_mouse)) +
  geom_point(size = .5, alpha = .25) + 
  geom_abline(intercept = 0, slope = 1, color = "red") +
  # stat_density2d(aes(fill = ..density..^0.25), geom = "tile", contour = FALSE, n = 200) +
  # scale_fill_continuous(low = "white", high = "dodgerblue4") +
  scale_x_log10() + scale_y_log10() +
  facet_grid(size_threshold ~ breakpoint) + 
  theme_bw()
p

fdrDF <- tidyDF %>% 
  # filter out gene pairs with larger distance than size threshold
  filter(dist <= size_threshold) %>% 
  group_by(size_threshold) %>%
  summarize(
    n = n(),
    n_adjacent = sum(adjacent_mouse, na.rm = TRUE),
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
    sensitivity = TP / (TP + FN),
    specificity = TN / (TN + FP)
  )


write_tsv(fdrDF, "results/syntenic_gene_pairs.hg38_mouse.fdrDF.tsv")  
write_rds(fdrDF, "results/syntenic_gene_pairs.hg38_mouse.fdrDF.rds")  
# fdrDF <- read_rds("results/syntenic_gene_pairs.hg38_mouse.fdrDF.rds")

# example testing -------------------------------------------------------------#

# false positiv breakpoints:
tidyDF %>% filter(size_threshold == 100000, dist <= size_threshold, syntenic, breakpoint) %>% 
  mutate(
    gene_name_1 = tssGR$gene_name[match(g1, tssGR$geneID)],
    gene_name_2 = tssGR$gene_name[match(g2, tssGR$geneID)]
  ) %>% 
  write_tsv("results/syntenic_gene_pairs.hg38_mouse.false_positive_pairs_100000.tsv")

# false positiv breakpoints:
tidyDF %>% filter(size_threshold == 1000000, dist <= size_threshold, syntenic, breakpoint) %>% 
  mutate(
    gene_name_1 = tssGR$gene_name[match(g1, tssGR$geneID)],
    gene_name_2 = tssGR$gene_name[match(g2, tssGR$geneID)]
  ) %>% 
  write_tsv("results/syntenic_gene_pairs.hg38_mouse.false_positive_pairs_1000000.tsv")


tssGR[match(c("ENSG00000103202", "ENSG00000242612"), tssGR$geneID)]
tssGR[match(c("ENSG00000103202", "ENSG00000242612"), tssGR$geneID)]


