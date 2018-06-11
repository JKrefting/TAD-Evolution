
# Refind breakpoint by searching for matching fills in adjacent regions ========


require(TxDb.Hsapiens.UCSC.hg38.knownGene)
source("R/functions.R")
require(tidyverse)

# Read metadata for analysis
SPECIES <- read_tsv("species_meta.tsv")
METADATA <- read_tsv("metadata.tsv")

# minimum size threasholds for chains (syntenic blocks) to be considered rearrangement blocks
THRESHOLDS <- unlist(METADATA %>% 
                       filter(!is.na(min_size_threshold)) %>%
                       dplyr::select(min_size_threshold)
)

# Load human seqinfo
hum_seqinfo <- seqinfo(TxDb.Hsapiens.UCSC.hg38.knownGene)

# -------------------------------------------------------------
# Fills
# -------------------------------------------------------------

# prepare data for plotting
# gather all date in one tibble

col_names = c("chr", "start", "end", "type", "chr_species", "strand_relative")
# , "start_species", "size_species")

# file = "data/fills/hg38.mm10.fill.bed.tsv"
# df <- read_tsv(file, col_names = col_names)

all_fill_df <- SPECIES %>% 
  select(genome_assembly, species = trivial_name) %>% 
  mutate(
    fill_tab_file = paste0("data/fills/hg38.", genome_assembly, ".fill.bed.tsv"),
    fill_df = map(fill_tab_file, read_tsv,
             col_names = col_names,
             col_types = list(chr = 'c', 
                              start = 'i',
                              end = 'i',
                              type = 'c',
                              chr_species = 'c',
                              strand_relative = 'c')
             ),
    fill_gr = map(fill_df, GRanges)
  )


# -------------------------------------------------------------
# Breakpoints
# -------------------------------------------------------------

bp_df <- crossing(
  genome_assembly = SPECIES$genome_assembly,
  size_threshold = THRESHOLDS
) %>% 
  mutate(
    bp_file = str_c("data/breakpoints/hg38.", genome_assembly, ".",
                    size_threshold,  ".bp.flt.bed"),
    bp = map(bp_file, import.bed, seqinfo = hum_seqinfo),
  )


# bpGR <- bp_df %>% filter(genome_assembly == "mm10", size_threshold == 100000) %>% pull(bp)
# bpGR <- bpGR[[1]]
# 
# fill_df <- all_fill_df %>% 
#   filter(genome_assembly == "mm10") %>% 
#   pull(df)
# fill_df <- fill_df[[1]]
# size_threshold = 100000

#' Test if differnt fills with same chr and strand are in flanking regions
flt_by_adjacent_fill <- function(bpGR, size_threshold, fill_df){
  
  
  stopifnot(all(width(bpGR) == 1))
  
  # filter by size
  fill_df <- fill_df %>% 
    mutate(size = end - start) %>% 
    filter(size >= size_threshold) %>% 
    mutate(id = row_number())
  
  # convert to GRanges object
  fillGR <- GRanges(fill_df)

  fills_left <- resize(bpGR, size_threshold, fix = "end") %>% 
    findOverlaps(fillGR) %>% 
    as.data.frame() %>% as.tibble()
  
  fills_right <- resize(bpGR, size_threshold, fix = "start") %>% 
    findOverlaps(fillGR) %>% 
    as.data.frame() %>% as.tibble()
  
  
  
  adj_fills <- full_join(fills_left, fills_right, 
                         by = "queryHits", suffix = c("_left", "_right")) %>% 
    
    # filter aut same fill in left and righ regions
    filter(subjectHits_left != subjectHits_right) %>% 
    
    # annotate fills with chrom and strand
    left_join(select(fill_df, id, chr, strand_relative), by = c("subjectHits_left" = "id")) %>% 
    left_join(select(fill_df, id, chr, strand_relative), by = c("subjectHits_right" = "id"),
              suffix = c("_left", "_right")) %>% 
    filter(chr_left == chr_right & strand_relative_left == strand_relative_right)
  
  # filter out false positive breakpoints
  fltGR <- bpGR[-unique(adj_fills$queryHits)]
  
  return(fltGR)
}


write_bed_wrapper <- function(gr, out_file, ...){
  if(length(gr) > 0){
    export.bed(gr, out_file, ...)
  }else{
    write("", out_file)
  }
}


# apply refinment to all breakpoint files
bp_df <- bp_df %>%  
  left_join(all_fill_df, by = "genome_assembly") %>% 
  mutate(
    flt_bp = pmap(list(bp, size_threshold, fill_df), flt_by_adjacent_fill),
    n_bp = map_int(bp, length),
    n_flt_bp = map_int(flt_bp, length),
    flt_ratio = n_flt_bp / n_bp * 100,
    out_file = str_replace(bp_file, ".bed$", ".flt_adj_fill.bed"),
    flt_bp = map2(flt_bp, out_file, write_bed_wrapper)
  )

# plot filtering rate
p <- ggplot(bp_df, aes(x = species, y = flt_ratio, fill = as.factor(size_threshold))) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = round(flt_ratio, 2), hjust = "inward"), position = position_dodge(.9)) +
  theme_bw() +
  coord_flip()
p


# analysis of size to nearest other breakpoint ---------------------------------
fpGR <- bpGR[unique(adj_fills$queryHits)]

fp_df <- fpGR %>%
  as.data.frame() %>% as.tibble() %>%
  mutate(
    str_loc = paste0(seqnames, ":", start-100000, "-", start + 100000)
  )

fltGR <- bpGR[-unique(adj_fills$queryHits)]
flt_df <- fltGR %>%
  as.data.frame() %>% as.tibble() %>%
  mutate(
    str_loc = paste0(seqnames, ":", start-100000, "-", start + 100000)
  )

d_flt <- distanceToNearest(fltGR) %>% as.data.frame() %>% as.tibble()
ggplot(d_flt, aes(x = distance)) +
  geom_histogram() + scale_x_log10()

d_fp <- distanceToNearest(fpGR) %>% as.data.frame() %>% as.tibble()
ggplot(d_fp, aes(x = distance)) +
  geom_histogram() + scale_x_log10()

