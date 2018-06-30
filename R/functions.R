require(rtracklayer)

# Read in fill BED file
readFillFile <- function(species, seqinfo){
 fi_file <- paste0("data/fills/hg38.", species, ".fill.bed")
  return(import(fi_file, seqinfo = seqinfo))
}

# Read in breakpoint BED file
readBPFile <- function(species, threshold, seqinfo){
  bp_file <- paste0("data/breakpoints/hg38.", species, ".", 
                    as.character(format(threshold, scientific = FALSE)), ".bp.flt.flt_adj_fill.bed")
return(import(bp_file, seqinfo = seqinfo))
}

# ' Calculates the number of breakpoints that fall into each bin for each region
# ' 
# ' @param breakpoints GRanges with breakpoint coordinates
# ' @param region_bins GRangesList where every region is subdivided into n_bins
# ' @param n_bins Integer indicating the number of bins each region is subdivided into
# ' @return The hits per bin summed up over all regions
calcHitsPerBin <- function(breakpoints, region_bins, n_bins){
  # check which breakpoint falls into which bin
  all_bins <- unlist(region_bins)
  hits_list_all <- countOverlaps(all_bins, breakpoints)
  # add up vectors as rows in a matrix
  hits_matrix <- matrix(hits_list_all, ncol = n_bins, byrow = TRUE)
  # sum up over all regions
  hits_per_bin <- colSums(hits_matrix)
  
  return(hits_per_bin)
}

# ' Sample random breakpoints for each actual breakpoint per chromosome
# ' 
# ' @param breakpoints_per_chr Dataframe providing the number of actual breakpoints per chromosome
# ' @param region_bins hum_seqinfo Seqinfo object for the construction of random breakpoint GRanges 
# ' @return GRanges containing the sampled random breakpoints
sampleBreakpoints <- function(breakpoints_per_chr, hum_seqinfo){
  
  breakpoints_per_chr <- breakpoints_per_chr %>% 
    mutate(
      chrom_len = seqlengths(hum_seqinfo)[seqnames]
      )
  
  # sample the coordinates of random breakpoints
  starts <- map2(breakpoints_per_chr$chrom_len, breakpoints_per_chr$count,
                 sample.int)
  
  seqnames <- rep(breakpoints_per_chr$seqnames, breakpoints_per_chr$count)  
  
  rdm_breakpoints <- GRanges(seqnames, IRanges(unlist(starts), width = 1), seqinfo = hum_seqinfo)
  
  return(rdm_breakpoints)
}

# ' A vector of p values is transformed into a vector of asterisks representing significance levels
# ' 
# ' @param p_vals Integer (vector) containing p values
# ' @return String (vector) containing asterisks for the respective significance level
# none = p_val > 0.05; * = p_val <= 0.05; ** = p_val <= 0.01; *** = p_val <= 0.001 
asterisks <- function(p_vals){
  blank <- p_vals > 0.05 
  sig <- p_vals <= 0.05
  hsig <- p_vals <= 0.01
  hhsig <- p_vals <= 0.001
  p_vals[blank] <- ""
  p_vals[sig] <- "*"
  p_vals[hsig] <- "**"
  p_vals[hhsig] <- "***"
  
  return(p_vals)
}

# ' Calculate the standard error of the mean
# ' source: https://stackoverflow.com/questions/2676554/in-r-how-to-find-the-standard-error-of-the-mean
# ' @param x Numeric vector
# ' @return Numeric indicating standard error of x
std <- function(x) sd(x)/sqrt(length(x))

# ' Capitalize first letters of all words in a string
# ' source: https://stackoverflow.com/questions/6364783/capitalize-the-first-letter-of-both-words-in-a-two-word-string
# ' @param x String
# ' @return String
simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}

#' Get gene coordinates as GRanges for largest transcript per gene from ensembl
get_species_geneGR <- function(species_str, ensembl_url = "aug2017.archive.ensembl.org"){
  
  # Load gene ensembl species gene ids and orthologs
  species_ensembl <- useMart(host = ensembl_url, 
                           biomart = "ENSEMBL_MART_ENSEMBL", 
                           dataset = str_c(species_str, "_gene_ensembl"))
  
  # Get species transcription start sites
  species_gene_df <- as.tibble(getBM(attributes=c(
    'ensembl_gene_id', 'external_gene_name', 'chromosome_name', 'strand',
    'start_position', 'end_position'), 
    mart = species_ensembl)
  )

  # Generate GRanges
  speciesGeneGR <- GRanges(seqnames = species_gene_df$chromosome_name,
                        strand = ifelse(species_gene_df$strand == 1, "+", "-"),
                        ranges = IRanges(start = species_gene_df$start_position,
                                         end = species_gene_df$end_position),
                        geneID = species_gene_df$ensembl_gene_id,
                        gene_name = species_gene_df$external_gene_name
  )
  
  return(speciesGeneGR)
}



#' get one-to-one orthologs to all human genes from target species via ensembl
get_orthologs <- function(species_str, mart){
  
  # Attributes that are returned
  orth_attr = c("ensembl_gene_id",  
                str_c(species_str, "_homolog_ensembl_gene"), 
                str_c(species_str, "_homolog_orthology_type"), 
                str_c(species_str, "_homolog_orthology_confidence"))
  
  # Query orthologs by human gene ID
  orthologsDF = getBM(attributes = orth_attr, mart = mart) 
  
  names(orthologsDF) <- names(orthologsDF) %>% 
    str_replace(species_str, "species")
  
  # filter all human genes to have an "one2one" ortholog in mouse
  orthologs <- as.tibble(orthologsDF) %>% 
    filter(
      species_homolog_ensembl_gene != "",
      species_homolog_orthology_type == "ortholog_one2one"
    ) %>% 
    select(ensembl_gene_id, species_homolog_ensembl_gene)
  
  return(orthologs)  
}

#' get all adjacent ranges  
getAdjacentPairs <- function(geneGR){
  
  # get the next tss for each gene along the chromsome
  # using the precede() function from GenomicRanges package
  nextGene = precede(geneGR, ignore.strand = TRUE)
  
  firstGR <- geneGR[!is.na(nextGene)]
  nextGR <- geneGR[nextGene[!is.na(nextGene)]]
  
  gP = tibble(
    g1 = firstGR$geneID, 
    g2 = nextGR$geneID,
    dist = start(nextGR) - end(firstGR),
    strand = case_when(
      as.logical(strand(firstGR) == strand(nextGR)) ~ "same",
      as.logical(strand(firstGR) == "+" & strand(nextGR) == "-") ~ "convergent",
      as.logical(strand(firstGR) == "-" & strand(nextGR) == "+") ~ "divergent"
    )
  )
  
}

  
#' For all human adjacent genes get orthologs in target speceis and whether they
#' are syntenic.
#' 
get_syntenic_pairs <- function(species_str, assembly_str, geneGR, 
                               ensembl_url = "aug2017.archive.ensembl.org",
                               size_thresholds = c(10000, 100000, 1000000)) {
  
  mart <- useMart(host = ENSEMBL_URL,
                         biomart = "ENSEMBL_MART_ENSEMBL", 
                         dataset = "hsapiens_gene_ensembl")
  
  speciesGeneGR <- get_species_geneGR(species_str, ensembl_url)
  orthologs <- get_orthologs(species_str, mart)
  
  human_ortholog_GR <- geneGR[geneGR$geneID %in% orthologs$ensembl_gene_id]
  species_ortholog_GR <- speciesGeneGR[speciesGeneGR$geneID %in% orthologs$species_homolog_ensembl_gene]
  
  human_adjacent <- getAdjacentPairs(human_ortholog_GR)
  species_adjacent <- getAdjacentPairs(species_ortholog_GR)

  # douplicate species adjacent pairs to have A-B and B-A pairs
  species_adjacent_all <- bind_rows(
    species_adjacent,
    rename(species_adjacent, g1 = g2, g2 = g1)
  ) %>% 
    rename(dist_species = dist, 
         strand_species = strand) %>% 
    mutate(adjacent_species = TRUE)
  
  df <- human_adjacent %>% 
    # add ortholog to first gene
    left_join(orthologs, by = c("g1" = "ensembl_gene_id")) %>% 
    rename(g1_species = species_homolog_ensembl_gene) %>% 
    # add ortholog to second gene
    left_join(orthologs, by = c("g2" = "ensembl_gene_id")) %>% 
    rename(g2_species = species_homolog_ensembl_gene) %>% 
    # test if orthologs are adjacent
    left_join(species_adjacent_all, by = c("g1_species" = "g1", "g2_species" = "g2")) %>% 
    mutate(
      adjacent_species = !is.na(adjacent_species),
      same_strand = strand == strand_species,
      syntenic = adjacent_species & same_strand
    )
  
  # build GRanges for span of adjacent paris
  gr1 = geneGR[match(df$g1, geneGR$geneID)]
  gr2 = geneGR[match(df$g2, geneGR$geneID)]
  
  adjacentGR <- GRanges(
    seqnames = seqnames(gr1),
    IRanges(start = pmin(end(gr1), end(gr2)),
            width = df$dist),
    seqinfo = seqinfo(gr1)
  )
  
  # read breakpoint data and compute overlap with adjacent gene regions
  rearranged = size_thresholds %>%
    set_names(str_c("breakpoints_", size_thresholds %>% format(scientific = FALSE) %>% str_trim)) %>%
    
    # build path to breakpoint file
    map(~ paste0("data/breakpoints/hg38.", assembly_str, ".",
                 as.character(format(.x, scientific = FALSE)),  ".bp.flt.bed")) %>%
    # read breakopoints
    map(import.bed) %>% 
    
    # get rearraged pairs by overlap with adjacent region
    map(~ overlapsAny(adjacentGR, .x)) %>% 
    as.tibble()
  
  # add rearraged columns to df
  df <- df %>% 
    bind_cols(rearranged)
  
}
