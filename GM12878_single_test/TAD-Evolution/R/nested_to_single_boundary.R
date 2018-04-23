require(tidyverse)
require(BSgenome.Hsapiens.UCSC.hg19)

source("R/functions.R")


hg19_seqinfo <- seqinfo(BSgenome.Hsapiens.UCSC.hg19)
tad_path <- "data/TADs/GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.bed"


tads <- import.bed(tad_path, seqinfo = hg19_seqinfo)

# extract boundaries
bdies <- c(
  GRanges(seqnames = seqnames(tads), 
                 ranges = IRanges(start = start(tads), width = 1), 
                 strand = strand(tads)),
  GRanges(seqnames = seqnames(tads), 
          ranges = IRanges(start = end(tads), width = 1), 
          strand = strand(tads))
  )

# 'assign' start and end boundaries to their corresponding TADs
start_hits <- hits <- findOverlaps(bdies, tads, type = "start")
end_hits <- hits <- findOverlaps(bdies, tads, type = "end")


# find boundaries which are within a TAD (includes the start and end boundaries 
# as hits)
all_hits <- findOverlaps(bdies, tads, type = "within")

# filter out start and end hits from all hits to identify boundaries which are 
# located in at least one other TAD
tmp_hits <- setdiff(all_hits, start_hits)
final_hits <- setdiff(tmp_hits, end_hits)

# get boundaries overlapping TADs (= nested)
bdy_idx <- unique(queryHits(final_hits))
nested_bdies <- bdies[bdy_idx]

# get boundaries NOT overlapping another TAD (= single)
bdy_idx <- setdiff(seq(1:length(bdies)), queryHits(final_hits))
single_bdies <- bdies[bdy_idx]

# save as bed files
export.bed(nested_bdies,
           "data/TADs/GSE63525_GM12878_primary+replicate_Arrowhead_domainlist_nested_boundaries.bed",
           format = "bed")

export.bed(single_bdies,
           "data/TADs/GSE63525_GM12878_primary+replicate_Arrowhead_domainlist_single_boundaries.bed",
           format = "bed")
