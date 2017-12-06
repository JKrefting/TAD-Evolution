#*******************************************************************************
# Associate TADs to GRBs by grouping TADs in GRB-TADs, and Non-GRB-TADs
#*******************************************************************************
require(BSgenome.Hsapiens.UCSC.hg19)
require(rtracklayer)
require(tidyverse)

# Read metadata for analysis
DOMAINS <- read_tsv("domains_meta.tsv")

# Load human seqinfo
genome <- BSgenome.Hsapiens.UCSC.hg19
hum_seqinfo <- seqinfo(genome)


#' Screen TADs for GRB overlp and report for each TAD one of three classes.
#'
#' This should be equivalent (but faster) than the screen.grb function from
#' Harmston 2016 et al. which can be found here:
#' https://github.com/ComputationalRegulatoryGenomicsICL/tad_cnes_harmston2017/blob/master/R/classify_dev_nondev.Rmd
#' 
tidy_screen_grb = function(tadGR, grbGR, top = 0.8, bottom = 0.2){
  
  # check that GRBs are non-overlapping, otherwise implementation does not work
  stopifnot(length(findOverlaps(grbGR, grbGR)) == length(grbGR))
  
  hits <- findOverlaps(tadGR, grbGR)
  
  hitsDF <- hits %>% 
    as.data.frame() %>% 
    as.tibble() %>% 
    mutate(
      ol_bp = width(pintersect(tadGR[queryHits], grbGR[subjectHits])),
      ol_frac = ol_bp / width(tadGR[queryHits])
    ) %>% 
    # since one TAD might overlap several GRB, summarize overlap bp and frac
    group_by(queryHits) %>% 
    summarize(
      overlap_count = n(),
      TAD_bp = sum(ol_bp),
      TAD_frac = sum(ol_frac)
    ) %>% 
    ungroup()
  
  # define class based on 
  GRB_class <- tibble(
    id = 1:length(tadGR)
  ) %>% 
    left_join(hitsDF, by = c("id" = "queryHits")) %>% 
    mutate(
      class = case_when(
        is.na(overlap_count) ~ "nonGRB",
        TAD_frac > top & overlap_count == 1 ~ "GRB",
        TAD_frac > bottom ~ "screened",
        TRUE ~ "nonGRB")
    )
  return(pull(GRB_class, class))
}


#*******************************************************************************
# Read GRB as GRanges ----
#*******************************************************************************
grb_path = DOMAINS %>%
  filter(domain_type == "GRB") %>% 
  pull(domain_path)

grbGR <- import.bed(grb_path, seqinfo = hum_seqinfo)

#*******************************************************************************
# calssify each TAD by GRB overlap ----
#*******************************************************************************

domain_to_GRB <- DOMAINS %>% 
  filter(domain_type != "GRB") %>% 
  mutate(
    tadGR = map(domain_path, import.bed, seqinfo = hum_seqinfo),
    GRB_class = map(tadGR, tidy_screen_grb, grbGR, top = 0.8, bottom = 0.2),
    domain_id = map(tadGR, ~1:length(.x))
  ) %>% 
  select(domain_id, domain_type, GRB_class) %>% 
  unnest(domain_id, GRB_class)

# write to output file
write_rds(domain_to_GRB, "results/domain_to_GRB.rds")
write_tsv(domain_to_GRB, "results/domain_to_GRB.tsv")

