################################################################################
# Get gene expression data in mouse and human from FANTOM5 via EBI Expression Atlas
################################################################################

require(tidyverse)
require(stringr)

#-------------------------------------------------------------------------------
# input parameters and data sources
#-------------------------------------------------------------------------------

human_exp_prefix <- "data/ExpressionAtlas/human_E-MTAB-3358"
mouse_exp_prefix <- "data/ExpressionAtlas/mouse_E-MTAB-3579"

expression_suffix <- "_baseline_expression.TPM.tsv"
metadata_suffix <- "_experimental_design.tsv"

# define file path prefix for output files and create folder
output_prefix = "data/ExpressionAtlas/formated"
dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)

#-------------------------------------------------------------------------------
# Read data
#-------------------------------------------------------------------------------

replaceNAwithZero <- function(x){
  x[is.na(x)] <- 0
  return(x)
}

# read human data
human_exp <- str_c(human_exp_prefix, expression_suffix) %>% 
  read_tsv(comment = "#", col_types = cols(
    .default = col_double(),
    `Gene ID` = col_character(),
    `Gene Name` = col_character()
  )) %>% 
  # treat missing expression values as zeros
  mutate_at(replaceNAwithZero, .vars = 3:ncol(.))


human_design <- str_c(human_exp_prefix, metadata_suffix) %>% 
  read_tsv(col_types = cols(.default = col_character()))

# read mouse data
mouse_exp <- str_c(mouse_exp_prefix, expression_suffix) %>% 
  read_tsv(comment = "#", col_types = cols(
    .default = col_double(),
    `Gene ID` = col_character(),
    `Gene Name` = col_character()
  )) %>% 
  # treat missing expression values as zeros
  mutate_at(replaceNAwithZero, .vars = 3:ncol(.))

mouse_design <- str_c(mouse_exp_prefix, metadata_suffix) %>% 
  read_tsv(col_types = cols(.default = col_character()))

#-------------------------------------------------------------------------------
# Map samples between human and mouse
#-------------------------------------------------------------------------------

# get mapping from human tissue name to UBERON term ID
human_design_sub <- human_design %>%
  select(name = `Sample Characteristic[organism part]`, 
         term_url = `Sample Characteristic Ontology Term[organism part]`) %>% 
  distinct()


# get metadat of each condition from column name and design table
human_meta <- tibble(
    human_col_name = names(human_exp)[-c(1,2)],
    human_col_idx = 3:ncol(human_exp),
    human_tissue = str_split_fixed(human_col_name, ", ", n = 2)[,1],
    human_dev_stage = str_split_fixed(human_col_name, ", ", n = 2)[,2]
  ) %>% 
  left_join(human_design_sub, by = c("human_tissue" = "name"))

# get mapping from mouse tissue name to UBERON term ID
mouse_design_sub <- mouse_design %>%
  select(name_mouse = `Sample Characteristic[organism part]`, 
         term_url = `Sample Characteristic Ontology Term[organism part]`) %>% 
  distinct()

# get metadat of each condition from column name and design table
mouse_meta <- tibble(
    mouse_col_name = names(mouse_exp)[-c(1,2)],
    mouse_col_idx = 3:ncol(mouse_exp),
    mouse_tissue = str_split_fixed(mouse_col_name, ", ", n = 2)[,1],
    mouse_dev_stage = str_split_fixed(mouse_col_name, ", ", n = 2)[,2]
  ) %>% 
  left_join(mouse_design_sub, by = c("mouse_tissue" = "name_mouse"))

#-------------------------------------------------------------------------------
# match conditions across species
#-------------------------------------------------------------------------------

cond_match <- inner_join(
  human_meta,
  mouse_meta,
  by = c("term_url", "human_dev_stage" = "mouse_dev_stage")
)

#-------------------------------------------------------------------------------
# select subset of matching columns from full expression data sets
#-------------------------------------------------------------------------------

human_exp_match <- human_exp %>%
  select(1, 2, cond_match$human_col_idx)

mouse_exp_match <- mouse_exp %>%
  select(1, 2, cond_match$mouse_col_idx)

#-------------------------------------------------------------------------------
# Write tables to output files
#-------------------------------------------------------------------------------

write_tsv(cond_match, str_c(output_prefix, ".condition_matching.mouse_human.tsv"))

write_tsv(human_exp_match, str_c(output_prefix, ".human_exp.matching_subset.tsv"))
write_tsv(mouse_exp_match, str_c(output_prefix, ".mouse_exp.matching_subset.tsv"))


#-------------------------------------------------------------------------------
# Filter human and mouse separatley for only adult tissue
#-------------------------------------------------------------------------------
human_adult_cols <- human_meta %>% 
  filter(human_dev_stage == "adult")

human_exp_adult <- human_exp %>% 
  select(1, 2, human_adult_cols$human_col_idx)

write_tsv(human_exp_adult, str_c(output_prefix, ".human_exp.adult_subset.tsv"))

# for mouse expression data
mouse_adult_cols <- mouse_meta %>% 
  filter(mouse_dev_stage == "adult")

mouse_exp_adult <- mouse_exp %>% 
  select(1, 2, mouse_adult_cols$mouse_col_idx)

write_tsv(mouse_exp_adult, str_c(output_prefix, ".mouse_exp.adult_subset.tsv"))

