
#*******************************************************************************
# Analysis of gene expression variance from Breschi et al. 2016
#*******************************************************************************

require(tidyverse)
require(ggsignif)

# my_palette = c(brewer.pal(5, "Dark2")[2], brewer.pal(5, "Set1")[c(2, 1, 3)])
GENE_CATEGORY_COLORS = brewer.pal(9, "Set1")[c(2, 1, 3, 9)]
GENE_VAR_GROUP_COL = brewer.pal(8, "Dark2")[c(3, 1, 2)]
  
#*******************************************************************************
# Parameters
#*******************************************************************************
exp_variance_url <- "http://public-docs.crg.es/rguigo/Papers/breschi_clustering/Supplementary_Tables.xlsx"
exp_variance_file <- "data/Breschi2016/Supplementary_Tables.xlsx"

outPrefix = "results/expression_variance.v02"
dir.create(dirname(outPrefix), showWarnings = FALSE)

SPECIES_VARIANCE = c("galGal5", "")

#*******************************************************************************
# Read genes with expression correlation and grouping by TADs
#*******************************************************************************

# read genes and with categories from ortholog_expression_analysis.R script
# expCorDF <- read_rds("results/ortholog_expression_correlation.rds")
tidy_domain_groups <- read_rds("results/tidy_domain_groups.rds")
genes_to_domains <- read_rds("results/genes_to_domains.rds")
genes_by_domain_categories <- read_rds("results/genes_by_domain_categories.rds")

# filter out mESCs because they does not make sense on human genome
expCorDF <- genes_by_domain_categories %>% 
  filter(species == "mm10" | is.na(species)) %>% 
  mutate(gene_category = factor(
    gene_category,
    levels = c("Conserved", "Rearranged", "Outside"),
    labels = c("Conserved TAD", "Rearranged TAD", "Outside TAD"))
  )

# count gens per category
countDF <- expCorDF %>% 
  count(domain_type, gene_category)

p <- ggplot(countDF, aes(x = gene_category, y = n, fill = gene_category)) +
  geom_bar(stat = "identity") + 
  geom_text(aes(label = n), vjust = "inward") +
  facet_grid(domain_type ~ .) + 
  scale_fill_manual(values = GENE_CATEGORY_COLORS) +
  theme_bw() + 
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste0(outPrefix, ".genes_by_gene_category_and_domain_type.barplot.pdf"), w = 3, h = 3)

#*******************************************************************************
# Download and read data from Breschi et al. 2016
#*******************************************************************************

# create directeory and download file
dir.create(dirname(exp_variance_file), , showWarnings = FALSE)
download.file(exp_variance_url, exp_variance_file)

expVarDF <- readxl::read_excel(exp_variance_file, sheet = 2)

#*******************************************************************************
# Add expression variance estimates to expression correlation data set
#*******************************************************************************
expCorDF <- expCorDF %>% 
  left_join(expVarDF, by = c("ensembl_gene_id" = "gene_id"))


#*******************************************************************************
# plot variance across species ----
#*******************************************************************************
p <- ggplot(expCorDF, aes(x = gene_category, y = var.Species, fill = gene_category)) + 
  geom_boxplot() + 
  geom_signif(comparisons = list(
    c("Conserved TAD", "Rearranged TAD"),
    c("Conserved TAD", "Outside TAD")
  ), y = log10(c(110, 150))) +
  facet_grid(. ~ domain_type) +
  scale_fill_manual(values = GENE_CATEGORY_COLORS) +
  scale_y_log10() +
  theme_bw() + 
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(
  paste0(outPrefix, ".var.Species_by_gene_category_and_domain_type.boxplot.pdf"), 
  w = 6, h = 6)

#*******************************************************************************
# plot percent of variance across species ----
#*******************************************************************************
p <- ggplot(expCorDF, aes(x = gene_category, y = percent.Species, fill = gene_category)) +
  geom_boxplot() + 
  geom_signif(comparisons = list(
    c("Conserved TAD", "Rearranged TAD"),
    c("Conserved TAD", "Outside TAD")
  ), y = c(1, 1.1)) +
  facet_grid(. ~ domain_type) +
  scale_fill_manual(values = GENE_CATEGORY_COLORS) +
  theme_bw() + 
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(
  paste0(outPrefix, ".percent.Species_by_gene_category_and_domain_type.boxplot.pdf"), 
  w = 6, h = 6)

#*******************************************************************************
# plot  variance across organse ----
#*******************************************************************************
p <- ggplot(expCorDF, aes(x = gene_category, y = var.Organ, fill = gene_category)) + 
  geom_boxplot() + 
  geom_signif(comparisons = list(
    c("Conserved TAD", "Rearranged TAD"),
    c("Conserved TAD", "Outside TAD")
  ), y = log10(c(110, 150))) +
  facet_grid(. ~ domain_type) +
  scale_fill_manual(values = GENE_CATEGORY_COLORS) +
  scale_y_log10() +
  theme_bw() + 
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(
  paste0(outPrefix, ".var.Organ_by_gene_category_and_domain_type.boxplot.pdf"), 
  w = 6, h = 6)


#*******************************************************************************
# plot grouping of genes by gene_category ----
#*******************************************************************************
classCountDF <- expCorDF %>% 
  group_by(domain_type, gene_category, percent.class) %>% 
  summarize(n = n()) %>% 
  filter(!is.na(percent.class)) %>% 
  mutate(
    percent = n / sum(n) * 100
  )

p <- ggplot(classCountDF, aes(x = gene_category, y = n, fill = gene_category)) +
  geom_bar(stat = "identity") + 
  geom_text(aes(label = n), vjust = "inward") +
  facet_grid(domain_type ~ percent.class) + 
  scale_fill_manual(values = GENE_CATEGORY_COLORS) +
  theme_bw() + 
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste0(outPrefix, ".n_of_percent.class_by_gene_category_and_domain_type.barplot.pdf"), w = 6, h = 6)

p <- ggplot(classCountDF, aes(x = gene_category, y = percent, fill = percent.class)) +
  geom_bar(stat = "identity") + 
  geom_text(aes(label = round(percent, 2)), position = position_stack(vjust = 0.5)) +
  facet_grid(domain_type ~ .) + 
  theme_bw() + 
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = GENE_VAR_GROUP_COL) 
  # scale_fill_brewer(type = "qual", palette = 2)
ggsave(paste0(outPrefix, ".percent_of_percent.class_by_gene_category_and_domain_type.barplot.pdf"), w = 3, h = 6)
  
