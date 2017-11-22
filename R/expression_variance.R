
#*******************************************************************************
# Analysis of gene expression variance from Breschi et al. 2016
#*******************************************************************************

require(tidyverse)
require(ggsignif)

#*******************************************************************************
# Parameters
#*******************************************************************************
exp_variance_url <- "http://public-docs.crg.es/rguigo/Papers/breschi_clustering/Supplementary_Tables.xlsx"
exp_variance_file <- "data/Breschi2016/Supplementary_Tables.xlsx"

outPrefix = "results/expression_variance.v01"
dir.create(dirname(outPrefix), showWarnings = FALSE)

#*******************************************************************************
# Read genes with expression correlation and grouping by TADs
#*******************************************************************************

# read genes and with categories from ortholog_expression_analysis.R script
expCorDF <- read_rds("results/ortholog_expression_correlation.rds")

# count gens per category
countDF <- expCorDF %>% 
  count(domain_type, category)

p <- ggplot(countDF, aes(x = category, y = n, fill = category)) +
  geom_bar(stat = "identity") + 
  geom_text(aes(label = n), vjust = "inward") +
  facet_grid(domain_type ~ .) + 
  theme_bw() + 
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste0(outPrefix, ".genes_by_category_and_domain_type.barplot.pdf"), w = 3, h = 3)

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
  left_join(expVarDF, by = c("human_gene_id" = "gene_id"))


#*******************************************************************************
# plot variance across species ----
#*******************************************************************************
p <- ggplot(expCorDF, aes(x = category, y = var.Species, fill = category)) + 
  geom_boxplot() + 
  geom_signif(comparisons = list(
    c("Conserved", "Rearranged")
  )) +
  facet_grid(. ~ domain_type) +
  theme_bw() + 
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(
  paste0(outPrefix, ".var.Species_by_category_and_domain_type.boxplot.pdf"), 
  w = 6, h = 6)

#*******************************************************************************
# plot percent of variance across species ----
#*******************************************************************************
p <- ggplot(expCorDF, aes(x = category, y = percent.Species, fill = category)) +
  geom_boxplot() + 
  geom_signif(comparisons = list(
    c("Conserved", "Rearranged")
  )) +
  facet_grid(. ~ domain_type) +
  theme_bw() + 
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(
  paste0(outPrefix, ".percent.Species_by_category_and_domain_type.boxplot.pdf"), 
  w = 6, h = 6)

#*******************************************************************************
# plot percent of variance across organse ----
#*******************************************************************************
p <- ggplot(expCorDF, aes(x = category, y = var.Organ, fill = category)) + 
  geom_boxplot() + 
  geom_signif(comparisons = list(
    c("Conserved", "Rearranged")
  )) +
  facet_grid(. ~ domain_type) +
  theme_bw() + 
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(
  paste0(outPrefix, ".var.Organ_by_category_and_domain_type.boxplot.pdf"), 
  w = 6, h = 6)


#*******************************************************************************
# plot grouping of genes by TAD category ----
#*******************************************************************************
classCountDF <- expCorDF %>% 
  group_by(domain_type, category, percent.class) %>% 
  summarize(n = n()) %>% 
  filter(!is.na(percent.class)) %>% 
  mutate(
    percent = n / sum(n) * 100
  )

p <- ggplot(classCountDF, aes(x = category, y = n, fill = category)) +
  geom_bar(stat = "identity") + 
  geom_text(aes(label = n), vjust = "inward") +
  facet_grid(domain_type ~ percent.class) + 
  theme_bw() + 
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste0(outPrefix, ".n_of_percent.class_by_category_and_domain_type.barplot.pdf"), w = 6, h = 6)

p <- ggplot(classCountDF, aes(x = category, y = percent, fill = percent.class)) +
  geom_bar(stat = "identity") + 
  geom_text(aes(label = round(percent, 2)), position = position_stack(vjust = 0.5)) +
  facet_grid(domain_type ~ .) + 
  theme_bw() + 
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(type = "qual", palette = 2)
ggsave(paste0(outPrefix, ".percent_of_percent.class_by_category_and_domain_type.barplot.pdf"), w = 3, h = 6)
  