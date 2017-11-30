require(tidyverse)
require(RColorBrewer)
require(ggsignif)

# parameters

# my_palette = c(brewer.pal(5, "Dark2")[2], brewer.pal(5, "Set1")[c(2, 1, 3)])
GENE_CATEGORY_COLORS = brewer.pal(5, "Set1")[c(2, 1, 3)]

#*******************************************************************************
# read table with genes, expression correlation, domains, and categories
#*******************************************************************************
genes_by_domain_categories <- read_rds("results/genes_by_domain_categories.rds")


#*******************************************************************************
# count categories per species and domain type
#*******************************************************************************
count_df <- tidy_domain_groups %>% 
  count(species, domain_type, conserved, rearranged, category) %>% 
  write_tsv("results/tidy_domain_groups.species_domain_rearraged_conserved.counts.tsv")

count_df <- tidy_domain_groups %>% 
  count(species, domain_type, category) %>% 
  write_tsv("results/tidy_domain_groups.species_domain_category.counts.tsv")

#*******************************************************************************
# Analyse correlation by TAD category for mouse in hESC
#*******************************************************************************

# filter for only mouse data and hESC TADs
data <- genes_by_domain_categories %>% 
  filter(species == "mm10" | is.na(species)) %>% 
  filter(domain_type == "hESC")

# summarize statistics on gene level
summaryDF <- data %>%
  group_by(gene_category) %>%
  summarise(
    n_TADs_with_gene = length(unique(domain_id)),
    n_genes = n(), 
    n_genes_with_cor = sum(!is.na(correlation)),
    cor_mean = mean(correlation, na.rm = TRUE),
    cor_median = median(correlation, na.rm = TRUE),
    cor_sd = sd(correlation, na.rm = TRUE)
    ) %>% 
  write_tsv("results/ortholog_expression.mm10.TAD_hESC.summaryDF.tsv")

plot_data <- data %>% 
  filter(
    !is.na(gene_category),
    !is.na(correlation)
    ) %>% 
  mutate(gene_category = factor(
    gene_category, 
    levels = c("Conserved", "Rearranged", "Outside"),
    labels = c("Conserved TAD", "Rearranged TAD", "Outside TAD"))
    )


p <- ggplot(plot_data, aes(x = gene_category, y = correlation, color = gene_category, fill = gene_category)) + 
  geom_violin(color = NA, adjust = 0.25, alpha = 0.3) +
  geom_boxplot(fill = "white", outlier.size = 0.3, width = 0.2) + 
  geom_signif(comparisons = list(
    c("Conserved TAD", "Rearranged TAD"),
    c("Conserved TAD", "Outside TAD")), 
    color = "black",
    map_signif_level = FALSE,
    test = wilcox.test,
    tip_length = 0,
    y_position = c(1.05, 1.2, 1.35, 1.5),
    step_increase = 0.13) +
  scale_color_manual(values = GENE_CATEGORY_COLORS) + 
  scale_fill_manual(values = GENE_CATEGORY_COLORS) +
  theme_bw() +
  theme(legend.position = "none", 
        # text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(-.5, 1.3) +
  labs(x = "", y = "Ortholog expression\ncorrelation")

ggsave("results/ortholog_expression.mm10.TAD_hESC.cor_by_gene_category.boxplot.pdf",
       h = 3, w = 3)

