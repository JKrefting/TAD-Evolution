require(tidyverse)
require(RColorBrewer)
require(ggsignif)

# parameters

# my_palette = c(brewer.pal(5, "Dark2")[2], brewer.pal(5, "Set1")[c(2, 1, 3)])
# GENE_CATEGORY_COLORS = brewer.pal(9, "Set1")[c(2, 1, 3, 9)]
TAD_CATEGORY_COLORS = brewer.pal(9, "Set1")[c(2, 1, 9)]

#*******************************************************************************
# read table with genes, expression correlation, domains, and categories
#*******************************************************************************
tidy_domain_groups <- read_rds("results/tidy_domain_groups.rds")


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
# Plot categorization of domains in classes
#*******************************************************************************
TADcounts <- tidy_domain_groups %>% 
  mutate(
    category = replace(category, is.na(category), "Other"),
    category = factor( category, 
                       levels = c("Conserved", "Rearranged", "Other"))
  ) %>% 
  filter(
    species == "mm10",
    domain_type %in% c("hESC", "GM12878")
  ) %>% 
  count(domain_type, category) %>% 
  write_tsv("results/domain_classes.mouse.counts.tsv")

p <- ggplot(TADcounts, aes(x = category, y = n, fill = category)) +
  geom_bar(stat = "identity") + 
  geom_text(aes(label = n), vjust = "inward") +
  facet_grid(. ~ domain_type) + 
  theme_bw() + 
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = TAD_CATEGORY_COLORS) + 
  labs(x = "", y = "Number of TADs")
ggsave("results/domain_classes.mouse.counts.barplot.pdf", w = 6, h = 3)



