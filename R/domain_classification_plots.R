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

domain_to_GRB <- read_rds("results/domain_to_GRB.rds")

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
# count gene categories and GRB-TADs
#*******************************************************************************
TADs_mm10_hESC <- tidy_domain_groups %>% 
  filter(
    species == "mm10" | is.na(species),
    domain_type == "hESC"
    ) %>% 
  left_join(domain_to_GRB, by = c("domain_id", "domain_type"))
  
TADs_mm10_hESC_category <- TADs_mm10_hESC %>% 
  count(category) %>%
  mutate(
    percent = n / sum(n) * 100
  ) %>% 
  write_tsv("results/tidy_domain_groups.mm10_hESC.category.counts.tsv")

TADs_mm10_hESC_GRB_class <- TADs_mm10_hESC %>%
  count(GRB_class) %>%
  mutate(
    percent = n / sum(n) * 100
  ) %>% 
  write_tsv("results/tidy_domain_groups.mm10_hESC.GRB_class.counts.tsv")

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



