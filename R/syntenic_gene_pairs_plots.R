#'##############################################################################
#' Plot results of syntenic gene pair analysis 
#'##############################################################################

require(tidyverse)
require(RColorBrewer)


# Read metadata for analysis
SPECIES <- read_tsv("species_meta.tsv", col_types = "ccc")

# define colors for species
COL_SPECIES = brewer.pal(12, "Set3") %>% 
  set_names(SPECIES$trivial_name)

COL_SYNTENIC <- RColorBrewer::brewer.pal(9, "Paired")[2:1]
COL_REARRANGED <- RColorBrewer::brewer.pal(9, "Paired")[5:6]


syntenic_performance_DF <- read_rds("results/syntenic_performance_DF.rds") %>% 
  mutate(trivial_name = factor(trivial_name, rev(SPECIES$trivial_name))) %>% 
  filter(!is.na(FDR))


# plot FDR ---------------------------------------------------------------------
p <- ggplot(syntenic_performance_DF, aes(x = trivial_name, 
                                         y = FDR, 
                                         fill = trivial_name,
                                         color = trivial_name)) + 
  geom_bar(aes(alpha = as.character(size_threshold / 1000)), stat = "identity", position = position_dodge(.9), color = "black") +
  geom_text(aes(label = round(FDR, 2), group = as.character(size_threshold)), position = position_dodge(.9), hjust = "inward", color = "black") +
  scale_color_manual(values = COL_SPECIES, guide = FALSE) +
  scale_fill_manual(values = COL_SPECIES, guide = FALSE) +
  scale_alpha_discrete(range = c(0.33, 1)) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "bottom") + 
  labs(x = "", y = "Breakpoint false discovery rate") + 
  guides(alpha = guide_legend(title = "Size threshold (kb)"))

ggsave("results/syntenic_genes.performance.FDR_by_species.pdf", w = 6, h = 6)


# plot FPR ---------------------------------------------------------------------
p <- ggplot(syntenic_performance_DF, aes(x = trivial_name, 
                                         y = FPR * 100, 
                                         fill = trivial_name,
                                         color = trivial_name)) + 
  geom_bar(aes(alpha = as.character(size_threshold / 1000)), stat = "identity", position = position_dodge(.9), color = "black") +
  geom_text(aes(label = round(FPR*100, 2), group = as.character(size_threshold)), position = position_dodge(.9), hjust = "inward", color = "black") +
  scale_color_manual(values = COL_SPECIES, guide = FALSE) +
  scale_fill_manual(values = COL_SPECIES, guide = FALSE) +
  scale_alpha_discrete(range = c(0.33, 1)) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "bottom") + 
  labs(x = "", y = "Breakpoint false positive rate (%)") + 
  guides(alpha = guide_legend(title = "Size threshold (kb)"))

ggsave("results/syntenic_genes.performance.FPR_by_species.pdf", w = 6, h = 6)
p
# plot PPV ---------------------------------------------------------------------
p <- ggplot(syntenic_performance_DF, aes(x = trivial_name, 
                                         y = PPV, 
                                         fill = trivial_name,
                                         color = trivial_name)) + 
  geom_bar(aes(alpha = as.character(size_threshold / 1000)), stat = "identity", position = position_dodge(.9), color = "black") +
  geom_text(aes(label = round(PPV, 2), group = as.character(size_threshold)), position = position_dodge(.9), hjust = "inward", color = "black") +
  scale_color_manual(values = COL_SPECIES, guide = FALSE) +
  scale_fill_manual(values = COL_SPECIES, guide = FALSE) +
  scale_alpha_discrete(range = c(0.33, 1)) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "bottom") + 
  labs(x = "", y = "Positive predictive value (PPV)") + 
  guides(alpha = guide_legend(title = "Size threshold (kb)"))
ggsave("results/syntenic_genes.performance.PPV_by_species.pdf", w = 6, h = 6)


mean_performance <- syntenic_performance_DF %>% 
  summarize_at(
    .vars = c("FDR", "FPR", "PPV", "sensitivity", "specificity"), 
    .funs = c("mean", "median")) %>% 
  write_tsv("results/syntenic_genes.mean_performance.tsv")


gene_pair_stats <- read_tsv("results/syntenic_genes.gene_pair_stats.tsv")


# RColorBrewer::display.brewer.all()
p <- ggplot(gene_pair_stats, aes(x = species_name, y = n, 
                            fill = syntenic,
                            color = rearranged)) +
  geom_bar(stat = "identity") +
  facet_grid(. ~ size_threshold, scales = "free_x") +
  coord_flip() +
  theme_bw() + 
  scale_color_manual(values = COL_REARRANGED) +
  scale_fill_manual(values = COL_SYNTENIC) +
  theme(legend.position = "bottom")

