require(tidyverse)
require(RColorBrewer)
require(ggsignif)

# data <- read_rds("results/ortholog_expression_correlation_old.rds")
# data <- read_rds("results/ortholog_expression_correlation.rds")
genes_by_domain_categories <- read_rds("results/genes_by_domain_categories.rds")

data <- genes_by_domain_categories %>% 
  filter(species == "mm10" | is.na(species)) %>% 
  filter(domain_type == "hESC")


data_stats <- data %>%
  group_by(domain_type, gene_category) %>%
  summarise(count = n(), 
            n_tads = length(unique(domain_id)))

plot_data <- data %>% 
  filter(
    domain_type == "hESC", 
    !is.na(gene_category), 
    !is.na(correlation)
    )

# my_palette = c(brewer.pal(5, "Dark2")[2], brewer.pal(5, "Set1")[c(2, 1, 3)])
my_palette = brewer.pal(5, "Set1")[c(2, 1, 3)]
# Plot correlations seperately depending on class
ggplot(plot_data, aes(x = gene_category, y = correlation)) + 
    geom_violin(aes(color = gene_category, fill = gene_category), adjust = 0.25, alpha = 0.3) +
    geom_boxplot(aes(color = gene_category), outlier.size = 0.3, width = 0.2) + 
    geom_signif(comparisons = list(
      c("Conserved", "Rearranged"),
      c("Conserved", "Outside")), 
      map_signif_level = FALSE,
      test = wilcox.test,
      textsize = 3,
      tip_length = 0,
      y_position = c(1.05, 1.2, 1.35, 1.5),
      step_increase = 0.13) + 
    scale_color_manual("", labels = c(
      "Conserved TAD", "Rearranged TAD", "Outside TAD"),
      values = my_palette
    ) +
    scale_fill_manual("", labels = c(
      "Conserved TAD", "Rearranged TAD", "Outside TAD"),
      values = my_palette
    ) +
    theme_minimal() + 
    theme(plot.title = element_text(size=12, face = "bold", colour = "black"),
          legend.text = element_text(size = 10, face = "bold"),
          legend.position = "bottom",
          # legend.position = c(0.5, -0.03),
          legend.direction = "horizontal",
          legend.key.width = unit(0.3, "cm"),
          legend.key.height = unit(0.3, "cm"),
          strip.background = element_blank(),
          strip.text = element_text(size = 10, face = "bold"),
          axis.title.x = element_text(face="bold", size=10),
          axis.title.y = element_text(face="bold", size=10),
          axis.text.x = element_blank(), 
          # axis.text.x = element_text(face="bold", size=10),
          axis.text.y = element_text(face="bold", size=10),
          axis.ticks.x = element_blank(),
          panel.grid.major = element_blank(),
          # panel.grid.minor = element_blank(),
          plot.margin=unit(c(0.1,0.5,0,0.5), "cm")) +
    ylim(c(-1, 1.35)) + 
    xlab("") + 
    ylab("Ortholog expression correlation")
