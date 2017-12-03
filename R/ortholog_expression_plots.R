require(tidyverse)
require(RColorBrewer)
require(ggsignif)

# parameters

# my_palette = c(brewer.pal(5, "Dark2")[2], brewer.pal(5, "Set1")[c(2, 1, 3)])
GENE_CATEGORY_COLORS = brewer.pal(9, "Set1")[c(2, 1, 3, 9)]
TAD_CATEGORY_COLORS = brewer.pal(9, "Set1")[c(2, 1, 9)]

#*******************************************************************************
# read table with genes, expression correlation, domains, and categories
#*******************************************************************************
genes_by_domain_categories <- read_rds("results/genes_by_domain_categories.rds")

#*******************************************************************************
# Plot number of gnes in gene_category 
#*******************************************************************************
geneCounts <- genes_by_domain_categories %>% 
  filter(
    is.na(species) | species == "mm10",
    domain_type %in% c("hESC", "GM12878"),
    # !is.na(gene_category)
  ) %>%
  mutate(gene_category = factor(
    gene_category,
    levels = c("Conserved", "Rearranged", "Outside"),
    labels = c("Conserved TAD", "Rearranged TAD", "Outside TAD"))
  ) %>%
  count(species, domain_type, gene_category) %>% 
  write_tsv("results/domain_classes.mouse.counts.gene_category.tsv")

p <- ggplot(geneCounts, aes(x = gene_category, y = n, fill = gene_category)) +
  geom_bar(stat = "identity", color = "black", lwd = .5) + 
  geom_text(aes(label = n), vjust = "inward") +
  facet_grid(. ~ domain_type) + 
  theme_bw() + 
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = GENE_CATEGORY_COLORS) + 
  labs(x = "", y = "Number of genes")
ggsave("results/domain_classes.mouse.counts.gene_category.barplot.pdf", w = 6, h = 3)


#*******************************************************************************
# Compare correlation for gene within and ouside TADs
#*******************************************************************************
# filter for only mouse data and hESC TADs
inOutTAD <- genes_by_domain_categories %>% 
  filter(species == "mm10" | is.na(species)) %>% 
  filter(domain_type == "hESC") %>% 
  mutate(
    inTAD = factor(
      !is.na(domain_id), 
      levels = c(TRUE, FALSE),
      labels = c("Inside TAD", "Outside TAD")
      ))
p <- ggplot(inOutTAD, aes(x = inTAD, y = correlation, color = inTAD, fill = inTAD)) + 
  geom_violin(color = NA, adjust = 0.25, alpha = 0.3) +
  geom_boxplot(fill = "white", outlier.size = 0.3, width = 0.2) + 
  geom_signif(comparisons = list(
    c("Inside TAD", "Outside TAD")), 
    color = "black",
    map_signif_level = FALSE,
    test = wilcox.test,
    tip_length = 0,
    y_position = 1.05,
    step_increase = 0.13) +
  scale_color_manual(values = GENE_CATEGORY_COLORS[3:4]) +
  scale_fill_manual(values = GENE_CATEGORY_COLORS[3:4]) +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(-.75, 1.2) +
  labs(x = "", y = "Ortholog expression\ncorrelation")
ggsave("results/ortholog_expression.mm10.TAD_hESC.cor_by_inTAD.boxplot.pdf",
       h = 3, w = 1.5)


#*******************************************************************************
# Compare correlation for gene within and ouside TADs
#*******************************************************************************
# filter for only mouse data and hESC TADs
conservedRearragedTAD <- genes_by_domain_categories %>% 
  filter(
    species == "mm10" | is.na(species),
    domain_type == "hESC", 
    # filter for only gene within TADs
    !is.na(domain_id),
    !is.na(gene_category)
  ) %>% 
  mutate(
    gene_category = factor(
      gene_category, 
      levels = c("Conserved", "Rearranged"),
      labels = c("Conserved TAD", "Rearranged TAD")
  ))
  
p <- ggplot(conservedRearragedTAD, 
            aes(x = gene_category, y = correlation,  color = gene_category, 
                fill = gene_category)) + 
  geom_violin(color = NA, adjust = 0.25, alpha = 0.3) +
  geom_boxplot(fill = "white", outlier.size = 0.3, width = 0.2) + 
  geom_signif(comparisons = list(
    c("Conserved TAD", "Rearranged TAD")), 
    color = "black",
    map_signif_level = FALSE,
    test = wilcox.test,
    tip_length = 0,
    y_position = 1.05,
    step_increase = 0.13) +
  scale_color_manual(values = GENE_CATEGORY_COLORS[1:2]) +
  scale_fill_manual(values = GENE_CATEGORY_COLORS[1:2]) +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(-.75, 1.2) +
  labs(x = "", y = "Ortholog expression\ncorrelation")
ggsave("results/ortholog_expression.mm10.TAD_hESC.cor_by_rearraged_vs_conserved.boxplot.pdf",
       h = 3, w = 1.5)


#*******************************************************************************
# Analyse correlation by TAD category for mouse in hESC
#*******************************************************************************

# filter for only mouse data and hESC TADs
mouse_hESC_data <- genes_by_domain_categories %>% 
  filter(species == "mm10" | is.na(species)) %>% 
  filter(domain_type == "hESC") %>% 
  mutate(
    gene_category = factor(
      gene_category, 
      levels = c("Conserved", "Rearranged", "Outside"),
      labels = c("Conserved TAD", "Rearranged TAD", "Outside TAD"))
  )
    

# summarize statistics on gene level
summaryDF <- mouse_hESC_data %>%
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

plot_data <- mouse_hESC_data %>% 
  filter(
    !is.na(gene_category),
    !is.na(correlation)
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

#*******************************************************************************
# Analyse number of highly correlated genes
#*******************************************************************************
COR_TH = 0.75
# summarize statistics on gene level
countHighCor <- mouse_hESC_data %>%
  filter(
    !is.na(correlation),
    !is.na(gene_category)
    ) %>% 
  group_by(gene_category) %>%
  summarise(
    n_genes = n(), 
    high_cor = sum(correlation >= COR_TH),
  ) %>% 
  mutate(
    percent_highCor = high_cor / n_genes * 100
  )

p <- ggplot(countHighCor, aes(x = gene_category, y = percent_highCor, fill = gene_category)) +
    geom_bar(stat = "identity", color = "black") +
    geom_text(aes(label = round(percent_highCor, 2), vjust = 1.2)) +
      scale_color_manual(values = GENE_CATEGORY_COLORS) + 
      scale_fill_manual(values = GENE_CATEGORY_COLORS) +
      theme_bw() +
      theme(legend.position = "none", 
            # text = element_text(size = 20),
            axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = paste0("Gens with high (R>", COR_TH, ")\n expression correaltion [%]"))

ggsave(paste0("results/ortholog_expression.mm10.TAD_hESC.cor_geq", COR_TH, "_by_gene_category.barplot.pdf"),
       h = 3, w = 3)
