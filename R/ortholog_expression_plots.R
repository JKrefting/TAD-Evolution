require(tidyverse)
require(RColorBrewer)
require(ggsignif)

# parameters

# my_palette = c(brewer.pal(5, "Dark2")[2], brewer.pal(5, "Set1")[c(2, 1, 3)])
GENE_CATEGORY_COLORS = brewer.pal(9, "Set1")[c(2, 1, 3, 9)]
TAD_CATEGORY_COLORS = brewer.pal(9, "Set1")[c(2, 1, 9)] # blue, red, gray
TAD_GRB_CLASS_COLORS = brewer.pal(9, "Set1")[c(5, 4, 9)] # orange, purple, gray
#*******************************************************************************
# read table with genes, expression correlation, domains, and categories
#*******************************************************************************
genes_by_domain_categories <- read_rds("results/genes_by_domain_categories.rds")

#*******************************************************************************
# Count genes with orthologs and expression data
#*******************************************************************************
geneCount <- genes_by_domain_categories %>% 
  distinct(ensembl_gene_id, .keep_all = TRUE) %>% 
  mutate(
    hasExp = !is.na(correlation),
    hasOrtholog = !is.na(mmusculus_homolog_ensembl_gene)
    ) %>% 
  count(hasExp, hasOrtholog) %>% 
  write_tsv("results/genes_numbers_by_orthologs_and_exp_data.tsv")

#*******************************************************************************
# Count number of genes by combinations of groups 
#*******************************************************************************
countTADs <- genes_by_domain_categories %>% 
  filter(
    is.na(species) | species == "mm10",
    domain_type %in% c("hESC"),
  ) %>%
  mutate(inTAD = !is.na(domain_id)) %>% 
  count(inTAD, category, GRB_class)
  
#*******************************************************************************
# Plot number of gnes in gene_category 
#*******************************************************************************

summaryDF <- genes_by_domain_categories %>% 
  filter(
    is.na(species) | species == "mm10",
    domain_type %in% c("hESC"),
    # !is.na(gene_category)
  ) %>%
  mutate(gene_category = factor(
    gene_category,
    levels = c("Conserved", "Rearranged", "Outside"),
    labels = c("Conserved TAD", "Rearranged TAD", "Outside TAD"))
  ) %>%
  group_by(species, domain_type, gene_category) %>%
  summarize(
    n = n(),
    n_TADs_with_gene = length(unique(domain_id)),
    n_genes_with_cor = sum(!is.na(correlation)),
    cor_mean = mean(correlation, na.rm = TRUE),
    cor_median = median(correlation, na.rm = TRUE),
    cor_sd = sd(correlation, na.rm = TRUE)
  ) %>% 
  write_tsv("results/domain_classes.mouse.cor_summary.gene_category.tsv")

p <- ggplot(summaryDF, aes(x = gene_category, y = n_genes_with_cor, fill = gene_category)) +
  geom_bar(stat = "identity", color = "black", lwd = .5) + 
  geom_text(aes(label = n_genes_with_cor), vjust = "inward") +
  # facet_grid(. ~ domain_type) + 
  theme_bw() + 
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = GENE_CATEGORY_COLORS) + 
  labs(x = "", y = "Genes")
ggsave("results/domain_classes.mouse_hESC.counts.gene_category.barplot.pdf", w = 3, h = 3)


#*******************************************************************************
# Compare correlation for gene within and ouside TADs
#*******************************************************************************

inTADsummaryDF <- genes_by_domain_categories %>% 
  filter(species == "mm10" | is.na(species)) %>% 
  filter(domain_type == "hESC") %>% 
  mutate(
    inTAD = factor(
      !is.na(domain_id), 
      levels = c(TRUE, FALSE),
      labels = c("Inside TAD", "Outside TAD")
    )) %>% 
  group_by(species, domain_type, inTAD) %>%
  summarize(
    n = n(),
    n_TADs_with_gene = length(unique(domain_id)),
    n_genes_with_cor = sum(!is.na(correlation)),
    cor_mean = mean(correlation, na.rm = TRUE),
    cor_median = median(correlation, na.rm = TRUE),
    cor_sd = sd(correlation, na.rm = TRUE)
  ) %>% 
  write_tsv("results/domain_classes.mouse.cor_summary.inTAD.tsv")


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
# Compare correlation for genes in conserved and rearranged TADs
#*******************************************************************************
conservedSummaryDF <- genes_by_domain_categories %>% 
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
    )) %>% 
  group_by(species, domain_type, gene_category) %>%
  summarize(
    n = n(),
    n_TADs_with_gene = length(unique(domain_id)),
    n_genes_with_cor = sum(!is.na(correlation)),
    cor_mean = mean(correlation, na.rm = TRUE),
    cor_median = median(correlation, na.rm = TRUE),
    cor_sd = sd(correlation, na.rm = TRUE)
  ) %>% 
  write_tsv("results/domain_classes.mouse.cor_summary.gene_category.tsv")


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
# count genes in GRB-TADs vs. nonGRB-TADs
#*******************************************************************************
GRBsummaryDF <- genes_by_domain_categories %>% 
  filter(
    is.na(species) | species == "mm10",
    domain_type %in% c("hESC"),
  ) %>%
  mutate(
    GRB_class = factor(
      GRB_class, 
      levels = c("GRB", "nonGRB", "screened"),
      labels = c("GRB-TAD", "Non-GRB-TAD", "other")
    )) %>% 
  group_by(species, domain_type, GRB_class) %>%
  summarize(
    n = n(),
    n_TADs_with_gene = length(unique(domain_id)),
    n_genes_with_cor = sum(!is.na(correlation)),
    cor_mean = mean(correlation, na.rm = TRUE),
    cor_median = median(correlation, na.rm = TRUE),
    cor_sd = sd(correlation, na.rm = TRUE)
  ) %>% 
  write_tsv("results/domain_classes.mouse.cor_summary.GRB_class.tsv")


p <- ggplot(GRBsummaryDF, aes(x = GRB_class, y = n_genes_with_cor, fill = GRB_class)) +
  geom_bar(stat = "identity", color = "black", lwd = .5) + 
  geom_text(aes(label = n_genes_with_cor), vjust = "inward") +
  # facet_grid(. ~ domain_type) + 
  theme_bw() + 
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = TAD_GRB_CLASS_COLORS) + 
  labs(x = "", y = "Genes")
ggsave("results/domain_classes.mouse_hESC.counts.GRB_class.barplot.pdf", w = 3, h = 3)

#*******************************************************************************
# Compare correlation for genes in GRB and non-GRB TADs
#*******************************************************************************

# filter for only mouse data and hESC TADs
grbTAD <- genes_by_domain_categories %>% 
  filter(
    species == "mm10" | is.na(species),
    domain_type == "hESC", 
    # filter for only gene within TADs
    !is.na(domain_id),
    GRB_class %in% c("GRB", "nonGRB")
  ) %>% 
  mutate(
    GRB_class = factor(
      GRB_class, 
      levels = c("GRB", "nonGRB", "screened"),
      labels = c("GRB-TAD", "Non-GRB-TAD", "other")
    ))

p <- ggplot(grbTAD, 
            aes(x = GRB_class, y = correlation,  color = GRB_class, 
                fill = GRB_class)) + 
  geom_violin(color = NA, adjust = 0.25, alpha = 0.3) +
  geom_boxplot(fill = "white", outlier.size = 0.3, width = 0.2) + 
  geom_signif(comparisons = list(
    c("GRB-TAD", "Non-GRB-TAD")), 
    color = "black",
    map_signif_level = FALSE,
    test = wilcox.test,
    tip_length = 0,
    y_position = 1.05,
    step_increase = 0.13) +
  scale_color_manual(values = TAD_GRB_CLASS_COLORS) +
  scale_fill_manual(values = TAD_GRB_CLASS_COLORS) +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(-.75, 1.2) +
  labs(x = "", y = "Ortholog expression\ncorrelation")

ggsave("results/ortholog_expression.mm10.TAD_hESC.cor_by_GRB-TAD_vs_nonGRB-TAD.boxplot.pdf",
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
