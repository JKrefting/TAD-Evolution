
# =============================================================================================================================
# Plots of data after preprocessing: fills and breakpoints extracted from whole-genome alignemnts (net files)
# =============================================================================================================================

require(TxDb.Hsapiens.UCSC.hg38.knownGene)
source("R/functions.R")

require(tidyverse)
require(ggplot2)
require(RColorBrewer)


# Read metadata for analysis
SPECIES <- read_tsv("species_meta.tsv", col_types = "ccc")
METADATA <- read_tsv("metadata.tsv", col_types = "iiiiii")

# minimum size threasholds for chains (syntenic blocks) to be considered rearrangement blocks
THRESHOLDS <- METADATA %>% pull(min_size_threshold)

# define colors for species
COL_SPECIES = brewer.pal(12, "Set3")

# Load human seqinfo
hum_seqinfo <- seqinfo(TxDb.Hsapiens.UCSC.hg38.knownGene)

# -------------------------------------------------------------
# Fills
# -------------------------------------------------------------

# prepare data for plotting
# gather all date in one tibble

fill_data <- SPECIES %>% 
  select(genome_assembly, species = trivial_name) %>% 
  mutate(
    fills = map(genome_assembly, readFillFile, seqinfo = hum_seqinfo),
    size = map(fills, width),
    type = map(fills, function(gr) mcols(gr)$name)
  ) %>% 
  select(-fills, -genome_assembly) %>% 
  unnest(size, type) %>% 
  # order species as factor and rename types
  mutate(
    species = factor(species, levels = SPECIES$trivial_name),
    type = case_when(
      type == "top" ~ "top-level",
      type == "nonSyn" ~ "non-syntenic",
      type == "syn" ~ "syn / inv",
      type == "inv" ~ "syn / inv"
    )
  )

# store data
write_rds(fill_data, "results/fill_data.rds")
write_tsv(fill_data, "results/fill_data.tsv")

# -------------------------------------------------------------
# Plot number of fills per species / type
# -------------------------------------------------------------

# ignore syn / inv fills
fill_data <- fill_data %>% 
  filter(type %in% c("top-level", "non-syntenic")) %>% 
  mutate(
    type = factor(type, levels = c("top-level", "non-syntenic"))
  )

plot_df <- fill_data %>%
  group_by(species, type) %>%
  summarise(n = n())

ggplot(plot_df, aes(x = species, y = n,fill = species)) + 
  geom_bar(stat = "identity",
           color = "black") + 
  geom_text(aes(label = n), 
            # hjust = -0.2, 
            vjust = -0.3,
            size = 3,
            angle = 0) +
  facet_grid(type ~ .) +
  scale_y_continuous(name = "Number of fills", 
                     limits = c(0, 150000),
                 breaks = c(0, 25000, 50000, 75000, 100000, 125000, 150000),
             labels = c("0", "25000", "", "75000", "", "125000", "")) +
  theme_bw() + 
  scale_fill_brewer("", palette = "Set3") +
  theme(plot.title = element_text(size = 14, face = "bold", colour = "black"),
        strip.background = element_blank(),
        strip.text.y = element_text(size = 10, face = "bold"),
        legend.title = element_text(size = 7, face = "bold"), 
        legend.text = element_text(size = 7, face = "bold"),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = 10),
        axis.text.x = element_text(face = "bold", size = 10, angle = 45, 
                                   vjust = 1, hjust = 1), 
        axis.text.y = element_text(face = "bold", size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  )

ggsave("results/fill_numbers.pdf", w = 6, h = 4)


# Horizontal plot

plotDF <- plot_df %>%
  ungroup() %>% 
  mutate(
    species = factor(species, levels = rev(levels(species)))
  )

ggplot(plotDF, aes(x = species, y = n, fill = species)) + 
  geom_bar(stat = "identity", color = "black") +
  coord_flip() +
  geom_text(aes(label = n), hjust = "inward") +
  facet_grid(. ~ type) +
  scale_y_continuous(name = "Number of fills",
                     limits = c(0, 150000),
                     breaks = c(0, 25000, 50000, 75000, 100000, 125000, 150000),
                     labels = c("0", "25000", "", "75000", "", "125000", "")) +
  theme_bw() + 
  scale_fill_manual(values = rev(COL_SPECIES)) +
  # scale_fill_brewer("", palette = "Set3") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_blank()
  ) +
  labs(x = "", y = "")

ggsave("results/fill_numbers_horiz.pdf", w = 3, h = 6)



# -------------------------------------------------------------
# Depict distributions of fill sizes
# -------------------------------------------------------------
# ignore syn / inv fills
plot_df <- filter(fill_data, type != "syn / inv")
plot_df$type <- ordered(factor(plot_df$type, levels = c("top-level", "non-syntenic")))

ggplot(plot_df) + 
  geom_violin(aes(x = species, y = size, fill = species), adjust = 0.5) +
  stat_summary(aes(x = species, y = size), fun.y = "median", geom = "point", size = 2, color = "black") +
  scale_y_log10(name = "Fill sizes [bp]", breaks = scales::trans_breaks("log10", function(x) 10^x),
                limit = c(100, 1e+9), labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  facet_grid(type ~ .) +
  theme_bw() + 
  scale_fill_brewer("", palette = "Set3") +
  theme(plot.title = element_text(size = 14, face = "bold", colour = "black"),
        strip.background = element_blank(),
        strip.text.y = element_text(size = 10, face = "bold"),
        legend.title = element_text(size = 7, face = "bold"), legend.text = element_text(size = 7, face = "bold"),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = 10),
        axis.text.x = element_text(face = "bold", size = 10, angle = 45, vjust = 1, hjust = 1), 
        axis.text.y = element_text(face = "bold", size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  )


ggsave("results/fill_size_distributions.pdf", w = 6, h = 4)

# Horizontal plot of fill sizes

plotDF <- fill_data %>% 
  filter(type != "syn / inv") %>% 
  ungroup() %>% 
  mutate(
    species = factor(species, levels = rev(levels(species)))
  )

ggplot(plotDF) + 
  geom_violin(aes(y = size, x = species, fill = species), adjust = 0.5) +
  stat_summary(aes(x = species, y = size), fun.y = "median", geom = "point", size = 2, color = "black") +
  scale_y_log10(name = "Fill sizes [bp]", breaks = scales::trans_breaks("log10", function(x) 10^x),
                limit = c(100,1e+9), labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  coord_flip() +
  facet_grid(. ~ type) +
  theme_bw() + 
  scale_fill_manual(values = rev(COL_SPECIES)) +
  theme(
    legend.position = "none",
    # axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_blank()
  ) +
  labs(x = "", y = "")
                    

ggsave("results/fill_size_distributions_horiz.pdf", w = 3, h = 6)

# -------------------------------------------------------------
# Breakpoints
# -------------------------------------------------------------

# horizontal plot --------------------------------------------------------------

breakpoint_count <- crossing(
  species = SPECIES,
  threshold = THRESHOLDS
) %>% 
  mutate(
    bp_path = str_c("data/breakpoints/hg38.", genome_assembly, ".", 
                    threshold, ".bp.flt.flt_adj_fill.bed"),
    bp_gr = map(bp_path, import.bed, seqinfo = hum_seqinfo),
    n = map_int(bp_gr, length),
    species = factor(trivial_name, levels = SPECIES$trivial_name),
    threshold = factor(threshold/1000, levels = rev(c("10", "100", "1000")))
    
  )


# plot number of breakpoints for species
ggplot(breakpoint_count, aes(x = threshold, y = n, fill = species)) + 
  geom_bar(stat = "identity", color = "black") + 
  geom_text(aes(label = n), hjust = "inward") +
  coord_flip() +
  facet_grid(species ~ .) +
  ylab("Number of breakpoints") +
  xlab("Size threshold [kb]") +
  scale_fill_manual(values = COL_SPECIES) +
  theme_bw() + 
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        legend.position = "none",
  )

ggsave("results/breakpoint_numbers_horiz.pdf", w = 3, h = 6)
