
# =============================================================================================================================
# Plots of data after preprocessing: fills and breakpoints extracted from whole-genome alignemnts (net files)
# =============================================================================================================================

require(BSgenome.Hsapiens.UCSC.hg19)
source("R/functions.R")

require(tidyverse)
require(ggplot2)
require(RColorBrewer)


# Read metadata for analysis
SPECIES <- read_tsv("species_meta.tsv")
METADATA <- read_tsv("metadata.tsv")

# minimum size threasholds for chains (syntenic blocks) to be considered rearrangement blocks
THRESHOLDS <- unlist(METADATA %>% 
                       filter(!is.na(min_size_threshold)) %>%
                       dplyr::select(min_size_threshold)
)

# define colors for species
COL_SPECIES = brewer.pal(12, "Set3")

# Load human seqinfo
genome <- BSgenome.Hsapiens.UCSC.hg19
hum_seqinfo <- seqinfo(genome)

# -------------------------------------------------------------
# Fills
# -------------------------------------------------------------

# prepare data for plotting
# gather all date in one tibble

data_df <- SPECIES %>% 
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
write_rds(data_df, "results/fill_data.rds")
write_tsv(data_df, "results/fill_data.tsv")
# data_df2 <- read_rds("results/fill_data.rds")

# -------------------------------------------------------------
# Plot number of fills per species / type
# -------------------------------------------------------------

# ignore syn / inv fills
data_df <- data_df %>% 
  filter(type %in% c("top-level", "non-syntenic")) %>% 
  mutate(
    type = factor(type, levels = c("top-level", "non-syntenic"))
  )

plot_df <- data_df %>%
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
  theme(plot.title = element_text(size=14, face = "bold", colour = "black"),
        strip.background = element_blank(),
        strip.text.y = element_text(size = 10, face="bold"),
        legend.title = element_text(size=7, face="bold"), legend.text = element_text(size = 7, face="bold"),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(face="bold", size=10),
        axis.text.x = element_text(face="bold", size=10, angle=45, vjust=1, hjust=1), 
        axis.text.y = element_text(face="bold", size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  )

ggsave("results/fill_numbers.pdf", w=6, h=4)


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

ggsave("results/fill_numbers_horiz.pdf", w=3, h=6)



# -------------------------------------------------------------
# Depict distributions of fill sizes
# -------------------------------------------------------------
# ignore syn / inv fills
plot_df <- filter(data_df, type != "syn / inv")
plot_df$type <- ordered(factor(plot_df$type, levels = c("top-level", "non-syntenic")))

ggplot(plot_df) + 
  geom_violin(aes(x = species, y = size, fill = species), adjust = 0.5) +
  stat_summary(aes(x = species, y = size), fun.y="median", geom="point", size=2, color="black") +
  scale_y_log10(name = "Fill sizes [bp]", breaks = scales::trans_breaks("log10", function(x) 10^x),
                limit=c(100,1e+9), labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  facet_grid(type ~ .) +
  theme_bw() + 
  scale_fill_brewer("", palette = "Set3") +
  theme(plot.title = element_text(size=14, face = "bold", colour = "black"),
        strip.background = element_blank(),
        strip.text.y = element_text(size = 10, face="bold"),
        legend.title = element_text(size=7, face="bold"), legend.text = element_text(size = 7, face="bold"),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(face="bold", size=10),
        axis.text.x = element_text(face="bold", size=10, angle=45, vjust = 1, hjust = 1), 
        axis.text.y = element_text(face="bold", size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  )


ggsave("results/fill_size_distributions.pdf", w=6, h=4)

# Horizontal plot of fill sizes

plotDF <- data_df %>% 
  filter(type != "syn / inv") %>% 
  ungroup() %>% 
  mutate(
    species = factor(species, levels = rev(levels(species)))
  )

ggplot(plotDF) + 
  geom_violin(aes(y = size, x = species, fill = species), adjust = 0.5) +
  stat_summary(aes(x = species, y = size), fun.y="median", geom="point", size=2, color="black") +
  scale_y_log10(name = "Fill sizes [bp]", breaks = scales::trans_breaks("log10", function(x) 10^x),
                limit=c(100,1e+9), labels = scales::trans_format("log10", scales::math_format(10^.x))) +
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
                    

ggsave("results/fill_size_distributions_horiz.pdf", w=3, h=6)

# -------------------------------------------------------------
# Breakpoints
# -------------------------------------------------------------

# join all breakpoints into one tibble
data_df <- tibble()
for (S in SPECIES$genome_assembly) {
  
  for (THR in THRESHOLDS){
    
    breakpoints <- readBPFile(S, THR)
    
    trivial_name <- unlist(SPECIES %>%
                             filter(genome_assembly == S) %>%
                             dplyr::select(trivial_name)
    )
    
    single_df <- tibble(breakpoint = start(breakpoints), 
                        species = factor(trivial_name), 
                        threshold = factor(THR)
                        )
    
    data_df <- rbind(data_df, single_df)
  }
}

# store data
write_rds(data_df, "results/breakpoint_data.rds")
write_tsv(data_df, "results/breakpoint_data.tsv")
# data_df2 <- read_rds("results/breakpoint_data.rds")

# prepare data for plotting
plot_df <- data_df %>%
  group_by(species, threshold) %>%
  summarise(n = n())

# plot number of breakpoints for species
ggplot(plot_df, aes(x = threshold, y = n)) + 
  geom_bar(aes(fill = species),
           stat = "identity",
           color = "black") + 
  geom_text(aes(label = n), 
             hjust = -0.2, 
             vjust = 0.5,
             size = 3,
             angle = 90,
             nudge_x = 0,
             nudge_y = 0) +
  facet_grid( ~ species) +
  ylab("Number of breakpoints") +
  ylim(c(0, 7000)) +
  scale_fill_brewer("", palette = "Set3") +
  scale_x_discrete(name ="", 
                   labels=c("10","100","1000")) +
  theme_bw() + 
  theme(plot.title = element_text(size=14, face = "bold", colour = "black"),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 10, face="bold", angle = 90),
        legend.title = element_text(size=7, face="bold"), legend.text = element_text(size = 7, face="bold"),
        legend.position = "none",
        axis.title.x = element_text(face="bold", size=10),
        axis.title.y = element_text(face="bold", size=10),
        axis.text.x = element_text(face="bold", size=10, angle=90, hjust = 1, vjust = 0.5), 
        axis.text.y = element_text(face="bold", size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  )

ggsave("results/breakpoint_numbers.pdf", w=6, h=4)

# horizontal plot --------------------------------------------------------------
breakpoint_count <- SPECIES %>% 
  select(species = trivial_name) %>% 
  tidyr::expand(species, threshold = THRESHOLDS) %>% 
  left_join(select(SPECIES, trivial_name, genome_assembly), by = c("species" = "trivial_name")) %>% 
  mutate(
    breakpoints = map2(genome_assembly, threshold, readBPFile),
    n = map_int(breakpoints, length),
    species = factor(species, levels = SPECIES$trivial_name),
    threshold = factor(threshold/1000, levels = rev(c("10", "100", "1000")))
  ) %>% 
  select(-breakpoints, -genome_assembly) 


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

ggsave("results/breakpoint_numbers_horiz.pdf", w=3, h=6)
