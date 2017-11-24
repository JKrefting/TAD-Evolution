
# =============================================================================================================================
# Plots of data after preprocessing: fills and breakpoints extracted from whole-genome alignemnts (net files)
# =============================================================================================================================

require(tidyverse)
require(ggplot2)
require(RColorBrewer)

source("R/functions.R")

# Read metadata for analysis
SPECIES <- read_tsv("species_meta.tsv")
METADATA <- read_tsv("metadata.tsv")

# minimum size threasholds for chains (syntenic blocks) to be considered rearrangement blocks
THRESHOLDS <- unlist(METADATA %>% 
                       filter(!is.na(min_size_threshold)) %>%
                       dplyr::select(min_size_threshold)
)


# -------------------------------------------------------------
# Fills
# -------------------------------------------------------------

# prepare data for plotting
plot_df <- tibble()
# gather all date in one tibble
for (S in SPECIES$genome_assembly) {
  
  fills <- readFillFile(S)
  
  trivial_name <- unlist(SPECIES %>%
                           filter(genome_assembly == S) %>%
                           dplyr::select(trivial_name)
                         )
  
  single_df <- tibble(size = width(fills), 
                      species = factor(trivial_name), 
                      type = mcols(fills)$name)
  
  plot_df <- rbind(plot_df, single_df)
}

# rename fill types
plot_df[plot_df$type == "top", ]$type <- "top-level"
plot_df[plot_df$type == "nonSyn", ]$type <- "non-syntenic"
plot_df[plot_df$type == "syn", ]$type <- "syn / inv"
plot_df[plot_df$type == "inv", ]$type <- "syn / inv"

# -------------------------------------------------------------
# Plot number of fills per species / type
# -------------------------------------------------------------

# ignore syn / inv fills
plot_df <- filter(plot_df, type != "syn / inv")
plot_df$type <- ordered(factor(plot_df$type, levels = c("top-level", "non-syntenic")))

ggplot(plot_df, aes(x = species, y = ..count.., )) + 
  geom_bar(aes(fill = species)) + 
  scale_y_continuous(name = "Number of fills", 
                 breaks = c(0, 25000, 50000, 75000, 100000, 125000, 150000),
             labels = c("0", "25000", "", "75000", "", "125000", "")) +
  theme_minimal() + 
  scale_fill_brewer("", palette = "Set3") +
  theme(plot.title = element_text(size=14, face = "bold", colour = "black"),
        strip.background = element_blank(),
        strip.text.y = element_text(size = 10, face="bold"),
        legend.title = element_text(size=7, face="bold"), legend.text = element_text(size = 7, face="bold"),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(face="bold", size=10),
        axis.text.x = element_text(face="bold", size=10, angle=0, vjust=0.5), 
        axis.text.y = element_text(face="bold", size=10),
        panel.grid.major = element_blank()
        # panel.grid.minor = element_blank()
  )

ggsave("results/fill_numbers.pdf", w=12, h=4)

# -------------------------------------------------------------
# Depict distributions of fill sizes
# -------------------------------------------------------------

ggplot(plot_df) + 
  geom_violin(aes(x = species, y = size, fill = species), adjust = 0.5) +
  stat_summary(aes(x = species, y = size), fun.y="median", geom="point", size=2, color="black") +
  scale_y_log10(name = "Chain sizes [bp]", breaks = scales::trans_breaks("log10", function(x) 10^x),
                limit=c(100,1e+9), labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  facet_grid(type ~ .) +
  theme_minimal() + 
  scale_fill_brewer("", palette = "Set3") +
  theme(plot.title = element_text(size=14, face = "bold", colour = "black"),
        strip.background = element_blank(),
        strip.text.y = element_text(size = 10, face="bold"),
        legend.title = element_text(size=7, face="bold"), legend.text = element_text(size = 7, face="bold"),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(face="bold", size=10),
        axis.text.x = element_text(face="bold", size=10, angle=0, vjust=0.5), 
        axis.text.y = element_text(face="bold", size=10),
        panel.grid.major = element_blank()
        # panel.grid.minor = element_blank()
  )

ggsave("results/fill_size_distributions.pdf", w=12, h=4)

# -------------------------------------------------------------
# Breakpoints
# -------------------------------------------------------------

# prepare data for plotting

plot_df <- tibble()

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
    
    plot_df <- rbind(plot_df, single_df)
  }
}

# plot number of breakpoints for species
ggplot(plot_df, aes(x = threshold, y = ..count..)) + 
  geom_bar(aes(fill = species)) + 
  facet_grid( ~ species) +
  theme_minimal() + 
  ylab("Number of breakpoints") +
  scale_fill_brewer("", palette = "Set3") +
  scale_x_discrete(name ="Chain size threshold [kb]", 
                   labels=c("10","100","1000")) +
  theme(plot.title = element_text(size=14, face = "bold", colour = "black"),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 10, face="bold"),
        legend.title = element_text(size=7, face="bold"), legend.text = element_text(size = 7, face="bold"),
        legend.position = "none",
        axis.title.x = element_text(face="bold", size=10),
        axis.title.y = element_text(face="bold", size=10),
        axis.text.x = element_text(face="bold", size=10, angle=45), 
        axis.text.y = element_text(face="bold", size=10),
        panel.grid.major = element_blank()
        # panel.grid.minor = element_blank()
  )

ggsave("results/breakpoint_numbers.pdf", w=12, h=4)