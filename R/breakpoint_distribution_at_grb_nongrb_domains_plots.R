
# =============================================================================================================================
# Plot breakpoint distributions with regard to stuctural domains of the human genome. 
# =============================================================================================================================

require(tidyverse)
require(ggplot2)
require(RColorBrewer)
require(ggsignif)
require(stringr)

source("R/functions.R")

# Read metadata for analysis
SPECIES <- read_tsv("species_meta.tsv")
DOMAINS <- read_tsv("domains_meta.tsv")
METADATA <- read_tsv("metadata.tsv")

# =============================================================================================================================
# Plot breakpoint distributions at whole domain
# =============================================================================================================================

# store results here
out_dir <- "results/whole_domain/"
dir.create(out_dir, showWarnings = FALSE)

data <- read_rds("results/breakpoints_at_domains_new.rds")

# combine hits by sample replicates
data_combined <- data %>% 
  filter(domain_subtype %in% c("GRB", "nonGRB")) %>%
  group_by(sample, replicate, species, threshold, domain_type, domain_subtype) %>% 
  mutate(
    percent = hits / sum(hits) * 100
  ) %>% 
  ungroup() %>% 
  group_by(bin, sample, species, threshold, domain_type, domain_subtype) %>% 
  summarise(
    n = n(),
    mean_hits = mean(hits),
    sd_hits = sd(hits),
    se_hits = std(hits),
    mean_percent = mean(percent),
    sd_percent = sd(percent),
    se_percent = std(percent)
  ) %>% 
  ungroup() 

# minimum size threasholds for chains (syntenic blocks) to be considered rearrangement blocks
THRESHOLDS <- unlist(METADATA %>% 
                       filter(!is.na(min_size_threshold)) %>%
                       dplyr::select(min_size_threshold)
)

# number of bins where breakpoint distribution is investigated
NBINS <- unlist(METADATA %>% 
                  filter(!is.na(n_bins)) %>%
                  dplyr::select(n_bins)
)

# -------------------------------------------------------------------------------------------------------------------
# Plot single species, all thresholds
# -------------------------------------------------------------------------------------------------------------------

# group sample and threshold together
group_levels <- str_c(rep(c("real", "random"), 2), "_", rep(c("GRB", "nonGRB"), each = 2))
group_labels <- str_c(rep(c("Breakpoints", "Background"), 2), " ", rep(c("GRB", "nonGRB"), each = 2))

plot_data <- data_combined %>% 
  left_join(SPECIES, by = c("species" = "genome_assembly")) %>% 
  mutate(sample_and_subtype = str_c(sample, "_", domain_subtype)) %>%
  mutate(sample_and_subtype = factor(sample_and_subtype, levels = group_levels, labels = group_labels)) %>%
  arrange(species, threshold, domain_type, sample_and_subtype)
  
# choose colors 
group_cols <- c(brewer.pal(10,"PuOr"))
group_cols <- group_cols[c(2, 4, 9, 7)]
grid::grid.raster(group_cols, int=F)
    
this_plot_data <- 
  filter(plot_data, species == "mm10", domain_type == "hESC", threshold == 10000)

trivial_name <- simpleCap(this_plot_data[1,] %>% pull(trivial_name))

ggplot(this_plot_data, 
       aes(x = bin + 0.5, y = mean_percent, colour = sample_and_subtype)) +
  
  geom_line(aes(linetype = sample)) + 
  
  # different linetypes for real and random breakpoints -> legend
  scale_color_manual(values = group_cols, guide = FALSE) +
  scale_linetype_manual(values=c("solid", "dashed"), 
                        limits = c("real", "random"),
                        labels = c("Breakpoints", "Background")) +
  
  # sd for background
  geom_ribbon(aes(ymin = mean_percent - se_percent,
                  ymax = mean_percent + se_percent,
                  group = sample_and_subtype,
                  alpha = 0.3),
              fill = "gray",
              color = NA,
              show.legend = FALSE) +
  
  # create invisible aes for subtype
  geom_boxplot(aes(fill = domain_subtype), color = NA) +
  scale_fill_manual(values = group_cols[c(1, 3)],
                    limits = c("real_GRB", "real_nonGRB"),
                    labels = c("GRB-TADs", "non-GRB-TADs")) +
  guides(fill = guide_legend(override.aes = list(fill = group_cols[c(1,3)]), order = 1)) +
  
  # scale_fill_manual(values = group_cols, guide = FALSE) + 
  scale_x_discrete(name="", 
                   limits=seq(1,NBINS + 1,NBINS/4), 
                   labels=c("-50%", "start", "TAD", "end", "+50%")) +  # adapt x axis
  geom_vline(xintercept = c(NBINS/4, (NBINS-(NBINS/4)))+1, linetype=3) + # TAD boundary lines
  ylab("Breakpoints [%]") +
  
  # from here all cosmetics like background, title, line colors...
  theme_bw() +
  theme(plot.title = element_text(size=12, face = "bold", colour = "black"),
        legend.title = element_blank(), legend.text = element_text(size = 11, face="bold"),
        legend.position = "bottom",
        axis.title.x = element_text(face="bold", size=11),
        axis.title.y = element_text(face="bold", size=11),
        axis.text.x = element_text(face="bold", size=11), 
        axis.text.y = element_text(face="bold", size=11),
        plot.margin=unit(c(0.1,0.5,0,0.5), "cm")
        # panel.grid.major = element_blank(),
        # panel.grid.minor = element_blank())
  )

ggsave(str_c(out_dir, trivial_name, "_", D ,"grb_subtypes_whole.pdf"), w=6, h=4)
