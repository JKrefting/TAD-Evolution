
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

data <- read_rds("results/breakpoints_at_domains.rds")

# combine hits by sample replicates
data_combined <- data %>% 
  group_by(sample, replicate, species, threshold, domains) %>% 
  mutate(
    percent = hits / sum(hits) * 100
  ) %>% 
  ungroup() %>% 
  group_by(bin, sample, species, threshold, domains) %>% 
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
group_levels <- str_c(rep(c("real", "random"), 3), "_", rep(THRESHOLDS, each = 2))
group_labels <- str_c(rep(c("Breakpoints", "Background"), 3), " ", rep(c("10 kb", "100 kb", "1000 kb"), each = 2))

plot_data <- data_combined %>% 
  left_join(SPECIES, by = c("species" = "genome_assembly")) %>% 
  mutate(sample_thr = str_c(sample, "_", threshold)) %>%
  mutate(sample_thr = factor(sample_thr, levels = group_levels, labels = group_labels)) %>%
  arrange(threshold, desc(sample))

# choose colors (sample is atlernating)
group_cols <- c(brewer.pal(3,"Blues"), brewer.pal(3,"Greys"))
group_cols <- group_cols[c(1, 4, 2, 5, 3, 6)]

for (D in DOMAINS$domain_type){
  
  for (S in SPECIES$genome_assembly) {
    
    this_plot_data <- filter(plot_data, species == S, domains == D)
    
    trivial_name <- simpleCap(this_plot_data[1,] %>% pull(trivial_name))

    ggplot(this_plot_data, 
           aes(x = bin, y = mean_percent, colour = sample_thr, fill = sample_thr)) +
      geom_line() + 
      geom_ribbon(aes(ymin = mean_percent - se_percent,
                      ymax = mean_percent + se_percent),
                  color = NA, alpha = .1) +
      scale_color_manual(values = group_cols) +
      scale_fill_manual(values = group_cols) + 
      scale_x_discrete(name="", 
                       limits=seq(1,NBINS + 1,NBINS/4), 
                       labels=c("-50%", "start", "TAD", "end", "+50%")) +  # adapt x axis
      geom_vline(xintercept = c(NBINS/4, (NBINS-(NBINS/4)))+1, linetype=3) + # TAD boundary lines
      ylab("Breakpoints [%]") +
      
      # from here all cosmetics like background, title, line colors...
      theme_classic() +
      theme(plot.title = element_text(size=12, face = "bold", colour = "black"),
            legend.title = element_blank(), legend.text = element_text(size = 11, face="bold"),
            legend.position = "bottom",
            axis.title.x = element_text(face="bold", size=11),
            axis.title.y = element_text(face="bold", size=11),
            axis.text.x = element_text(face="bold", size=11), 
            axis.text.y = element_text(face="bold", size=11),
            plot.margin=unit(c(0.1,0.5,0,0.5), "cm"))
    
    ggsave(str_c(out_dir, trivial_name, "_", D, "_whole.pdf"), w=6, h=4)
    
  }
}

# -------------------------------------------------------------------------------------------------------------------
# Plot several species, single threshold
# -------------------------------------------------------------------------------------------------------------------

# determine species, threshold, domain to be plotted
SELECTED_SPECIES <- c("panTro5", "bosTau8", "monDom5", "danRer10")
SPECIES_NAMES <- map(unname(unlist(SPECIES %>%
                                     filter(genome_assembly %in% SELECTED_SPECIES) %>%
                                     dplyr::select(trivial_name)
)
), simpleCap)

THR <- THRESHOLDS[1]
D <- DOMAINS$domain_type[1]

NSPECIES <- length(SELECTED_SPECIES)
NSAMPLES <- length(unique(data_combined$sample))
NTHR <- length(THRESHOLDS)

# sample names
SAMPLES <- unique(data_combined$sample)

# group species, sample and threshold together
group_levels <- str_c(rep(SELECTED_SPECIES, each = NSAMPLES),  
                      rep(SAMPLES, NSPECIES), 
                      rep(THR, NSPECIES * NSAMPLES), 
                      sep = "_")
group_labels <- str_c(rep(SPECIES_NAMES, each = NSAMPLES),
                      rep(c("Breakpoints", "Background"), NSPECIES),
                      rep(c("10 kb"), NSPECIES * NSAMPLES),
                      sep = " ")

plot_data <- data_combined %>% 
  left_join(SPECIES, by = c("species" = "genome_assembly")) %>% 
  filter(species %in% SELECTED_SPECIES, threshold == THR, domains == D) %>%
  mutate(species_sample_thr = str_c(species, "_", sample, "_", threshold)) %>%
  mutate(species_sample_thr = factor(species_sample_thr, levels = group_levels, labels = group_labels))

# choose colors (sample is atlernating)
group_cols <- c(brewer.pal(3,"Reds"), brewer.pal(3,"Blues"), brewer.pal(3,"Greens"), brewer.pal(3,"Purples"))
group_cols <- group_cols[c(2, 3, 5, 6, 8, 9, 11, 12)]

ggplot(plot_data, 
       aes(x = bin, y = mean_percent, colour = species_sample_thr, fill = species_sample_thr)) +
  geom_line() + 
  # geom_ribbon(aes(ymin = mean_percent - se_percent,
                 # ymax = mean_percent + se_percent),
              # color = NA, alpha = .1) +
  scale_color_manual(values = group_cols) +
  scale_fill_manual(values = group_cols) +
  scale_x_discrete(name="", 
                   limits=seq(1,NBINS + 1,NBINS/4), 
                   labels=c("-50%", "start", "TAD", "end", "+50%")) +  # adapt x axis
  geom_vline(xintercept = c(NBINS/4, (NBINS-(NBINS/4)))+1, linetype=3) + # TAD boundary lines
  
  # from here all cosmetics like background, title, line colors...
  ylab("Breakpoints [%]") +
  theme_classic() +
  theme(plot.title = element_text(size=12, face = "bold", colour = "black"),
        legend.title = element_text(size=11, face="bold"), legend.text = element_text(size = 11, face="bold"),
        legend.position = "bottom",
        axis.title.x = element_text(face="bold", size=11),
        axis.title.y = element_text(face="bold", size=11),
        axis.text.x = element_text(face="bold", size=11), 
        axis.text.y = element_text(face="bold", size=11),
        plot.margin=unit(c(0.1,0.5,0,0.5), "cm"))


ggsave(str_c(out_dir, str_c(SPECIES_NAMES, collapse = "_"), "_whole.pdf"), w=6, h=4)


# =============================================================================================================================
# Plot breakpoint distributions at boundaries
# =============================================================================================================================

# store results here
out_dir <- "results/boundaries/"
dir.create(out_dir, showWarnings = FALSE)

data <- read_rds("results/breakpoints_at_boundaries.rds")

# combine hits by sample replicates
data_combined <- data %>% 
  group_by(sample, replicate, species, threshold, domains) %>% 
  mutate(
    percent = hits / sum(hits) * 100
  ) %>% 
  ungroup() %>% 
  group_by(bin, sample, species, threshold, domains) %>% 
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

# area around a boundary where breakpoint distribution is investigated 
BOUNDARY_PLUS <-  unlist(METADATA %>% 
                           filter(!is.na(boundary_plus_adjacence)) %>%
                           dplyr::select(boundary_plus_adjacence)
)

# -------------------------------------------------------------------------------------------------------------------
# Plot single species, all thresholds
# -------------------------------------------------------------------------------------------------------------------

# group sample and threshold together
group_levels <- str_c(rep(c("real", "random"), 3), "_", rep(THRESHOLDS, each = 2))
group_labels <- str_c(rep(c("Breakpoints", "Background"), 3), " ", rep(c("10 kb", "100 kb", "1000 kb"), each = 2))

plot_data <- data_combined %>% 
  left_join(SPECIES, by = c("species" = "genome_assembly")) %>% 
  mutate(sample_thr = str_c(sample, "_", threshold)) %>%
  mutate(sample_thr = factor(sample_thr, levels = group_levels, labels = group_labels)) %>%
  arrange(threshold, desc(sample))

# choose colors (sample is atlernating)
group_cols <- c(brewer.pal(3,"Blues"), brewer.pal(3,"Greys"))
group_cols <- group_cols[c(1, 4, 2, 5, 3, 6)]

for (D in DOMAINS$domain_type){
  
  for (S in SPECIES$genome_assembly) {
    
    this_plot_data <- filter(plot_data, species == S, domains == D)
    
    trivial_name <- simpleCap(this_plot_data[1,] %>% pull(trivial_name))
    
    ggplot(this_plot_data, 
           aes(x = bin, y = mean_percent, colour = sample_thr, fill = sample_thr)) +
      geom_line() + 
      geom_ribbon(aes(ymin = mean_percent - se_percent,
                      ymax = mean_percent + se_percent),
                  color = NA, alpha = .1) +
      scale_color_manual(values = group_cols) +
      scale_fill_manual(values = group_cols) + 
      scale_x_discrete(name = "Distance to boundary [kb]", 
                       limits=seq(1,NBINS+1,NBINS/4), 
                       labels=c(as.character(-(BOUNDARY_PLUS/1000)), as.character(-(BOUNDARY_PLUS/2000)), 0, 
                                as.character((BOUNDARY_PLUS/2000)), as.character((BOUNDARY_PLUS/1000)))) +
      geom_vline(xintercept = (NBINS/2) + 1, linetype=3) +
      # annotate region where breakpoint enrichment is tested
      # annotate("rect", xmin=9, xmax=13, ymin=-Inf, ymax=Inf, alpha=0.6, fill="grey90") + 
      
      # from here all cosmetics like background, title, line colors...
      ylab("Breakpoints [%]") +
      theme_classic() +
      theme(plot.title = element_text(size=12, face = "bold", colour = "black"),
            legend.title = element_text(size=11, face="bold"), legend.text = element_text(size = 11, face="bold"),
            legend.position = "bottom",
            axis.title.x = element_text(face="bold", size=11),
            axis.title.y = element_text(face="bold", size=11),
            axis.text.x = element_text(face="bold", size=11), 
            axis.text.y = element_text(face="bold", size=11),
            plot.margin=unit(c(0.1,0.5,0,0.5), "cm"))
    
      ggsave(str_c(out_dir, trivial_name,"_", D, "_boundary.pdf"), w=6, h=4)
    
  }
}

# -------------------------------------------------------------------------------------------------------------------
# Plot several species, single threshold
# -------------------------------------------------------------------------------------------------------------------

# determine species, threshold, domain to be plotted
SELECTED_SPECIES <- c("panTro5", "bosTau8", "monDom5", "danRer10")
SPECIES_NAMES <- map(unname(unlist(SPECIES %>%
                                     filter(genome_assembly %in% SELECTED_SPECIES) %>%
                                     dplyr::select(trivial_name)
)
), simpleCap)

THR <- THRESHOLDS[1]
D <- DOMAINS$domain_type[1]

NSPECIES <- length(SELECTED_SPECIES)
NSAMPLES <- length(unique(data_combined$sample))
NTHR <- length(THRESHOLDS)

# sample names
SAMPLES <- unique(data_combined$sample)

# group species, sample and threshold together
group_levels <- str_c(rep(SELECTED_SPECIES, each = NSAMPLES),  
                      rep(SAMPLES, NSPECIES), 
                      rep(THR, NSPECIES * NSAMPLES), 
                      sep = "_")
group_labels <- str_c(rep(SPECIES_NAMES, each = NSAMPLES),
                      rep(c("Breakpoints", "Background"), NSPECIES),
                      rep(c("10 kb"), NSPECIES * NSAMPLES),
                      sep = " ")

plot_data <- data_combined %>% 
  left_join(SPECIES, by = c("species" = "genome_assembly")) %>% 
  filter(species %in% SELECTED_SPECIES, threshold == THR, domains == D) %>%
  mutate(species_sample_thr = str_c(species, "_", sample, "_", threshold)) %>%
  mutate(species_sample_thr = factor(species_sample_thr, levels = group_levels, labels = group_labels))

# choose colors (sample is atlernating)
group_cols <- c(brewer.pal(3,"Reds"), brewer.pal(3,"Blues"), brewer.pal(3,"Greens"), brewer.pal(3,"Purples"))
group_cols <- group_cols[c(2, 3, 5, 6, 8, 9, 11, 12)]

ggplot(plot_data, 
       aes(x = bin, y = mean_percent, colour = species_sample_thr, fill = species_sample_thr)) +
  geom_line() + 
  # geom_ribbon(aes(ymin = mean_percent - se_percent,
  # ymax = mean_percent + se_percent),
  # color = NA, alpha = .1) +
  scale_color_manual(values = group_cols) +
  scale_fill_manual(values = group_cols) +
  scale_x_discrete(name = "Distance to TAD boundary [kb]", 
                   limits=seq(1,NBINS+1,NBINS/4), 
                   labels=c(as.character(-(BOUNDARY_PLUS/1000)), as.character(-(BOUNDARY_PLUS/2000)), 0, 
                            as.character((BOUNDARY_PLUS/2000)), as.character((BOUNDARY_PLUS/1000)))) +
  geom_vline(xintercept = (NBINS/2) + 1, linetype=3) +
  
  # from here all cosmetics like background, title, line colors...
  ylab("Breakpoints [%]") +
  theme_classic() +
  theme(plot.title = element_text(size=12, face = "bold", colour = "black"),
        legend.title = element_text(size=11, face="bold"), legend.text = element_text(size = 11, face="bold"),
        legend.position = "bottom",
        axis.title.x = element_text(face="bold", size=11),
        axis.title.y = element_text(face="bold", size=11),
        axis.text.x = element_text(face="bold", size=11), 
        axis.text.y = element_text(face="bold", size=11),
        plot.margin=unit(c(0.1,0.5,0,0.5), "cm"))


ggsave(str_c(out_dir, str_c(SPECIES_NAMES, collapse = "_"), "_boundary.pdf"), w=6, h=4)

# =============================================================================================================================
# Plot fisher odds ratios of breakpoint enrichments at domain boundaries
# =============================================================================================================================

data <- read_rds("results/breakpoints_at_boundaries.rds")

# choose data for analysis
test_data <- data %>% 
  filter(domains %in% c("hESC", "GM12878"))

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

# bins (around boundaries) that are investigated for breakpoint enrichment (only relevant for 2nd analysis)
BOUNDARY_SBIN <-  unlist(METADATA %>% 
                           filter(!is.na(boundary_area_start_bin)) %>%
                           dplyr::select(boundary_area_start_bin)
)
BOUNDARY_EBIN <-  unlist(METADATA %>% 
                           filter(!is.na(boundary_area_end_bin)) %>%
                           dplyr::select(boundary_area_end_bin)
)

# -------------------------------------------------------------------------------------------------------------------
# Compute a fisher test for comparing the amount of breakpoints inside the range BOUNDARY_SBIN, BOUNDARY_EBIN to
# the proportion of breakpoints outside that range to the ratio of random breakpoints.
# -------------------------------------------------------------------------------------------------------------------

# bin ranges where breakpoint ratios are compared
TEST_BINS <- seq(BOUNDARY_SBIN, BOUNDARY_EBIN)
OTHER_BINS <- seq(1, NBINS)[-TEST_BINS]

# sum up breakpoints and random breakpoints for test and other bins
fisher_input <- test_data %>%
  group_by(domains, threshold, species, sample) %>%
  summarise(test_sum = sum(hits[bin %in% TEST_BINS]),
            other_sum = sum(hits[bin %in% OTHER_BINS])
            )

# loop over groups, calculate fisher test
domains_list <- map(unique(test_data$domains), function(D){
    
  species_list <- map(SPECIES$genome_assembly, function(S){
  
    thr_list <- map(THRESHOLDS, function(THR){   
      
      # filter 
      this_fisher_input <- fisher_input %>%
        filter(domains == D, species == S, threshold == THR) %>%
        ungroup()
      
      # handle no breakpoints
      if (this_fisher_input[1, 5] == 0){
        return(tibble(domains = D, 
                      species = S, 
                      threshold = THR, 
                      p_value = NA, 
                      odds_ratio = NA))
      }
      
      # create matrix for fisher test function: #1 col = real, #2 = random
      fisher_matrix <- t(matrix(
        unlist(
          this_fisher_input %>%
            arrange(desc(sample)) %>%
            select(test_sum, other_sum)
        ),
        nrow = 2
        )
      )
        
      # apply fisher test
      test_result <- fisher.test(fisher_matrix)
      
      # gather results in data frame
      return(tibble(domains = D, 
                    species = S, 
                    threshold = THR, 
                    p_value = test_result$p.value, 
                    odds_ratio = test_result$estimate)
      )
    }) # end thr
    
    return(bind_rows(thr_list))

  }) # end species
  
  return(bind_rows(species_list))

}) # end domains

fisher_results <- bind_rows(domains_list)

# -------------------------------------------------------------------------------------------------------------------
# Plot fisher odds ratios (log_2)
# -------------------------------------------------------------------------------------------------------------------

# add trivial names
fisher_results <- fisher_results %>%
  left_join(SPECIES, by = c("species" = "genome_assembly"))

ggplot(fisher_results, 
       aes(x = trivial_name, 
           y = log2(odds_ratio), 
           group=factor(threshold), 
           fill=factor(threshold))) + 
  geom_bar(stat="identity", 
           position="dodge", 
           width=0.5) +
  geom_text(aes(label=asterisks(p_value), 
                y = ifelse(odds_ratio < 1, 0, log2(odds_ratio))), 
                   position = position_dodge(width=0.5), 
            vjust=0.8, 
            hjust=-0.1, 
            angle=90,
            inherit.aes = TRUE) +
  facet_grid(domains ~ .) +

# cosmetics
  geom_hline(yintercept=0, linetype="solid", color="#666666", alpha=0.5) +
  scale_fill_brewer("Chain sizes >", palette = "Blues", labels = c("10 kb", "100 kb", "1000 kb")) +
  xlab("") + ylab("log2(odds ratio)") + 
  theme_minimal() +
  theme(plot.title = element_text(size=14, face = "bold", colour = "black"),
               strip.background = element_blank(),
               strip.text.y = element_text(size = 10, face="bold"),
               legend.title = element_text(size=7, face="bold"), legend.text = element_text(size = 7, face="bold"),
               legend.position = "none",
               axis.title.x = element_text(face="bold", size=10),
               axis.title.y = element_text(face="bold", size=10),
               axis.text.x = element_text(face="bold", size=10, angle=0, vjust=0.5), 
               panel.grid.major = element_blank(),
               # panel.grid.minor = element_blank(),
               axis.text.y = element_text(face="bold", size=10))

ggsave("results/breakpoint_enrichment_odds_ratios.pdf", w=12, h=4)

# =============================================================================================================================
# Display distance distributions of breakpoints to domain boundaries and random breakpoints to domain boundaries
# =============================================================================================================================

data <- read_tsv("results/distances_to_boundaries.rds")

# choose data for analysis
plot_data <- data %>% 
  filter(boundaries %in% c("hESC", "GM12878"))

# rename species
plot_data$species <- factor(
  SPECIES[match(plot_data$species, SPECIES$genome_assembly), ]$trivial_name,
  levels = unlist(SPECIES$trivial_name))

ggplot(plot_data, aes(x = type, y = distance)) + 
  geom_boxplot(aes(color = type), outlier.size = 0.01, outlier.alpha = 0.3, width = 0.2) +
  facet_grid(boundaries ~ species) +
  geom_signif(comparisons=list(c("breakpoint", "random")), 
              map_signif_level = FALSE,
              test = wilcox.test,
              textsize = 3,
              tip_length = 0,
              y_position = 0,
              step_increase = 0.13,
              na.rm = TRUE) +
  scale_y_log10(name = "Distance to boundary [bp]", limit=c(1,1e+8), 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_color_discrete(name = "", labels=c("Breakpoints", "Control")) +
  xlab("") +
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
        axis.text.y = element_text(face="bold", size=10),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        # panel.grid.minor = element_blank(),
        plot.margin=unit(c(0.1,0.5,0,0.5), "cm")
        )

ggsave("results/distances_to_boundaries.pdf", w=12, h=4)
