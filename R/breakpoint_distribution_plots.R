
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

# area around a boundary where breakpoint distribution is investigated (only relevant for 2nd analysis)
BOUNDARY_PLUS <-  unlist(METADATA %>% 
                           filter(!is.na(boundary_plus_adjacence)) %>%
                           dplyr::select(boundary_plus_adjacence)
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

# =============================================================================================================================
# Plot breakpoint distributions (whole domain)
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
group_cols_transp = adjustcolor( group_cols, alpha.f = 0.2)

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

SELECTED_SPECIES <- c("panTro5", "bosTau8", "monDom5", "danRer10")
SPECIES_NAMES <- map(unname(unlist(SPECIES %>%
                          filter(genome_assembly %in% SELECTED_SPECIES) %>%
                          dplyr::select(trivial_name)
                          )
                        ), simpleCap)

NSPECIES <- length(SELECTED_SPECIES)
NSAMPLES <- length(unique(data_combined$sample))
NTHR <- length(THRESHOLDS)

# reverse sample name positions
SAMPLES <- unique(data_combined$sample)[length(unique(data_combined$sample)):1]

# group species, sample and threshold together
group_levels <- str_c(rep(SELECTED_SPECIES, NSAMPLES * NTHR),  
                      rep(SAMPLES, each = NSPECIES * NTHR), 
                      rep(THRESHOLDS, NSPECIES * NSAMPLES), 
                      sep = "_")
group_labels <- str_c(rep(SPECIES_NAMES, NSAMPLES * NTHR),
                      rep(c("Breakpoints", "Background"), each = NSPECIES * NTHR),
                      rep(c("10 kb", "100 kb", "1000 kb"), NSPECIES * NSAMPLES),
                      sep = " ")

plot_data <- data_combined %>% 
  left_join(SPECIES, by = c("species" = "genome_assembly")) %>% 
  filter(species %in% SELECTED_SPECIES, threshold == 10000, domains == "hESC") %>%
  mutate(species_sample_thr = str_c(species, "_", sample, "_", threshold)) %>%
  mutate(sample_thr = factor(species_sample_thr, levels = group_levels, labels = group_labels))

# choose colors (sample is atlernating)
group_cols <- c(brewer.pal(3,"Reds"), brewer.pal(3,"Blues"), brewer.pal(3,"Greens"), brewer.pal(3,"Purples"))
group_cols <- group_cols[c(3, 1, 6, 4, 9, 7, 12, 10)]
group_cols_transp = adjustcolor( group_cols, alpha.f = 0.2)
  

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
# Plot breakpoint distributions of all thresholds for one species (at boundaries)
# =============================================================================================================================

# store results here
out_dir <- "results/boundaries/"
dir.create(out_dir, showWarnings = FALSE)

data <- read_rds("results/breakpoints_at_boundaries.rds")

# prepare data for plotting
plot_data <- data %>% 
  group_by(boundaries, species, threshold) %>%
  # normalise breakpoint hits for each species and threshold 
  mutate(hits = hits / sum(hits) * 100, 
         rdm_hits = rdm_hits / sum(rdm_hits) * 100,
         rdm_hits_sd = sd(rdm_hits),
         bin = bin + 0.5)

for (D in DOMAINS$domain_type){
  
  for (S in SPECIES$genome_assembly) {
    
    this_plot_data <- filter(plot_data, species == S, boundaries == D)
    
    trivial_name <- unlist(SPECIES %>%
                             filter(genome_assembly == S) %>%
                             dplyr::select(trivial_name)
    )
    
    plot_title <- trivial_name
    
    ggplot(this_plot_data, 
           aes(x = bin, y = hits, colour = threshold)) + # breakpoint distribution
      geom_point() + 
      geom_line() + 
      geom_line(data = this_plot_data %>% 
                  group_by(bin) %>%
                  summarise(rdm_hits = mean(rdm_hits),
                            rdm_hits_sd = mean(rdm_hits_sd)), # random distribution (was normalised over all thresholds)
                aes(x = bin, y = rdm_hits),
                inherit.aes = FALSE) +
      geom_ribbon(data = this_plot_data %>% 
                    group_by(bin) %>%
                    summarise(rdm_hits = mean(rdm_hits),
                              rdm_hits_sd = mean(rdm_hits_sd)), # add sd for random breakpoints
                  aes(x = bin, ymin = rdm_hits - rdm_hits_sd, ymax = rdm_hits + rdm_hits_sd),
                  fill = "grey70", 
                  alpha = 0.5, 
                  inherit.aes = FALSE) + 
      scale_x_discrete(name = "Distance to TAD boundary [kb]", 
                       limits=seq(1,NBINS+1,NBINS/4), 
                       labels=c(as.character(-(BOUNDARY_PLUS/1000)), as.character(-(BOUNDARY_PLUS/2000)), 0, 
                                as.character((BOUNDARY_PLUS/2000)), as.character((BOUNDARY_PLUS/1000)))) +
      geom_vline(xintercept = (NBINS/2) + 1, linetype=3) +
      # annotate region where breakpoint enrichment is tested
      # annotate("rect", xmin=9, xmax=13, ymin=-Inf, ymax=Inf, alpha=0.6, fill="grey90") + 
      
      # from here all cosmetics like background, title, line colors...
      ggtitle(plot_title) +
      scale_colour_brewer("Chain sizes larger", 
                          palette = "Blues", 
                          labels = c("10 kb", "100 kb", "1000 kb")) +
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
# =============================================================================================================================
# Plot breakpoint distributions for several species (at boundaries)
# =============================================================================================================================

# store results here
out_dir <- "results/boundaries/"
dir.create(out_dir, showWarnings = FALSE)

data <- read_rds("results/breakpoints_at_boundaries.rds")

SELECTED_SPECIES <- c("panTro5", "bosTau8", "monDom5", "danRer10")
SPECIES_NAMES <- unname(unlist(SPECIES %>%
                                 filter(genome_assembly %in% SELECTED_SPECIES) %>%
                                 dplyr::select(trivial_name)
                               )
                        )

# prepare data for plotting
plot_data <- data %>% 
  # <<< select species, domain and threshold here >>>
  filter(species %in% SELECTED_SPECIES,
         boundaries == "hESC",
         threshold == "10000") %>%
  # normalise breakpoint hits for each species and threshold 
  group_by(boundaries, species, threshold) %>%
  mutate(hits = hits / sum(hits) * 100, 
         rdm_hits = rdm_hits / sum(rdm_hits) * 100,
         rdm_hits_sd = sd(rdm_hits),
         # add 0.5 to each bin for cosmetic reasons
         bin = bin + 0.5) 

ggplot(plot_data, 
       aes(x = bin, y = hits, colour = species)) +
  geom_point() + 
  geom_line() + 
  geom_line(data = plot_data %>% 
              group_by(bin) %>%
              summarise(rdm_hits = mean(rdm_hits),
                        rdm_hits_sd = mean(rdm_hits_sd)),
            aes(x = bin, y = rdm_hits),
            inherit.aes = FALSE) +
  geom_ribbon(data = plot_data %>% 
                group_by(bin) %>%
                summarise(rdm_hits = mean(rdm_hits),
                          rdm_hits_sd = mean(rdm_hits_sd)),
              aes(x = bin, ymin = rdm_hits - rdm_hits_sd, ymax = rdm_hits + rdm_hits_sd),
              fill = "grey70", 
              alpha = 0.5, 
              inherit.aes = FALSE) + 
  scale_color_brewer("", palette = "Set1", labels = SPECIES_NAMES) +
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
  filter(boundaries %in% c("hESC", "GM12878"))

# -------------------------------------------------------------------------------------------------------------------
# Compute a fisher test for comparing the amount of breakpoints inside the range BOUNDARY_SBIN, BOUNDARY_EBIN to
# the proportion of breakpoints outside that range to the ratio of random breakpoints.
# -------------------------------------------------------------------------------------------------------------------

# bin ranges where breakpoint ratios are compared
TEST_BINS <- seq(BOUNDARY_SBIN, BOUNDARY_EBIN)
OTHER_BINS <- seq(1, NBINS)[-TEST_BINS]

# sum up breakpoints and random breakpoints for test and other bins
fisher_input <- test_data %>%
  group_by(boundaries, threshold, species) %>%
  summarise(bp_test_sum = sum(hits[bin %in% TEST_BINS]),
            rdm_test_sum = sum(rdm_hits[bin %in% TEST_BINS]),
            bp_other_sum = sum(hits[bin %in% OTHER_BINS]),
            rdm_other_sum = sum(rdm_hits[bin %in% OTHER_BINS]))

results <- tibble()

for (D in unique(test_data$boundaries)){
    
  for (S in SPECIES$genome_assembly){
  
    for (THR in THRESHOLDS){    
      
      # create matrix for fisher test function
        fisher_matrix <- matrix(
          unlist(
            fisher_input %>%
            filter(boundaries == D, species == S, threshold == THR) %>%
            ungroup() %>%
            select(bp_test_sum, rdm_test_sum, bp_other_sum, rdm_other_sum)
          ),
          nrow = 2
        )
      
      # apply fisher test
      test_result <- fisher.test(fisher_matrix)
      
      # get trivial name
      trivial_name <- unlist(SPECIES %>%
        filter(genome_assembly == S) %>%
        select(trivial_name)
      )
      
      # gather results in data frame
      tmp <- tibble(boundaries = factor(D), species = factor(trivial_name), threshold = factor(THR), 
                    p_value = test_result$p.value, odds_ratio = test_result$estimate)
      results <- rbind.data.frame(results, tmp)
    }
  }
}

# -------------------------------------------------------------------------------------------------------------------
# Plot fisher odds ratios (log_2)
# -------------------------------------------------------------------------------------------------------------------

ggplot(results, aes(x = species, y = log2(odds_ratio), group=threshold))+ 
  geom_bar(aes(fill=threshold),
                  stat="identity", position="dodge", width=0.5) +
  geom_text(aes(label=asterisks(p_value), y = ifelse(odds_ratio < 1, 0, log2(odds_ratio))), 
                   position = position_dodge(width=0.5), vjust=0.8, hjust=-0.1, angle=90,
            inherit.aes = TRUE) +
  facet_grid(boundaries ~ .) +

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
