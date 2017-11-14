# =============================================================================================================================
# Plot breakpoint distributions with regard to stuctural domains of the human genome. 
# =============================================================================================================================

require(tidyverse)

# read metadata for analysis
METADATA <- read_tsv("metadata.csv")
SPECIES <- select(METADATA, genome_assembly, trivial_name)
DOMAINS <- METADATA %>% 
  select(genomic_domain_type, genomic_domain_path) %>%
  filter(!is.na(genomic_domain_type))
THRESHOLDS <- unlist(METADATA %>% 
                       filter(!is.na(min_size_threshold)) %>%
                       select(min_size_threshold)
)

NBINS <- 20

data <- readRDS("results/breakpoints_at_domains.rds")

# normalise random breakpoint hits for each species
# normalise breakpoint hits for each species and threshold 
# add 0.5 to each bin for cosmetic reasons
data <- data %>% 
  group_by(domains, species) %>%
  mutate(rdm_hits = rdm_hits / sum(rdm_hits) * 100,
         rdm_hits_sd = sd(rdm_hits / sum(rdm_hits) * 100)) %>%
  group_by(domains, species, threshold) %>%
  mutate(hits = hits / sum(hits) * 100, 
         bin = bin + 0.5)
  
# -------------------------------------------------------------------------------------------------------------------
# Plot breakpoint distributions of all thresholds for one species
# -------------------------------------------------------------------------------------------------------------------

# example for mouse breakpoints at hESC TADs:
# SP <- factor("mm10", levels = levels(data$species))
# D <- factor("hESC", levels = levels(data$domains))
# then run plot code without the loops

plot_data <- filter(data, species == SP, domains == D)

plot_title <- unlist(SPECIES %>% 
  filter(genome_assembly == SP) %>%
  select(trivial_name)
)

ggplot(plot_data, 
       aes(x = bin, y = hits, colour = threshold)) + # breakpoint distribution
  geom_point() + 
  geom_line() + 
  geom_line(data = filter(plot_data, threshold == "10000"), # random distribution (was normalised over all thresholds)
                  aes(x = bin, y = rdm_hits),
           inherit.aes = FALSE) +
  geom_ribbon(data = filter(plot_data, threshold == "10000"), # add sd for random breakpoints
             aes(x = bin, ymin = rdm_hits - rdm_hits_sd, ymax = rdm_hits + rdm_hits_sd),
             fill = "grey70", 
             alpha = 0.5, 
             inherit.aes = FALSE) + 
  scale_x_discrete(name="", 
                   limits=seq(1,NBINS + 1,NBINS/4), 
                   labels=c("-50%", "start", "TAD", "end", "+50%")) +  # adapt x axis
  geom_vline(xintercept = c(NBINS/4, (NBINS-(NBINS/4)))+1, linetype=3) + # TAD boundary lines
  
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
