# TAD-Evolution

This repository contains the source code for all analysis in the study:

Krefting J, Andrade-Navarro MA, Ibn-Salem J. Evolutionary stability of topologically associating domains is associated with conserved gene regulation. 

The following documentation guides through all steps including downloading of external source data, retrieving ortholog genes from ENSEMBL, filtering, annotation, and running all analysis.

### Abstract 
The human genome is highly organized in the three-dimensional nucleus. Chromosomes fold locally into topologically associating domains (TADs) defined by increased intra-domain chromatin contacts. TADs contribute to gene regulation by restricting chromatin interactions of regulatory sequences, such as enhancers, with their target genes. Disruption of TADs can result in altered gene expression and is associated to genetic diseases and cancers. However, it is not clear to which extent TAD regions are conserved in evolution and whether disruption of TADs by evolutionary rearrangements can alter gene expression.
Here, we hypothesize, that TADs represent essential functional units of genomes, which are selected against rearrangements during evolution. We investigate this using whole-genome alignments to identify evolutionary rearrangement breakpoints of different vertebrate species. Rearrangement breakpoints are strongly enriched at TAD boundaries and depleted within TADs across species. Furthermore, using gene expression data across many tissues in mouse and human, we show that genes within TADs have more conserved expression patterns. Disruption of TADs by evolutionary rearrangements is associated with changes in gene expression profiles, consistent with a functional role of TADs in gene expression regulation.
Together, these results indicate that TADs are conserved building blocks of genomes with regulatory functions that are rather often reshuffled as a whole instead of being disrupted by rearrangements. 

### Requirements:

The analysis is mainly implemented in R by using several additional packages:

 - R version 3.3.3 (2017-03-06)
 - R packages:
	- BSgenome.Hsapiens.UCSC.hg19_1.4.0 
	- stringr_1.2.0
	- ggsignif_0.4.0
	- RColorBrewer_1.1-2
	- dplyr_0.7.4 
	- purrr_0.2.4
	- readr_1.1.1
	- tidyr_0.7.2
	- tibble_1.3.4 
	- ggplot2_2.2.1 
	- tidyverse_1.1.1 
	- BSgenome_1.42.0 
	- rtracklayer_1.34.2
	- Biostrings_2.42.1
	- XVector_0.14.1
	- biomaRt_2.30.0
	- GenomicRanges_1.26.4
	- GenomeInfoDb_1.10.3 
	- IRanges_2.8.2
	- S4Vectors_0.12.2
	- BiocGenerics_0.20.0


 - Python 2.7 

# Workflow

## Download and preprocessing
For downloading all external data execute the [bash/download.sh](bash/download.sh) bash script:
```{bash}
sh bash/download.sh
```
Than, rearrangement breakpoints can be extracted from net-files for all species.

```{bash}
sh bash/preprocess.sh
```

## Fill number and size distribution

The basic data can be evaluated in [fills_and_breakpoints_plots.R](R/fills_and_breakpoints_plots.R).

## Breakpoint distribution around TADs
Analyse breakpoint distributions at domains and domain boundaries in [breakpoint_distribution_analysis.R](R/breakpoint_distribution_analysis.R), the results can be visualised in [breakpoint_distribution_plots.R](R/breakpoint_distribution_plots.R).

## Classification of TADs
 - [domain_classification.R](R/domain_classification.R)
 - [domain_to_GRBs.R](R/domain_to_GRBs.R)
 
## Ortholog Expression correlation
 - [ortholog_expression_analysis.R](R/ortholog_expression_analysis.R)
 - [ortholog_expression_plots.R](R/ortholog_expression_plots.R)
