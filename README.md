# TAD-Evolution

This repository contains the source code for all analysis in the study:

Krefting J, Andrade-Navarro MA, Ibn-Salem J. (2017) **Evolutionary stability of topologically associating domains is associated with conserved gene regulation.** 
bioRxiv 231431; doi: https://doi.org/10.1101/231431 

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
	- rtracklayer_1.34.2
	- biomaRt_2.30.0
	- ggsignif_0.4.0
	- RColorBrewer_1.1-2
	- tidyverse_1.1.1 
 - Python 2.7 

# Workflow
The scripts have to be executed in the order as listed bellow.

## Download and preprocessing
For downloading all external data execute the [bash/download.sh](bash/download.sh) bash script:
```{bash}
sh bash/download.sh
```
Then, rearrangement breakpoints can be extracted from net-files for all species:

```{bash}
sh bash/preprocess.sh
```

## Fill number and size distribution
Syntenic regions and breakpoitns are analysed in:

 - [fills_and_breakpoints_plots.R](R/fills_and_breakpoints_plots.R).

## Breakpoint distribution around TADs
Analysis of breakpoint distributions at domains and domain boundaries and visualisation of results:

 - [breakpoint_distribution_analysis.R](R/breakpoint_distribution_analysis.R)
 - [breakpoint_distribution_at_all_domains_plots.R](R/breakpoint_distribution_at_all_domains_plots.R)
 - [breakpoint_distribution_at_grb_nongrb_domains_plots.R](R/breakpoint_distribution_at_grb_nongrb_domains_plots.R)


## Classification of TADs
TADs are classified into conserved or rearranged TADs, as well as, GRB-TADs and non-GRB-TADs:

 - [domain_classification.R](R/domain_classification.R)
 - [domains_to_GRBs.R](R/domains_to_GRBs.R)
 - [domain_classification_plots.R](R/domain_classification_plots.R)
 
## Ortholog Expression correlation
Ortholog expression correlation across matching tissues in human and mouse:

 - [ortholog_expression_data.R](R/ortholog_expression_data.R)
 - [ortholog_expression_analysis.R](R/ortholog_expression_analysis.R)
 - [ortholog_expression_plots.R](R/ortholog_expression_plots.R)


