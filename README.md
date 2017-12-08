# TAD-Evolution

### Background 
The human genome is highly organized in the three-dimensional nucleus. Chromosomes fold locally into topologically associating domains (TADs) defined by increased intra-domain chromatin contacts. TADs contribute to gene regulation by restricting chromatin interactions of regulatory sequences, such as enhancers, with their target genes. Disruption of TADs can result in altered gene expression and is associated to genetic diseases and cancers. However, it is not clear to which extent TAD regions are conserved in evolution and whether disruption of TADs by evolutionary rearrangements can alter gene expression.

### Results
Here, we hypothesize, that TADs represent essential functional units of genomes, which are selected against rearrangements during evolution. We investigate this using whole-genome alignments to identify evolutionary rearrangement breakpoints of different vertebrate species. Rearrangement breakpoints are strongly enriched at TAD boundaries and depleted within TADs across species. Furthermore, using gene expression data across many tissues in mouse and human, we show that genes within TADs have more conserved expression patterns. Disruption of TADs by evolutionary rearrangements is associated with changes in gene expression profiles, consistent with a functional role of TADs in gene expression regulation.

### Conclusions
Together, these results indicate that TADs are conserved building blocks of genomes with regulatory functions that are rather often reshuffled as a whole instead of being disrupted by rearrangements. 

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
