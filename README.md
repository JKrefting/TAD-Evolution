# TAD-Evolution

Genome organisation is crucial for a variety of genomic functions.
Topologically associating domains (TADs) are genomic regions characterised by frequent intra-regional contacts. These contacts are essential for gene regulation and the disruption of TADs should thus be negatively selected during evolution. In this project, we compare human TADs to rearrangements in whole-genome alignments of different vertebrate species to find out whether TADs can be considered conserved building blocks of genomes. Furthermore, gene expression analysis is conducted to investigate expression patterns of orthologs with regard to TAD rearrangements. Together, this tests the hypothesis that TADs constitute the structural basis for gene regulation whose conservation can be explained by maintaining important regulatory environments.   

# Workflow

## Download and preprocessing
For downloading all external data execute the [bash/download.sh](bash/download.sh) bash script:
```{bash}
sh bash/download.sh
```
Than, rearragement breakpoints can be extracted from net-files for all species.

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
