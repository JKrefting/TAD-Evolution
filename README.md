# TAD-Evolution

Genome organisation is crucial for a variety of genomic functions.
Topologically associating domains (TADs) are genomic regions characterised by frequent intra-regional contacts. These contacts are essential for gene regulation and the disruption of TADs should thus be negatively selected during evolution. In this project, we compare human TADs to rearrangements in whole-genome alignments of different vertebrate species to find out whether TADs can be considered conserved building blocks of genomes. Furthermore, gene expression analysis is conducted to investigate expression patterns of orthologs with regard to TAD rearrangements. Together, this tests the hypothesis that TADs constitute the structural basis for gene regulation whose conservation can be explained by maintaining important regulatory environments.   

# Workflow

First download and preprocess data:
./bash/download.sh
./bash/preprocess.sh

The basic data can be evaluated in fills_and_breakpoints_plots.R.

Then start analysis:

Analyse breakpoint distributions at domains and domain boundaries in breakpoint_distribution_analysis.R, the results can be visualised in breakpoint_distribution_plots.R.

Analyse the expression correlation of orthologs regarding conserved and rearranged domains by first classifying the domains in domain_classification.R. The correlation is conducted in
ortholog_expression_analysis.R and the results visualised in ortholog_expression_plots.R.
