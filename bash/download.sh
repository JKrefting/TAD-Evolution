#!/bin/bash

# ================================================================================================
# Download all data needed for the analysis to 'data' folder in project directory
# 1: LiftOver data and tool
# 2: TAD calls of Dixon et al. 2012 and Rao et al. 2014
# 3: Whole-genome alignments and chains for species listed in 'species.csv'
# ================================================================================================

# stop script on first error
set -o errexit

#--------------------------------------------------------------------------------------------------
# UCSC liftover
#--------------------------------------------------------------------------------------------------

# UCSC liftover chains
mkdir -p 'data/liftover/'
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz -P 'data/liftover/' 
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz -P 'data/liftover/' 
wget http://hgdownload.cse.ucsc.edu/goldenPath/mm9/liftOver/mm9ToMm10.over.chain.gz -P 'data/liftover/'
gunzip 'data/liftover/'*.gz

# download liftOver tool from UCSC:
mkdir -p 'bin/'
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver -P 'bin/'
chmod u+x 'bin/'liftOver

# -------------------------------------------------------------------------------------------------
# Genomic Domains
# -------------------------------------------------------------------------------------------------

mkdir -p 'data/TADs/'

# TAD calls of Rao et al 2014, GM12878 cells
wget -P 'data/TADs/' ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt.gz
gunzip 'data/TADs/'GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt.gz
# reformat into bed file
tail -n +2 'data/TADs/'GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt \
			|cut -f 1-3 \
			| sed -e 's/^/chr/' \
			> 'data/TADs/'GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.bed 

	
# TAD calls of Dixon et al 2012, hESC
wget -P 'data/TADs/' http://chromosome.sdsc.edu/mouse/hi-c/hESC.domain.tar.gz  
# extract and rename
tar xvfz 'data/TADs/'hESC.domain.tar.gz -C 'data/TADs/'
cp 'data/TADs/hESC/combined/'total.combined.domain 'data/TADs/'hESC.hg18.bed
rm -r 'data/TADs/hESC/'
# liftover to hg19
	'bin/'liftOver \
        'data/TADs/'hESC.hg18.bed \
        'data/liftover/'hg18ToHg19.over.chain \
        'data/TADs/'hESC.hg19.bed \
        'data/TADs/'hESC.hg18Tohg19_unmapped.bed

# TAD calls of Dixon et al 2012, mouse ESC
wget -P 'data/TADs/' http://chromosome.sdsc.edu/mouse/hi-c/mESC.domain.tar.gz
tar xvfz 'data/TADs/'mESC.domain.tar.gz -C 'data/TADs/'
cp 'data/TADs/'mESC/*_combined/total.*.combined.domain 'data/TADs/'mESC.mm9.bed
rm -r 'data/TADs/mESC/'
# liftover to mm10
	'bin/'liftOver \
        'data/TADs/'mESC.mm9.bed \
        'data/liftover/'mm9ToMm10.over.chain \
        'data/TADs/'mESC.mm10.bed \
        'data/TADs/'mESC.mm9Tomm10_unmapped.bed

# GRBs of Harmston et al 2016
wget -P 'data/TADs/' https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-017-00524-5/MediaObjects/41467_2017_524_MOESM2_ESM.txt
cut -f1-3 41467_2017_524_MOESM2_ESM.txt | tail -n +2 > 'data/TADs/'41467_2017_524_MOESM2_ESM.bed

# ------------------------------------------------------------------------------------------------
# Alignments and chains
# ------------------------------------------------------------------------------------------------

# create folders
mkdir -p 'data/alignments/'
mkdir -p 'data/chains/'

# init variables
DLPATH='http://hgdownload.soe.ucsc.edu/goldenPath/hg38/'
SPECIES="$(cut -f1 metadata.csv | tail -n +2)"

# download alignment / chains for every species
for S in ${SPECIES[@]}; do

	UPPER_S="$(tr '[:lower:]' '[:upper:]' <<< ${S:0:1})${S:1}" # uppercase species

	# alignment
	alnfile='hg38.'$S'.net.gz'
	wget $DLPATH'vs'$UPPER_S'/'$alnfile -P 'data/alignments/'
	gzip -d 'data/alignments/'$alnfile 

	# chains  
	chainfile='hg38.'$S'.all.chain.gz'
	wget $DLPATH'vs'$UPPER_S'/'$chainfile -P 'data/chains/'
	gzip -d 'data/chains/'$chainfile 

done
