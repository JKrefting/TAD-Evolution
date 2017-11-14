#!/bin/bash

# ================================================================================================
# Get the coordinates of syntenic alignment blocks (from chains) and rearrangement breakpoints (from whole-genome alignments) and them in BED format using python scripts.
# ================================================================================================

set -o errexit # stop script on first error

mkdir -p 'data/breakpoints/'

# species for the analyis
SPECIES="$(cut -f1 metadata.csv | tail -n +2)"

# min size thresholds for syntenic blocks
THRESHOLDS="$(cut -f3 metadata.csv | tail -n +2)"

for S in ${SPECIES[@]}; do
	
	# Extract rearrangement breakpoints from alignments
	for T in ${THRESHOLDS[@]}; do
	
			python Python/find_breakpoints.py \
			-i 'data/alignments/'hg38.$S.net \
			-o 'data/breakpoints/'hg38.$S.$T.bp.bed \
			-t $T
			
			'bin/'liftOver \
        	'data/breakpoints/'hg38.$S.$T.bp.bed \
        	'data/liftover/'hg38ToHg19.over.chain \
        	'data/breakpoints/'hg19.$S.$T.bp.bed \
        	'data/breakpoints/'hg38Tohg19.$S.$T.unmapped.bp.bed
	done
	
	# get the coordinates of syntenic alignment blocks (chains)
	python Python/chain2bed.py \
			-i 'data/chains/'hg38.$S.all.chain \
			-o 'data/chains/'hg38.$S.chain.bed
			
			'bin/'liftOver \
        	'data/chains/'hg38.$S.chain.bed \
        	'data/liftover/'hg38ToHg19.over.chain \
        	'data/chains/'hg19.$S.chain.bed \
        	'data/chains/'hg38Tohg19.$S.unmapped.chain.bed
done
