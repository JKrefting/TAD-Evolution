#!/bin/bash

# ================================================================================================
# Get the coordinates of rearrangement breakpoints (from whole-genome alignments) and save in BED format using python scripts.
# ================================================================================================

set -o errexit # stop script on first error

mkdir -p 'data/breakpoints/'
mkdir -p 'data/fills/'

# species for the analyis
SPECIES="$(cut -f1 species_meta.tsv | tail -n +2)"

# min size thresholds for syntenic blocks
THRESHOLDS="$(cut -f1 metadata.tsv | tail -n +2)"

for S in $SPECIES ; do
	
	# Extract rearrangement breakpoints from alignments
	for T in $THRESHOLDS ; do
	
			python Python/find_breakpoints.py \
			-i 'data/alignments/'hg38.$S.net \
			-o 'data/breakpoints/'hg38.$S.$T.bp.bed \
			-t $T
			
			# filter for plain chromosomes
      cat data/breakpoints/hg38.${S}.${T}.bp.bed \
        | egrep -v "(KI|GL|JH|KB)" \
        > data/breakpoints/hg38.$S.$T.bp.flt.bed

	done
	
	# get the coordinates of syntenic alignment blocks (fills)
	python Python/get_fills.py \
			-i 'data/alignments/'hg38.$S.net \
			-o 'data/fills/'hg38.$S.fill.bed
			
done
