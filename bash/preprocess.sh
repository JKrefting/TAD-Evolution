#!/bin/bash

# ==============================================================================
# Get the coordinates of rearrangement breakpoints (from whole-genome alignments) 
# and save in BED format using python scripts.
# ==============================================================================

set -o errexit # stop script on first error

mkdir -p 'data/breakpoints/'
mkdir -p 'data/fills/'

# species for the analyis
SPECIES="$(cut -f1 species_meta.tsv | tail -n +2)"

# min size thresholds for syntenic blocks
THRESHOLDS="$(cut -f1 metadata.tsv | tail -n +2)"

for S in $SPECIES ; do

    echo "INFO: "$S" -----------------------------------------------------------#"
  
    # get the coordinates of syntenic alignment blocks (fills)
    python Python/get_fills.py  \
        -i data/alignments/hg38.$S.net \
        -o data/fills/hg38.$S.fill.bed
    
    for T in $THRESHOLDS ; do
    
        # Extract rearrangement breakpoints from alignments
        python Python/find_breakpoints.py \
            -i 'data/alignments/'hg38.${S}.net \
            -o 'data/breakpoints/'hg38.${S}.${T}.bp.bed \
            -t ${T}
        
        # check if file is present and non-empty    
        if [ -s data/breakpoints/hg38.${S}.${T}.bp.be ]
        then
            # filter for plain chromosomes
            cat data/breakpoints/hg38.${S}.${T}.bp.bed \
              | egrep -v "(KI|GL|JH|KB)" \
              > data/breakpoints/hg38.${S}.${T}.bp.flt.bed
        else
            # create emty file
            touch data/breakpoints/hg38.${S}.${T}.bp.flt.bed
        fi
            
    done
done
