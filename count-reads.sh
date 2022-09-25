#!/bin/bash

humanbam=$1
sample=$2
bed=$3
binsize=$4
qual=$5
outdir=$6
refdir=$7

tranchenum=`basename "$bed" | sed 's/tmpbed//g'`

binbed="$refdir"/hg19_${binsize}bp_bins.gc.bed

if [[ "$qual" != "0" ]]; then
	bedtools multicov -bed $bed -bams $humanbam -q $qual | cut -f1-5 > $outdir/${sample}.tmp.${tranchenum}.bed
else    
	bedtools multicov -bed $bed -bams $humanbam | cut -f1-5 > $outdir/${sample}.tmp.${tranchenum}.bed
fi
