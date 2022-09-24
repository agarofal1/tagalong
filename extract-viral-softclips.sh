#!/bin/bash

sigbed=$1
sample=$2
bam=$3
pad=$4
mapq_thresh=$5
clip_length=$6
clust_support=$7
outdir=$8

if [[ ! -e "$outdir"/${sample}.softclip.q"$mapq_thresh".reads.sam ]]; then 
	bedtools slop -g /drive3/staging/agarofal/indexes/hg19.genome -i $sigbed -b $pad > $outdir/${sample}.sig.add${pad}bp.bed

	bedtools merge -i $outdir/${sample}.sig.add${pad}bp.bed > $outdir/${sample}.sig.add${pad}bp.merged.bed

	samtools view -F 2048 -q "$mapq_thresh" -L $outdir/${sample}.sig.add${pad}bp.merged.bed $bam | awk '$6 ~ "S"' > $outdir/${sample}.softclip.q"$mapq_thresh".reads.sam
fi

cat $outdir/${sample}.softclip.q"$mapq_thresh".reads.sam | awk 'substr($6, 0, 4) ~ "S" && $6 !~ /S$/' | cut -f1-6,10 > $outdir/${sample}.leftclip.q"$mapq_thresh".reads.txt

cat $outdir/${sample}.softclip.q"$mapq_thresh".reads.sam | awk '$6 ~ /S$/ && substr($6, 0, 4) !~ "S"' | cut -f1-6,10 > $outdir/${sample}.rightclip.q"$mapq_thresh".reads.txt

Rscript /drive3/staging/agarofal/CNV_integration_method/scripts/process-softclips.R $sample $outdir/${sample}.leftclip.q"$mapq_thresh".reads.txt $outdir/${sample}.rightclip.q"$mapq_thresh".reads.txt $clip_length $mapq_thresh $clust_support $outdir 

