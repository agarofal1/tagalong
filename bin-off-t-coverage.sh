#!/bin/bash

#humanbam=$1
sample=$2
refidx=$3
ref=$4
binsize=$5
qual=$6
selector=$7
outdir=$8

if [[ ! -e /drive3/staging/agarofal/indexes/hg19_${binsize}bp_bins.gc.bed ]]; then
	echo "Splitting genome into ${binsize} bp bins..."
	bedtools makewindows -w $binsize -g $refidx > /drive3/staging/agarofal/indexes/hg19_${binsize}bp_bins.bed
	echo "Computing GC content for each bin..."
	bedtools nuc -fi $ref -bed /drive3/staging/agarofal/indexes/hg19_${binsize}bp_bins.bed > /drive3/staging/agarofal/indexes/hg19_${binsize}bp_bins.info.bed
	bedtools nuc -fi $ref -bed /drive3/staging/agarofal/indexes/hg19_${binsize}bp_bins.bed | cut -f1-3,5 > /drive3/staging/agarofal/indexes/hg19_${binsize}bp_bins.gc.bed
fi

echo "Counting reads in each bin..."

if [[ ! -e "$outdir"/${sample}.depth.gc.bed ]]; then
	if [[ "$qual" != "0" ]]; then
		bedtools multicov -bed /drive3/staging/agarofal/indexes/hg19_${binsize}bp_bins.gc.bed -bams $humanbam -q $qual | cut -f1-5 > $outdir/${sample}.depth.gc.bed
	else
		bedtools multicov -bed /drive3/staging/agarofal/indexes/hg19_${binsize}bp_bins.gc.bed -bams $humanbam | cut -f1-5 > $outdir/${sample}.depth.gc.bed
	fi
fi

if [[ "$selector" != "NA" ]]; then
	echo "Subtracting on-target regions from coverage file..."
	bedtools intersect -v -a $outdir/${sample}.depth.gc.bed -b $selector > $outdir/${sample}.tmp
else
	cp $outdir/${sample}.depth.gc.bed $outdir/${sample}.tmp
fi

cat /drive3/staging/agarofal/HL-plasma-virome-capture/human_coverage/WG/head $outdir/${sample}.tmp > $outdir/${sample}.tmp2
echo "sample" > $outdir/${sample}.tmp3
numlines=`cat $outdir/${sample}.tmp | wc -l`
yes $sample | head -${numlines} >> $outdir/${sample}.tmp3
paste $outdir/${sample}.tmp3 $outdir/${sample}.tmp2 > $outdir/${sample}.depth.offtarget.gc.bed

rm $outdir/${sample}.tmp*

echo "Done."
