#!/bin/bash

humanbam=$1
sample=$2
refidx=$3
ref=$4
binsize=$5
qual=$6
selector=$7
outdir=$8
refdir=$9
threads=${10}

parallel=/drive3/staging/agarofal/parallel/src/parallel
scriptDir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

echo "Counting reads in each bin..."

binbed="$refdir"/hg19_${binsize}bp_bins.gc.bed

if [[ ! -e "$outdir"/${sample}.depth.gc.bed ]]; then

	split -l$((`wc -l < "$binbed"`/$threads)) $binbed "$outdir"/tmpbed -da 4

	readlink -f "$outdir"/tmpbed* | $parallel --bar "sh $scriptDir/count-reads.sh $humanbam $sample {} $binsize $qual $outdir $refdir"
	
	cat $outdir/*tmp*bed > $outdir/${sample}.depth.gc.bed
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

rm $outdir/*tmp*

