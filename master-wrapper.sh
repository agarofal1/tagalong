#!/bin/bash
usage() { echo "Usage: $0 [-b <bam file>] [-s <sample>] [-r <reference genome index>] [-h <human reference genome>] [-v <viral blast database>] [-l <bin size>] [-w <posterior window size>] [-Q <mapping quality threshold>] [-a <bin significance threshold>] [-L <human selector>] [-c <path to control data>] [-p <number of parallel threads>] [-S <seed for simulations>] [-o <output directory>]" 1>&2; exit 1; } 
 
while getopts ":b:s:r:h:v:l:w:Q:a:L:c:p:S:o:" o; do 
        case "${o}" in 
                b) 
                        bam=${OPTARG} 
                        ;;
		s)
			sample=${OPTARG}
			;;	 
                r) 
                        refidx=${OPTARG} 
                        ;; 
                h) 
                        humanRef=${OPTARG} 
                        ;; 
		v)
			viralDB=${OPTARG}
			;;
                l)
                        binSize=${OPTARG}
                        ;;
		w)
			windowSize=${OPTARG}
			;;
		Q)
			mapQ=${OPTARG}
			;;
		a)
			alpha=${OPTARG}
			;;
                L)
                        humanRegion=${OPTARG}
                        ;;
		c)
			controlData=${OPTARG}
			;;
		p)
			threads=${OPTARG}
			;;
		S)
			seed=${OPTARG}
			;;
                o)
                        outDir=${OPTARG} 
                        ;; 
                *)
                        usage
                        ;;

        esac
done
shift $((OPTIND-1))

if [[ -z "${bam}" || -z "${sample}" || -z "${outDir}" || -z "${controlData}" || -z "${viralDB}" ]]; then
	usage 
        exit                                                                                                                     
fi

if [[ -z "${binSize}" ]]; then
	binSize="10000"
fi   

if [[ -z "${mapQ}" ]]; then
	mapQ="30"
fi

if [[ -z "${windowSize}" ]]; then
	windowSize="100"
fi

if [[ -z "${threads}" ]]; then
	threads=12
fi

if [[ -z "${seed}" ]]; then
	seed=1
fi

if [[ -z "${alpha}" ]]; then
	alpha="0.01"
fi

if [[ ! -e /drive3/staging/agarofal/indexes/hg19_${binSize}bp_bins.gc.bed && -z $refidx ]] ; then
	usage
	exit
elif [[ -e /drive3/staging/agarofal/indexes/hg19_${binSize}bp_bins.gc.bed && -z $refidx ]]; then
	refidx="NA"
	humanRef="NA"
fi

if [[ -z "${humanRegion}" ]]; then
	humanRegion="NA"
fi

scriptDir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

echo "BAM file : $bam"
echo "Sample : $sample"
echo "Human reference genome : $humanRef"
echo "Reference index : $refidx"
echo "Viral BLAST DB : $viralDB"
echo "Bin size : $binSize"
echo "Window size : $windowSize"
echo "MapQ threshold : $mapQ"
echo "Sig alpha : $alpha"
echo "Num threads : $threads"
echo "Seed : $seed"
echo "Human selector : $humanRegion"
echo "Control data : $controlData"
echo "Output directory : $outDir"

refDir="$outDir"/ref
if [[ ! -d "$refDir" ]]; then
	mkdir $refDir
fi

if [[ ! -e "$refDir"/hg19_${binSize}bp_bins.gc.bed ]]; then
        echo "Splitting genome into ${binSize} bp bins..."
        bedtools makewindows -w $binSize -g $refidx > "$refDir"/hg19_${binSize}bp_bins.bed
        echo "Computing GC content for each bin..."
        bedtools nuc -fi $humanRef -bed "$refDir"/hg19_${binSize}bp_bins.bed > "$refDir"/hg19_${binSize}bp_bins.info.bed
        bedtools nuc -fi $humanRef -bed "$refDir"/hg19_${binSize}bp_bins.bed | cut -f1-3,5 > "$refDir"/hg19_${binSize}bp_bins.gc.bed
fi

if [[ ! -d "$outDir"/"$sample" ]]; then
	mkdir $outDir/$sample
fi

cnvDir=$outDir/$sample/cnv
breakDir=$outDir/$sample/breaks

if [[ ! -d "$cnvDir" ]]; then
	mkdir $cnvDir 
fi

if [[ ! -d "$breakDir" ]]; then
	mkdir $breakDir
fi

if [[ ! -e "$cnvDir"/${sample}.depth.offtarget.gc.bed ]]; then
	echo "Generating read counts in $binSize bp bins..."
	
	sh $scriptDir/bin-off-t-coverage.sh $bam $sample $refidx $humanRef $binSize $mapQ $humanRegion $cnvDir $refDir $threads

	echo "Done."
fi

echo "Normalizing read count data and identifying outlier bins..."

infoBed=/drive3/staging/agarofal/indexes/hg19_${binSize}bp_bins.info.bed

if [[ ! -e "$infoBed" ]]; then
	echo "Genomic bin metadata bed file is missing. Exiting."
	exit
fi

export LD_LIBRARY_PATH=$HOME/usr/lib64:$HOME/bin/lib:$HOME/bin/lib64

if [[ ! -e "$cnvDir"/${sample}.sig_bins.bed ]]; then
	Rscript $scriptDir/bayes-CNV.R $infoBed $controlData $cnvDir/${sample}.depth.offtarget.gc.bed $binSize $windowSize $threads $alpha $cnvDir $seed
fi

echo "Extracting and filtering soft-clipped reads from outlier bins..."

sh $scriptDir/extract-viral-softclips.sh $cnvDir/${sample}.sig_bins.bed $sample $bam 1000 60 15 10 $breakDir

rm $breakDir/*leftclip* $breakDir/*rightclip* $breakDir/${sample}.softclip*sam 

echo "BLASTing soft-clipped reads to viral DB..."

sh $scriptDir/run-blast.sh $breakDir $viralDB 90

cat $breakDir/*blast_results.txt > $breakDir/${sample}.all.clusters.blast.results.txt

rm $breakDir/*fa $breakDir/*blast_results.txt

echo "Interpreting BLAST hits and generating breakpoints..."

Rscript $scriptDir/generate-outputs.R $sample $breakDir/*blast.results.txt $breakDir/*passed.softclip.reads.txt 0.9 0.01 85 0.9 $breakDir

echo "Done." 

