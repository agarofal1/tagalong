#!/bin/bash

sample=$1
fasta=$2
outdir=$3
clustPerc=$4

cluster=`basename $fasta | cut -d '_' -f3-5`
species=`basename $fasta | cut -d '_' -f6`

if (( $(echo "$clustPerc >= 0.95" | bc -l) )) && (( $(echo "$clustPerc < 1.0" | bc -l) ))
then
        wordLength=10
elif (( $(echo "$clustPerc >= 0.9" | bc -l) )) && (( $(echo "$clustPerc < 0.95" | bc -l) ))
then
        wordLength=8
elif (( $(echo "$clustPerc >= 0.88" | bc -l) )) && (( $(echo "$clustPerc < 0.9" | bc -l) ))
then
        wordLength=7
elif (( $(echo "$clustPerc >= 0.85" | bc -l) )) && (( $(echo "$clustPerc < 0.88" | bc -l) ))
then
        wordLength=6
elif (( $(echo "$clustPerc >= 0.8" | bc -l) )) && (( $(echo "$clustPerc < 0.85" | bc -l) ))
then
        wordLength=5
else
        echo "ERROR: Please select a cluster identity percentage between 0.8 and 1.0."
        exit
fi

cd-hit-est -i $fasta -o $outdir/${sample}_${cluster}_${species}.cdhit.out.txt -c $clustPerc -n $wordLength -d 0 -M 10000
