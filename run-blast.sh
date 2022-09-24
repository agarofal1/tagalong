#!/bin/bash

fastaDir=$1
viralDB=$2
pidentThresh=$3

export blastdb="$viralDB"
export pthresh="$pidentThresh"

ls $fastaDir/*fa | /drive3/staging/agarofal/parallel/src/parallel --jobs "10%" --bar 'fasta={1}; outpath=`echo "$fasta" | sed 's/nonhuman_clips.fa/blast_results.txt/g'`; blastn -db "$blastdb" -query {1} -out "$outpath" -perc_identity "$pthresh" -word_size 11 -outfmt "6 qseqid sseqid pident qlen length mismatch gapopen qstart qend sstart send evalue bitscore"'
