#!/bin/bash
barcode=$1
date=$2
ml minimap2
##reads="data/rawdata/vanRsoilactino_${date}/basecalled_fastq/demultiplexed/$barcode.fastq"
reads="data/rawdata/${date}/$barcode.fastq"
output_folder="data/assembly/$date/$barcode/racon"
mkdir $output_folder
canu_assembly="data/assembly/$date/$barcode/canu/$barcode.contigs.fasta"
cp $canu_assembly $output_folder/${barcode}.0.fasta # copying as 0
###################################
#run three rounds of racon
for time in {1..3}
do
#align_reads
minimap2 -ax map-ont -a $output_folder/${barcode}.$((time-1)).fasta $reads > $output_folder/alignment_${time}.sam
#run racon
~/racon/build/bin/racon -m 8 -x -6 --gap -8 --window-length 500 --threads 16 \
$reads $output_folder/alignment_${time}.sam ${output_folder}/${barcode}.$((time-1)).fasta > $output_folder/${barcode}.$((time)).fasta
done 