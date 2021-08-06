#!/bin/bash
barcode=$1
date='09172020'
input_file="data/rawdata/vanRsoilactino_${date}/basecalled_fastq/demultiplexed/$barcode.fastq"
output_folder="data/assembly/$date/$barcode/wtdbg2/"
mkdir -p $output_folder
#first stage
~/wtdbg2/wtdbg2 -x ont -g 8m -t 16 \
 -i $input_file \
 -L 1000 -fo $output_folder/assm --edge-min 2 --rescue-low-cov-edges --node-max 6000 -S1 -X 6000 -K 10000
#second stage
~/wtdbg2/wtpoa-cns -t 16 -i $output_folder/assm.ctg.lay.gz -fo $output_folder/assm.ctg.fa