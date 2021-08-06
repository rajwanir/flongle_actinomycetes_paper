#!/bin/bash
barcode=$1
date='09172020'
##input_file="data/assembly/${date}/$barcode/canu/$barcode.contigs.fasta"
input_file="data/assembly/${date}/$barcode/wtdbg2/assm.ctg.fa"
module load kraken/1.1
kraken --db /fdb/kraken/20200428_bacteria_64GB \
--output data/kraken/$barcode.krakenout \
--threads ${SLURM_CPUS_PER_TASK} \
--fasta-input \
$input_file
##--report data/kraken/$barcode.krakenreport \
##--use-names \
##--use-mpa-style \
#############
#translate
kraken-translate --db /fdb/kraken/20200428_bacteria_64GB \
data/kraken/$barcode.krakenout > data/kraken/$barcode.translated