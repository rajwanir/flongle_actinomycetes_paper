#!/bin/bash
barcode=$1
date=$2
input_folder="data/assembly/$date"
assembler="canu"
reads="data/rawdata/vanRsoilactino_${date}/basecalled_fastq/demultiplexed/$barcode.fastq"
module load medaka || exit 1
medaka_consensus -i $reads \
-d $input_folder/$barcode/$assembler/$barcode.contigs.fasta \
-o $input_folder/$barcode/medeka/ \
-t $SLURM_CPUS_PER_TASK \
-b 150
mv $input_folder/$barcode/medeka/consensus.fasta $input_folder/$barcode/medeka/$barcode.contigs.fasta 