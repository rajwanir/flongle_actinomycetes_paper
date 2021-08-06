#!/bin/bash
barcode=$1
date='readlen'
#sample=$2
assmbler_from='racon'
input_folder="data/assembly/$date"
module load medaka || exit 1
medaka_consensus -i data/rawdata/$date/$barcode.fastq \
-d $input_folder/$barcode/$assmbler_from/$barcode.3.fasta \
-o $input_folder/$barcode/medeka/ \
-t $SLURM_CPUS_PER_TASK \
-b 150
mv $input_folder/$sample/$barcode/medeka/consensus.fasta $input_folder/$sample/$barcode/medeka/$barcode.contigs.fasta
##for x in 2 x7 x15 x30 x60;do sbatch --cpus-per-task=16 --mem=16g src/antismash.sh $x;done
##sbatch --cpus-per-task=16 --mem=16g src/antismash.sh x60