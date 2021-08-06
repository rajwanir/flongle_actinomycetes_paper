#! /bin/bash
barcode=$1
date=$2
assembler='medeka'
input_file="data/assembly/${date}/$barcode/$assembler/$barcode.contigs.fasta"
source /gpfs/gsfs11/users/rajwanir2/conda/etc/profile.d/conda.sh
conda activate Antismash_latest
antismash $input_file \
--cpus ${SLURM_CPUS_PER_TASK} \
--logfile data/antismash/${date}/$barcode/$assembler/antismash.log \
--output-dir data/antismash/${date}/$barcode/$assembler/ \
--cb-general \
--cb-subclusters \
--cb-knownclusters \
--cb-knownclusters \
--genefinding-tool prodigal-m \
--asf --pfam2go --smcog-trees