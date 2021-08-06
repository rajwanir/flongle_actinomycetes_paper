#! /bin/bash
barcode=$1
input_folder="aminocoumarin/input_cluster_sequences"
output_gcf="aminocoumarin"
source /gpfs/gsfs11/users/rajwanir2/conda/etc/profile.d/conda.sh
conda activate Antismash_latest
antismash $input_folder/$barcode.gbk \
--cpus ${SLURM_CPUS_PER_TASK} \
--logfile data/antismash/$output_gcf/$barcode/antismash.log \
--output-dir data/antismash/$output_gcf/$barcode/ \
--cb-general \
--cb-subclusters \
--cb-knownclusters \
--cb-knownclusters \
--genefinding-tool prodigal-m
## submit as a swarm jobs forexample
## swarm -t 8 -g 16 -f src/antismash_aminocoumarin_jobsub.txt --logdir data/antismash/aminocoumarin/