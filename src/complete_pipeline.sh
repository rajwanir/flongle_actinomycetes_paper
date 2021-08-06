#!/bin/bash
barcodes=(barcode04 barcode06 barcode07 barcode09)
##canujobs=(4228013 4227317 4227850 4227848)
date='12112020'
##fast5_folder=$1
##job_id=$(sbatch --partition=gpu --cpus-per-task=14 --mem=16g --gres=lscratch:200,gpu:p100:1 src/guppy.sh $fast5_folder)
##job_id=$(sbatch --cpus-per-task=16 --mem=16g --depend=afterany:$job_id src/qcat.sh $fast5_folder/basecalled_fastq)
##mv $fast5_folder/basecalled_fastq data/rawdata/vanRsoilactino_${date}/basecalled_fastq
for i in {0..3};do
##job_id=$(sbatch --cpus-per-task=16 --mem=16g src/canu.sh ${barcodes[$i]]} $date) ## incorrect
##job_id=$(sbatch --cpus-per-task=16 --mem=64g --partition=quick src/racon.sh ${barcodes[$i]]} $date)
job_id=$(sbatch --cpus-per-task=16 --mem=64g --partition=quick src/medaka.sh ${barcodes[$i]]} $date)
job_id=$(sbatch --cpus-per-task=16 --mem=16g --partition=quick --depend=afterok:$job_id src/antismash.sh ${barcodes[$i]]} $date)
done
