#!/bin/bash
barcodes='barcode10 barcode08 barcode07 barcode01'
for barcode in $barcodes
do
sbatch --cpus-per-task=16 --output=data/assm/$barcode_slurm.out --error=data/assm/$barcode_slurm.err \
--mem=8g src/wtdbg2.sh $barcode
done