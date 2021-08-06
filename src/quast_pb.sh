#! /bin/bash
barcode=$1
date='10132020'
assembler="medeka"
input_folder="data/assembly/$date/$barcode/$assembler"
assembly_file="$barcode.contigs.fasta"
orignal_name="GB4-014"
module load quast
quast.py --threads 16 \
$input_folder/$assembly_file  \
-R data/pacbio_reference/$orignal_name.fasta \
--glimmer \
 --output-dir data/quast/$orignal_name