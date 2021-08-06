#! /bin/bash
date='readlen'
assembler="medeka"
input_folder="data/assembly/${date}/"
##input_file="data/assembly/${date}/$barcode/wtdbg2/assm.ctg.fa"
module load quast
quast.py --threads 16 \
$input_folder/500/$assembler/500.contigs.fasta $input_folder/1000/$assembler/1000.contigs.fasta \
$input_folder/2000/$assembler/2000.contigs.fasta $input_folder/4000/$assembler/4000.contigs.fasta \
-R data/pacbio_reference/GB4-014.fasta \
--glimmer \
 --output-dir data/quast/$date/all/$assembler