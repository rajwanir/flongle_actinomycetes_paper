#! /bin/bash
barcode=$1
date='simulation'
assembler="medeka"
input_folder="data/assembly/${date}/$barcode/"
##input_file="data/assembly/${date}/$barcode/wtdbg2/assm.ctg.fa"
module load quast
quast.py --threads 16 \
$input_folder/x7/$assembler/x7.contigs.fasta $input_folder/x15/$assembler/x15.contigs.fasta \
$input_folder/x30/$assembler/x30.contigs.fasta $input_folder/x60/$assembler/x60.contigs.fasta \
-R $input_folder/$barcode/$assembler/$barcode.contigs.fasta \
--glimmer \
 --output-dir data/quast/$date/$barcode/$assembler