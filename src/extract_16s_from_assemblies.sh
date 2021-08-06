#! /bin/bash
assm_file=$1
file_path=$(dirname $assm_file|sed 's/\/medeka//g')
assm=$(basename --suffix '.contigs.fasta' $assm_file)
source /gpfs/gsfs11/users/rajwanir2/conda/etc/profile.d/conda.sh
conda activate Prokka
mkdir $file_path/16s
barrnap --quiet --threads 16 $assm_file --outseq $file_path/16s/16s.fasta > $file_path/16s/16s.gff
##assemblies=$(find data/ -wholename 'data/assembly/*/*/medeka/*.contigs.fasta' | grep -v 'simulation')
##for g in $assemblies; do bash src/extract_16s_from_assemblies.sh $g;done
##blastn -num_threads 16 -max_target_seqs 1 -db /fdb/blastdb/v4/16SMicrobial -query data/assembly/09182020/barcode10/16s/16s.fasta -outfmt '6 delim=@ staxids stitle pident' | head -n1