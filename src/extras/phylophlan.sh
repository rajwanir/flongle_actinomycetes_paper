#! /bin/bash
source /gpfs/gsfs11/users/rajwanir2/conda/etc/profile.d/conda.sh
conda activate phylophlan
input_genomes=$(find . -path "*data/assembly/*/*/medeka/bar*.contigs.fasta")
cp $input_genomes data/phylophylan/input
phylophlan -i data/phylophylan/input \
-d phylophlan \
-f data/phylophylan/02_tol.cfg \
--diversity high \
--fast \
-o output_tol \
--nproc 16 \
--verbose 2>&1 | tee logs/phylophlan.log