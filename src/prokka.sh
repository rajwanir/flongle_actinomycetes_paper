#! /bin/bash
input=$1
source /gpfs/gsfs11/users/rajwanir2/conda/etc/profile.d/conda.sh
conda activate Prokka
prokka --cpus ${SLURM_CPUS_PER_TASK} \
GPA_evolution/input_cluster_sequences_fasta/$input \
--outdir GPA_evolution/cluster_sequences_relabelled/$input \
--proteins GPA_evolution/multi.fa.centroids.tidynames.fa \