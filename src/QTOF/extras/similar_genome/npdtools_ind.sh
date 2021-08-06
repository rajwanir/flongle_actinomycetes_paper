#! /bin/bash
strain=$1
ref_acc=$2
input_data="data/QTOF/01082021/*$strain*"
program='moldiscovery' ##dereplicator+ or varquest or metaminer or moldiscovery
db="prism_similargenome_ind/$ref_acc" 
dbtools_options="--reuse --fdr --fdr-limit 1 --db-path npdtools/$db"
metaminer_options='--class all --correspondence data/QTOF/metaminer_prism/corresp.tsv --sequence data/extracted_seqs/ripps_prism/'
varquest_options='--accurate_p_value --extract-best --adduct_Na --adduct_K'
npdtoolsBin='npdtools2/bin'
general_options="--mode HH --fdr"
$npdtoolsBin/$program.py \
--debug $general_options $input_data -o data/QTOF/${program}_prismSimInd/${strain}_${ref_acc}/ $dbtools_options --threads 1
##sbatch --cpus-per-task=16 --ntasks=10 --partition=multinode --mem=200g src/npdtools.sh
##--threads $SLURM_CPUS_PER_TASK --preprocess