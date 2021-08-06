#! /bin/bash
input_data='data/QTOF/GB4-14/*.mzXML'
program='moldiscovery' ##dereplicator+ or varquest or metaminer or moldiscovery
db='npatlas' 
dbtools_options="--reuse --fdr --fdr-limit 1 --db-path npdtools/$db"
metaminer_options='--class all --keep-ga-files --correspondence data/QTOF/metaminer_prism/corresp.tsv --sequence data/extracted_seqs/ripps_prism/'
varquest_options='--accurate_p_value --extract-best --adduct_Na --adduct_K'
npdtoolsBin='npdtools2/bin'
general_options="--mode HH --fdr"
$npdtoolsBin/$program.py \
--debug $general_options $input_data -o data/QTOF/${program}_${db}_GB4-14/ $dbtools_options --threads 1 \
--preprocessed_targets data/QTOF/moldiscovery_npatlas/db_preproc/library.info.prob.mz.target.bin \
--preprocessed_decoys data/QTOF/moldiscovery_npatlas/db_preproc/library.info.prob.mz.decoy.bin
##sbatch --cpus-per-task=16 --ntasks=10 --partition=multinode --mem=200g src/npdtools.sh
##--threads $SLURM_CPUS_PER_TASK --preprocess