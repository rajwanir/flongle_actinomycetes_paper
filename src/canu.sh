#!/bin/bash
barcode=$1
date=$2
##input_file="data/rawdata/vanRsoilactino_${date}/basecalled_fastq/demultiplexed/$barcode.fastq"
input_file="data/rawdata/${date}/$barcode.fastq"
output_folder="data/assembly/$date/$barcode/canu"
module load canu
canu -p $barcode \
-d $output_folder \
-genomeSize=8m \
-nanopore-raw $input_file \
minReadLength=1000 minOverlapLength=400 \
corMhapSensitivity=high  \
usegrid=1 gridOptions="--partition norm" gridOptionsJobName=$barcode \
stopOnLowCoverage=2 \
minInputCoverage=2
#stopOnReadQuality=false  
## stageDirectory=/lscratch/$SLURM_JOBID gridEngineStageOption="--gres=lscratch:100"
##corMinCoverage=0