#!/bin/bash
set -e
module load guppy || exit 1
fast5_folder=$1
guppy_basecaller --input_path $fast5_folder/ --flowcell FLO-FLG001 --kit SQK-LSK109 \
--save_path $fast5_folder/basecalled_fastq/ \
--min_qscore 7 \
-x cuda:all \
--compress_fastq \
--records_per_fastq 0 \
--recursive