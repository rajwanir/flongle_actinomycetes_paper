#!/bin/bash
set -e
module load guppy || exit 1
fastq_folder=$1
guppy_barcoder --input_path $fastq_folder/ --barcode_kits EXP-NBD104 \
--save_path $fastq_folder/demultiplexed \
--compress_fastq \
--records_per_fastq 0 \
--trim_barcodes \
--detect_mid_strand_barcodes \
-x cuda:all