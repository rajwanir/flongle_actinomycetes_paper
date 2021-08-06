#!/bin/bash
module load qcat
FASTQ_DIR=$1
zcat $FASTQ_DIR/*.fastq.gz |qcat -b $FASTQ_DIR/demultiplexed2/ \
-k NBD103/NBD104 \
--detect-middle --trim --tsv > $FASTQ_DIR/demultiplexed.tsv
#--min-read-length 1000 