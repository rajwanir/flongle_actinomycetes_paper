#! /bin/bash
for i in 500 1000 2000 4000 8000; do
j=$(expr $i + $i)
badread simulate --reference data/pacbio_reference/GB4-014.fasta --quantity 10x --length $i,$j 2> data/rawdata/readlen/$i.log 1> data/rawdata/readlen/$i.fastq & done
