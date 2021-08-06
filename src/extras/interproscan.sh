#!/bin/bash
folder='data/rawdata/simulation/3'
input_file="$folder/3.fastq.gz"
orignal_reads=203702
orginal_bases=887484358
orignal_coverage=110
seqtk sample -s100 $input_file $(echo $(echo 60/$orignal_coverage|bc -l)*$orignal_reads|bc -l) > $folder/x60.fq
seqtk sample -s100 $input_file $(echo $(echo 30/$orignal_coverage|bc -l)*$orignal_reads|bc -l) > $folder/x30.fq
seqtk sample -s100 $input_file $(echo $(echo 15/$orignal_coverage|bc -l)*$orignal_reads|bc -l) > $folder/x15.fq
seqtk sample -s100 $input_file $(echo $(echo 7/$orignal_coverage|bc -l)*$orignal_reads|bc -l) > $folder/x7.fq
seqtk sample -s100 $input_file $(echo $(echo 3/$orignal_coverage|bc -l)*$orignal_reads|bc -l) > $folder/x3.fq