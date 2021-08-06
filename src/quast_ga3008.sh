#! /bin/bash
module load quast
np_folder='data/quast/GA3-008/input/minion'
fl_folder='data/quast/GA3-008/input/flongle'
quast.py --threads 16 \
$fl_folder/*.fasta \
-R $np_folder/*.fasta \
--output-dir data/quast/GA3-008/output \
--features $np_folder/*.gff \
--plots-format png \
--circos
