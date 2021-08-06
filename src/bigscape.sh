#! /bin/bash
target_cluster="gpa"
input_folder="data/antismash/$target_cluster"
source /gpfs/gsfs11/users/rajwanir2/conda/etc/profile.d/conda.sh
conda activate bigscape
python bigscape/bigscape.py \
--inputdir $input_folder \
--outputdir data/bigscape/$target_cluster \
--cutoffs 1 \
--clans-off \
--exclude_gbk_str region ## remove this line if aminocoumarin