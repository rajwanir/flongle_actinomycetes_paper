#! /bin/bash
reads_folder=$1
source /gpfs/gsfs11/users/rajwanir2/conda/etc/profile.d/conda.sh
conda activate bonito
mkdir -p $reads_folder/bonito
bonito basecaller dna_r9.4.1 $reads_folder > $reads_folder/bonito/bonito.fasta 
module load seqtk
seqtk seq -F '#' $reads_folder/bonito/bonito.fasta | gzip -c > $reads_folder/bonito/bonito.fastq.gz ## adding fake quality scores for conversion
bash src/qcat.sh $reads_folder/bonito
# bonito has the  --fastq flag but currently not complete. It will outputs constant quality scores.
#sbatch --partition=gpu --cpus-per-task=14 --mem=16g --gres=lscratch:200,gpu:p100:1 src/bonito.s /gpfs/gsfs11/users/rajwanir2/soil_metagenomics/data/rawdata/vanRsoilactino_10132020/20201013_1905_MN31218_AEZ324_3a6e0404/fast5