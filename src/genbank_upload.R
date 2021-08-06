library(tidyverse)

assm_stats = read_csv("data/metadata_assembled.csv")

assm_stats = assm_stats %>% mutate(fastq_path = sprintf("data/rawdata/vanRsoilactino_%s/basecalled_fastq/demultiplexed/barcode%02d.fastq",seq_date,barcode),
                                   paste_path = sprintf("data/genbank_upload/%s_%02d_%s.fastq",seq_date,barcode,sample))

parallel::mclapply(1:nrow(assm_stats),function(idx) 
  file.copy(assm_stats$fastq_path[idx],assm_stats$paste_path[idx]),
  mc.cores = 30)