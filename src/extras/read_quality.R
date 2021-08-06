library(tidyverse)
# library(nanoR)

guppy_stats = read.csv("data/rawdata/vanRsoilactino_10222020/basecalled_fastq/sequencing_summary.txt", 
         sep = '\t')

ggplot(data=guppy_stats %>% filter(mean_qscore_template<7),
       aes(x=sequence_length_template))+
  geom_histogram(binwidth = 1) +
  scale_x_continuous(breaks = 1:15)


guppy_stats %>% group_by(channel) %>% summarize(output = sum(sequence_length_template)) 