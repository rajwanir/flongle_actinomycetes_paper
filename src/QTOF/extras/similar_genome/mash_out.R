library(tidyverse)
library(data.table)
library(biomartr)

mash_out = fread("Y:/soil_metagenomics/data/similar_genome/similar_genomes_mash.tsv")
colnames(mash_out) = c("ref","query","ident","pval","sketches")
mash_out = separate(mash_out,col = "sketches", 
                    into = c("sketch_of","sketch_t"),
                    sep = '/',
                    convert =T)
mash_out = separate(mash_out,col = "query", 
                      into = c("query","q_comment"),
                    sep = ':')

mash_out = separate(mash_out,col = "ref", 
                    into = c("ref","r_comment"),
                    sep = ':')
mash_out = mash_out %>% group_by(query) %>%
  arrange(., -sketch_of) %>% slice(1)
mash_out$ref_acc = str_remove(mash_out$ref,
                              "_ASM(.)+")

mash_out$ref_path = lapply(mash_out$ref_acc, function(x) getGenome(db = "refseq", 
                                           organism = x, 
                                           path = "data/similar_genome/ref_genomes",
                                           reference = FALSE,
                                           gunzip = F))
mash_out$ref_path = as.character(mash_out$ref_path)
write.csv(mash_out,"mash_out.csv",row.names = F)


#confirm ani
mash_out = read.csv("data/similar_genome/mash_out.csv")

mash_out$cmds = sprintf("./fastANI -t 8 -q %s -r data/similar_genome/ref_genomes/fasta/%s* -o /dev/stdout | cut -f3",
        mash_out$query,mash_out$ref_acc)

mash_out$ani = sapply(mash_out$cmds,function(x) as.numeric(system(x,intern = T)),
       USE.NAMES = F)

write.csv(mash_out,"data/similar_genome/mash_out.csv",row.names = F)