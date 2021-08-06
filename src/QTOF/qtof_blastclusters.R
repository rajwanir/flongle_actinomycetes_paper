library(tidyverse)


## blast ##############
  date='12022020'
  barcode='barcode05'
  compound='HerbicolinB'
  blast_cmd = sprintf('tblastn -subject data/assembly/%s/%s/medeka/%s.contigs.fasta -query data/QTOF/dereplicator+_npatlas/known_bgcs/%s/protein.fasta -outfmt \'6 delim=@ qseqid sseqid pident length sstart\' -max_target_seqs 1',
                      date,barcode,barcode,compound)
  blast_out = read.delim(text = system(blast_cmd, intern = T), header = F,
                         col.names = c("qseqid", "sseqid", "pident", "length", "sstart"),
                         sep = '@')
  blast_out
  
  p = ggplot(blast_out) +
    geom_histogram(mapping = aes(x=sstart, fill = qseqid)) + theme_classic() +
    theme(legend.position = "bottom")
  
  p + scale_x_continuous(breaks=seq(1,max(blast_out$sstart),100000)) + theme(axis.text.x = element_text(angle=90))

plotly::ggplotly(p)

### nrps A domains ####

# json_filenames = "data/QTOF/dereplicator+_npatlas/known_bgcs/aureobasidin/antismash_QuerySubset/query_subset.json"
json_filenames = "data/antismash/12082020/barcode03/medeka/barcode03.contigs.json"
nrps_data = jsonlite::stream_in(file(json_filenames),,pagesize = 100000)$records[[1]]$modules$antismash.modules.nrps_pks$domain_predictions
mod_genes = str_subset(colnames(nrps_data),"nMT") %>% str_extract("ctg\\d+_\\d+")
nrps_data = nrps_data %>% select(contains("AMP-binding"))
adomain_genes=colnames(nrps_data) %>% str_extract("ctg\\d+_\\d+")
nrps_data = lapply(1:ncol(nrps_data),function(idx) nrps_data[idx]) %>% sapply(.,'[[',1) %>% bind_rows() %>% filter(!is.na(uncertain))
nrps_data = mutate(nrps_data,across(where(is.list), as.character))
nrps_data$adomain_genes = adomain_genes
nrps_data = mutate(nrps_data,methylation = adomain_genes %in% mod_genes)


query_aa = c("gln","phe","ile","ahp","tyr","thr")
nrps_data$query_aa = nrps_data$stachelhaus_predictions %in% query_aa

p =  ggplot(nrps_data,aes(x=adomain_genes,fill=stachelhaus_predictions) ) +
  geom_bar(stat = "count",color = "black") +
  geom_text(aes(label = stachelhaus_predictions, color = query_aa),stat="count",
            position = position_stack(vjust = .5), fontface = "bold") +
  labs(y= NULL, x= NULL) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"))+
  theme_classic() + theme(legend.position = "none") +
  coord_flip()


