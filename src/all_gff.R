library(tidyverse)
library(gggenes)
library(RColorBrewer)

target_type = "amglyccycl"

#all data import
samples = Sys.glob(file.path(getwd(), "data/antismash/*/*/medeka/*.gff.1"))
names(samples) = samples
tags_to_load = c("product","type","gene","specificity","aSDomain","gene_kind","locus_tag","gene_functions")
samples = lapply(samples,function(x)
  rtracklayer::readGFF(x,tags=tags_to_load))
samples = bind_rows(samples,.id="path")
samples$gene_functions = samples$gene_functions %>% str_extract("(.*):(.*)\\(") %>% str_remove_all(pattern = "(.*)smcogs\\) |\\(") 

# select cluster
target_type = samples %>% filter(product == target_type) %>% select(path,seqid,start,end) %>% distinct(path,.keep_all = T)

target_type = lapply(seq(1,nrow(target_type)),
       function(x)
         samples %>% filter(path == target_type$path[x] &
                              start > target_type$start[x] &
                              end < target_type$end[x] &
                              seqid == target_type$seqid[x])
       ) %>% bind_rows()

target_type = target_type %>% mutate(id_sim = target_type$path %>% str_extract("\\d+/barcode\\d+"))


## homologous clusters

homolgous = Sys.glob(file.path(getwd(), "data/antismash/aminoglycoside/*/*.gff.1"))
names(homolgous) = homolgous
homolgous = lapply(homolgous,function(x)
  rtracklayer::readGFF(x,tags=tags_to_load))
homolgous = bind_rows(homolgous,.id="path")
homolgous$gene_functions = homolgous$gene_functions %>% str_extract("(.*):(.*)\\(") %>% str_remove_all(pattern = "(.*)smcogs\\) |\\(") 


##merge with oringal
target_type = bind_rows(target_type,homolgous)  
target_type = target_type %>% mutate(id_sim = case_when(
  is.na(id_sim) ~ as.character(seqid), TRUE ~ as.character(id_sim)
  ))  

##specific filter
target_type = target_type %>% filter((seqid == "BGC0000283" & start > 17799 & end < 39144)| (seqid != "BGC0000283" ))
target_type = target_type %>% mutate(sample_name = case_when(
  id_sim == "10222020/barcode10" ~ "GB15-009",
  id_sim == "12022020/barcode08" ~ "GA10-003",
  id_sim == "12082020/barcode01" ~ "GA10-004",
  id_sim == "BGC0000283" ~ "cetoniacytoneA",
  id_sim == "BGC0001538" ~ "caboxamycin "))

# #start all sequences with 1
# lapply(unique(target_type$id_sim), function(cluster_id) {
# GenomicRanges::makeGRangesFromDataFrame(target_type %>% filter(id_sim == cluster_id)
#                                       , keep.extra.columns = T) %>%
#     GenomicRanges::shift(x = .,shift = -start(.)[1] ) %>% .[1]
# })

p = ggplot(data = target_type %>% filter(type...4=="CDS"),
           aes(xmin=start,xmax=end, y=sample_name,
               fill = gene_functions)) +
  geom_gene_arrow(size = 0.8) + facet_wrap(~id_sim, scales = "free", ncol = 1) +
  geom_gene_label(aes(label = gene), align = "left", fontface="bold", min.size =6) +
  scale_fill_brewer(palette = "Set3", na.value = "dimgray" ) +
  theme_genes() + labs(y=NULL) +
  theme(axis.text.y = element_text(size = 14,color = "black"),legend.position = "bottom", legend.background = element_rect(color = "black",size =1), panel.grid.major.y = element_line(color="black")) +
  guides(fill = guide_legend(nrow=4,title = "Gene function", title.position = "top", title.hjust = 0.5, title.theme = element_text(face="bold")))


ggsave(p,filename = "figures/aminoglycoside.png",
       height = 4, width =14, dpi = 300)