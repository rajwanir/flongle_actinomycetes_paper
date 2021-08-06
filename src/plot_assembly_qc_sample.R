library(tidyverse)
library(data.table)
library(RColorBrewer)

## functions #####

read_quast_data <- function(tsv_report,sample_name){
  tsv_report = fread(tsv_report)
  tsv_report = tsv_report %>% pivot_longer(starts_with("x"), names_to ="depth", values_to ="value")
  tsv_report$value = as.numeric(tsv_report$value)
  tsv_report$depth = tsv_report$depth %>% str_remove("x") %>% str_remove(".contigs") %>% as.numeric()
  tsv_report %>% mutate(sample = sample_name)
}

get_cluster_num <- function(antismash_json) {
  # ncol(jsonlite::stream_in(file(antismash_json), , pagesize = 100000)$timings)
  file = jsonlite::stream_in(file(antismash_json), , pagesize = 100000)
  clust_num  = lapply(file$records[[1]]$features, function(x) {nrow (x %>% filter(type == "region"))}) %>% unlist() %>% sum()
  return(clust_num)
  }


get_cluster_num_all <- function(sample, assembler = "medeka", quast_data = quast_data) {
  lapply(paste("x", unique(quast_data$depth), sep = ""), function(each_depth) {
    bgc_num = get_cluster_num(
      sprintf(
        "data/antismash/simulation/%s/%s/%s/%s.contigs.json",
        sample,
        each_depth,
        assembler,
        each_depth
      )
    )
    data.frame(
      Assembly = "# BGC",
      depth = each_depth %>% str_remove("x") %>% as.numeric(),
      value = as.numeric(bgc_num),
      sample = as.character(sample),
      assembler = assembler
    )

  })
}



#### metadata import ####

sample_metadata = read_csv("data/assm_stats.csv")
assm_stats = sample_metadata
flowcell_colors = colorRampPalette(brewer.pal(9,"Set1"))(length(unique(sample_metadata$flow_cell_id)))
names(flowcell_colors) = unique(sample_metadata$flow_cell_id)


# #### assm stats #####
# 
# # need to compute once 
# 
# sample_metadata = read_csv("data/metadata_assembled.csv")
# 
# # genomes_path = Sys.glob(file.path(
# #   #getwd(),
# #   "data/assembly/*/*/medeka/*.contigs.fasta"))
# # genomes_path = bind_cols(s = str_extract(genomes_path, pattern = "\\w+/barcode\\d+"),
# #                          genomes_path = genomes_path) %>% separate(col = "s",
# #                                                                    into = c("seq_date", "barcode"),
# #                                                                    sep = '/')
# # genomes_path$barcode = str_remove(genomes_path$barcode, "barcode0|barcode")
# # sample_metadata = merge(
# #   sample_metadata,
# #   genomes_path,
# #   by.x = c("barcode", "seq_date"),
# #   by.y = c("barcode", "seq_date"),
# #   all.x = F,
# #   all.y = T
# # )
# 
# 
# assm_stats_prog = "/home/rajwanir2/assembly-stats/build/assembly-stats -t"
# genomes_path_expr = "data/assembly/*/*/medeka/*.contigs.fasta"
# assm_stats = system(paste(assm_stats_prog,genomes_path_expr),intern = T)
# assm_stats = fread(text = assm_stats)
# 
# assm_stats = merge(sample_metadata,assm_stats,
#       by.x="genomes_path",by.y = "filename",
#       all.x = F, all.y = T)
# 
# assm_stats =assm_stats[order(assm_stats$number),]
# 
# 
# 
# 
# 

#### coverage ####

coverage = Sys.glob(file.path("data/assembly/*/*/canu/*.contigs.layout.tigInfo"))
names(coverage) = coverage
coverage = lapply(coverage,data.table::fread)
coverage = bind_rows(coverage,.id="path")
coverage$barcode = str_extract(coverage$path, pattern = "barcode\\d+") %>% str_remove(pattern = "barcode") %>% as.integer()
coverage$date = str_extract(coverage$path, pattern = "(\\d{8})")
coverage = coverage %>% filter(tigClass != "unassm")
coverage = coverage %>% group_by(path) %>% summarize(coverage = mean(coverage),
                                          barcode = unique(barcode),
                                          date = unique(date)) 

assm_stats = merge(assm_stats,coverage,all.x = T,all.y = F,
                   by.x = c("barcode","seq_date"), by.y = c("barcode","date"))


#### BGC count ####
bgcs = Sys.glob(file.path(getwd(), "data/antismash/*/barcode*/medeka/*.gff.1"))
names(bgcs) = bgcs
tags_to_load = c("product","type","gene","specificity","aSDomain","gene_kind","locus_tag","gene_functions")
bgcs = lapply(bgcs,function(x)
  rtracklayer::readGFF(x,tags=tags_to_load))
bgcs = bind_rows(bgcs,.id="path")
# bgcs$gene_functions = bgcs$gene_functions %>% str_extract("(.*):(.*)\\(") %>% str_remove_all(pattern = "(.*)smcogs\\) |\\(") 

bgcs = bgcs %>% filter(type...3 == "region") %>% group_by(path) %>% count(name = "bgc_count")
bgcs$barcode = str_extract(bgcs$path, pattern = "barcode\\d+") %>% str_remove(pattern = "barcode") %>% as.integer()
bgcs$date = str_extract(bgcs$path, pattern = "(\\d{8})")

assm_stats = merge(assm_stats,bgcs,all.x = T,all.y = F,
                   by.x = c("barcode","seq_date"), by.y = c("barcode","date"))



### ploting #############

assm_stats = assm_stats %>% pivot_longer(cols = c("coverage","total_length","number","bgc_count"), names_to = "variable", values_to = "value") 

# p_assm_stats = ggplot(data = assm_stats,
#        aes(x = coverage,y = value, color = flow_cell_id)) + 
#   geom_point() +
#   ggrepel::geom_text_repel(aes(label = sample)) +
#   facet_wrap(~variable, scales = "free_y", ncol =1,
#              labeller = as_labeller(c(bgc_count = "Number of BGCs", number = "Number of contigs", total_length = "Assembly size"))) +
#   theme_bw() + theme(axis.title.x = element_text(face="bold"), legend.position = "bottom", strip.background = element_blank(), strip.text = element_text(face = "bold", size = 10), legend.background = element_rect(color = "black")) +
#   labs(y=NULL, x = "Coverage") +
#   guides(color = guide_legend(title = "flow cell", title.position = "left", title.hjust = 0.5, nrow = 1)) +
#   scale_color_manual(values = flowcell_colors)


assm_stats = assm_stats %>% arrange(flow_cell_id) %>% 
              mutate(sample_id = factor(paste(sample,flow_cell_id), 
                                        levels = unique(paste(sample,flow_cell_id)),
                                        ordered = T ),
                     variable = factor(variable,
                                       levels = sort(unique(assm_stats$variable),decreasing = T),
                                       ordered = T))

p_assm_stats = ggplot(data = assm_stats,
       aes(x = sample_id, y= value, fill = flow_cell_id)) + 
  geom_bar(stat="identity",color = "black") +
  facet_wrap(~variable, scales = "free_x", ncol =4,
             labeller = as_labeller(c(bgc_count = "Number of BGCs", number = "Number of contigs", total_length = "Assembly size", coverage = "Coverage"))) +
  scale_x_discrete(labels = str_extract(levels(assm_stats$sample_id),pattern = '.* ') %>% str_remove(pattern = ' ') )+
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90), legend.position = "top", strip.background = element_blank(), strip.text = element_text(face = "bold", size = 10), legend.background = element_rect(color = "black")) +
  labs(y=NULL, x = NULL) +
  guides(fill = guide_legend(title = "Flow cell :", title.position = "left", title.hjust = 0.5, nrow = 1)) +
  scale_fill_manual(values = flowcell_colors) +
  coord_flip()



ggsave("figures/assembly_stats_sample.svg",p_assm_stats,
       height = 9, width = 8)


# bgc_data = data.frame(Assembly = as.character(),depth = as.numeric(),value=as.numeric(),sample = as.numeric())


