library(tidyverse)
library(gggenes)
library(ggsci)


## samples and paths input ####
# date = "09182020"
# genomes_to_check = c("barcode08", "barcode07", "barcode10")
# genomes_path = lapply(genomes_to_check, function (x)
#   sprintf("data/assembly/09182020/%s/medeka/%s.contigs.fasta", x, x))

sample_metadata = read.csv(
  "data/metadata.csv",
  colClasses = c(
    "integer",
    "character",
    rep("numeric", 5),
    "integer",
    "character",
    "character",
    "logical",
    "character",
    "character"
  )
)
genomes_path = Sys.glob(file.path(getwd(), "data/assembly/*/*/medeka/*.contigs.fasta"))
genomes_path = bind_cols(s = str_extract(genomes_path, pattern = "\\w+/barcode\\d+"),
                         genomes_path = genomes_path) %>% separate(col = "s",
                                                                   into = c("seq_date", "barcode"),
                                                                   sep = '/')
genomes_path$barcode = str_remove(genomes_path$barcode, "barcode0|barcode")
sample_metadata = merge(
  sample_metadata,
  genomes_path,
  by.x = c("barcode", "seq_date"),
  by.y = c("barcode", "seq_date"),
  all.x = F,
  all.y = T
)


## blast search ######
query_path = Sys.glob(file.path(getwd(), "GPA_evolution/select_genes/split/*.fasta"))

get_blast_hits = function(q){
blast_cmds = lapply(sample_metadata$genomes_path, function(x)
  sprintf(
    "tblastn -query %s -subject %s -max_target_seqs 1 -outfmt \'6 delim=@  qseqid pident qlen length sseqid sstart send evalue score\' | head -n1",
    q,x
  ))
blast_results = lapply(blast_cmds, function(x)
  system(x, intern = T))
  
blast_results = lapply(blast_results, function(x)
  read.table(text = x, sep = '@')) %>%  bind_rows() %>% 
  rename(target = V1, pident = V2,query_len = V3, alignment_length = V4,
         contig = V5, start = V6, end = V7, evalue = V8, score = V9) %>% 
  bind_cols(.,sample_metadata)

return(blast_results)
}

blast_results = lapply(query_path,get_blast_hits)

blast_results = blast_results %>% bind_rows() %>% separate(col = "target", into = c("gene","source_bgc"), sep = "_")
blast_results = blast_results %>% mutate(
  contig_clean = str_remove(blast_results$contig,pattern = "tig(0)\\1+") %>% str_remove(pattern =  "_segment\\d+")  %>%  str_pad(., max(nchar(.)), side = "left", pad = 0)
)
blast_results = blast_results %>% mutate(sample_clean = blast_results$sample %>% str_pad(., max(nchar(.)), side = "right", pad = ' '))
  
blast_results$contig_full = paste(blast_results$sample_clean,blast_results$contig_clean, sep = '\t') 
blast_results = blast_results[order(blast_results$contig_full),] 

# blast_results_filtered = blast_results %>% filter(evalue < 1)

p_blast_hits = ggplot(blast_results, aes(xmin = start, xmax = end, 
                                   y = contig_full, alpha = pident, fill = gene)) +
  geom_gene_arrow(arrowhead_height = unit(20, "mm"), arrowhead_width = unit(20, "mm"),  arrow_body_height = unit(15, "mm")) +
  facet_wrap(drop = T,~contig_full, scales = "free",ncol = 2, dir = "v") +
  theme_genes() + 
  labs(title = "Targeted BLAST screen - Glycopeptide") +
  theme(#axis.text.y = element_blank(),
        #axis.text.x = element_blank(),
       # strip.text.y.right = element_text(angle = 0),
        #strip.text = element_blank(), 
        axis.text.y = element_text(face = "bold", size = 30, hjust = 0),
        axis.title.y = element_blank(),
        legend.position="top",legend.text=element_text(face="bold",size=30),
        text = element_text(size = 30),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  guides(fill = guide_legend(nrow = 1),override.aes = list(size = 5)) +
  scale_y_discrete(expand = c(0.01,0.01)) +
  scale_fill_jco() + scale_alpha(range =c(0.5,1))

ggsave(p_blast_hits, filename =  "figures/p_blast_hits.png", height = 25, width = 49)


## create tables ####
create_table <- function(blast_result){
table = condformat(blast_result) %>% rule_fill_discrete(sample,
                                                 expression = percent_identity > pident_highight & coverage > cov_highlight,
                                                 colours  = c("TRUE" = "lightgreen")) %>% condformat2grob()
return(table)}

t_blast_results = lapply(genomes_to_check,
       function(x){
         create_table(blast_results %>% filter(sample == x))
         })

png("figures/blast_hits.png", width = 13, height = 5, units = "in", res = 300)
grid.arrange(gtable_combine(t_blast_results[[1]]),
             gtable_combine(t_blast_results[[2]]),
             gtable_combine(t_blast_results[[3]]),
             ncol = length(genomes_to_check)) 
dev.off()


# t_coverage = blast_results %>% select(-percent_identity) %>% pivot_wider(names_from = sample, values_from =
#                                                                            coverage)
# t_percent_identity = blast_results %>% select(-coverage) %>% pivot_wider(names_from = sample, values_from =
#                                                                            percent_identity)
# ## sepearate table for coverage
# p_coverage = condformat(t_coverage) %>% rule_fill_discrete(barcode08,
#                                                  expression = barcode08 > 80,
#                                                  colours  = c("TRUE" = "lightgreen")) %>%
#   rule_fill_discrete(barcode07,
#                      expression = barcode07 > 80,
#                      colours  = c("TRUE" = "lightgreen")) %>%
#   rule_fill_discrete(barcode10,
#                      expression = barcode10 > 80,
#                      colours  = c("TRUE" = "lightgreen")) %>%  condformat2grob()
# 
# ## separate table for identity
# p_percent_identity = condformat(t_percent_identity) %>% rule_fill_discrete(barcode08,
#                                                                            expression = barcode08 > 50,
#                                                                            colours  = c("TRUE" = "lightgreen")) %>%
#   rule_fill_discrete(barcode07,
#                      expression = barcode07 > 50,
#                      colours  = c("TRUE" = "lightgreen")) %>%
#   rule_fill_discrete(barcode10,
#                      expression = barcode10 > 50,
#                      colours  = c("TRUE" = "lightgreen")) %>%  condformat2grob()
# 
# grid.arrange(gtable_combine(p_percent_identity,p_coverage))
