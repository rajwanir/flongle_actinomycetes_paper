library(tidyverse)
library(condformat)
library(grid)
library(gridExtra)


pident_highight = 40
cov_highlight = 70

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
blast_cmds = lapply(sample_metadata$genomes_path, function(x)
  sprintf(
    "tblastn -query GPA_evolution/select_genes/selected_genes.fasta -subject %s -max_target_seqs 1 -outfmt \'6 delim=@  qseqid pident qlen length\'",
    x
  ))
blast_results = lapply(blast_cmds, function(x)
  system(x, intern = T) %>% read.table(sep = '@'))
  
  %>%  bind_cols(.,))  
# blast_results = lapply(blast_results, function(x)
#   read.table(text = x, sep = '@'))

for (s in 1:length(sample_metadata$sample)) {
  blast_results[[s]] <-
    blast_results[[s]] %>% mutate(sample = sample_metadata$sample[s])
}

blast_results = blast_results %>% bind_rows()
colnames(blast_results) <-
  c("gene",
    "percent_identity",
    "query_length",
    "alignment_length",
    "sample")
blast_results$gene = blast_results$gene %>% str_remove("_\\w+")
blast_results = blast_results %>%  distinct(gene, sample, .keep_all = T) # only keeps the first hit
blast_results = blast_results %>% mutate(coverage = alignment_length / query_length * 100) %>% select(-c(alignment_length, query_length))
blast_results = blast_results %>% mutate(across(is.numeric, ~ round(., 1)))

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
