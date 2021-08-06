library(tidyverse)

sample_metadata = read.csv("data/metadata.csv", colClasses = rep("character",13) )
sample_metadata$barcode = sprintf("barcode%02d",as.numeric(sample_metadata$barcode))
assm_stats = read.csv("data/quast/assm_stats_samples.csv")
assm_stats = select(assm_stats,genomes_path = filename, n_contigs = number,assm_size=total_length)

#### 16s ######################
genomes_path = Sys.glob(file.path(
  #getwd(),
  "data/assembly/*/*/medeka/*.contigs.fasta"))
blast_16s = lapply(assm_stats$genomes_path,function(x)
            read.table(text = system(
              sprintf("blastn -num_threads 16 -max_target_seqs 1 -db /fdb/blastdb/v4/16SMicrobial -query %s -outfmt '6 delim=@ staxids stitle pident' | head -n1",x),
              intern = T)
              ,sep = '@', col.names = c("taxid","title","p_ident")
              ))

blast_16s = bind_rows(blast_16s)
blast_16s = separate(blast_16s, col = title, into = c("genus","species", "strain"), sep = ' ')

blast_16s = cbind(assm_stats,blast_16s)

#alternate read pre-processed
# blast_16s = read_csv("tables/16s_taxonomy.csv")

##### Alignments ####

#convert cluster gbk to fasta - single file per genome
genomes = data.frame(
  date = assm_stats$genomes_path %>% str_extract('\\d{8}'),
  barcode = assm_stats$genomes_path %>% str_extract('barcode\\d+')
)
genomes = genomes %>% mutate(gb2fasta_cmd = sprintf('sed -f src/gb2fasta.sed data/antismash/%s/%s/medeka/*region*.gbk > data/bgcs_fasta/%s.fasta',
                   date,barcode,paste(date,barcode,sep = '-')),
                   date_bc = paste(date,barcode,sep = '-'))


# genomes = genomes %>% mutate(copy_gbk_cmd = sprintf('cp data/antismash/%s/%s/medeka/*.contigs.gbk data/bgcs_gbk/%s.gbk',
#                                                     date,barcode,date_bc_sample))
# 
# 
# lapply(genomes$copy_gbk_cmd,system)

genomes=cbind(genomes,blast_16s)

genomes = merge(genomes,sample_metadata,
      by.x = c("date","barcode"),
      by.y = c("seq_date","barcode"),
      all.x = T, all.y = F)


#alignments

minimap_cols = c("q_seq_name", "q_seq_len", "q_start_coord", "q_end_coord",
                "strand",
                "t_seq_name", "t_seq_len", "t_start_coord", "t_end_coord",
                "n_matches","n_bases","map_q", "other")
minimap_colclasses = c(rep("character",13))

perms = gtools::permutations(n = length(genomes$date_bc), v = genomes$date_bc, r =2) %>% 
        as_tibble()
colnames(perms) = c("target","query")

alignments = parallel::mclapply(1:nrow(perms),function(idx){
  system(sprintf("minimap2 -v1 -x asm5 data/bgcs_fasta/%s.fasta data/bgcs_fasta/%s.fasta --paf-no-hit ",
                 perms$target[idx],perms$query[idx]),
       intern = T) %>% 
  read.table(text = ., sep='\t', header = F, fill = T, 
             #col.names = minimap_cols, colClasses = minimap_colclasses
             ) %>% 
  mutate(target = perms$target[idx], query = perms$query[idx]) %>% 
  mutate(t_genus = genomes[genomes$date_bc == perms$target[idx],]$genus,
         t_species = genomes[genomes$date_bc == perms$target[idx],]$species,
         q_genus = genomes[genomes$date_bc == perms$query[idx],]$genus,
         q_species = genomes[genomes$date_bc == perms$query[idx],]$species,
         t_sample_name = genomes[genomes$date_bc == perms$target[idx],]$sample,
         q_sample_name = genomes[genomes$date_bc == perms$query[idx],]$sample,
         q_n_contigs = genomes[genomes$date_bc == perms$query[idx],]$n_contigs,
         t_n_contigs = genomes[genomes$date_bc == perms$target[idx],]$n_contigs)

  },mc.cores=30) 

alignments = lapply(alignments,function(x) x %>% mutate_if(~ !is.character(.),as.character))
alignments = alignments %>% bind_rows()
colnames(alignments) = c(minimap_cols,sprintf("V%d",14:18),colnames(alignments)[19:length(colnames(alignments))])

#remove same sample comparision
alignments = alignments %>% filter(!q_sample_name == t_sample_name)
#remove contigs > 50
alignments = alignments %>% filter(q_n_contigs < 50 & t_n_contigs < 50)



#create an alignment summary by counting number of unique BGCs
alignments_summary = filter(alignments,strand == "*") %>% group_by(target,query) %>% 
  summarize(n_unique_bgcs = n(),
            t_genus = unique(t_genus),q_genus = unique(q_genus),
            t_species = unique(t_species),q_species = unique(q_species))

#plot same genus and same species  
p_same_species = ggplot(alignments_summary %>% filter(t_species == q_species),
       mapping = aes(y=n_unique_bgcs, x = paste(t_genus,t_species,sep = '\n' ) )) + 
  geom_boxplot(fill="burlywood") + geom_point(position = position_dodge(width=0.75)) + ggprism::theme_prism() + 
  labs(title = "Within species",x= "\nSpecies", y = "# of unique BGCs") +
  theme(plot.title = element_text(hjust=0.5, face = "bold"), panel.grid.minor = element_blank())



p_same_genus = ggplot(alignments_summary %>% filter( (t_genus == q_genus) & (!t_species == q_species) ),
                      mapping = aes(y=n_unique_bgcs, x = t_genus )) + 
  geom_boxplot(fill="brown4") + geom_point(position = position_dodge(width=0.75)) + ggprism::theme_prism() + 
  labs(title = "Between species",x= "\nGenus", y = "# of unique BGCs") +
  theme(plot.title  = element_text(hjust=0.5, face = "bold"), panel.grid.minor = element_blank())

p_biosyn_diversity = cowplot::plot_grid(p_same_species,p_same_genus, 
                                        ncol = 2, labels = c("A","B"))

ggsave(p_biosyn_diversity,filename = "figures/biosyn_diversity.svg", height=6,width = 14)

## BGCs ####

# read all gffs

read_all_gffs = function(){
samples = Sys.glob(file.path(getwd(), "data/antismash/*/*/medeka/*.gff.1"))
names(samples) = samples
tags_to_load = c("product","type")
samples = lapply(samples,function(x)
  rtracklayer::readGFF(x,tags=tags_to_load))
samples = bind_rows(samples,.id="path") %>% filter(type...3 == "region") 
return(samples)
}

bgcs = read_all_gffs() 
bgcs = mutate(bgcs,
                 date = path %>% str_extract('\\d{8}'),
                 barcode = path %>% str_extract('barcode\\d+'),
                 date_bc = paste(date,barcode,sep = '-')
                 )

bgcs = merge(bgcs,genomes,
      by.x = "date_bc", by.y = "date_bc",
      all.x = T, all.y = T)

p_bgc__type_tiles = ggplot(bgcs) +
  geom_tile(mapping  = aes(x = product, y = paste(sample, flow_cell_id)),
            fill = "orange", color = "black") +
  scale_fill_discrete(na.value = "white") +
  labs(y=NULL, x = "BGC type") +
  theme_linedraw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90),
        plot.margin = margin(t=0, unit = "pt")
        ) 

p_bgc_type_bar = ggplot(bgcs, mapping  = aes(x = product)) + 
  geom_bar(stat = "count", fill = "orange", color = "black") +
  geom_text(stat='count', aes(label=..count..), vjust=-1) +
  scale_y_continuous(expand = c(0,30)) +
  theme_void() +
  theme(plot.margin = margin(l = 3, b = -0.4, unit = "cm"))

p_bgc_types = gridExtra::grid.arrange(p_bgc_type_bar,p_bgc__type_tiles, 
                        nrow =2,ncol = 1,
                        heights = c(20,60),
                        as.table =F)

ggsave(p_bgc_types, filename = "figures/bgc_types.png",
       height = 5.5, width = 8, dpi = 300)