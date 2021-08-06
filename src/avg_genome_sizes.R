library(tidyverse)
library(data.table)
library(cowplot)

img_actino_data = fread("Y:/soil_metagenomics/data/misc_survey/taxontable95359_29-dec-2020.xls",
                           sep = '\t', header = T) %>% janitor::clean_names()

#keep only isolate genomes
img_actino_data =  filter(img_actino_data,
            gold_analysis_project_type == "Genome Analysis (Isolate)")

#Technology, number of contigs and number of genomes

img_actino_data = img_actino_data %>% mutate(sequencing_technology = case_when(
  sequencing_method %like% "PacBio" ~ "PacBio",
  sequencing_method %like% "Nanopore" ~ "Nanopore",
  sequencing_method %like% "Illumina" ~ "Illumina",
  TRUE ~ "Others (454/Sanger)"
))

#excluding Nanopore genomes
#only 10 and BGC count unknown
img_actino_data =  filter(img_actino_data,
                          sequencing_technology != "Nanopore")

#convert genome scale to Mb
img_actino_data = img_actino_data %>% filter(genome_size_assembled > 1000000 )


#taxa filtering
  taxa_to_filter=c("Streptomycetales", "Pseudonocardiales", "Corynebacteriales")
  taxa_to_filter=c("Streptomyces", "Amycolatopsis", "Nocardia")
  
  img_actino_data = img_actino_data %>% filter(genus %in% taxa_to_filter)

  img_actino_data = img_actino_data %>% filter(!is.na(biosynthetic_cluster_count))



p_genome_size = ggplot(img_actino_data,mapping = aes(x=fct_reorder(genus, genome_size_assembled),y=genome_size_assembled/1000000)) +
  geom_boxplot(outlier.shape = NA) +
  #geom_jitter(alpha=0.5,size=0.1) +
  scale_y_continuous(breaks = seq(1,15,by = 2)) +
  labs(y="Genome size (Mb)", x="Genus") +
 # facet_wrap(~order, scales = "free_y", drop = T, ncol = 3) +
  theme_classic() +
  theme(strip.background =element_blank(),
        text = element_text(size=20)) 
  #coord_flip()

# ggsave("figures/avg_genome_sizes.png",p_genome_size, limitsize = F, width = 25, height = 10)




p_bgc_count = ggplot(img_actino_data %>% filter(biosynthetic_cluster_count>5),mapping = aes(x=fct_reorder(genus, biosynthetic_cluster_count),y= biosynthetic_cluster_count)) +
  geom_boxplot(outlier.shape = NA) +
 # geom_jitter(alpha=0.5,size=0.1) +
#  scale_y_continuous(breaks = seq(1,50,by = 1)) +
  labs(y="# of Secondary \n metabolite gene clusters", x="Genus") +
 # theme(axis.text.x = element_text(angle = 90, hjust = 1),panel.grid.minor.y =element_blank()) +
  ylim(0,50) +
  theme_classic() +
  theme(strip.background =element_blank(),
        text = element_text(size=20)) 
  stat_summary(fun.y = max, fun.ymax = length,
               geom = "text", aes(label = ..ymax..), vjust = -1) 



#summarize number of genomes by technology
# img_actino_data %>% group_by(genus) %>% filter(scaffold_count_assembled < 30) %>% summarize(number_of_genomes = n())
# 
# 
# ggplot(data = img_actino_data) +
#   geom_boxplot(aes(x=sequencing_technology, y = scaffold_count_assembled),
#                outlier.shape = NA) + 
#   scale_y_continuous(limits = quantile(img_actino_data$scaffold_count_assembled, c(0.1, 0.9)))


p_sequencing_tech = ggplot(data = img_actino_data %>% filter(biosynthetic_cluster_count <100)) +
  geom_point(aes(x=scaffold_count_assembled, y = biosynthetic_cluster_count)) + 
  facet_wrap(~sequencing_technology, ncol = 3, scales = "fixed" ) +
  theme_classic() + theme(strip.background = element_blank(),
                          text = element_text(size=20)) +
  labs(x = "# of scaffolds", y= "# of BGC")


  
## combined plot genome sizes  
p_combined = plot_grid(plot_grid(p_genome_size,p_bgc_count, labels = c("A","B")),
                       p_sequencing_tech,labels = c("","C"), ncol = 1)

ggsave("figures/genome_sizes.svg",height = 8, width = 11)
