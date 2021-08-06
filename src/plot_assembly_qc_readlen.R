library(tidyverse)
library(data.table)
library(cowplot)


## functions #####


# fastqstats = tibble(sample = 1:3)
# fastqstats$fqpath = sapply(1:3,function(x) sprintf("data/rawdata/readlen/%s/%s.fastq.gz",x,x))
# fastqstats$readN50 = sapply(fastqstats$fqpath,function(y) get_read_N50(y)) 

get_read_N50 = function(fastqFilepath){
  read_lengths = tibble(readlen = sort(yield(FastqStreamer(fastqFilepath,n=5000))@sread@ranges@width) ,
                        cumulative = cumsum(readlen))
  readN50 = read_lengths$readlen[read_lengths$cumulative > (sum(read_lengths$readlen)/2)][1]
  return(readN50)
}


get_canu_coverage = function(){
  coverage = Sys.glob(file.path("data/assembly/readlen/*/canu/*.contigs.layout.tigInfo"))
  names(coverage) = coverage
  coverage = lapply(coverage,data.table::fread)
  coverage = bind_rows(coverage,.id="path")
  coverage = coverage %>% filter(tigClass == "contig")
  coverage$readlen = str_extract(coverage$path, pattern = "/\\d+/") %>% str_remove_all(pattern = "/") 
  coverage = coverage %>% group_by(path) %>% summarize(canu_coverage = mean(coverage),
                                                       readlen = unique(readlen)
                                                       )
  coverage = coverage %>% select(-c(path))
  return(coverage)
}



### main ####

# downsampling = c("1","2","3")
# readlen = c("500","1000","2000","4000")
# samples = readlen
# main = "readlen"


attributes_to_plot = c("# contigs","# mismatches per 100 kbp","# predicted genes (unique)","Total length", "# BGC (antiSMASH)")

#load quast data ######

# canu      =  bind_rows(read_quast_data("data/quast/readlen/2/canu/report.tsv","2"),
#               read_quast_data("data/quast/readlen/3/canu/report.tsv","3"),
#               read_quast_data("data/quast/readlen/1/canu/report.tsv","1")) %>% mutate(assembler = "canu")
# 


medeka    =  read_tsv("data/quast/readlen/all/medeka/report.tsv") %>% 
  pivot_longer(!Assembly,values_to = "value",names_to = "readlen") %>% 
  mutate(readlen = str_extract(readlen,"\\d+"),
         value = as.numeric(value),
         assembler = "medeka")
  
  
#antismash BGCs #####
#calcluating BGC num based on number of region gbk files
bgcs = tibble(path = Sys.glob("data/antismash/readlen/*/medeka/*region*.gbk")) %>%
  separate(
    path,
    into = c(
      "data",
      "antismash",
      "folder",
      "readlen",
      "assembler",
      "BGC"
    ),
    sep = "/"
  ) %>%
  count(readlen, assembler, name = "value") %>%
  mutate(Assembly = "# BGC (antiSMASH)")

assembly_stats = bind_rows(medeka,
                           bgcs
)


#calculating BGC number by reading antismash json files
# canu_bgc = lapply(1:3, function(x) {get_cluster_num_all(sample = x, assembler = "canu", quast_data = canu) %>% bind_rows() } ) %>% bind_rows()
# medeka_bgc = lapply(1:3, function(x) {get_cluster_num_all(sample = x, assembler = "medeka", quast_data = medeka) %>% bind_rows() } ) %>% bind_rows()
# assembly_stats = bind_rows(canu,
#                            medeka,
#                            canu_bgc,
#                            medeka_bgc
#                            )

## filtering assmbly stats ####
assembly_stats  = assembly_stats  %>% filter(Assembly %in% attributes_to_plot)
assembly_stats$Assembly = factor(assembly_stats$Assembly,
                                 levels = as.character(c("Total length","# contigs", "# mismatches per 100 kbp", "# predicted genes (unique)","# BGC (antiSMASH)" )))
assembly_stats = assembly_stats %>% filter(assembler == "medeka")


# # prism BGCs #####
# json_files = Sys.glob("data/prism_downsampling/*.json")
# prism_bgcs = sapply(json_files,function(x)
#   length(jsonlite::fromJSON(x,flatten = T,simplifyDataFrame =T)$prism_results$clusters$start)
# )
# 
# prism_bgcs = tibble(Assembly = "# BGC (PRISM)",
#                     assembler= "medeka",
#                     value = prism_bgcs,
#                     depth = as.numeric(sapply(str_extract_all(names(prism_bgcs),pattern = "\\d+"),'[[',2)),
#                     sample = sapply(str_extract_all(names(prism_bgcs),pattern = "\\d+"),'[[',1))
# 
# assembly_stats = bind_rows(assembly_stats,prism_bgcs)
# 
# 
# 
#get canu coverage ####
canu_coverage = get_canu_coverage()
canu_coverage = mutate(canu_coverage,
                       value = canu_coverage,
                       assembler = "medeka",
                       Assembly = "coverage")

assembly_stats = bind_rows(assembly_stats,canu_coverage) %>%
  select(-canu_coverage)


#final editing dataframe
assembly_stats = assembly_stats %>% mutate(Assembly = ifelse(Assembly == "# BGC", "# BGC (antiSMASH)",Assembly))


### final combined plot #######
p_assembly_stats = lapply(unique(assembly_stats$Assembly), function(attr)
  ggplot(data = assembly_stats %>% filter(Assembly == attr),
         aes(x = as.numeric(readlen), y = value)) +
    geom_line(size = 1.5) +
    geom_point(size = 5) +
    expand_limits(y = 0,x=0) +
   theme_classic() +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.placement = "outside",
      text = element_text(size = 15)
    ) +
    xlab("Average read length") + ylab(attr)
)

#edit the scale for mismatches graph
# names(p_assembly_stats) = unique(assembly_stats$Assembly)
# p_assembly_stats[["# mismatches per 100 kbp"]] = p_assembly_stats[["# mismatches per 100 kbp"]] + scale_y_continuous(limits=c(0,100))


p_assembly_stats = plot_grid(get_legend(p_assembly_stats[[1]] + theme(legend.position = "bottom")),
                             plot_grid(plotlist = p_assembly_stats, labels = LETTERS[1:6]),
                             ncol = 1,
                             rel_heights = c(0.1,1))

ggsave("figures/assembly_stats_readlen.svg",p_assembly_stats, height = 10, width = 10)


# bgc_data = data.frame(Assembly = as.character(),depth = as.numeric(),value=as.numeric(),sample = as.numeric())


#### contig edge ##########
# calculating BGC number by reading gff files and filtering those on contig edge
gff_paths=Sys.glob("data/antismash/readlen/*/medeka/*.gff.1")
gffs = lapply(gff_paths, function(g) rtracklayer::readGFF(g[1], filter = list(type=c("region"))))
names(gffs) = gff_paths
gffs = bind_rows(gffs,.id="path")
gffs = gffs %>% separate(
  path,
  into = c(
    "data",
    "antismash",
    "folder",
    "readlen",
    "assembler",
    "BGC"
  ),
  sep = "/"
)

gffs = mutate(gffs,readlen = as.numeric(readlen))




p_contig_edge = ggplot(gffs, aes(x=fct_reorder(as.character(readlen),
                                               as.numeric(readlen)),
                                 fill=contig_edge)) +
  geom_bar() +
  scale_y_continuous(limits = c(0,40), expand = c(0,0)) +
  ylab("# BGC") + xlab("\nAverage read length (nucleotides)") +
  theme_classic() +
  theme(legend.position = "top",
        strip.background = element_blank(),
        strip.text = element_text(size = 12,face = "bold")) +
  scale_fill_manual(values = c("True" = "brown","False"="grey"),
                    labels = c("On Contig edge","Not on contig edge"))+
  guides(fill = guide_legend(title = "# of BGCs on contig edge",
                             title.position = "top",
                             title.theme = element_text(face = "bold",hjust = 0.5)))

ggsave(p_contig_edge,
       filename = "figures/readlen_contigedge.svg",
       height = 4,
       width = 7)

