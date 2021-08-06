library(tidyverse)
library(data.table)
library(cowplot)


## functions #####


# fastqstats = tibble(sample = 1:3)
# fastqstats$fqpath = sapply(1:3,function(x) sprintf("data/rawdata/simulation/%s/%s.fastq.gz",x,x))
# fastqstats$readN50 = sapply(fastqstats$fqpath,function(y) get_read_N50(y)) 

get_read_N50 = function(fastqFilepath){
  read_lengths = tibble(readlen = sort(yield(FastqStreamer(fastqFilepath,n=5000))@sread@ranges@width) ,
                        cumulative = cumsum(readlen))
  readN50 = read_lengths$readlen[read_lengths$cumulative > (sum(read_lengths$readlen)/2)][1]
  return(readN50)
}


read_quast_data <- function(tsv_report,sample_name){
  tsv_report = fread(tsv_report)
  tsv_report = tsv_report %>% pivot_longer(starts_with("x"), names_to ="depth", values_to ="value")
  tsv_report$value = as.numeric(tsv_report$value)
  tsv_report$depth = tsv_report$depth %>% str_remove("x") %>% str_remove(".contigs") %>% as.numeric()
  tsv_report %>% mutate(sample = sample_name)
}

get_cluster_num <- function(antismash_json, type = "antismash") {
  
  #reads the antismash json
  # ncol(jsonlite::stream_in(file(antismash_json), , pagesize = 100000)$timings)
  #file = jsonlite::stream_in(file(antismash_json), , pagesize = 100000)
  
  clust_num  = lapply(jsonlite::stream_in(file(antismash_json), pagesize = 100000)$records[[1]]$features, function(x) 
    {nrow (x %>% filter(type == "region"))}) %>%
    unlist() %>% sum()
  
  #reads gff and check number of regions not on contig edge
  rtracklayer::readGFF(gff_paths[1], filter = list(type=c("region"),contig_edge=c("FALSE")))
  
  
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


get_canu_coverage = function(){
  coverage = Sys.glob(file.path("data/assembly/simulation/*/*/canu/*.contigs.layout.tigInfo"))
  names(coverage) = coverage
  coverage = lapply(coverage,data.table::fread)
  coverage = bind_rows(coverage,.id="path")
  coverage = coverage %>% filter(tigClass == "contig")
  coverage$sample = str_extract(coverage$path, pattern = "/\\d/") %>% str_remove_all(pattern = "/") 
  coverage$depth = str_extract(coverage$path, pattern = "x\\d+|\\d/\\d") %>% str_remove_all(pattern = "\\d/|x") %>% as.integer()
  coverage = coverage %>% group_by(path) %>% summarize(canu_coverage = mean(coverage),
                                                       sample = unique(sample),
                                                       depth = unique(depth),
                                                       )
  coverage = coverage %>% select(-c(path))
  return(coverage)
}



### main ####


attributes_to_plot = c("# contigs","# mismatches per 100 kbp","# predicted genes (unique)","Total length", "# BGC")

#load quast data ######

# canu      =  bind_rows(read_quast_data("data/quast/simulation/2/canu/report.tsv","2"),
#               read_quast_data("data/quast/simulation/3/canu/report.tsv","3"),
#               read_quast_data("data/quast/simulation/1/canu/report.tsv","1")) %>% mutate(assembler = "canu")
# 

medeka    =  bind_rows(read_quast_data("data/quast/simulation/2/medeka/report.tsv","2"),
                        read_quast_data("data/quast/simulation/3/medeka/report.tsv","3"),
                        read_quast_data("data/quast/simulation/1/medeka/report.tsv","1")) %>% mutate(assembler = "medeka")

#antismash BGCs #####
#calcluating BGC num based on number of region gbk files
bgcs = tibble(path = Sys.glob("data/antismash/simulation/*/*/medeka/*region*.gbk")) %>%
  separate(
    path,
    into = c(
      "data",
      "antismash",
      "simulation",
      "sample",
      "depth",
      "assembler",
      "BGC"
    ),
    sep = "/"
  ) %>%
  count(sample, depth, assembler, name = "value") %>%
  mutate(Assembly = "# BGC",
         depth = as.numeric(str_extract(depth, pattern = "\\d+")))

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
                                    levels = as.character(c("Total length","# contigs", "# mismatches per 100 kbp", "# predicted genes (unique)","# BGC" )))
assembly_stats = assembly_stats %>% filter(assembler == "medeka")


# prism BGCs #####
json_files = Sys.glob("data/prism_downsampling/*.json")
prism_bgcs = sapply(json_files,function(x)
  length(jsonlite::fromJSON(x,flatten = T,simplifyDataFrame =T)$prism_results$clusters$start)
)

prism_bgcs = tibble(Assembly = "# BGC (PRISM)",
                    assembler= "medeka",
                    value = prism_bgcs,
                    depth = as.numeric(sapply(str_extract_all(names(prism_bgcs),pattern = "\\d+"),'[[',2)),
                    sample = sapply(str_extract_all(names(prism_bgcs),pattern = "\\d+"),'[[',1))

assembly_stats = bind_rows(assembly_stats,prism_bgcs)



#get canu coverage ####
canu_coverage = get_canu_coverage()
assembly_stats = merge(assembly_stats,canu_coverage,
                  by = c("sample","depth"))


#final editing dataframe
assembly_stats = assembly_stats %>% mutate(Assembly = ifelse(Assembly == "# BGC", "# BGC (antiSMASH)",Assembly))
assembly_stats = assembly_stats %>% dplyr::rename(genome = sample)


### final combined plot #######
p_assembly_stats = lapply(sort(unique(assembly_stats$Assembly),decreasing = T), function(attr)
  ggplot(data = assembly_stats %>% filter(Assembly == attr),
                   aes(x = canu_coverage, y = value, 
                       color = genome)) +
    geom_line(
      size = 1.5,
      #aes(linetype =assembler)
      ) +
      geom_point(size = 5) +
        #scale_linetype_manual(values = c("canu" = "dotted", "medeka" = "solid")) +
        #scale_alpha_discrete(range=c(0.5,1))+
        #facet_wrap(vars(Assembly),scales = "free", nrow=2,strip.position="left") +
        theme_classic() +
        theme(
          legend.position = "none",
          strip.background = element_blank(),
          strip.placement = "outside",
          text = element_text(size = 15)
        ) +
        xlab("Sequencing coverage") + ylab(attr)
)

p_assembly_stats = plot_grid(get_legend(p_assembly_stats[[1]] + theme(legend.position = "bottom")),
          plot_grid(plotlist = p_assembly_stats, labels = LETTERS[1:6]),
          ncol = 1,
          rel_heights = c(0.1,1))
    
ggsave("figures/assembly_stats_simulation.svg",p_assembly_stats, height = 10, width = 10)


# bgc_data = data.frame(Assembly = as.character(),depth = as.numeric(),value=as.numeric(),sample = as.numeric())


#### contig edge ##########
#calculating BGC number by reading gff files and filtering those on contig edge
gff_paths=Sys.glob("data/antismash/simulation/*/*/medeka/*.gff.1")
gffs = lapply(gff_paths, function(g) rtracklayer::readGFF(g[1], filter = list(type=c("region"))))
names(gffs) = gff_paths
gffs = bind_rows(gffs,.id="path")
gffs = gffs %>% separate(
  path,
  into = c(
    "data",
    "antismash",
    "simulation",
    "sample",
    "depth",
    "assembler",
    "BGC"
  ),
  sep = "/"
)

gffs = mutate(gffs,depth = as.numeric(str_extract(depth,pattern = "\\d+"))) %>%
  merge(.,canu_coverage,
        by=c("sample","depth")) %>% 
  mutate(coverage_label = paste(depth,
                                              sprintf('%.1f',canu_coverage),
                                              sep = '\n\n'))


gffs = mutate(gffs,
              sample = paste0("Genome-",sample))


p_contig_edge = ggplot(gffs, aes(x=factor(depth, levels = c(7,15,30,60)),
                                   #fct_reorder(coverage_label,depth),
                                 fill=contig_edge)) +
  geom_bar() + facet_wrap(~sample, scales = "free") +
  scale_y_continuous(limits = c(0,40), expand = c(0,0)) +
  ylab("# BGC") + xlab("\nSequencing coverage") +
  theme_classic() +
  theme(legend.position = "top",
        strip.background = element_blank(),
        strip.text = element_text(size = 12,face = "bold")) +
  scale_fill_manual(values = c("True" = "brown","False"="grey"), 
                    labels = c("On contig edge","Not on contig edge"))+
  guides(fill = guide_legend(title = "# of BGCs on contig edge",
                             title.position = "top",
                             title.theme = element_text(face = "bold",hjust = 0.5)))

ggsave(p_contig_edge,
       filename = "figures/downsampling_contigedge.svg",
       height = 4,
       width = 7)

