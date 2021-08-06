library(tidyverse)
library(MSnbase)
library(cowplot)

prism_matches = read_tsv(file = "data/QTOF/moldiscovery_prism/significant_matches.tsv") 

# filter any thing seen in media
prism_matches = prism_matches[!prism_matches$Name %in%
                                prism_matches[grepl(prism_matches$SpecFile, pattern = "MC-00"),]$Name
                                  ,]

#create columns for media and solvent
prism_matches$sample = str_extract(prism_matches$SpecFile,pattern = "(R2A|ISP1)_(.)+_(.)+\\.") %>% str_remove(pattern = "\\.")
prism_matches = separate(prism_matches,col = "sample",into = c("medium","solvent","m_strain"), sep = '_')
prism_matches = mutate(prism_matches,
                       condition=paste(medium,solvent,m_strain,sep = ' '))
#get info from prism json path
prism_matches = separate(prism_matches,col = "Name", into = c("date","barcode","p_strain","bgc_n","smile_n"),sep = '_', remove = F)
prism_matches$json = str_extract(prism_matches$Name,pattern = "\\d{8}_barcode\\d+_(.)+-\\d+")
prism_matches$json = paste0("data/prism/",prism_matches$json,".gbk.json")
prism_matches$bgc_n = str_extract(prism_matches$Name, pattern = "_c\\d+_") %>% str_remove_all(pattern = '_|c')


#same strain and metabolite #
  # filter wrong strain and metabolite matches
  prism_matches = prism_matches %>% filter(m_strain == p_strain)
  
  # remove duplicated matches - same strain and compound
  prism_matches = mutate(prism_matches,
         condition=paste(medium,solvent,m_strain,sep = ' '),
         name_simple=case_when(Name == "12112020_barcode09_GA6-002_c3_s2" ~ "cWW-Me2",
                               Name == "10222020_barcode11_GA6-002_c10_s1" ~ "cWW")
         )
  
  # one row - one compound
  # prism_matches = prism_matches[!duplicated(prism_matches[c('m_strain','SMILES')]),] 

  
# different strain and metabolite
  prism_matches = prism_matches %>% filter(m_strain != p_strain)
  

  #GB4-14
  prism_matches = prism_matches %>% filter(m_strain == "GB4-14" & p_strain == "GB14-14")
  prism_matches
  

# one row - one condition & compound 
prism_matches = prism_matches[!duplicated(prism_matches[c('condition','name_simple')]),]








### BGC #################
# retrieve genomic bgc coordinates
prism_bgcs = lapply(1:nrow(prism_matches),function(idx) {
  x = prism_matches$json[idx]
  #list of dfs, each df corresponds to cluster
  prism_json = jsonlite::fromJSON(x,flatten = T,simplifyDataFrame =T)
  prism_df_int = tibble(seqid = prism_json$prism_results$clusters$contig_name,
                        start = prism_json$prism_results$clusters$start,
                        end = prism_json$prism_results$clusters$end,
                        type = prism_json$prism_results$clusters$family %>% as.character(),
                        path_json = x,
                        smiles_n = sapply(prism_json$prism_results$clusters$predicted_molecule_masses, 
                                          function(y) dim(y)[1])
  )
  # remmove clusters that do not have smiles
  #prism_df_int = filter(prism_df_int, smiles_n > 0)
  # get specific cluster only
  prism_df_int = prism_df_int[as.numeric(prism_matches[idx,'bgc_n']),]
  
  return(prism_df_int)}) %>% bind_rows()

prism_matches = cbind(prism_matches,prism_bgcs)
rm(prism_bgcs)


# get antismash genes for prism bgc coordinates
prism_matches = mutate(
  prism_matches,
  gff_path = sprintf(
    "data/antismash/%s/%s/medeka/%s.contigs.gff.1",
    date,
    barcode,
    barcode
  )
)
tags_to_load = c("product", "type", "gene", "locus_tag", "gene_functions","translation")
prism_matches_gff = lapply(prism_matches$gff_path, function(x)
  rtracklayer::readGFF(x, tags = tags_to_load))
names(prism_matches_gff) = prism_matches$gff_path
prism_matches_gff = bind_rows(prism_matches_gff, .id = "gff_path")

#filter only bgc region
prism_matches_gff = lapply(1:nrow(prism_matches),function(idx)
  prism_matches_gff %>% filter(
    gff_path == prism_matches$gff_path[idx] &
      seqid == prism_matches$seqid[idx] &
      start > prism_matches$start[idx] &
      end <= prism_matches$end[idx] & 
    type...3 == "CDS" 
    ))


prism_matches_gff[[1]]$gene_functions = simplify_genefunction(prism_matches_gff[[1]]$gene_functions,
                                                              type = "cluster")

## ploting #####################

#genes
prism_matches_gff = retrieve_prism_genes()
p_genes = ggplot(data = prism_matches_gff,
       aes(xmin=start_gene,xmax=stop_gene,fill=clean_name,y=1, label = name)) +
  gggenes::geom_gene_arrow(arrowhead_height = unit(10, "mm"), 
                           arrowhead_width = unit(3, "mm"),
                           arrow_body_height = unit(8, "mm"),
                           color = "black") +
  geom_text(aes(y=1,x=(start_gene+stop_gene)/2))+
  scale_fill_brewer(palette = "Dark2") +
  theme_void() +
  ggtitle("Biosynthetic gene cluster") +
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5,face = "bold")) + 
  guides(fill = guide_legend(title = "Gene function", 
                             title.position = "top", 
                             title.theme = element_text(hjust = 0.5))) 

#spectra
# prism_spectra = retrieveSpectra(prism_matches[1,])

#eic
prism_eic = bind_rows(lapply(1:nrow(prism_matches),
                             function(idx) retrieveEIC(idx = idx)))
  
p_eic = ggplot(prism_eic,
               aes(x=rtime,y=intensity, 
                   color = condition
               )) +
  geom_line(size=1,
            aes(linetype=name)
  ) +
  labs(x="retention time (min)")+
  theme_classic() +
  theme(legend.position = "right") + 
  scale_x_continuous(limits = c(4,6)) +
  scale_y_continuous(labels=scales::scientific_format()) +
  scale_color_brewer(palette = "Set1") +
  ggtitle("Extracted ion chromatogram") +
  theme(plot.title = element_text(hjust = 0.5,face = "bold")) +
  ggprism::theme_prism()
  

#structure
prism_unique_compounds = unique(select(prism_matches,c(name_simple,SMILES)))
p_struc = cowplot::plot_grid(plotlist = lapply(prism_unique_compounds$SMILES,
                             function(smile) {
                               cowplot::ggdraw() +
                                 cowplot::draw_image(sprintf(
                                   "https://gnps-structure.ucsd.edu/structureimg?smiles=%s",
                                   smile
                                 ))}),
labels = prism_unique_compounds$name_simple,
nrow = 1, hjust =-0.7,vjust =2,label_fontface = "plain"
)


p_struc= plot_grid(ggdraw() + draw_label("Predicted structure", fontface = "bold"),
                    p_struc,
                    ncol=1,nrow=2,
                    rel_heights = c(0.1,0.9)
                    )
  
#combined
p_prism_combined = cowplot::plot_grid(cowplot::plot_grid(p_genes,p_struc,nrow = 2,
                                      rel_heights = c(0.4,0.9),
                                      labels=c("A","B")),
p_eic,ncol = 1,
labels = c("","C"))


svg("figures/prism.svg",height = 11,width=8)
p_prism_combined
dev.off()
