library(tidyverse)
library(gggenes)
library(cowplot)

prism_file="data/similar_genome/prism_out/GCF_001509505.1_ASM150950v1_genomic.gbff.json"
cluster_num=6
cluster_name="GCF_001509505"

# prism_file="data/prism/09182020_barcode08_B06-02.gbk.json"
# cluster_num=26
# cluster_name="B06-02"


prism_file = jsonlite::fromJSON(prism_file)
prism_genes = prism_file$prism_results$clusters$orfs[[cluster_num]] %>% select(-sequence)
prism_genes = mutate(prism_genes,
       name_mod = case_when(is.na(type) ~ NA_character_,
                           !is.na(type ) ~ str_extract(name,pattern = "_\\w+") %>% str_remove("_")
                           ))


p_genes = ggplot(prism_genes, aes(xmin = start, xmax= stop, fill = type, y= 1)) +
  geom_gene_arrow() +
  ggrepel::geom_text_repel(aes(label=name_mod, x = (start+stop)/2, y=1))+
  scale_fill_discrete(na.value="white")  +
  theme_void() +
  ggtitle(sprintf("BGC structure \n%s cluster %d",
                  cluster_name,cluster_num)) +
  guides(fill = guide_legend(title = "Gene type",
                             title.theme = element_text(hjust = 0.5),
                             title.position = "top",
                             nrow = 1)) +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust=0.5),
        legend.background = element_rect(color="black"))


prism_domains = prism_file$prism_results$clusters$orfs[[cluster_num]]$domains
names(prism_domains) = prism_genes$name
prism_domains = prism_domains %>% bind_rows(.id = "gene_id") 
prism_domains = separate(prism_domains,col = "gene_id",sep = "_",into = c("locus_extra","gene_id"))


names(prism_domains$substrates) = paste(prism_domains$gene_id,prism_domains$module_number,sep = "_")
prism_domains_substrates = prism_domains$substrates %>% 
  bind_rows(.id = "gene_moduleno") %>% 
  distinct(gene_moduleno,.keep_all = T) %>% 
  separate(gene_moduleno,sep = "_", into = c("gene_id","module_no"))

prism_domains = merge(prism_domains,prism_domains_substrates,
      all.x = T,
      by.x = c("gene_id","module_number"),
      by.y = c("gene_id","module_no"),
      suffixes = c("domain","substrate"))



prism_domains = mutate(prism_domains,
       name_domain_specific = ifelse(!namedomain %in% c("A","AT"),namedomain,namesubstrate),
       name_domain_a = ifelse(namedomain %in% c("A","AT"),"yes",NA_character_))



#sort by domain start
prism_domains = prism_domains %>% arrange(startdomain)
prism_domains$row_name = row.names(prism_domains)


p_domains = ggplot(data = prism_domains,
       aes(x=fct_reorder(row_name,as.numeric(row_name)),
                         label = name_domain_specific,fill=name_domain_a)) +
         geom_point(shape=21,size=15,aes(y=1, color = active)) +
         geom_text(aes(y=1)) +
         facet_wrap(~gene_id, scales="free",
                    strip.position = "top",nrow = 2
                    ) +
         scale_y_continuous(expand=c(0,0), breaks = c(0,1,2))+
         theme_void() +
         ggtitle("Domain organization")+
         theme(legend.position = "none",
               plot.title = element_text(hjust=0.5),
               panel.grid.major.y = element_line(color="grey", size = 2),
              panel.spacing.x = unit(15,"mm")
               ) +
  scale_fill_discrete(na.value="white")


p_genes_domains = plot_grid(p_genes,p_domains,
          ncol=1,
          rel_heights = c(0.3,0.5),
          scale = 0.9)

ggsave(p_genes_domains,filename = sprintf("figures/%s_%d.svg",
                          cluster_name,cluster_num),
       width = 11, height=5)



##spectra ###
spec_details = read_tsv("data/QTOF/moldiscovery_prism_similargenome/significant_matches.tsv")
spec_details = spec_details %>% filter(str_detect(SpecFile,"B06-02"))
spectra = retrieveSpectra(spec_details)
spectra$df = spec2df(spectra$spectra[[1]])

spectra$plot = ggplot(spectra$df,
       aes(x=mz,y=intensity)) +
  geom_bar(stat = "identity",width=0.2, fill = "black") +
  theme_linedraw() + theme(panel.grid = element_blank()) +
  scale_y_continuous(expand = c(0,10)) +
  scale_x_continuous(expand = c(0,0), limits = c(25,spec_details$SpectrumMass[1]))


spectra$plot = label_peaks(spectra$plot, by = "top_intensity")
spectra$plot = spectra$plot + ggtitle(sprintf('precursor m/z: %.3f, collision energy: %d, charge: M+2H',
                               spectra$df$precursor[1],
                               spectra$ce[[1]]),
                       ) 
  


#structure ###
p_structure = ggdraw()+draw_image("figures/GCF_001509505_rban_structure.svg",
                                   #scale = 2, 
                                  # halign =0.8, valign = -0.5,
                                  # width  =2, height = 2
                                  )



### combine all ###
p_combined =  plot_grid(plot_grid(p_simgenome,p_B0602,nrow = 1,labels = c("A","B")),
                        plot_grid(spectra$plot,
                                 p_structure,
                                  nrow = 1, labels = c("C","D")),
               nrow=2, rel_widths = c(1,0.6))


ggsave(filename = "figures/GCF_001509505_combined.svg",
       p_combined, width = 20, height = 8)

