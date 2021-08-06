library(tidyverse)
library(data.table)

mibig_npatlas_matches = read_tsv(file = "data/QTOF/moldiscovery_npatlas/significant_matches_filtered.tsv") %>% 
  filter(mibig_ids != "[]") 
mibig_npatlas_matches = mibig_npatlas_matches[!duplicated(mibig_npatlas_matches[c('strain','mibig_ids')]),]  
mibig_npatlas_matches$mibig_ids = str_extract(mibig_npatlas_matches$mibig_ids,
                                              pattern = 'BGC\\d+')

mibig_npatlas_matches = mibig_npatlas_matches %>%
  filter(!mibig_ids %in% c("BGC0001180"))

# mibig_npatlas_matches$tax_division = sapply(mibig_npatlas_matches$genus,function(genus){
# tax_id =  rentrez::entrez_search(db = "taxonomy", 
#                                  term = genus)$ids 
# 
# tax_division = rentrez::entrez_summary(db="taxonomy",
#                                        id = tax_id)$division
# 
# return(tax_division)
# })
# 
# mibig_npatlas_matches$tax_division = as.character(mibig_npatlas_matches$tax_division)


## pre-processing ####

#mibig
mibig_npatlas_matches$exit_stat = download_mibig(mibig_ids = mibig_npatlas_matches$mibig_ids,
                                                 dest_dir = "data/QTOF/moldiscovery_npatlas/")

gbk2gff(file_pat = "data/QTOF/moldiscovery_npatlas/mibig/gbk/*.gbk")


sapply(Sys.glob("data/QTOF/moldiscovery_npatlas/mibig/gbk/*.gff"),
       extract_proteins_from_gff)


# samples pre-process

gbk2gff(file_pat = "data/bgcs_gbk/*.gbk")
sapply(Sys.glob("data/bgcs_gbk/*.gff.1"),
       extract_proteins_from_gff)

# mibig_npatlas_matches = mibig_npatlas_matches %>%
#   filter(genus %in% c("Streptomyces","Amycolatopsis"))
# 
# mibig_npatlas_matches = mibig_npatlas_matches %>% 
#   filter(tax_division != "cyanobacteria")
# 
# mibig_npatlas_matches = mibig_npatlas_matches %>% 
#   filter(!strain %in% "GA3-001A")


## blasting mibig to sample ####
blast_results = lapply(1:nrow(mibig_npatlas_matches),function(idx){
blast_cmd = sprintf('blastp -subject %s -query data/QTOF/moldiscovery_npatlas/mibig/gbk/%s_protein.fasta -outfmt \'6 delim=@ qseqid sseqid pident length sstart\' -max_target_seqs 1',
                    Sys.glob(sprintf("data/bgcs_gbk/*_%s_protein.fasta.1",mibig_npatlas_matches$strain[idx])),
                    mibig_npatlas_matches$mibig_ids[idx])
blast_out = read.delim(text = system(blast_cmd, intern = T), header = F,
                       col.names = c("qseqid", "sseqid", "pident", "length", "sstart"),
                      sep = '@')
return(blast_out)
}) 


# antismash known cluster blast #####

sample_metadata = read.csv("data/metadata_assembled.csv") %>% 
  select(sample, barcode, seq_date) %>%
  mutate(barcode = sprintf("barcode%02d", barcode),
         seq_date = sprintf("%08d", seq_date))
sample_metadata = sample_metadata[!duplicated(sample_metadata$sample),] #GA3-008 is duplicated



mibig_npatlas_matches = merge(mibig_npatlas_matches,sample_metadata,
                              by.x = "strain", by.y = "sample",
                              all.x = T,all.y = F)


rm(sample_metadata)

antismash_clusterblast = lapply(1:nrow(mibig_npatlas_matches),function(idx){

  cmd = sprintf("grep %s data/antismash/%s/%s/medeka/knownclusterblastoutput.txt",
           mibig_npatlas_matches$mibig_ids[idx],
           mibig_npatlas_matches$seq_date[idx],
           mibig_npatlas_matches$barcode[idx]
           )
  
  result = system(cmd,intern = T)
  
  result = read_antismash_clusterblast(result)
  
  return(result)
}) 

antismash_clusterblast = lapply(antismash_clusterblast,keep_best_clusterblast)


#single cluster

# no longer needed due to keep_best_clusterblast earliar
singlecluster_blast = sapply(antismash_clusterblast,function(x)
  select(x,"sample_bgc_num") %>% unique() %>% nrow())

#need loading mibig gff
# coverage_clusterblast = sapply(1:length(antismash_clusterblast),function(idx){
#   sum(antismash_clusterblast[[idx]]$alignment_length) /
#     sum(abs(gffs$mibig[[idx]]$end - gffs$mibig[[idx]]$start)/3) 
# })

#sample:mibig cds aligment ratio
cds_aln_ratio = sapply(antismash_clusterblast,function(blastdf){
  blastdf = blastdf %>% 
    pivot_longer(contains('_locus_id'),names_to="locus_id_type",values_to="locus_id") %>% 
    pivot_longer(contains('cds_length'),names_to="cds_length_type",values_to="cds_length")  
  
  cds_sum = function(blastdf_i,type) {
    return(blastdf_i %>% filter(locus_id_type %like% type & 
                                  cds_length_type %like% type) %>% 
             distinct(locus_id,.keep_all = T) %>% 
             summarize(cds_length_sum = sum(cds_length)))
  }
  
  cds_sum(blastdf,type = "sample")/
    cds_sum(blastdf,type = "mibig")
  
}) %>% as.numeric()


antismash_clusterblast = antismash_clusterblast[singlecluster_blast == 1 & 
                                                  (cds_aln_ratio > 0.7 & cds_aln_ratio <1.1)]
mibig_npatlas_matches = mibig_npatlas_matches[singlecluster_blast == 1 & 
                                                (cds_aln_ratio > 0.7 & cds_aln_ratio <1.1),]




gffs = list()

gffs$sample = lapply(1:nrow(mibig_npatlas_matches),function(idx){
  gff_path = Sys.glob(sprintf("data/bgcs_gbk/%s_%s_*.gff.1",
                         mibig_npatlas_matches$seq_date[idx],
                         mibig_npatlas_matches$barcode[idx]
                         ))
  
  gff = read_gff_clusterblast(gff_path, idx = idx, type = "sample")
  return(gff)
})

gffs$mibig = lapply(1:nrow(mibig_npatlas_matches),function(idx){
  gff_path = sprintf("data/QTOF/moldiscovery_npatlas/mibig/gbk/%s.gff.1",
                     mibig_npatlas_matches$mibig_ids[idx])
  gff = rtracklayer::readGFF(gff_path,
                             filter=list(type = c("CDS"))
  )
  return(gff)
})

gffs$merged = lapply(1:nrow(mibig_npatlas_matches),function(idx){
  print(idx)
  bind_rows(as.tibble(gffs$mibig[[idx]]),
            as.tibble(gffs$sample[[idx]])
                      )
})

antismash_clusterblast = lapply(antismash_clusterblast, function (blast_df){
  blast_df = mutate(blast_df,
                    sample_cds_start = as.numeric(sample_cds_start),
                    sample_cds_end = as.numeric(sample_cds_end),
                    mibig_cds_start = as.numeric(mibig_cds_start),
                    mibig_cds_end = as.numeric(mibig_cds_end)
                    )
  
  blast_df = mutate(blast_df,
                    sample_cds_start = abs(sample_cds_start-min(sample_cds_start)),
                    sample_cds_end = abs(sample_cds_end-min(sample_cds_start))
                                           )
  return(blast_df)
})





# antismash_clusterblast = lapply(antismash_clusterblast, function (blast_df){
#   
#   #not all mibig gffs contain locus_id, needs solution
#   
#   # blast_df = mutate(blast_df,
#   #                   tag2 = sample_cds_locus_id)
#   
#   blast_df= pivot_longer(blast_df,
#                cols=contains("locus_id"),
#                names_to="locus_type",
#                values_to="blast_locus_id")
#   
#   # blast_df= pivot_longer(blast_df,
#   #                        cols=c("mibig_add_id","sample_cds_locus_id"),
#   #                        names_to="add_id_type",
#   #                        values_to="add_id")
#   
# }
# )

# gffs$merged = lapply(1:nrow(mibig_npatlas_matches), function(idx) {
#   if ("locus_tag" %in% colnames(gffs$merged[[idx]]))
#   {
#     merge(
#       gffs$merged[[idx]],
#       antismash_clusterblast[[idx]],
#       all.x = T,
#       all.y = F,
#       by.x = "locus_tag",
#       by.y = "blast_locus_id"
#     )
#   }
#   
#   else{
#     merge(
#       gffs$merged[[idx]],
#       antismash_clusterblast[[idx]],
#       all.x = T,
#       all.y = F,
#       by.x = "protein_id",
#       by.y = "mibig_add_id"
#     )
#   }
#   
# })


# zfill_names = sapply(gffs$merged,'[[','gene_kind') %>% unlist() %>% unique()
# fill_colors = RColorBrewer::brewer.pal(length(fill_names),"Dark2")
# names(fill_colors) = fill_names

gffs$plot = lapply(1:nrow(mibig_npatlas_matches),function(idx){
  ggplot()+
   gggenes::geom_gene_arrow(data = gffs$merged[[idx]],
                            mapping = aes(xmin=start,xmax=end,y=seqid, 
                                          fill = gene_kind
                                          )) +
  geom_segment(data=antismash_clusterblast[[idx]],
    aes(x=mibig_cds_start,
        xend=sample_cds_start,
        y=1.2,yend=1.8,
        color=p_id,
        group=sample_cds_locus_id
        ),
    inherit.aes = F,size = 1) +
  scale_color_continuous(na.value=NA, limits = c(0,100),low="white",high="navy") +
    scale_fill_brewer(palette = "Dark2",na.value=NA) +
    theme_void() +
    scale_y_discrete(labels = c(sprintf("%s\n%s %s",
                                        mibig_npatlas_matches$compound_names[idx],
                                      mibig_npatlas_matches$mibig_ids[idx],
                                      mibig_npatlas_matches$genus[idx]
                                      ),
                                mibig_npatlas_matches$strain[idx])) +
    theme(axis.text.y = element_text(face="italic")) +
    theme(legend.position = "none") 
})


p_mibig_legend = cowplot::get_legend(gffs$plot[[3]] + 
                                       theme(legend.position = "bottom")+
                                       guides(color = guide_colourbar(title = "Amino acid identity",
                                                                       title.position = "top"),
                                              fill = guide_legend(title = "Gene type",
                                                                  nrow = 2,
                                                                  title.position = "top",
                                                                  title.theme = element_text(hjust = 0.5))))


# p_mibig = cowplot::plot_grid(plotlist = gffs$plot,
#                              p_mibig_legend,
#                              ncol=1)


svg(file = "figures/npatlas_matches_mibig.svg", width = 8, height = 11)
cowplot::plot_grid(plotlist = gffs$plot,
                   p_mibig_legend,
                   ncol=1)
dev.off()