library(Biostrings)
library(tidyverse)
require(xcms)
library(stringr)
library(cowplot)
library(gggenes)
library(RColorBrewer)


ripps = c("LAP","lanthipeptide","lassopeptide","linaridin","thiopeptide",
          "lipolanthine","microviridin","proteusin","sactipeptide",
          "glycocin","cyanobactin")


##  ripps prism ####

json_files = paste0("data/prism/",list.files("data/prism/",))
prismdf = lapply(json_files,function(x) {
  #list of dfs, each df corresponds to cluster
  prism_json = jsonlite::fromJSON(x,flatten = T,simplifyDataFrame =T)
  prism_df_int = tibble(seqid = prism_json$prism_results$clusters$contig_name,
                        start = prism_json$prism_results$clusters$start,
                        end = prism_json$prism_results$clusters$end,
                        type = prism_json$prism_results$clusters$family %>% as.character(),
                        path_json = x
  )
  return(prism_df_int)}) %>% bind_rows() %>% filter(type == "RIBOSOMAL")


prismdf = mutate(prismdf,
                   width = abs(start - end),
                   date = str_extract(prismdf$path_json, pattern = "\\d{8}"),
                   barcode = str_extract(prismdf$path_json, pattern = "barcode\\d+"),
                   path = sprintf("/gpfs/gsfs11/users/rajwanir2/soil_metagenomics/data/antismash/%s/%s/medeka/%s.contigs.gff.1",
                                  date,barcode,barcode),
                   )
regions = prismdf

### extract precursors ####

samples = Sys.glob(file.path(getwd(), "data/antismash/*/*/medeka/*.gff.1"))
names(samples) = samples
tags_to_load = c("product","type","gene","locus_tag","translation","gene_functions")
samples = lapply(samples,function(x)
  rtracklayer::readGFF(x,tags=tags_to_load))
samples = bind_rows(samples,.id="path")
samples = mutate(samples, width = abs(start-end))

product_type = ripps

#if not prism
regions = filter(samples,product %in% product_type &
                   type...3 == "region") 


# regions = filter(regions, width < 50000) #exluding large BGC boundries

precursorRegions = lapply(1:nrow(regions),function(idx)
  samples %>% filter(
    path == regions$path[idx] &
    seqid == regions$seqid[idx] &
      start > regions$start[idx] &
       end < regions$end[idx] & 
      type...3 == "CDS" &
      width < 600
  )
  )

precursorSeqs = lapply(precursorRegions,function(r) {
  precSeq = AAStringSet(x = r$translation,use.names=TRUE)
  names(precSeq) = as.character(r$locus_tag)
  return(precSeq)
}
)




ripps_folder = sprintf("data/extracted_seqs/%s/",
                       "ripps_prism")
system(sprintf('rm data/extracted_seqs/%s/*',ripps_folder))
sapply(1:nrow(regions),function(idx){
if( nrow(precursorRegions[[idx]]) > 0 ) {

date_bc = paste(str_extract(precursorRegions[[idx]]$path, pattern = "\\d{8}"),
str_extract(precursorRegions[[idx]]$path, pattern = "barcode\\d+"),
sep = '_')
precursorSeqsFasta = sprintf("data/extracted_seqs/ripps_prism/%s.fasta",
                             unique(date_bc))
writeXStringSet(x=precursorSeqs[[idx]],
               append=T,
                filepath=precursorSeqsFasta)

}
# return(unique(date_bc))  
})


### create correspondence file ####

regions = mutate(regions, date_bc = paste(str_extract(regions$path, pattern = "\\d{8}"),
                                          str_extract(regions$path, pattern = "barcode\\d+"),
                                          sep = '_')
)

sample_names = read_csv("data/metadata.csv") %>% mutate(date_bc = paste(seq_date,
                                                         sprintf("barcode%02d",barcode),
                                                         sep="_")) %>% select(date_bc,sample)

correspondance = merge(regions,sample_names,
      all.x=T,all.y=F,
      by.x="date_bc",by.y="date_bc") %>% select(sample,date_bc) %>%
  unique() %>% 
  mutate(SequenceFile = paste0(ripps_folder,date_bc,".fasta"),
         spec_isp_ea = sprintf("data/QTOF/01082021/ISP1_EA_%s.mzXML",sample),
         spec_isp_buoh = sprintf("data/QTOF/01082021/ISP1_BuOH_%s.mzXML",sample),
         spec_r2a_ea = sprintf("data/QTOF/01082021/R2A_EA_%s.mzXML",sample),
         spec_r2a_buoh = sprintf("data/QTOF/01082021/R2A_BuOH_%s.mzXML",sample)) %>% 
  pivot_longer(cols = starts_with("spec_"), names_to = "spec_condition", values_to ="SpectraFile") %>% 
  select(SpectraFile,SequenceFile) %>% unique()


write_tsv(correspondance,
          file = "data/QTOF/metaminer_prism/corresp.tsv")


### plot ###############

sig_matches = read_tsv(file = "data/QTOF/metaminer_/significant_matches.tsv")


# draw spectra 
sample_plots = lapply(1:nrow(sig_matches),function(idx){
  spectra = list()
  spectra$spectra = spectra(filterPrecursorScan(readMSData(sig_matches$SpecFile[idx],mode = "onDisk"),
                                                sig_matches$Scan[idx])
  )[[2]]
  spectra$ce = collisionEnergy(spectra$spectra)
  #spectra$plot = plot_spectra(spectra$spectra)
  return(spectra)
})

spectra_df = tibble(
mz = as.numeric(sprintf("%.1f",mz(sample_plots[[1]]$spectra))),
intensity = intensity(sample_plots[[1]]$spectra))

# par(mar=c(5,5,1,1))
# plot(sample_plots[[1]]$spectra, 
#      y = "GYPWWQNRDIFGGRTFL",
#      modifications=c(D=-18.010, Q =+13.0269) , type = c("y","b"),
#      neutralLoss=NULL,
#     # neutralLoss=defaultNeutralLoss(disableAmmoniaLoss=c("K", "N", "Q", "R")),
#      main = sprintf("MS/MS spectra precursor m/z %f",
#                     precursorMz(sample_plots[[1]]$spectra)),
#      col="black"
#      ) 
# 
# p_spectra = recordPlot()

calc_fragments = calculateFragments(sig_matches$FragmentSeq, modifications=c(D=-18,Q =+13.0269))
calc_fragments$mz = as.numeric(sprintf("%.1f",calc_fragments$mz))
spectra_df = merge(spectra_df,calc_fragments,by.x="mz",by.y="mz",all.x=T,all.y=F)

ggplot(spectra_df, aes(x=mz,y=intensity)) +
  geom_bar(stat = "identity",width=2, fill = "black") +
  geom_text(aes(y=1500,label = ion), color = "grey") +
  geom_segment(data = spectra_df %>% filter(!is.na(ion)),
               aes(x=mz,xend=mz, y=intensity, yend = 1500),
               linetype = "dashed",
               color = "grey") +
  #ggrepel::geom_text_repel(aes(label=ion,x=mz,y=intensity),
  #                        min.segment.length = 0, segment.colour = "darkblue", color ="darkblue") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(limits = c(0,2000), expand = c(0,0)) +
  theme_linedraw() + theme(panel.grid = element_blank())




#draw precursor peptide sequence
precMShit = tibble(
aa=str_extract_all("MKKIAYIAPVIEKIGGFREVTNGYPWWDNRDIFGGRTFL",boundary ("character"))[[1]],
pos=1:length(aa))
precMShit = mutate(precMShit,prec = case_when(pos > 22 ~ "core",
                                  pos < 22|pos == 22 ~ "leader",
                                  ))

p_precSeq = ggplot(precMShit, aes(x=pos,y= 1, label = aa, fill = prec),color = "black") +
  geom_curve(y=1,yend=1,x=23,xend=31,
             curvature=0.5,angle=270,linetype="dashed") +
  geom_point(shape=21, size = 8) + 
  geom_text() +
  scale_y_discrete() + scale_x_discrete() +
  scale_fill_manual(values = c("leadeer" = "lightgrey",
                               "core" = "darkorange")) +
  ggtitle("Precursor peptide sequence") +
  theme_void() +
  guides(fill=guide_legend(title=NULL,nrow = 2)) + 
    theme(legend.position = "right",plot.title = element_text(hjust=0.5,face = "bold"))



# draw genes

# precLocustag= "ctg1_5502"
# hit_rip_genes = samples %>% filter(path==regions$path[26] ) 
# hit_rip_idx = which(hit_rip_genes$locus_tag ==precLocustag)[1]
# hit_rip_genes = hit_rip_genes[(hit_rip_idx-10):(hit_rip_idx+10),] %>% filter(type...3 == "CDS")
# hit_rip_genes$gene_functions = hit_rip_genes$gene_functions %>% str_extract("(.*):(.*)\\(") %>% str_remove_all(pattern = "(.*)smcogs\\) |\\(")
# hit_rip_genes = mutate(hit_rip_genes,gene_functions = case_when(locus_tag != precLocustag ~ gene_functions,
#                                                 locus_tag == precLocustag ~ "precursor peptide"))
# 


gff = "data/bgcs_gbk/12112020_barcode09_GA6-002.gff.1"
gene_start = 5116341 #
gff = rtracklayer::readGFF(gff,
                           tags = c("locus_tag","gene_functions"),
                           filter = list(type = c("CDS"),
                                         seqid = c("tig00000001")
                           )
) %>% filter(start > gene_start-(4*1000) & 
               end < gene_staWrt+(2*1000))

gff = mutate(gff,
             gene_func2 = simplify_genefunction(gene_functions,type = "smcogs")) %>% 
  mutate(gene_func2 = ifelse(is.na(gene_func2),str_remove(gene_functions,".*\\) "),gene_func2)) 




p_ripgenes = ggplot(gff, aes(xmin=start,xmax=end,y=1,fill=gene_func2, 
                          label=str_remove(locus_tag,pattern ="ctg\\d+_")
                          )) +
  geom_gene_arrow() + geom_gene_label(min.size = 5) +
  scale_y_discrete()+
  theme_void() + 
  ggtitle("Lassopeptide BGC in GA6-002") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"),legend.position = "bottom", 
        #panel.grid.major.y = element_line(size=2,color="black")
        ) + 
  scale_fill_brewer(palette = "Dark2") +
  guides(fill=guide_legend(nrow = 3, title.position = "top", title.hjust = 0.5, title = "Gene function")) 


# combine plots

p_combined= cowplot::plot_grid(p_ripgenes,
                   p_precSeq,
                   ggdraw(p_spectra),
                   ncol=1, rel_heights = c(0.4,0.7,0.7),
                   labels = "auto")


svg("figures/lassopeptide.svg",width=11,height=8)
p_combined
dev.off()

