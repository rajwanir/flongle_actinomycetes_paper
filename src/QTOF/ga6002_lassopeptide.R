library(tidyverse)
library(MSnbase)
library(gggenes)
#load simplify_genefunction from bgc_functions.R

# draw precursor ####

lassoseqs = read_csv("data/extracted_seqs/lassopeptide_ga6002.csv",
                     na = "-")
lassoseqs = lassoseqs %>% pivot_longer(c(GA6_002,propeptin,propeptin2), names_to = "peptide", values_to = "residue")
lassoseqs = mutate(lassoseqs,
                    residue_color = case_when(section == "core" ~ residue,
                                              section == "leader" ~ NA_character_
                                              ))


p_lassoseq = ggplot(lassoseqs,aes(x=position, y=peptide)) +
  geom_point(aes(fill = residue_color, color = residue_color), shape =21, size = 8) +
  geom_text(aes(label=residue)) 


for (i in 1:3){
    p_lassoseq = p_lassoseq + geom_segment(x = 1, xend = 9,
               y = i+0.5, yend = i+0.5, linetype = "dashed") +
  geom_segment(x = 1, xend = 1,
               y = i+0.5, yend = i+0.1, linetype = "dashed") +
  geom_segment(x = 9, xend = 9,
               y = i+0.5, yend = i+0.1, linetype = "dashed") } 

p_lassoseq = p_lassoseq +
  theme_classic() +
  ggtitle("Lasso peptide sequences") +
  theme(axis.line.y = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, face="bold"),
        axis.ticks.y= element_blank(),
        axis.text.y= element_text(size=12),
        axis.text = element_text(color="black")) +
  labs(y=NULL) +
  scale_fill_discrete(na.value = "white") +
  scale_color_discrete(na.value = "white") +
  scale_x_continuous(breaks = c(seq(-1,-20,by=-3),seq(1,20,by=3)))

# draw spectra ####

sig_matches = read_tsv(file = "data/QTOF/metaminer_/significant_matches.tsv")



# draw spectra 
# sample_plots = lapply(1:nrow(sig_matches),function(idx){
#   spectra = list()
#   spectra$spectra = spectra(filterPrecursorScan(readMSData(sig_matches$SpecFile[idx],mode = "onDisk"),
#                                                 sig_matches$Scan[idx])
#   )[[2]]
#   spectra$ce = collisionEnergy(spectra$spectra)
#   #spectra$plot = plot_spectra(spectra$spectra)
#   return(spectra)
# })
# 
# spectra_df = tibble(
#   mz = as.numeric(sprintf("%.1f",mz(sample_plots[[1]]$spectra))),
#   intensity = intensity(sample_plots[[1]]$spectra))


calc_fragments = bind_rows(calculateFragments(sig_matches$FragmentSeq,
                                              neutralLoss=defaultNeutralLoss(disableAmmoniaLoss=c("K", "N", "Q", "R")),
                                              z=1),
                           calculateFragments(sig_matches$FragmentSeq,
                                    neutralLoss=defaultNeutralLoss(disableAmmoniaLoss=c("K", "N", "Q", "R")),
                                    z=2))
                           

spectra_df = read_csv("data/QTOF/extracted_spectrum/lassopeptide 3plus ion 10V MSMS peak list.csv",comment = "#") %>% janitor::clean_names() %>%
  select(mz = center_x, intensity = height, charge = z)


spectra_df = fuzzyjoin::difference_left_join(spectra_df,calc_fragments,max_dist = 0.001, distance_col = "mz_diff") 
spectra_df = mutate(spectra_df,mz_observed = mz.x, mz_theoretical = mz.y,
                    ppm = sprintf("%.1f",(mz_diff/mz_theoretical)*1000000)) %>% select(-c(mz.x,mz.y))

spectra_df = mutate(spectra_df,ion_label = ifelse(!is.na(ion), sprintf("%s (%s ppm)\n[M+%dH]", ion,ppm, z), NA_real_))
spectra_df = mutate(spectra_df,ion_label = ifelse(mz_observed == 694.6700,"*",ion_label))



# spectra_df = merge(spectra_df,calc_fragments,by.x="mz",by.y="mz",all.x=T,all.y=F)

p_spectra = ggplot(spectra_df, aes(x = mz_observed, y = intensity)) +
  geom_bar(stat = "identity", width = 1, fill = "black") +
  geom_text(aes(
    y = intensity + 5000,
    label = ion_label
  ),
  color = "black",
  size = 4) +
  geom_segment(
    data = spectra_df %>% filter(!is.na(ion)),
    aes(
      x = mz_observed,
      xend = mz_observed,
      y = intensity,
      yend = intensity + 5000
    ),
    linetype = "dashed",
    color = "black"
  ) +
  #ggrepel::geom_text_repel(aes(label=ion,x=mz,y=intensity),
  #                        min.segment.length = 0, segment.colour = "darkblue", color ="darkblue") +
  ggtitle("MS/MS spectra GA6-002 \nprecursor m/z 694.672, charge +3, collision energy 10") +
  labs(x="m/z") +
  theme_classic(base_size=14) + 
  theme(panel.grid = element_blank(),
                          plot.title = element_text(hjust = 0.5, face =
                                                      "bold"),
        axis.text = element_text(color="black")) 
  


#adjust scale
p_spectra = p_spectra + expand_limits(x=0,y=0) +
  scale_x_continuous(breaks = seq(0,1500,by=200),expand = c(0,0)) + 
  scale_y_continuous(limits = c(0,22000),expand = c(0,0)) 
  

# draw genes ####
gff = "data/bgcs_gbk/12112020_barcode09_GA6-002.gff.1"
gene_start = 5116341 #
gff = rtracklayer::readGFF(gff,
                           tags = c("locus_tag","gene_functions"),
                           filter = list(type = c("CDS"),
                                         seqid = c("tig00000001")
                           )
) %>% filter(start > gene_start-(4*1000) & 
               end < gene_start+(2*1000))

gff = mutate(gff,
             gene_func2 = simplify_genefunction(gene_functions,type = "smcogs")) %>% 
  mutate(gene_func2 = ifelse(is.na(gene_func2),str_remove(gene_functions,".*\\) "),gene_func2)) 




p_ripgenes = ggplot(gff, aes(xmin=start,xmax=end,y=1,fill=gene_func2, 
                             label=str_remove(locus_tag,pattern ="ctg\\d+_")
)) +
  geom_gene_arrow() + geom_gene_label(min.size = 5) +
  scale_y_discrete()+
  theme_genes() + 
  ggtitle("Lasso peptide BGC in GA6-002") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"),
        legend.position = "bottom", 
        axis.title.y = element_blank(),
        axis.text = element_text(color="black")
        #panel.grid.minor.y = element_line(size=2,color="black")
  ) + 
  scale_fill_brewer(palette = "Dark2") +
  guides(fill=guide_legend(nrow = 2, title.position = "top", title.hjust = 0.5, title = "Gene function")) 



# eic ####
eic_df = bind_rows(sig_matches,sig_matches,sig_matches,sig_matches)
eic_df$SpecFile = Sys.glob("data/QTOF/01082021/*GA6-002*")
eic_df = mutate(eic_df,
                condition = str_extract(SpecFile,"(R2A|ISP1)_(EA|BuOH)"),
                name_simple = FragmentSeq)
eic_df = lapply(1:nrow(eic_df),function(i) retrieveEIC(df=eic_df,idx=i))
eic_df = bind_rows(eic_df)  
plot_chromatogram(eic_df)
  
##combine_plots ####

p_combined= cowplot::plot_grid(p_ripgenes,
                               p_lassoseq,
                               p_spectra,
                               ncol=1, 
                               labels = c("A","B","C"))


ggsave(p_combined,
       filename = "figures/lassopeptide.svg",
       width=11,height=9)