require(tidyverse)
require(MSnbase)
source("src/QTOF/qtof_functions.R")
library(parallel)
library(cowplot)


npatlas_matches = annotateSigmatches("moldiscovery_npatlas")

#simplify gnps ids
npatlas_matches_gnps = npatlas_matches %>% filter(gnps_ids!="[]") 
npatlas_matches_gnps$gnps_ids_clean = npatlas_matches_gnps %>% select(gnps_ids) %>%
  sapply(.,function(x) str_extract(x,pattern = "CCMSLIB\\d+")) %>% 
  as.character()

#populate sample spectra
npatlas_matches_gnpsSpectra = retrieveSpectraParallel(npatlas_matches_gnps)
# npatlas_matches_gnps = npatlas_matches_gnps[!duplicated(npatlas_matches_gnps[c("compound_names")] ),]


#for polting selected:
npatlas_matches_gnps = read_tsv("data/QTOF/moldiscovery_npatlas/significant_matches_filtered_gnps.tsv")
npatlas_matches_gnps = npatlas_matches_gnps %>% filter(gnps_common > 4 & gnps_dotprod > 0.6) %>% arrange(gnps_dotprod) %>%  distinct(compound_names,.keep_all=T)
npatlas_matches_gnpsSpectra = retrieveSpectra(npatlas_matches_gnps)

#populate gnps spectra
gnps = list()
gnps$spectra = lapply(npatlas_matches_gnps$gnps_ids_clean,loadGNPSSpectra)
gnps$plot = lapply(gnps$spectra,plot_spectra)

#compare observed spectra wtih gnps 
npatlas_matches_gnps$gnps_common = sapply(1:nrow(npatlas_matches_gnps),function(idx){
  compareSpectra(npatlas_matches_gnpsSpectra$spectra[[idx]],
                 gnps$spectra[[idx]], fun = "common")
})


npatlas_matches_gnps$gnps_dotprod = sapply(1:nrow(npatlas_matches_gnps),function(idx){
  compareSpectra(npatlas_matches_gnpsSpectra$spectra[[idx]],
                 gnps$spectra[[idx]], fun = "dotproduct")
})


commonSpectra = lapply(1:length(gnps$spectra), function(idx){
  consensusSpectrum(list(gnps$spectra[[idx]],
                         npatlas_matches_gnpsSpectra$spectra[[idx]]),
                    mzd = 0.1,ppm = 50,minProp = 1)
})



# label common peaks

gnps$plot = lapply(1:length(gnps$plot),function(idx){
  label_peaks(gnps$plot[[idx]], common_idx = idx, by = "gnps_common")
})

npatlas_matches_gnpsSpectra$plot = lapply(1:length(gnps$plot),function(idx){
  label_peaks(npatlas_matches_gnpsSpectra$plot[[idx]], common_idx = idx, by = "gnps_common")
})

#add ce annotation
npatlas_matches_gnps$ce = as.numeric(npatlas_matches_gnpsSpectra$ce)

# add spectra title
gnps$plot = lapply(1:length(gnps$plot),function(idx){
  add_spectra_title(gnps$plot[[idx]], common_idx = idx, 
                    type = "gnps")
})

npatlas_matches_gnpsSpectra$plot = lapply(1:length(gnps$plot),function(idx){
  add_spectra_title(npatlas_matches_gnpsSpectra$plot[[idx]], common_idx = idx,
                    type = "observed")
})




#plot the top matches
#add scanids to plot names
names(gnps$plot) = npatlas_matches_gnps$Scan
names(npatlas_matches_gnpsSpectra$plot) = npatlas_matches_gnps$Scan

npatlas_matches_gnpsBest = npatlas_matches_gnps %>% arrange(compound_names,-gnps_common) %>% distinct(compound_names, .keep_all = T)

# plot all combined
p_gnps_combine_theme = theme(axis.text.y = element_blank(),
                             axis.title.y = element_blank(),
                             axis.ticks.y = element_blank())


combined_plots = lapply(1:nrow(npatlas_matches_gnps),function(scanid)
  plot_grid(gnps$plot[[scanid]]+p_gnps_combine_theme,
            npatlas_matches_gnpsSpectra$plot[[scanid]]+ p_gnps_combine_theme)
)


plot_grid(plotlist = combined_plots,ncol = 1,labels = toupper(letters[1:4]),hjust = 0.5)



#msnbase plot

combined_plots = lapply(1:nrow(npatlas_matches_gnps), function(i) {
  par(mar = c(4, 4.8, 2, 0.1))
  gnps$spectra[[i]]@centroided = T
  npatlas_matches_gnpsSpectra$spectra[[i]]@centroided = T
  plot(gnps$spectra[[i]],npatlas_matches_gnpsSpectra$spectra[[i]],
       #tolerance = 0.00050
       cex=0.5)
  title(sprintf("%s detected in %s",npatlas_matches_gnps$compound_names[i],npatlas_matches_gnps$strain[i]),
        adj=0.5,line=1,
        font.main = 2,cex.main=1)
  x = recordPlot() 
  plot.new()
  return(x)
  })

combined_plots = plot_grid(plotlist = combined_plots,
                           labels = toupper(letters[1:4]),
                           ncol = 1,
                           scale = 0.9)
ggsave(combined_plots,filename = "figures/gnps_validated-v.svg",
       width=7,height=10)


for (i in seq(1,length(combined_plots),by=6)){
combined_plot_split = plot_grid(plotlist = combined_plots[i:(i+5)], ncol = 1)
ggsave(filename = sprintf("figures/gnps_validated_%d_%d.svg",i,(i+5)),
       combined_plot_split,
       width = 11,
       height = 8)
  }


# ploting all spectra not just best
# combined_plots = lapply(seq(1, length(npatlas_matches_gnpsSpectra$plot),
#                             6), function(idx)
#                               cowplot::plot_grid(
#                                 cowplot::plot_grid(plotlist = npatlas_matches_gnpsSpectra$plot[idx:(idx +
#                                                                                                       5)],
#                                                    ncol = 1),
#                                 cowplot::plot_grid(plotlist = gnps$plot[idx:(idx + 5)],
#                                                    ncol = 1)
#                               ))
# 
# 
# pdf(file="figures/gnps_validated_matches.pdf", onefile = TRUE, height = 8, width = 11)
# 
# for (i in 1:length(combined_plots)) {
#   print(combined_plots[[i]])
#  # grid::grid.newpage()
# }
# 
# dev.off()


# export write tables

write.table(npatlas_matches, 
            file = "data/QTOF/moldiscovery_npatlas/significant_matches_filtered.tsv",
            quote = F,row.names = F, sep = '\t')

write.table(npatlas_matches_gnps, 
            file = "data/QTOF/moldiscovery_npatlas/significant_matches_filtered_gnps.tsv",
            quote = F,row.names = F, sep = '\t')