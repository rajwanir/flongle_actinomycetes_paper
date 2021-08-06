require(xcms)
require(magrittr)
require(data.table)
require(tidyverse)
library(ggplot2)
library(cowplot)

##https://bioconductor.org/packages/devel/bioc/vignettes/xcms/inst/doc/xcms-lcms-ms.html
## add 0 at last line: sed -i '$ s/.*/\n0/g' data/QTOF/cfmid_compounds/actinorhodin/*.log
## convert dos to linux: dos2unix data/QTOF/cfmid_compounds/actinorhodin/*.log


get_ms_hits = function(file_path){
  ### import MS MS data ####
  msdata = readMSData(files = file_path, mode = "onDisk")
  # to check ms level spectra
  # table(msLevel(msdata))
  # msdata %>% filterMsLevel(2L) %>% precursorMz()
  cwp <- CentWaveParam(snthresh = 5, noise = 100, ppm = 10,
                       peakwidth = c(3, 30))
  msdata <- findChromPeaks(msdata, param = cwp)
  msspectra = chromPeakSpectra(msdata)
  # msdata = chromPeaks(msdata,ppm = 20)
  # msdata = as_tibble(msdata, rownames = "peak_id") 
  
  #roun of mz digits and get hits
  # msdata$mz_2 = format(round(msdata$mz, 4), nsmall = 4)
  # msdata_hits = msdata %>% filter(mz_2 %in% substructures$ExactMass_2)
  
  return(list(msspectra,msdata))
}

plot_spectra = function(input_spec){
  spect_df = data.frame(mz = as.numeric(mz(input_spec)),intensity = as.numeric(intensity(input_spec)))
  peaks_to_label = spect_df %>% filter(signif(mz,digits=3) %in% common_peaks) %>% arrange(-intensity) %>% head(10)
  spectra_plot = ggplot(spect_df, aes(x=mz,y=intensity)) +
    geom_bar(stat = "identity",width=3) +
    ggrepel::geom_text_repel(data=peaks_to_label,
              inherit.aes = T,
              aes(label=mz,x=mz,y=intensity)) +
    scale_x_continuous(breaks = seq(0,1000,by = 100), limits = c(0,1000)) +
    theme_linedraw() + theme(panel.grid = element_blank())
  return(spectra_plot)
}


### Get mz values for query functional group ####
functional_group = 74217  #3-Aminocoumarin
substructures = data.table::fread(paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsubstructure/cid/",
                                         functional_group,"/property/MolecularFormula,ExactMass,MonoisotopicMass,CanonicalSMILES/CSV?StripHydrogen=true&fullsearch=true"),
                                  header = T)
substructures$ExactMass_2 = format(round(substructures$ExactMass, 4), nsmall = 4)

#### data import ####
msFiles = Sys.glob("data/QTOF/B06-03/*.mzXML")

ms_complete = lapply(Sys.glob(msFiles),get_ms_hits)


cfmid_paths = Sys.glob("data/QTOF/pchid121491875_insilico/*.txt")
q_pred_specs = lapply(cfmid_paths,function(cfmid_path){
q_pred_spec = cfmR::cfm_predict_readfile(cfmid_path)
q_pred_spec$spec = lapply(q_pred_spec$spec, function(x)x %>% rename(intensity = int))
q_pred_spec$spec = lapply(q_pred_spec$spec, function(x) new("Spectrum2", mz = as.numeric(x$mz), intensity = as.numeric(x$intensity), centroided =T))
})


spectra = lapply(ms_complete,'[[',1) 
msdata = lapply(ms_complete,'[[',2) 

sp_comparision = 
lapply(q_pred_specs, function(q_pred_spec){
lapply(c("energy0","energy1","energy2"), function(energyLevel){
lapply(1:length(msFiles),function(y)
  lapply(spectra[[y]],function(x) compareSpectra(x,q_pred_spec[[energyLevel]],fun = "dotproduct")) %>% 
    as_tibble() %>% 
    pivot_longer(col = everything(), names_to = "peak_id") %>% 
    mutate(energy = energyLevel, sample = msFiles[y])
)
})
})


names(sp_comparision) = cfmid_paths
sp_comparision = lapply(sp_comparision,bind_rows)
lapply(sp_comparision,function(x)x[order(-x$value) %>% head(1),])


#ploting spectra ####
ex_spectra  = spectra[[2]][mcols(spectra[[2]])$peak_id == 'CP4970']
plot(ex_spectra$CP4970.F1.S3685,q_pred_specs[[3]]$energy2)
common_peaks = intersect(signif(mz(ex_spectra$CP4970.F1.S3685),digits = 3),signif(mz(q_pred_specs[[3]]$energy2),digits = 3))


hit_spectra = ex_spectra$CP4970.F1.S3685
hit_spectra = removePeaks(hit_spectra,t=25)
query_spectra = q_pred_specs[[3]]$energy2
common_peaks = intersect(signif(mz(hit_spectra),digits = 3),signif(mz(query_spectra),digits = 3))

# ggdraw() +
#   draw_image(image ="data/QTOF/compound_structures/102115546.png",  x = 1500, y = 10) +
#   draw_plot(plot_spectra(query_spectra))

spectra_plot=cowplot::plot_grid(plot_spectra(hit_spectra),plot_spectra(query_spectra),
                   ncol = 1, labels = c("A","B")) 

ggsave(spectra_plot,width=15,height = 6,dpi = 500,filename = "figures/spectra.png")

# p_spectra = ggplotify::as.ggplot(~plot(hit_spectra))
# p_spectra + labs(x=NULL)

##./print_structure --structure ../npatlas/mols/6038.mol --configs_dir ../../data/QTOF/dereplicator+_npatlas/work/configs --print_mgf_spectrum



## Drawing chemical structures ####
library("ChemmineR")
smiles = c(Rubradirin = 'CC1C=C(C(=O)C2=C(C(=CC3=C2C(=O)C4=C(C3=O)NCC(O4)(C(=O)C(C1OC(=O)C5=C(C(=CC(=N5)C(=O)NC6=C(C7=C(C=C(C=C7)OC)OC6=O)O)OC8CC(C(C(O8)C)OC)(C)[N+](=O)[O-])O)O)C)C)O)C')
system(sprintf("echo '%s' | obabel -ismi -omol2 > %s/smiles.sdf",smiles,tempdir()))
sdfset = read.SDFset(sprintf("%s/smiles.sdf",tempdir()))