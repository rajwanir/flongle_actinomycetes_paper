require(xcms)
require(magrittr)
require(data.table)
require(tidyverse)
library(ggplot2)
library(cowplot)
# library(parallel)
# library(foreach)
# library(doParallel)


##https://bioconductor.org/packages/devel/bioc/vignettes/xcms/inst/doc/xcms-lcms-ms.html
## add 0 at last line: sed -i '$ s/.*/\n0/g' data/QTOF/cfmid_compounds/actinorhodin/*.log
## convert dos to linux: dos2unix data/QTOF/cfmid_compounds/actinorhodin/*.log

## functions #################
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
  spect_df = data.frame(mz = as.numeric(mz(input_spec)),
                        intensity = as.numeric(intensity(input_spec) * 1000),
                        precursor = as.numeric(precursorMz(input_spec) )) 
  spect_df = spect_df %>% mutate(mz_label = case_when(intensity > sort(intensity,decreasing = T)[10] ~ signif(mz,digits=4), TRUE ~ NA_real_))
  # peaks_to_label = spect_df %>% filter(signif(mz,digits=3) %in% common_peaks) %>% arrange(-intensity) %>% head(10)
  
  spectra_plot = ggplot(spect_df, aes(x=mz,y=intensity, label = mz_label)) +
    geom_bar(stat = "identity",width=2, fill = "black") +
    ggtitle(label = "Observed", 
            subtitle = paste('precursor m/z',spect_df$precursor[1],sep=' ')) +
    geom_text(aes(y=intensity+200)) +
    # ggrepel::geom_text_repel(data=peaks_to_label,
    #           inherit.aes = T,
    #           aes(label=mz,x=mz,y=intensity)) +
    scale_x_continuous(breaks = seq(0,spect_df$precursor[1]+100,by = 100), limits = c(0,spect_df$precursor[1]+100)) +
    theme_linedraw() + theme(panel.grid = element_blank())
  return(spectra_plot)
}

to_cannonicalSmiles = function(smile){
  #smile = "C1C[C@@H](C(=O)N(C1)O)NC(=O)CNC(=O)[C@@H](CCCN=C(N)NC(=O)C2=C(C(=CC=C2)O)O)NC(=O)C3=C(C(=CC=C3)O)O"
  smile = tryCatch(readLines(url(sprintf("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/%s/property/CanonicalSMILES/txt",
                                         smile))),
                   error=function(e) NULL)
                   return(smile)}                   


smiles2monomers = function(smile){
  #smile = "C1C[C@@H](C(=O)N(C1)O)NC(=O)CNC(=O)[C@@H](CCCN=C(N)NC(=O)C2=C(C(=CC=C2)O)O)NC(=O)C3=C(C(=CC=C3)O)O"
  smile = tryCatch(readLines(url(sprintf("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/%s/property/CanonicalSMILES/txt",
                        smile))),
                   error=function(e) NULL)
  monomers = tryCatch(readLines(url(sprintf("https://bioinfo.cristal.univ-lille.fr/norine/rest/s2m/graph/%s",
                                   smile))),
                      error=function(e) NULL)
return(monomers)}

smiles_compare = function(smile1, smile2){
  # smile1='Cn1c(=O)c2c(ncn2C)n(C)c1=O'
  # smile2='Cn1c(=O)c2c(ncn2C)n(C)c1=O'
  comparision = tryCatch(readLines(url(sprintf("https://gnps-structure.ucsd.edu/structuresimilarity?smiles1=%s&smiles2=%s",
                smile1,smile2))),
                error=function(e) NULL)
  
  return(comparision)
}


### Get mz values for query functional group ####
functional_group = 74217  #3-Aminocoumarin
substructures = data.table::fread(paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsubstructure/cid/",
                                         functional_group,"/property/MolecularFormula,ExactMass,MonoisotopicMass,CanonicalSMILES/CSV?StripHydrogen=true&fullsearch=true"),
                                  header = T)
substructures$ExactMass_2 = format(round(substructures$ExactMass, 4), nsmall = 4)

#### data import ####

#read and plot sample spectra
sample_plots = lapply(1:nrow(npatlas_matches),function(idx){
  spectra = list()
  spectra$spectra = spectra(filterPrecursorScan(readMSData(npatlas_matches$SpecFile[idx],mode = "onDisk"),
                                         npatlas_matches$Scan[idx])
                     )[[2]]
  spectra$ce = collisionEnergy(spectra$spectra)
  spectra$plot = plot_spectra(spectra$spectra)
  return(spectra)
})

# read and plot predicted npatlas spectra
npatlas_plots = lapply(1:nrow(npatlas_matches),function(idx){
  spectra = list()
  spectra$spectra = spectra(readMgfData(paste0("data/QTOF/dereplicator+_npatlas/work/sig_mgfs/",
                                        npatlas_matches$Name[idx],".mgf")))
  names(spectra$spectra) = c("10","20","40")
  spectra$plot = lapply(spectra$spectra,plot_spectra)
  return(spectra)
  })


# comparision
comparision = sapply(1:nrow(npatlas_matches),function(idx){
  compareSpectra(sample_plots[[idx]]$spectra,
                 npatlas_plots[[idx]]$spectra[
                   names(npatlas_plots[[idx]]$spectra) %in% sample_plots[[idx]]$ce
                 ][[1]]
                 , fun = "dotproduct")
})


combine_spectra_plot = function(idx){
  p = cowplot::plot_grid(sample_plots[[idx]]$plot,npatlas_plots[[idx]]$plot$`10`,
                           ncol = 1, labels = c("A","B"))
  p = cowplot::plot_grid(ggdraw()+draw_label(paste(npatlas_matches$compound_names[idx],
                                               paste("Strain:",npatlas_matches$strain[idx],"Solvent:",npatlas_matches$solvent[idx]),sep = '\n')
  ),
  p, nrow = 2, rel_heights = c(0.1,1))
return(p)
  }



pdf(file = "figures/spectra_comparision_01142021.pdf")
combine_spectra_plot(43)
combine_spectra_plot(97)
combine_spectra_plot(34)
dev.off()

# read and plot experimental spectra
experimental_plots = lapply(c("NPA016182","NPA008538"),function(name){
  plot_spectra(spectra(readMgfData(paste0("data/QTOF/dereplicator+_npatlas/work/sig_mgfs_experimental//",
                                          name,".mgf")))[[1]])
})

# sample_npatlas_cplots = lapply(1:nrow(npatlas_matches),function(idx){
#   cowplot::plot_grid(sample_plots[[idx]],npatlas_plots[[idx]],
#                    ncol = 1, labels = c("A","B")) 
# })

names(sample_plots) = npatlas_matches$Name
names(npatlas_plots) = npatlas_matches$Name
names(experimental_plots) = c("NPA016182","NPA008538")

sample_npatlas_cplots = lapply(names(experimental_plots),function(name){
  cowplot::plot_grid(sample_plots[[name]],
                     experimental_plots[[name]],
                     npatlas_plots[[name]],
                   ncol = 1, labels = c("A","B","C"))
})

# peak intesities differ significantly based on collision energy
predicted_energy_plots  = lapply(1:3,function(idx){
  (spectra(readMgfData(paste0("data/QTOF/dereplicator+_npatlas/work/sig_mgfs/",
                                          "NPA008538",".mgf")))[[idx]])
})



ms_complete = lapply(Sys.glob(msFiles),get_ms_hits)



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


### smiles2monomers ####

monomers = sapply(npatlas_matches$SMILES,smiles2monomers)

npatlas_matches$cannonicalSmiles = sapply(npatlas_matches$SMILES,to_cannonicalSmiles) %>% as.character()
prism = separate(prism,col = "name",into = c("date","barcode","strain","cluster"), sep = '_')
prism$cannonicalSmiles = sapply(prism$smiles,to_cannonicalSmiles) %>% as.character()
prism_copy = prism # saving a copy
prism = prism %>% filter(!cannonicalSmiles %in% c("NULL",""))

comparision = sapply(npatlas_matches$cannonicalSmiles,function(x){
  np_strain = npatlas_matches %>% filter(cannonicalSmiles == x) %>% select(strain) %>% as.character()
  prism_strain_subset = prism %>% filter(strain == np_strain)
  combination = crossing(x,prism_strain_subset$cannonicalSmiles)
  
  if(nrow(combination) > 0){ # check if strain data is available
  combination$comparision = sapply(1:nrow(combination),function(idx) 
    smiles_compare(combination[idx,1],combination[idx,2]))
  combination$comparision = as.character(combination$comparision)
  combination = arrange(combination,desc(comparision)) %>% filter(!comparision == "NULL") %>% .[1,]
  combination$strain = np_strain
  #cbind(combination,prism_strain_subset %>% filter(smiles == combination$`prism$smiles`))
  }
  return(combination)})


## Drawing chemical structures ####
library("ChemmineR")
smiles = c(Rubradirin = 'CC1C=C(C(=O)C2=C(C(=CC3=C2C(=O)C4=C(C3=O)NCC(O4)(C(=O)C(C1OC(=O)C5=C(C(=CC(=N5)C(=O)NC6=C(C7=C(C=C(C=C7)OC)OC6=O)O)OC8CC(C(C(O8)C)OC)(C)[N+](=O)[O-])O)O)C)C)O)C')
system(sprintf("echo '%s' | obabel -ismi -omol2 > %s/smiles.sdf",smiles,tempdir()))
sdfset = read.SDFset(sprintf("%s/smiles.sdf",tempdir()))