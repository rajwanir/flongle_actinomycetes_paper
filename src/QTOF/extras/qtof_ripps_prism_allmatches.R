require(tidyverse)
require(MSnbase)
all_matches = read_tsv(file = "data/QTOF/metaminer_prism/all_matches.tsv")
idxseq = 1

all_matches = all_matches %>% arrange(`P-Value`) %>% head(idxseq) # taking only one lowest p-value hit
spectra = retrieveSpectra(all_matches)
eic = retrieveEIC(all_matches,idx = idxseq)
p_eic = plot_chromatogram(eic)
all_matches$FragmentSeq[idxseq]

ripseq = "TGRPQQPRSPYSCTRTTARS" 
  #change log for ripseq 
  #changed 11S to Y
  #changed TT to Q
plot(spectra$spectra[[1]], 
     y = ripseq,
     modifications=c(T =-18.010, 
                     Y =-56.03,#convert back Y to S (-76.0313) and loss of -20
                     Q =-7.01 # convert back Q to T (-27.01) and loss of -20
                     ), 
     type = c("y","b"),
     neutralLoss=NULL,
     # neutralLoss=defaultNeutralLoss(disableAmmoniaLoss=c("K", "N", "Q", "R")),
     # main = sprintf("MS/MS spectra precursor m/z %f",
     #                precursorMz(sample_plots[[1]]$spectra)),
      col="black"
) 
