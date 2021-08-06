# spec_name = '00028'
# sapply (30:40, function(spec_name)
# sprintf('/gpfs/gsfs11/users/rajwanir2/soil_metagenomics/npdtools2/bin/dereplicate "data/QTOF/moldiscovery_npatlas/work/spectra-000%s.mzXML" "npdtools2/npatlas" "npdtools2/npatlas/library.info" --configs_dir "data/QTOF/moldiscovery_npatlas/work/configs" --product_ion_thresh 0.02 --pm_thresh 0.02 --max_charge 2 --decoy_type 2 --decoy_fraction 1  --min_score 9.0 --fragmentation_tree --no_pvalue --min_num_comp 0 --fragment "data/QTOF/moldiscovery_npatlas/work/configs/configs/common/fragment_list_prob.txt" --precomputed_target_db_path data/QTOF/moldiscovery_npatlas/db_preproc/library.info.prob.mz.target.bin --precomputed_decoy_db_path data/QTOF/moldiscovery_npatlas/db_preproc/library.info.prob.mz.decoy.bin --skip_overcharged --online_null_model --max_rank 7 --molDiscovery data/QTOF/moldiscovery_npatlas/work/configs/configs/molDiscovery/model.online.charge2.rank7  > "data/QTOF/moldiscovery_npatlas/work/spectra-000%s_dereplicator.log"',
#         spec_name,spec_name)
# )



loadGNPSSpectra = function(gnps_id){
#retreives spectra from GNPS and returns as msnbase object
  #example: gnps_id = "CCMSLIB00004722223"
json_file = sprintf("https://gnps.ucsd.edu/ProteoSAFe/SpectrumCommentServlet?SpectrumID=%s",
                    gnps_id)
json_file = json_file %>% jsonlite::fromJSON() 
peak_list = json_file$spectruminfo$peaks_json %>% jsonlite::fromJSON(flatten = T) 
GNPSspectra = new("Spectrum2",
    mz = as.numeric(peak_list[,1]),
    intensity = as.numeric(peak_list[,2]),
    precursorMz = as.numeric(sort(json_file$annotations$Precursor_MZ,decreasing = T)[1]),
    precursorCharge = as.integer(json_file$annotations$Charge)
    )
return(GNPSspectra)
}

annotateSigmatches = function(folder_name) {
  #example: folder_name = "moldiscovery_npatlas"
  #read NPatlas database
  npatlas_metadata = read.csv(file = "npdtools/npatlas/NPatlas_01132020.tsv", sep = '\t')
  
  # data filtering after matches
  npatlas_matches = read_tsv(sprintf("data/QTOF/%s/significant_matches.tsv",folder_name))
  #add metadata from database
  npatlas_matches= merge(npatlas_matches,npatlas_metadata,
                       by.x = "Name",by.y = "npaid",
                       all.x = T,all.y = F)
  #filter compounds detected in media
  npatlas_matches = npatlas_matches[!npatlas_matches$Name %in%
                                    npatlas_matches[grepl(npatlas_matches$SpecFile, pattern = "MC-00"),]$Name
                                  ,]
  #create columns for media and solvent
  npatlas_matches$sample = str_extract(npatlas_matches$SpecFile,pattern = "(R2A|ISP1)_(.)+_(.)+\\.") %>% str_remove(pattern = "\\.")
  npatlas_matches = separate(npatlas_matches,col = "sample",into = c("medium","solvent","strain"), sep = '_')

  
  #npatlas_matches = npatlas_matches[!duplicated(npatlas_matches[c("Name","SpecFile")] ),]
  
return(npatlas_matches)
}

retrieveSpectra = function(df) {
  #given the dataframe of type npatlas matches returns a list
  #1. spectra 2. plot 3. ce
  dfSpectra = list()
  dfSpectra$spectra = lapply(1:nrow(df), function(idx) {
    spectra(filterPrecursorScan(
      readMSData(df$SpecFile[idx], mode = "onDisk"),
      df$Scan[idx])
    )[[2]]
  })
  dfSpectra$plot = lapply(dfSpectra$spectra,plot_spectra)
  dfSpectra$ce = lapply(dfSpectra$spectra,collisionEnergy) 
return(dfSpectra)
}

retrieveSpectraParallel = function(df) {
  #given the dataframe of type npatlas matches returns a list
  #1. spectra 2. plot 3. ce
  dfSpectra = list()
  dfSpectra$spectra = mclapply(1:nrow(df), function(idx) {
    spectra(filterPrecursorScan(
      readMSData(df$SpecFile[idx], mode = "onDisk"),
      df$Scan[idx])
    )[[2]]
  },
  mc.cores = 60)
  dfSpectra$plot = lapply(dfSpectra$spectra,plot_spectra)
  dfSpectra$ce = lapply(dfSpectra$spectra,collisionEnergy) 
  return(dfSpectra)
}

retrieveEIC = function(df=prism_matches,idx){
  
  mz_target= (df$PeptideMass[idx] + df$Charge[idx])
  rt_target=as.numeric(df$Retention)
  eic = MSnbase::chromatogram(MSnbase::readMSData(df$SpecFile[idx], mode = "onDisk"),
               mz = c( mz_target - (mz_target * 20 / 10^6),
                       mz_target + (mz_target * 20 / 10^6)),
              # rt = c(rt_target-120,rt_target+120),
              aggregationFun = "sum",
               rt = c(0,10*60),
               missing = 0
               )
  
  
  eic = tibble(rtime = rtime(eic[1,1])/60, 
               intensity = intensity(eic[1,1]),
               # name = df$name_simple[idx],
               # condition = df$condition[idx]
               )

  return(eic)
}


plot_chromatogram = function(chrom_df){
  p_chrom = ggplot(chrom_df,
         aes(x=rtime,y=intensity, 
             #color = condition
             )) +
    geom_line(size=0.5,
             # aes(linetype=name)
              ) +
    labs(x="retention time (min)")+
    theme_classic() +
    theme(legend.position = "right")
  
  return(p_chrom)
}

spec2df =  function(input_spec){
  spect_df = data.frame(mz = as.numeric(mz(input_spec)),
                        intensity = as.numeric(intensity(input_spec) * 1000),
                        precursor = as.numeric(precursorMz(input_spec) )) 
  return(spect_df)}

plot_spectra = function(input_spec){
  spect_df = spec2df(input_spec)
  spectra_plot = ggplot(spect_df, aes(x=mz,y=intensity)) +
    geom_bar(stat = "identity",width=2, fill = "black") +
    ggtitle(label = "Observed", 
            subtitle = paste('precursor m/z',spect_df$precursor[1],sep=' ')) +
    #geom_text(aes(y=intensity+200)) +
    # ggrepel::geom_text_repel(data=peaks_to_label,
    #           inherit.aes = T,
    #           aes(label=mz,x=mz,y=intensity)) +
    scale_x_continuous(breaks = seq(100,spect_df$precursor[1]+100,by = 100), limits = c(50,spect_df$precursor[1]+25)) +
    theme_linedraw() + theme(panel.grid = element_blank())
  return(spectra_plot)
}

label_peaks = function(spectra_plot, by = "gnps_common", common_idx = 1){
  
  if (by == "top_intensity") {
    spectra_plot$data = spectra_plot$data %>% mutate(mz_label = case_when(
      intensity > sort(intensity, decreasing = T)[10] ~ signif(mz, digits = 4),
      TRUE ~ NA_real_
    ))
    # peaks_to_label = spect_df %>% filter(signif(mz,digits=3) %in% common_peaks) %>% arrange(-intensity) %>% head(10)
  }
  
  else if (by == "none") {
    spectra_plot$data = spectra_plot$data %>% mutate(mz_label = NA_real_)
  }
  
  else if (by == "gnps_common") {
    spectra_plot$data = mutate(spectra_plot$data, rel_int = intensity/max(intensity))
    spectra_plot$data = spectra_plot$data %>% mutate(mz_label = case_when(
      signif(mz, digits = 4) %in% signif(mz(commonSpectra[[common_idx]]), digits = 4) & 
        rel_int > 0.01
        ~ signif(mz, digits = 4),
      TRUE ~ NA_real_
    ))
  }
  
  
  spectra_plot = spectra_plot + ggrepel::geom_text_repel(aes(label = mz_label,y=intensity+200),
                                                         color = "brown",
                                                         min.segment.length = 0,
                                                         force = 1.2,
                                                         max.overlaps=100)
  
  return(spectra_plot)
}

add_spectra_title = function(spectra_plot, type = "observed", common_idx = idx){
  
  precursor_mz = sprintf('precursor m/z:%.3f',
                        spectra_plot$data$precursor[1])
  
  
  if (type == "observed"){
    spectra_plot = spectra_plot + ggtitle(subtitle = NULL, 
                                          label = paste(precursor_mz,
                                                           "CE:",npatlas_matches_gnps[common_idx,'ce'],
                                                           "strain:",npatlas_matches_gnps[common_idx,'strain'],
                                                           "\nmedium:",npatlas_matches_gnps[common_idx,'medium'],
                                                           "solvent:",npatlas_matches_gnps[common_idx,'solvent'],
                                                           sep=' ')) +
      theme(panel.border = element_rect(color="DarkGreen", size = 2))
  }
  
  if (type == "gnps"){
    spectra_plot = spectra_plot + ggtitle(subtitle = NULL, 
                                          label = paste(precursor_mz,
                                                           "gnps_id:",npatlas_matches_gnps[common_idx,'gnps_ids_clean'],
                                                           "\ncompound:",npatlas_matches_gnps[common_idx,'compound_names'],
                                                           sep=' '))
  }
  
  
  spectra_plot = spectra_plot + theme(title = element_text(size = 7))
  
  return(spectra_plot)
}

smiles2monomers = function(smile, mode = "rban"){
  # smile = tryCatch(readLines(url(sprintf("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/%s/property/CanonicalSMILES/txt",
  #                                        smile))),
  #                  error=function(e) NULL)
  
  if (mode == "web"){
  monomers = tryCatch(readLines(url(sprintf("https://bioinfo.cristal.univ-lille.fr/norine/rest/s2m/graph/%s",
                                            smile))),
                      error=function(e) NULL)
  }
  
  if (mode == "rban"){
  dir = tempdir()
  file = basename(tempfile(pattern=as.character(Sys.getpid()) ))
  system(sprintf("java -jar rBAN-1.0.jar -inputSmiles \"%s\" -outputFolder %s/ -inputId a -outputFileName %s",
                 smile,dir,file),
         ignore.stdout = T) 
  monomers = paste(jsonify::from_json(sprintf("%s/%s",dir,file))$monomersNames,collapse = ",")  
  
  }
  
  return(monomers)}

to_cannonicalSmiles = function(smile, mode = "openbabel"){
 
   if (mode == "pubchem"){
  smile = tryCatch(readLines(url(sprintf("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/%s/property/CanonicalSMILES/txt",
                                         smile))),
                   error=function(e) NULL)
  }
  
  else if (mode == "openbabel"){
    smile = system(sprintf("obabel -:\"%s\" -ocan",smile),
                   intern = T)
  }
    
  return(smile)}   