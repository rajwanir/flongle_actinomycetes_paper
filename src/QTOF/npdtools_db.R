  database_name = 'prism_similargenome'
  
  construct_db = function(database){
  library_metadata = read.csv(sprintf('npdtools/%s/library_metadata.csv',database))
  library_metadata$cmds = sapply(1:nrow(library_metadata), function(idx) {
  sprintf("obabel --gen3D -:\"%s\" -omol -x3 > mols/%s.mol",
                 library_metadata$smiles[idx],
                 library_metadata$name[idx]
                 )
  })
  
  write.table(library_metadata$cmds,
              file = sprintf("npdtools/%s/openbabel_cmds.swarm",database),
              quote = F,row.names = F,col.names = F)
  
  
  library_metadata$molfile = paste0("mols/",library_metadata$name,".mol")
  library_metadata$metadata = "NA"
  
  write.table(library_metadata[c('molfile','name','mass',
                                 'n_aminoacids','metadata')],
              file = sprintf('npdtools/%s/library.info',database),
              col.names = F,
              row.names = F,
              quote = F)
  
  write.table(library_metadata['smiles'],
              file = sprintf('npdtools/%s/library.smiles',database),
              col.names = F,
              row.names = F,
              quote = F)
  
  }
  
  
  construct_db(database_name)


# prism ################
require(tidyverse)

json_files = paste0("data/prism/",list.files("data/prism/",))
prism = lapply(json_files,function(x) {
  #list of dfs, each df corresponds to cluster
  prismdf1 = jsonlite::fromJSON(x,flatten = T,simplifyDataFrame =T)$prism_results$clusters$predicted_molecule_masses
  # #filter empty dfs/clusters with no smiles
  # prismdf1 = prismdf1[sapply(prismdf1,function(x) dim(x)[1] > 0)]
  # #1:1 for smiles and monoisotopic masses
  # prismdf1 = lapply(prismdf1,function(prismdf){
  #   prismdf$smiles = sapply(prismdf$smiles,"[[",1)
  #   prismdf$monoisotopic_mass = sapply(prismdf$monoisotopic_mass,"[[",1)
  # return(prismdf)})
  # append the rownames to include cluster and sample idx
  for (idx in seq(1, length(prismdf1), 1)) {
    
    if (nrow(prismdf1[[idx]]) > 0){
      
      #smiles is a list
      names(prismdf1[[idx]]$smiles) = prismdf1[[idx]]$monoisotopic_mass
      prismdf1[[idx]] = tibble(
        smiles = unlist(prismdf1[[idx]]$smiles),
        monoisotopic_mass = names(smiles))
      
      
      #adding cluster and smiles number to rowname
      rownames(prismdf1[[idx]]) = paste(x,
                                      paste0("c",idx),
                                      paste0("s",rownames(prismdf1[[idx]]) ),
                                      sep = "_")

    }
  }
  #single df per sample including all clusters
  prismdf1 = bind_rows(prismdf1, .id = "name")
  
  return(prismdf1)
  })

# names(prism) = 


prism = bind_rows(prism, .id = "name")
prism$name = stringr::str_remove_all(rownames(prism),pattern = ".gbk.json|data/prism/")
prism$mass = as.numeric(prism$monoisotopic_mass)
prism$n_aminoacids = 0


prism = prism %>% filter(as.numeric(mass) > 180) %>% select(-c(monoisotopic_mass))

write.table(prism,
            file = sprintf("npdtools/%s/library_metadata.csv",
                           database_name),
            quote = F, row.names = F,sep = ',')



# prism varquest ###

prism_db = read.csv("npdtools/prism/library_metadata.csv")

# smiles2monomores
# ...


prism_db = prism_db %>%
  filter(!monomors %in% c("NULL","character(0)"))

prism_db$n_monomers = str_count(prism_db$monomors,",")

prism_db = prism_db %>%
  filter(n_monomers > 3)

#get amino acids from noraine db
noraine_monomers = jsonify::from_json(url('https://bioinfo.cristal.univ-lille.fr/norine/rest/monomers/flat/json'),
                                      simplify = T)
noraine_monomers = bind_rows(noraine_monomers)
noraine_monomers = filter(noraine_monomers,
                          type == "NRPS")
aa = paste(as.character(noraine_monomers$code),collapse = '|')


#count and filter n_aminoacids < 3
prism_db$n_aminoacids = str_count(prism_db$monomors,aa)
prism_db = filter(prism_db,
                  n_aminoacids > 3)


write.csv(prism_db,file = "npdtools/prism_varquest/library_metadata.csv",
          row.names = F)

## prism varquest with cannonical ####
prism_db = read.csv('npdtools/prism/library_monomors.csv')
prism_db = prism_db %>% filter(monomors == "NULL")

prism_db$canonical_smiles = mclapply(prism_db$smiles,to_cannonicalSmiles,mc.cores = 71)
prism_db$canonical_smiles = as.character(prism_db$canonical_smiles)   
prism_db$canonical_smiles = str_remove(prism_db$canonical_smiles,'\t')


prism_db$monomors = mclapply(prism_db$canonical_smiles,smiles2monomers,
                             mc.cores = 40)
prism_db$monomors = as.character(prism_db$monomors)

prism_db$smiles = prism_db$canonical_smiles

#repeat steps as for orignal
#...

prism_exist = read.csv('npdtools/prism_varquest/library_metadata.csv')
prism_db = bind_rows(prism_exist,prism_db)





## npatlas ####

npatlas_metadata = read.csv(file = "npdtools/npatlas/NPatlas_01132020.tsv", sep = '\t')

#retriving by indexing in sdf file is faster, but save by compound id.
npatlas_metadata$npatlas_cmds = sapply(1:nrow(npatlas_metadata),function(idx){
sprintf("obabel NPAtlas_download.sdf --gen3D -omol -x3 -f%s -l%s > %s",
        idx,
        idx,
        paste0("mols/",npatlas_metadata$compound_id[idx],".mol") )
})

# write.table(npatlas_metadata$compound_smiles,
#             file = "npdtools/npatlas/library.smiles",
#             quote = F,row.names = F,col.names = F)
#swarm --maxrunning=10000 --merge-output --partition=quick --time=2:00 -g 1 --module=OpenBabel/3.0.0 -f openbabel_cmds.swarm


npatlas_metadata = npatlas_metadata %>% mutate(molfile =  paste0("mols/",compound_id,".mol"), n_aminoacids = 0, metadata = NA)
npatlas_libraryinfo = npatlas_metadata %>% select(molfile,npaid,compound_accurate_mass,n_aminoacids,metadata)

# write.table(npatlas_libraryinfo,
#             file = "npdtools/npatlas/library.info",
#             quote = F,row.names = F,col.names = F)

## data filtering after matches
npatlas_matches = read_tsv("data/QTOF/moldiscovery_npatlas/significant_matches.tsv")
#add metadata from database
npatlas_matches= merge(npatlas_matches,npatlas_metadata,
      by.x = "Name",by.y = "npaid",
      all.x = T,all.y = F)
#filter compounds detected in media
npatlas_matches = npatlas_matches[!npatlas_matches$Name %in%
npatlas_matches[grepl(npatlas_matches$SpecFile, pattern = "MC-00"),]$Name
,]

npatlas_matches$sample = str_extract(npatlas_matches$SpecFile,pattern = "(R2A|ISP1)_(.)+_(.)+\\.") %>% str_remove(pattern = "\\.")
npatlas_matches = separate(npatlas_matches,col = "sample",into = c("medium","solvent","strain"), sep = '_')

#npatlas_matches = npatlas_matches[!duplicated(npatlas_matches[c("Name","SpecFile")] ),]



write.table(npatlas_matches, 
                 file = "data/QTOF/moldiscovery_npatlas/significant_matches_filtered.tsv",
                 quote = F,row.names = F, sep = '\t')
#slurm jobs to predict spectra
npatlas_matches$cmd_cfmid = sprintf("source /gpfs/gsfs11/users/rajwanir2/conda/etc/profile.d/conda.sh && conda activate cfm && cfm-predict \'%s\' 0.001 src/metab_ce_cfm/param_output0.log src/metab_ce_cfm/param_config.txt 1 sig_mgfs/%s.mgf",
        npatlas_matches$compound_inchi,
        npatlas_matches$Name)

write.table(npatlas_matches$cmd_cfmid, 
            file = "cfmid_sig.swarm",
            quote = F,row.names = F,col.names = F)

# puchchemid to exact mass############

# metadata = read.csv(file = "npdtools/mibig/library_metadata.csv")
# metadata = metadata[!is.na(metadata$pchids),]
# metadata$mass = sapply(metadata$pchids,
# FUN = function(id){
# id_mass = readLines(url(sprintf("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/property/ExactMass/txt",
#         id)))
#           Sys.sleep(1)
#           return(id_mass)},
# simplify = T
# )
# 
# write.table(metadata, file = "npdtools/mibig/library_metadata.csv",
#             quote = F, row.names = F)