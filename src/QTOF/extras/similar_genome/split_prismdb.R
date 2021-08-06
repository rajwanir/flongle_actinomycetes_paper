database='prism_similargenome_ind'
copy_mols_from = 'prism_similargenome'
library_metadata = read.csv(sprintf('npdtools/%s/library_metadata.csv',database))
library_metadata$group = str_extract(library_metadata$name,'GCF_\\d+')
groups = unique(library_metadata$group)
library_metadata = library_metadata %>% group_by(group) %>% group_split()
names(library_metadata) = groups

#create group folders
sapply(groups,function(x) dir.create(sprintf('npdtools/%s/%s',
                                        database,x)))

#write group metadata file
sapply(groups,function(x)
write.table(select(library_metadata[x][[1]],-group),
            file = sprintf("npdtools/%s/%s/library_metadata.csv",
                           database,x),
            quote = F, row.names = F,sep = ',')
)

#create mols folder
sapply(groups,function(x) dir.create(sprintf('npdtools/%s/%s/mols',
                                             database,x)))
#copy molecules
sapply(groups,function(eachgroup)
sapply(library_metadata[eachgroup][[1]]$name, function(molecule)
       file.copy(from = sprintf("npdtools/%s/mols/%s.mol",
                                copy_mols_from,molecule),
                 to = sprintf("npdtools/%s/%s/mols/%s.mol",
                             database,eachgroup,molecule)
       )
)
)

#write additional info files
sapply(groups,function(eachgroup) {
  
  library_metadata_sbset = library_metadata[eachgroup][[1]]
  
  library_metadata_sbset$molfile = paste0("mols/",library_metadata_sbset$name,".mol")
  library_metadata_sbset$metadata = "NA"
  
  write.table(library_metadata_sbset[c('molfile','name','mass',
                                       'n_aminoacids','metadata')],
              file = sprintf('npdtools/%s/%s/library.info',database,eachgroup),
              col.names = F,
              row.names = F,
              quote = F)
  
  write.table(library_metadata_sbset['smiles'],
              file = sprintf('npdtools/%s/%s/library.smiles',database,eachgroup),
              col.names = F,
              row.names = F,
              quote = F)
  
}
)

#
mash_out$ref_acc = str_remove_all(mash_out$ref_acc,pattern = "\\.\\d")
sapply(1:nrow(mash_out),function(idx)
system(sprintf("bash src/QTOF/similar_genome/npdtools_ind.sh %s %s",
               mash_out$ref_acc[idx],mash_out$strain[idx]))
)
