library(tidyverse)
library(parallel)
library(doParallel)


## construct sdf db for structure searching #####
prism_db = read.csv("library_metadata.csv")
prism_db$foropenbabel = sprintf("-:\"%s %s\"",
        prism_db$smiles,prism_db$name)
# write.csv(prism_db[,5],file = "for_openbabel.csv",
#           row.names = F, quote = F)
# obabel --gen3D -ismi for_openbabel.smi -O prism_combined.sd
# obabel prism_combined.sdf -ofs

#read npatlas matches unique
npatlas_metadata = read.csv(file = "npdtools/npatlas/NPatlas_01132020.tsv", sep = '\t')
npatlas_matches = read_tsv(sprintf("data/QTOF/%s/significant_unique_matches.tsv","moldiscovery_npatlas_GB4-14"))
npatlas_matches= merge(npatlas_matches,npatlas_metadata,
                       by.x = "Name",by.y = "npaid",
                       all.x = T,all.y = F)
npatlas_matches$sample = str_extract(npatlas_matches$SpecFile,pattern = "(R2A|ISP1)_(.)+_(.)+\\.") %>% str_remove(pattern = "\\.")
npatlas_matches = separate(npatlas_matches,col = "sample",into = c("medium","solvent","strain"), sep = '_')


## structure search ####
require(parallel)
npatlas_matches$prismMatch = mclapply(npatlas_matches$SMILES, function(smile) {
  system(
    sprintf(
      "obabel npdtools/prism_GB4-14_pb/prism_combined.fs -osmi -s\"%s\" -at0.85",
      smile
    ),
    intern = T,
    ignore.stderr = T
  )
}, mc.cores = 70)

#number of matches
npatlas_matches$n_prismSMILES = sapply(npatlas_matches$prismMatch,length) 
npatlas_matches = npatlas_matches %>% filter(n_prismSMILES >0)

npatlas_matches$prismMatch = unlist(npatlas_matches$prismMatch) %>% .[1:5]
write.csv(npatlas_matches,
          file = "data/QTOF/moldiscovery_npatlas_GB4-14/significant_unique_matches_85sim_prism.csv")

#need to match if same strain
npatlas_matches$prismMatch = sapply(npatlas_matches$prismMatch, function(x)
  paste(x, collapse = ','),USE.NAMES =F)
npatlas_matches$same_strain = sapply(1:nrow(npatlas_matches), function(idx)
  str_detect(npatlas_matches$prismMatch[idx],
             npatlas_matches$strain[idx]))

# filter only same strain matches
npatlas_matches = npatlas_matches %>% filter(same_strain ==T)



npatlas_matches$n_prismStrains = sapply(npatlas_matches$prismMatch,function(match) { 
  str_extract_all(match,pattern = "[A-Z]+[0-9]+-[0-9]+") %>%
  .[[1]] %>% 
  unique() %>% length()
},USE.NAMES =F)


npatlas_StrainSummary = npatlas_matches %>% group_by(compound_names) %>% distinct(strain) %>% summarise(
  n_npatlasStrains = n(),
  npatlas_OtherStrains = paste(strain, collapse = ',')
)
npatlas_matches = merge(npatlas_matches,npatlas_StrainSummary,
      all.x = T,all.y = F,
      by.x = "compound_names",by.y = "compound_names")

# npatlas_matches$n_npatlasStrains = sapply(npatlas_matches$compound_names,function(compound){
#   sum(str_count(unique(paste(npatlas_matches$compound_names,
#                   npatlas_matches$strain)),
#             pattern = compound))
# },USE.NAMES =F)

npatlas_matches= mutate(npatlas_matches,
            perc_expressed = n_npatlasStrains/n_prismStrains)

npatlas_matches = npatlas_matches %>% 
  select(perc_expressed,n_prismStrains,n_npatlasStrains,n_prismSMILES,strain,compound_names,Name,SMILES,prismMatch) %>%
  distinct() %>% 
  filter(perc_expressed >0.2) 


##retain only the strain match 
npatlas_matches$smiles_SameStrain = sapply(1:nrow(npatlas_matches),function(idx)
str_split(npatlas_matches$prismMatch[idx],',') %>% .[[1]] %>% 
  .[str_detect(.,
               sprintf("(.)*%s(.)*",npatlas_matches$strain[idx]))] %>% 
  str_split(.,'\t')
,USE.NAMES =F)


##calculate score for smiles from same strain
npatlas_matches$scores_SameStrain = sapply(1:nrow(npatlas_matches),function(idx){
  
    sapply(1:length(npatlas_matches$smiles_SameStrain[[idx]]),
           function(idx2) 
           system(sprintf("obabel -:\"%s\" -:\"%s\" -ofpt",
                          npatlas_matches$SMILES[idx],
                          npatlas_matches$smiles_SameStrain[[idx]][[idx2]][1]),
                  intern = T)
           ,USE.NAMES =F) %>% 
  str_extract("\\d+\\.\\d+") %>% 
    na.omit() %>% as.numeric()
  
})

# need to retain only the maximum



npatlas_matches = separate(npatlas_matches,
                                col = "smiles_SameStrain",
                             into = c("smiles_SameStrain","prism_id"),
                             sep = '\t')


bind_rows(
select(npatlas_matches,c(id=prism_id,smiles=smiles_SameStrain)),
select(npatlas_matches,c(id=Name,smiles=SMILES))
) %>% jsonlite::toJSON() %>% 
  write(file="smiles.json")
  
#java -jar rBAN-1.0.jar -inputFile smiles.json -outputFolder smiles_out -imgs -imgsFolderName img_out

npatlas_matches$similarity_score = sapply(1:nrow(npatlas_matches),function(idx){
  system(sprintf("obabel -:\"%s\" -:\"%s\" -ofpt",
                 npatlas_matches$SMILES[idx],
                 npatlas_matches$smiles_SameStrain[idx]),
         intern = T) %>%  
   str_extract("\\d+\\.\\d+") %>% 
   na.omit() %>% as.numeric()
})


#convert imgs to grayscale
imgs_path  = c(sprintf("smiles_outimg_out/compounds/%s.png",
                       npatlas_matches$prism_id),
               sprintf("smiles_outimg_out/compounds/%s.png",
                       npatlas_matches$Name))

greyscale_image =function(img_path){
  img = magick::image_read(img_path)
  img = magick::image_convert(img,colorspace = 'gray')
  magick::image_write(img,img_path)
}

lapply(imgs_path,function(path) tryCatch(greyscale_image(path),
                                         error=function(e) NULL))

#combine all imgs
combine_structure_imgs =  lapply(1:nrow(npatlas_matches),function(idx){
  tryCatch(plot_grid(
    ggdraw() + draw_image(sprintf("smiles_outimg_out/compounds/%s.png",
                                  npatlas_matches$Name[idx])),
    ggdraw() + draw_image(sprintf("smiles_outimg_out/compounds/%s.png",
                                  npatlas_matches$prism_id[idx])))
    ,error=function(e) NULL)
  
})

combine_structure_imgs = compact(combine_structure_imgs)

pdf("figures/BGC_structurematches.pdf")  
for (i in seq(1,10,5)){
  cowplot::plot_grid(plotlist = combine_structure_imgs[i:(i+4)],
                     ncol=1)
}
dev.off()

# p_eic =  p_eic + scale_x_continuous(breaks = seq(0,6,by=0.5)) + 
#   ggtitle("EIC m/z 218.14") +
#   theme(plot.title = element_text(hjust=0.5),legend.position = "none") + 
#   geom_vline(xintercept = 3.471,linetype = "dashed") 
# 
# 
# ggsave(p_eic,filename = "figu")

##obabel -:"O=CCC(=O)[C@@H](NC(=O)[C@H]1N(C(=O)CCCCCCCCC)CCC1)[C@@H](O)C" -:"CCC(C)C(CO)NC(=O)C1CCCN1C(=O)C(C)(C)NC(=O)C(CC(C)C)NC(=O)C(CC(C)C)NC(=O)C2CCCN2C(=O)C(C)(C)NC(=O)C(C(C)CC)NC(=O)C(C(C)CC)NC(=O)C(CC(=O)N)NC(=O)C(C)(C)NC(=O)C" -ofpt