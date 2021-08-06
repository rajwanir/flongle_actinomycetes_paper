library(tidyverse)

#read prism unique hits and tidy it
prism_matches_u = read_tsv("data/QTOF/moldiscovery_prism_similargenome/significant_matches.tsv")
prism_matches_u$prism_filepath = str_extract(prism_matches_u$Name,pattern = "(.)+.json")
prism_matches_u$prism_filepath = sprintf("data/similar_genome/prism_out/%s",prism_matches_u$prism_filepath)
prism_matches_u$cid = str_extract(prism_matches_u$Name,"_c[0-9]+") %>% str_remove("_c")
prism_matches_u$rgenome_filepath = str_extract(prism_matches_u$prism_filepath, pattern = "GCF_[0-9]+.[0-9]+")
prism_matches_u$rgenome_filepath = sprintf("data/similar_genome/ref_genomes/fasta/%s_genomic_refseq.fna",prism_matches_u$rgenome_filepath)
prism_matches_u$cname = sprintf("%s_c%s",
                                str_extract(prism_matches_u$prism_filepath, pattern = "GCF_[0-9]+.[0-9]+"),
                                prism_matches_u$cid)
prism_matches_u$ref_acc = str_extract(prism_matches_u$prism_filepath, pattern = "GCF_[0-9]+.[0-9]+")
prism_matches_u$strain = str_extract(prism_matches_u$SpecFile, "[A-Z]+[0-9]+-[0-9]+")

#check correct strain-reference matches
mash_out = read_csv("data/similar_genome/mash_out.csv") %>%select(genomes_path = query,ref_acc,ani,r_comment)
mash_out$org = str_extract(mash_out$r_comment,"[A-Z][a-z]+ [a-z]+")
asm_metadata = read_csv("data/metadata_assembled.csv") %>% select(genomes_path,strain = sample)
mash_out = merge(mash_out,asm_metadata,all.x = T,all.y = F)
  
  # # #plot ani
  # p_ani = ggplot(mash_out,
  #                aes(
  #                  x = fct_reorder(paste(strain, ref_acc, org),
  #                                  ani),
  #                  y = ani,
  #                  label = sprintf("%.1f",ani)
  #                )) +
  #   stat_identity(geom = "bar", color = "black", fill = "white") +
  #   geom_text(fontface ="bold", color = "black",
  #             #size = 8
  #             ) +
  #   scale_y_continuous(limits = c(80, 100), oob = scales::rescale_none) +
  #   labs(x=NULL, y= "Average Nucletide Identity") +
  #   theme_classic() +
  #   theme(#text = element_text(size = 18),
  #         axis.text = element_text(color="black")) +
  #   coord_flip()
  # # 
  # ggsave(p_ani,
  #        filename = "figures/similar_genomes_ani.png",
  #        width = 8, height = 7, dpi = 300)
  # ###


prism_matches_u = merge(prism_matches_u,mash_out,
      by.x = "strain",by.y = "strain",
      all.x = T,all.y = F,
      suffixes = c("spectra","mash"))

#consistent matches
prism_matches_u = prism_matches_u %>% filter(ref_accmash == ref_accspectra)

#inconsistent strain genome matches
#prism_matches_u = prism_matches_u %>% filter(ref_accmash != ref_accspectra)


#blast searches to confirm BGC presence
# blast_out = lapply(1:nrow(prism_matches_u),function(idx){
#   
#   q=sprintf("data/similar_genome/prism_extracted_bgcs/%s.fasta",prism_matches_u$cname[idx])
#   
#   system(sprintf("blastn -query %s -subject %s -outfmt \'6 sseqid length pident\'",
#          q,
#          prism_matches_u$genomes_path[idx]),
#          intern = T)
# })


#for B0602 molecule
require(MSnbase)
B0602 = prism_matches_u %>% filter(strain == "B06-02")
B0602Spectra = retrieveSpectra(B0602)
B0602Spectra$plot[[1]] = label_peaks(B0602Spectra$plot[[1]],by = "top_intensity")
B0602Spectra$eic = plot_chromatogram(retrieveEIC(B0602,idx =1))
B0602Spectra$eic = B0602Spectra$eic + 
  geom_vline(xintercept = B0602$Retention/60, linetype = "dotted") 
p_B0602 = cowplot::plot_grid(B0602Spectra$eic,B0602Spectra$plot[[1]],nrow = 2) 


B0602_predicted_fragments = read_csv("data/QTOF/moldiscovery_prism_similargenome/cfidspec_GCF_001509505.1_ASM150950v1_genomic.gbff.json_c6_s2.csv")

## npdtools2/bin/print_structure --print_mgf_spectrum --configs_dir data/QTOF/moldiscovery_prism_similargenome/work/configs/ --structure npdtools/prism_similargenome/mols/GCF_001509505.1_ASM150950v1_genomic.gbff.json_c6_s2.mol --fragmentation_tree --mrank 7 --charge 2 --fragment GCF_001509505.1_ASM150950v1_genomic.gbff.json_c6_s2.mol.cereal.mz.target
# B0602_predicted_fragments = read.delim(text = c("0 C4H5NO2 99.032
# 1 C3H5NO2 87.032
# 2 C5H9NO 99.0684
# 3 C9H9NO 147.068"),sep = ' ', header=F, col.names = c("frag_no","formula","mass"))
# 

#mz 98.9591 matches the closest
# B0602_predicted_fragments$mass_low = B0602_predicted_fragments$mass - 0.1
# B0602_predicted_fragments$mass_high = B0602_predicted_fragments$mass + 0.1
# sapply(1:nrow(B0602_predicted_fragments),function(idx)
# subset(mz(B0602Spectra$spectra[[1]]),mz(B0602Spectra$spectra[[1]]) > B0602_predicted_fragments$mass_low[idx] &
#          mz(B0602Spectra$spectra[[1]]) <B0602_predicted_fragments$mass_high[idx])
# )

predicted_spectra = new("Spectrum2",
                        mz=B0602_predicted_fragments$mz,
                        intensity = B0602_predicted_fragments$intensity)


compareSpectra(B0602Spectra$spectra[[1]],predicted_spectra, fun="common")

intersect(sprintf("%.1f",mz(B0602Spectra$spectra[[1]])),
sprintf("%.1f",B0602_predicted_fragments$mz)
)


#load prism data
prism_paths = Sys.glob("data/similar_genome/prism_out/*")

prism_json = lapply(prism_paths,function(prism_file)
  jsonlite::fromJSON(prism_file)$prism_results$clusters %>%
    select(where(is.atomic)) %>% 
    # mutate(prism_filepath = prism_matches_u$prism_filepath[idx]) %>% 
    rownames_to_column(.,"cid")
    )
prism_json = bind_rows(prism_json,.id = "prism_filepath")
prism_json = unique(prism_json)


#add the information
prism_matches_u = merge(prism_matches_u,prism_json,
      by.x = c("prism_filepath","cid"),by.y = c("prism_filepath","cid"),
      all.x = T,all.y = F) 

write_tsv(prism_matches_u,"data/QTOF/moldiscovery_prism_similargenome/significant_unique_matches_coords.tsv")

#subset sequences
library(biomaRt)
library(Biostrings)


cluster_sequeneces = lapply(1:nrow(prism_matches_u),function(idx){
g = read_genome(prism_matches_u$rgenome_filepath[idx])
g = subseq(g[str_detect(names(g), pattern = prism_matches_u$contig_name[idx])][[1]],
       start = (prism_matches_u$start[idx]+1), end = prism_matches_u$end[idx])
g = DNAStringSet(g)
names(g)= prism_matches_u$cname[idx]
writeXStringSet(g, 
                filepath = sprintf("data/similar_genome/prism_extracted_bgcs/%s.fasta",
                                   prism_matches_u$cname[idx])
                )
})



#explore prismdb
prismdb = read.csv("npdtools/prism/library_metadata.csv")
prismdb = separate(prismdb,
                   col = "name",
                   into = c("date","barcode","strain","cid","smile_n"),
                   sep = '_')
