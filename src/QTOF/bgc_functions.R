


download_mibig = function(mibig_ids = NULL,dest_dir = NULL){
  # given a list of mibig ids downloads gbk files
  mibig_web = sprintf("https://mibig.secondarymetabolites.org/repository/%s/%s.1.region001.gbk",
          mibig_ids,mibig_ids)
  
  exit_stat = sapply(1:length(mibig_ids),function(idx){
    tryCatch(download.file(mibig_web[idx],
                           destfile = sprintf("%s/mibig/gbk/%s.gbk",
                                              dest_dir,mibig_ids[idx])),
             error=function(e) 1)
  })
  
  return(exit_stat)
}


gbk2gff = function(file_pat = NULL){
  
  gbk_files = Sys.glob(file_pat)
  
  sapply(gbk_files, function(gbkfile) system(paste0("bash src/gb2gff.sh ",
                                  gbkfile)))
  
  
}


extract_proteins_from_gff = function(gff_path){
  #read gff
  gff = rtracklayer::readGFF(gff_path,tags=c("translation","locus_tag","type","ID"))
  #filter CDS
  gff = gff[gff[,3]=="CDS",]
  
  #get Biostring object
  protSeq = Biostrings::AAStringSet(x = gff$translation,use.names=TRUE)
  names(protSeq) = as.character(gff$ID)
  #write
  exit_stat = tryCatch(Biostrings::writeXStringSet(x=protSeq,
                  append=F,
                  filepath=str_replace(gff_path,
                                       pattern = ".gff",
                                       replacement = "_protein.fasta")),
           error=function(e) 1)
  
  
  return(exit_stat)
}

simplify_genefunction = function(gene_funtion_col, type = "cluster"){
  
  if (type == "smcogs"){
  simplified_gene_function = str_extract(gene_funtion_col,"(.*):(.*)\\(") %>% 
    str_remove_all(pattern = "(.*)smcogs\\) |\\(")
  }
  else if (type == "cluster"){
    simplified_gene_function = str_extract(gene_funtion_col,"\\w+=[a-z]+")  
  }
  
  return(simplified_gene_function)
}

extract_proteins_prism = function(){
  
  prism_matches_gff = lapply(1:length(prism_matches_gff),function(idx) 
    mutate(prism_matches_gff[[idx]],bgc_id = prism_matches$Name[idx])
  )
  
  
  sapply(prism_matches_gff, function(gff) {
    #get Biostring object
    protSeq = Biostrings::AAStringSet(x = gff$translation, use.names = TRUE)
    names(protSeq) = as.character(gff$locus_tag)
    #write
    Biostrings::writeXStringSet(
      x = protSeq,
      append = F,
      filepath = sprintf(
        "data/QTOF/moldiscovery_prism/extracted_seqs/%s_protein.fasta",
        gff$bgc_id[1]
      ))
  })
}

blast_prism = function(){
  lapply(1:nrow(prism_matches),function(idx){
    blast_cmd = sprintf("blastp -query %s -subject %s -outfmt \'6 delim=@ qseqid sseqid pident length sstart\' -max_target_seqs 1",
                        sprintf("data/QTOF/moldiscovery_prism/extracted_seqs/%s_protein.fasta",prism_matches$Name[idx]), 
                        Sys.glob(sprintf("data/bgcs_gbk/*%s_protein.fasta.1", 
                                         prism_matches$m_strain[idx])))
    
    blast_out = read.delim(text = system(blast_cmd, intern = T), header = F,
                           col.names = c("qseqid", "sseqid", "pident", "length", "sstart"),
                           sep = '@')
    
    return(blast_out)
  })
  
}


retrieve_prism_genes = function(cluster_num = 3,
                                json_path = "data/prism/12112020_barcode09_GA6-002.gbk.json") {
  prism_json = jsonlite::fromJSON(json_path,
                                  flatten = T,
                                  simplifyDataFrame = T)
  prism_json = prism_json$prism_results$clusters$orfs[[cluster_num]] %>%
    select(start, stop, name, domains, frame, type)
  prism_json$n_domains = sapply(prism_json$domains, length)
  names(prism_json$domains) = prism_json$name
  domains = bind_rows(prism_json$domains, .id = "locus_tag")
  prism_json = merge(
    prism_json,
    domains,
    by.x = "name",
    by.y = "locus_tag",
    suffixes = c("_gene", "_domain")
  )
  return(prism_json)
}



read_antismash_clusterblast = function(result) {
  blast_col_fmt8 = c(
    "query",
    "subject",
    "p_id",
    "alignment_length",
    "mismatches",
    "gap_openings",
    "query_start",
    "query_end",
    "subject_start",
    "subject_end",
    "Evalue",
    "bitscore"
  )
  antismash_blast_seq_cols = c("input",
                               "bgc_num",
                               "cds_range",
                               "cds_strand",
                               "cds_locus_id",
                               "gene",
                               "add_id")
  
  result = read.table(text = result,
                      sep = '\t',
                      col.names = blast_col_fmt8)
  
  result = result %>% separate(col = "query",
                               sep = "\\|",
                               into = c(paste(
                                 "sample", antismash_blast_seq_cols, sep = '_'
                               )))
  
  
  result = result %>% separate(col = "subject",
                               sep = "\\|",
                               into = c(paste(
                                 "mibig", antismash_blast_seq_cols, sep = '_'
                               )))
 
  result = mutate(result,
                  contig = as.numeric(unique(str_extract(sample_cds_locus_id,"\\d+")) )
                    )
  
  
  result = result %>% separate(col = "sample_cds_range",
                               sep = "-",
                               into = c("sample_cds_start","sample_cds_end"))
  
  result = result %>% separate(col = "mibig_cds_range",
                               sep = "-",
                               into = c("mibig_cds_start","mibig_cds_end"))
  
  
  result = mutate_at(result,vars(contains(c("_start","_end"))),
                                 as.numeric)
  
  result = mutate(result, 
                  sample_cds_length = abs(sample_cds_end-sample_cds_start),
                  mibig_cds_length = abs(mibig_cds_end-mibig_cds_start))
  
  
  return(result) 
}


read_gff_clusterblast = function(gff,idx = idx, type = "sample"){
  
  
  all_start = antismash_clusterblast[[idx]][paste0(type,'_cds_start')] %>% unlist() %>% as.numeric() %>% min()
  all_end = antismash_clusterblast[[idx]][paste0(type,'_cds_end')] %>% unlist() %>% as.numeric() %>% max()
  contig = system(sprintf("grep %s %s |cut -f1 | uniq",
                          paste0("ctg",unique(antismash_clusterblast[[idx]]$contig)),
                          gff),
                 intern = T)
  
  
  gff = rtracklayer::readGFF(gff,
                             filter = list(
                               type = c("CDS"),
                               start = as.character(all_start:all_end),
                               seqid = contig
                             ))
  
  
  
  gff = mutate(gff,
               start = abs((start - all_start)),
               end = abs((end - all_start)))
  
  
  return(gff)
}

keep_best_clusterblast = function(clusterblast_df){
  bestcluster= clusterblast_df %>% 
    group_by(sample_bgc_num) %>% 
    summarize(sum(bitscore),
              n_cds = length(unique(mibig_cds_locus_id)),
              pid = mean(p_id),
              coverage = mean(abs(as.numeric(mibig_cds_end)-as.numeric(mibig_cds_start))/as.numeric(alignment_length)),
              eval = mean(Evalue)
    ) %>%
    arrange(-n_cds,-coverage,-pid) %>% 
    .[1,"sample_bgc_num"] %>% 
    as.character()
  
  clusterblast_df = filter(clusterblast_df,
                           sample_bgc_num == bestcluster)
  
  return(clusterblast_df)
}