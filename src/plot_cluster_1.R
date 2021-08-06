library(gggenes)
library(tidyverse)
library(stringr)
library(data.table)
library(RColorBrewer)
library(ggsci)
library(ggfittext)


#### prerequisite functions ####

# to get nested list from json file into a clean dataframe
# reference: https://stackoverflow.com/questions/35444968/read-json-file-into-a-data-frame-without-nested-lists
col_fixer <- function(x, vec2col = FALSE) {
  if (!is.list(x[[1]])) {
    if (isTRUE(vec2col)) {
      as.data.table(data.table::transpose(x))
    } else {
      vapply(x, toString, character(1L))
    }
  } else {
    temp <- rbindlist(x, use.names = TRUE, fill = TRUE, idcol = TRUE)
    temp[, .time := sequence(.N), by = .id]
    value_vars <- setdiff(names(temp), c(".id", ".time"))
    dcast(temp, .id ~ .time, value.var = value_vars)[, .id := NULL]
  }
}
Flattener <- function(indf, vec2col = FALSE) {
  require(data.table)
  require(jsonlite)
  indf <- flatten(indf)
  listcolumns <- sapply(indf, is.list)
  newcols <- do.call(cbind, lapply(indf[listcolumns], col_fixer, vec2col))
  indf[listcolumns] <- list(NULL)
  cbind(indf, newcols)
}


load_antismash <- function(antismash_json,cluster_num) {
  input_cluster = jsonlite::stream_in(file(antismash_json),,pagesize = 100000)
  input_cluster = input_cluster$records[[1]][cluster_num,] # get the specific cluster
  input_cluster$features[[1]]$qualifiers$note<-NULL # sometimes an additional note section may be present which do not work well.
  input_cluster = Flattener(input_cluster$features[[1]]) # convert nested list into a flattened dataframe where each row describes a feature (could be gene or a domain)
  input_cluster = input_cluster %>% separate(col = "location", sep = ":", into = c("start","end")) %>% 
    separate(col = "end", sep = "]", into = c("end","strand")) # extract start,end and strand from location column
  input_cluster$start = input_cluster$start %>% str_remove(.,pattern = "\\[") %>% as.numeric()
  input_cluster$end = as.numeric(input_cluster$end)
  input_cluster$direction = input_cluster$strand %>% str_replace_all(.,pattern = "\\(\\+\\)", "1") %>% str_replace_all(.,pattern = "\\(-\\)", "-1") %>% as.numeric()
  input_cluster$strand = input_cluster$strand %>% str_replace_all(.,pattern = "\\(\\+\\)", "forward") %>% str_replace_all(.,pattern = "\\(-\\)", "reverse")
  input_cluster = input_cluster %>% mutate(cluster =  basename(antismash_json))
  return(input_cluster)
}



get_nrps_polymer <- function(antismash_json, cluster_num) {
  input_cluster = jsonlite::stream_in(file(antismash_json), , pagesize = 100000)
  polymer = input_cluster$records[[1]][cluster_num,]$modules$antismash.modules.nrps_pks$region_predictions[[1]][[1]][1,]
  polymer = polymer %>% mutate(cluster = basename(antismash_json))
  return(polymer)
}
  

get_nrps_polymer_sample <- function(antismash_json, cluster_num) {
  input_cluster = jsonlite::stream_in(file(antismash_json), , pagesize = 100000)
  polymer = input_cluster$records[[1]]$modules$antismash.modules.nrps_pks$region_predictions$`2`[2][[1]][1,]
  polymer = polymer %>% mutate(cluster = basename(antismash_json))
  return(polymer)
}

cleanup_antismash <- function(clusters, multiple_clusters = T){
  if(multiple_clusters == T){clusters = do.call(bind_rows,clusters)}
  clusters = clusters %>% mutate_all(na_if,"")
  clusters$cluster = clusters$cluster %>% str_remove(pattern="\\d+_cluster_") %>% str_remove(".json")
  return(clusters)
}


cleanup_nrps_polymer <- function(polymer_df,modules_to_show){
  polymer_df = polymer_df %>% separate(polymer, into = modules_to_show, sep = "(\\) \\+ \\()| - ")
  polymer_df = lapply(polymer_df, gsub, pattern="\\(|\\)", replacement="") %>% bind_cols() 
  polymer_df  = polymer_df %>% pivot_longer(!sc_number&!docking_used&!smiles&!cluster,names_to = "module_no", values_to = "amino_acid")
  return(polymer_df)
}

trim_cluster = function(df = df, bgcs_to_trim = bgcs_to_trim){
  # bgc_to_trim = "WAC1416cluster068"
  for (bgc_to_trim in bgcs_to_trim) {
  marker = "hydrolase"
  marker_start = df %>% filter(cluster == bgc_to_trim & type == "CDS" & qualifiers.sec_met_domain_simplified == marker) %>% select("start") %>% .[1,]
  trim_end = marker_start + 50000
  trim_start = marker_start - 50000
  df = df[!(df$cluster == bgc_to_trim & (df$start <= trim_start | df$start >= trim_end)),]
  }
  
  return(df)
}


c25_colpal <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)



## genes ####



# known glycopeptide clusters
gpa_characterized = fread("data/antismash/gpa/characterized_clusters.csv") %>% filter(characterized == T)
gpa_clusters = list.files(path = "data/antismash/gpa/", pattern = "json", recursive = T, full.names = T)
gpa_clusters = lapply(gpa_clusters,function(x) load_antismash(x,1))
gpa_clusters = do.call(bind_rows,gpa_clusters)
#gp_clusters_select = 'Vancomycin|Teicoplanin|feglymycin|complestatin'
#gpa_clusters = gpa_clusters %>% filter(str_detect(cluster,gp_clusters_select))
gpa_clusters = gpa_clusters %>% mutate_all(na_if,"")
gpa_clusters$cluster = gpa_clusters$cluster %>% str_remove(pattern="\\d+_cluster_") %>% str_remove(".json")



# query cluster
# barcode08 = load_antismash("data/antismash/09182020/barcode08/medeka/consensus.json",2)
# barcode08 = barcode08 %>% mutate(cluster = "B06-02")
# barcode08 = barcode08 %>%  filter(start > 230000 & end < 300000) %>% mutate_all(na_if,"")
# gpa_clusters = bind_rows(barcode08,gpa_clusters)

#simplifying domain/gene annotation
gpa_clusters$qualifiers.sec_met_domain = str_remove_all(gpa_clusters$qualifiers.sec_met_domain, pattern = "\\((.*?)\\)") %>% trimws(which="both")

domain_labels =read.csv("data/antismash/gpa/sec_met_domains.csv")
domain_labels$Domain = domain_labels$Domain %>% trimws(which="both")
gpa_clusters = gpa_clusters %>%  mutate(qualifiers.sec_met_domain_simplified = case_when(
  (qualifiers.sec_met_domain %in% filter(domain_labels, label == "NRPS")$Domain ) ~ "NRPS",
  (qualifiers.sec_met_domain %in% filter(domain_labels, label == "epimerase")$Domain) ~ "epimerase",
  (qualifiers.sec_met_domain %in% filter(domain_labels, label == "glycotransferase")$Domain) ~ "glycotransferase",
  (qualifiers.sec_met_domain %in% filter(domain_labels, label == "p450")$Domain) ~ "p450",
  (qualifiers.sec_met_domain %in% filter(domain_labels, label == "polyketide synthase")$Domain) ~ "polyketide synthase",
  (qualifiers.sec_met_domain %in% filter(domain_labels, label == "sugar isomerase")$Domain) ~ "sugar isomerase",
  (qualifiers.sec_met_domain %in% filter(domain_labels, label == "aminotransferase")$Domain) ~ "aminotransferase",
  (qualifiers.sec_met_domain %in% filter(domain_labels, label == "Acyl-CoA dehydrogenase")$Domain) ~ "Acyl-CoA dehydrogenase",
  (qualifiers.sec_met_domain %in% filter(domain_labels, label == "3,5-dihydroxyphenylacetyl-CoA synthase")$Domain) ~ "3,5-dihydroxyphenylacetyl-CoA synthase",
  (qualifiers.sec_met_domain %in% filter(domain_labels, label == "hydrolase")$Domain) ~ "hydrolase",
  (qualifiers.sec_met_domain %in% filter(domain_labels, label == "halogenase")$Domain) ~ "halogenase"
#  TRUE ~ "other"
))


# trim large bgcs
# bgcs_to_trim = c("WAC1416cluster068","WAC1376cluster143","WAC4182cluster089")
# gpa_clusters = trim_cluster(df = gpa_clusters, bgcs_to_trim = bgcs_to_trim)


# ploting genes and coloring by domain/function
p_gpa_domains = ggplot(
  gpa_clusters %>% filter(type == 'CDS'),
  aes(
    xmin = start,
    xmax = end,
    y = cluster,
    forward = direction,
    fill = str_remove_all(qualifiers.sec_met_domain_simplified, pattern = "\\((.*?)\\)"),
  )
) +
  geom_gene_arrow(arrowhead_height = unit(9, "mm"), arrowhead_width = unit(3, "mm"),arrow_body_height = unit(9, "mm"),size =1.1) + 
  theme_genes() +
  ggtitle("BGC\nstructure") +
  #scale_fill_brewer(palette = "Set3") +
  scale_fill_manual(values = colorRampPalette(brewer.pal(11, "Paired"))(11), name = "Gene function") +
  theme(
    legend.position = "bottom",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_blank(),
    text = element_text(size=42,face = "bold"),
    panel.spacing.y = unit(0.5, "lines"),
    plot.title = element_text(hjust = 0.5, face = "bold")
    #axis.text.y = element_blank(),
  ) +
  scale_x_continuous(expand = c(0,0))+
  ylab(NULL) + guides(fill = guide_legend(nrow = 2, title.position = "top", title.hjust = 0.5)) 

##ggsave("figures/gpa_domains.png", p_gpa_domains, height = 25, width = 45)

## Showing domains
# ggplot(gpa_clusters %>% filter(type == 'CDS' ), aes(xmin = start,xmax = end, y = cluster,forward = direction)) +
#   facet_wrap(~ cluster, scales = "free", ncol = 1) +
#   geom_gene_arrow(fill  = "white") +
#   geom_subgene_arrow(data = gpa_clusters %>% filter(type == 'aSDomain' ),
#                      aes(xsubmin = start, xsubmax = end, y = cluster, fill = qualifiers.aSDomain ) , color = "black" )+
#   theme_genes() +
#   scale_fill_brewer(palette = display.brewer.pal(20, "Set3"))







### polymer ####

gpa_clusters_files = list.files(
  path = "data/antismash/gpa/",
  pattern = "json",
  recursive = T,
  full.names = T
)

gpa_nrps_polymers = lapply(gpa_clusters_files, function(x)
                      get_nrps_polymer(x, 1))
gpa_nrps_polymers = cleanup_antismash(gpa_nrps_polymers)
modules_to_show = sprintf("%d", 01:10)
gpa_nrps_polymers = gpa_nrps_polymers %>% separate(polymer, into = modules_to_show, sep = "(\\) \\+ \\()| - ")
gpa_nrps_polymers = lapply(gpa_nrps_polymers, gsub, pattern="\\(|\\)", replacement="") %>% bind_cols() 

#consistent coloring with filtering or subsetting
gpa_nrps_polymers = gpa_nrps_polymers %>% pivot_longer(!sc_number&!docking_used&!smiles&!cluster,names_to = "module_no", values_to = "amino_acid")
d_gpa_nrps_polymers = gpa_nrps_polymers %>% mutate(amino_acid_col = str_remove(pattern = 'D-', amino_acid))

amino_acid_col_pal = colorRampPalette(brewer.pal(9, "Set1"),bias = 0.5)(length(unique( na.omit(d_gpa_nrps_polymers$amino_acid_col) )))
#amino_acid_col_pal = colorRampPalette(c25_colpal)(length(unique( na.omit(d_gpa_nrps_polymers$amino_acid_col) )))
names(amino_acid_col_pal) = unique(na.omit(d_gpa_nrps_polymers$amino_acid_col))
amino_acid_col_pal[['X']] = "#D3D3D3"
amino_acid_col_pal[['tyr']] = "#D2F8D2"
amino_acid_col_pal[['asn']] = "#9FD7FB"


d_gpa_nrps_polymers$module_no = as.numeric(d_gpa_nrps_polymers$module_no)

p_gpa_nrps_polymers = ggplot(data = d_gpa_nrps_polymers, 
       aes(x=module_no,
           y=cluster,
           label=amino_acid
           )) +
  geom_point(aes(fill=amino_acid_col ),size = 16,shape=21, color = "black", stroke = 1.2) +
  scale_x_continuous(breaks = 1:10) +
  scale_fill_manual(values = amino_acid_col_pal, na.value = "white") +
  theme_minimal() +
  ggtitle("A-domain\nspecificity") +
  scale_y_discrete(position = "right") +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),axis.ticks.y = element_blank(),
        legend.position = "bottom", 
        panel.grid = element_blank(), panel.spacing.x = unit(2, "lines"),
        text = element_text(size=42,face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  xlab(NULL) +  guides(fill = guide_legend(nrow = 2, title.position = "top", title.hjust = 0.5,
                                           title = "Amino acid")) 


##ggsave("figures/gpa_nrps_polymers_trp.png",p_gpa_nrps_polymers, height = 30, width = 15)
