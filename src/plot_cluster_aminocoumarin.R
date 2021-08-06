library(tidyverse)
library(gggenes)
library(data.table)
library(RColorBrewer)

target_cluster = "aminocoumarin"

clusters = list.files(path = sprintf("data/antismash/%s//",target_cluster),
                       pattern = ".gff.1", recursive = T, full.names = T)
clusters = lapply(clusters,
            rtracklayer::readGFF)
clusters = bind_rows(clusters, .id = "id")

clusters = clusters %>% mutate(seqid = case_when(seqid == "c00001_tig0000.." ~ "B06-3",
                                     TRUE ~ as.character(seqid)))

clusters = clusters %>% mutate(molecule = case_when(seqid == "BGC0001287" ~ "chaxamycin",
                                                    seqid == "BGC0000833" ~ "coumermycinA1",
                                                    seqid == "BGC0000832" ~ "clorobiocin",
                                                    seqid == "BGC0000141" ~ "rubradirin",
                                                    seqid == "AF170880" ~ "novobiocin",
                                                    seqid == "B06-3" ~ "B06-3"))


clusters = clusters %>% filter(!(seqid == "B06-3" & start > 100244 ))
clusters$gene_functions = clusters$gene_functions %>% str_extract("(.*):(.*)\\(") %>% str_remove_all(pattern = "(.*)smcogs\\) |\\(")
#get high frequency gene function 
freq_g_func = clusters %>% group_by(gene_functions) %>% count(gene_functions) %>% filter(n>3)
clusters = clusters %>% mutate(g_func_simple = case_when(
                        gene_functions %in% freq_g_func$gene_functions ~ as.character(gene_functions)))
                                                
g_func_colors = colorRampPalette(brewer.pal(12, "Set3"))(nrow(freq_g_func))

# clusters$g_func_simple = clusters$gene_functions %>% str_extract(pattern = '\\)(.*)') %>% str_remove_all(pattern = 'SMC\\w+:|\\) |\\((.*)')
# clusters = clusters %>% mutate(g_func_simple = case_when(
#   g_func_simple %like% "regulator" ~ "regulator",
#   g_func_simple %like% "other" ~ "other",
#   g_func_simple %like% "transporter" ~ "transporter",
#   g_func_simple %like% "dehydrogenase" ~ "dehydrogenase",
#   g_func_simple %like% "p450|P450" ~ "p450",
#   g_func_simple %like% "tra_KS|hyb_KS|PKS_AT" ~ "PKS",
#   g_func_simple %like% "halogenase" ~ "halogenase",
#   g_func_simple %like% "methyltransferase" ~ "methyltransferase",
#   g_func_simple %like% "DegT/DnrJ/EryC1/StrS aminotransferase" ~ "aminotransferase",
#   g_func_simple %like% "malonyl CoA-acyl carrier protein transacylase" ~ "malonyl CoA-ACP transacylase",
#   g_func_simple %like% "Beta-ketoacyl synthase" ~ "Beta-ketoacyl synthase",
#   TRUE ~ "other"
#   ))


# clusters = clusters %>% mutate(g_func_simple = na_if(g_func_simple,"other"))
  

p = ggplot(clusters %>% filter(type...3 == "CDS")) +
  geom_gene_arrow(aes(xmin=start, xmax = end, y = seqid,
                      fill = g_func_simple), size = 0.6) + theme_genes() +
  theme(legend.position = "bottom") +
  # scale_x_discrete(expand = c(0,0),
  #                  breaks = scales::trans_breaks(identity, identity, n = 5)) +
  scale_fill_discrete(type = g_func_colors,na.value = "dimgrey") + 
  guides(fill = guide_legend(ncol=3, title = "Gene function", 
                             title.position = "top", title.hjust = 0.5, 
                             title.theme = element_text(face = "bold")))



##add AHBA synthesis genes
aminocoumarin_regions = clusters %>% filter(product == "aminocoumarin" & type...3 == "region" & seqid == "B06-3") %>% select(seqid,start,end,id,product)
aminocoumarin_regions$molecule = factor("B06-3",levels = sort(unique(clusters$molecule)))
aminocoumarin_regions$id = as.numeric(aminocoumarin_regions$id)
p = p + geom_segment(data = aminocoumarin_regions, aes(y=id+0.2, x = start, xend = end, yend = id+0.2),color = "darkgreen", size = 0.7) +
  geom_text(data = aminocoumarin_regions, aes(y=id+0.4, x = start+end/2, label = "AHBA biosynthesis genes"), color = "darkgreen")


#AHBA domain
CAL_AHBA = clusters %>% filter(specificity == "Minowa=AHBA") %>% select(seqid,start,end,locus_tag,aSDomain,specificity,id)
CAL_AHBA_cds = clusters %>% filter((locus_tag %in% as.character(CAL_AHBA$locus_tag) | (gene %in% as.character(CAL_AHBA$locus_tag))) & type...3 == "CDS") %>% select(seqid,gstart=start,gend=end,id)
CAL_AHBA = merge(CAL_AHBA,CAL_AHBA_cds) %>% unique()
CAL_AHBA$id = as.numeric(CAL_AHBA$id)
CAL_AHBA = CAL_AHBA %>% mutate(aSDomain = "CAL domain\nAHBA specific")

p = p + geom_subgene_arrow(data = CAL_AHBA, aes(xmin = gstart, xmax = gend, 
                                            xsubmin = start, xsubmax = end,
                                            y = id),fill = "darkgreen") +
  ggrepel::geom_text_repel(data = CAL_AHBA, aes(label = " ", x = (start+end)/2, y =id +0.2), color = "darkgreen", fontface = "bold", size = 3.5, min.segment.length = 0) +
  geom_text(data = CAL_AHBA, aes(label = aSDomain, x = (start+end)/2, y =id +0.5), color = "darkgreen") 


p = p + scale_y_discrete(labels = unique(clusters$molecule), expand = c(0.1,0.3))+
        theme(axis.text.y = element_text(size = 14,face="bold", color = "black"),
              panel.grid.major.y = element_line(color = "black"),
              legend.background = element_rect(color = "black")) +
    labs(x="Scale (bp)", y = NULL) 
    

ggsave(p, filename = "figures/aminocoumarin_cluster.png",
       width = 15, height = 6, dpi = 300)