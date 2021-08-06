library(tidyverse)
library(ggtree)
library(aplot)
library(cowplot)
# https://yulab-smu.top/treedata-book/chapter13.html?q=align#genome-locus


list2dist = function(dat){
  #from https://rdrr.io/cran/spaa/src/R/list2dist.R
  dat.name1 <- as.character(dat[,1])
  dat.name2 <- as.character(dat[,2])
  dat.value <- dat[,3]
  names1 <- sort(unique(as.character(dat[,1])))
  names2 <- sort(unique(as.character(dat[,2])))
  total.names <- unique(c(names1, names2))
  elements <- rep(NA, length(total.names)^2)
  dim(elements) <- c(length(total.names),length(total.names))
  rownames(elements) <- total.names
  colnames(elements) <- total.names
  
  for(i in 1:length(total.names)){
    for(j in 1:length(total.names)){
      for(k in 1:length(dat.name1)){
        if((total.names[i] == dat.name1[k])&(total.names[j] == dat.name2[k])){
          elements[i,j] <- dat.value[k]
        }
      }
    }
  }
  res <- as.dist(t(elements))
  return(res)
}



# construct tree
distances = read.csv("/gpfs/gsfs11/users/rajwanir2/soil_metagenomics/data/bigscape/gpa/network_files/2020-12-23_14-03-25_hybrids_glocal/NRPS/NRPS_c1.00.network",
         sep = '\t')
distances = distances %>% select(Clustername.1,Clustername.2,raw.DSS.anchor)
distances = list2dist(distances)
tree = ape::nj(distances)
tree$tip.label = tree$tip.label %>%  str_remove(pattern="\\d+_cluster_",.)
p_tree = ggtree(tree)


this_study = list.files(path = "data/antismash/gpa/", pattern = "json", recursive = T, full.names = T) %>% grep(pattern = "00_cluster_(.*)\\/", value = T) %>% str_extract(pattern = "00_cluster_(.*)\\/") %>% str_remove_all(pattern = "00_cluster_|\\/")
tip_label_color = data.frame(taxa = tree$tip.label,
                             tipcolor = if_else(tree$tip.label %in% this_study,"red","black"))  

p_tree = p_tree %<+% tip_label_color + 
  geom_tiplab(align = T, aes(color = tipcolor), size = 10) + 
  scale_color_manual(values = c(black = "black", red ="red")) + 
  theme(legend.position = "none") + scale_x_continuous(expand = c(0,5))

# Get legends separate
p_gpa_nrps_polymers_leg = cowplot::get_legend(p_gpa_nrps_polymers)
p_gpa_domains_leg = cowplot::get_legend(p_gpa_domains)
p_gpa_nrps_polymers = p_gpa_nrps_polymers + theme(legend.position = "none")
p_gpa_domains = p_gpa_domains + theme(legend.position = "none")

# combine tree with additional plots
p_tree_plus = p_gpa_domains %>%
  insert_left(p_gpa_nrps_polymers,width = 0.2) %>%
  insert_left(p_tree,width = 0.42) 

# save figure without legend then add legend and resave
ggsave(p_tree_plus,file="figures/gpa_tree.png",width = 48,height=30)
p_tree_plus = plot_grid(ggdraw() + draw_image("figures/gpa_tree.png"),
                     p_gpa_nrps_polymers_leg,
                     p_gpa_domains_leg,
                     ncol = 1,nrow = 3,
                     rel_heights = c(0.9,0.1,0.1)
                   )

ggsave(p_tree_plus,file="figures/gpa_tree.png",width = 48,height=30)
