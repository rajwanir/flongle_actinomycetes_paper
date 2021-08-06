library(gggenes)
library(tidyverse)
library(stringr)
library(data.table)

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








## main ####

barcode08 = jsonlite::stream_in(file("data/antismash/barcode08/barcode08.contigs.json"),,pagesize = 100000)
barcode08 = barcode08$records[[1]][2,] # get the specific cluster
barcode08 = Flattener(barcode08$features[[1]]) # convert nested list into a flattened dataframe where each row describes a feature (could be gene or a domain)
barcode08 = barcode08 %>% separate(col = "location", sep = ":", into = c("start","end")) %>% 
              separate(col = "end", sep = "]", into = c("end","strand")) # extract start,end and strand from location column
barcode08$start = barcode08$start %>% str_remove(.,pattern = "\\[") %>% as.numeric()
barcode08$end = as.numeric(barcode08$end)
barcode08$direction = barcode08$strand %>% str_replace_all(.,pattern = "\\(\\+\\)", "1") %>% str_replace_all(.,pattern = "\\(-\\)", "-1") %>% as.numeric()
barcode08$strand = barcode08$strand %>% str_replace_all(.,pattern = "\\(\\+\\)", "forward") %>% str_replace_all(.,pattern = "\\(-\\)", "reverse")


# ####
# 
# 
features_to_plot = barcode08 %>%  filter(start > 235000   & end < 275877 & type == "CDS")


ggplot(features_to_plot, aes(xmin = start,xmax = end, y = "strand",forward = direction,
                             color = str_extract(features_to_plot$qualifiers.sec_met_domain, pattern = "\\w+"))) +
  geom_gene_arrow() +
#  geom_gene_label(align = "right") +
  theme_genes() +
  scale_fill_brewer(palette = "Set3")





