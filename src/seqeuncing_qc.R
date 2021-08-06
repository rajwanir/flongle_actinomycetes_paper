library(tidyverse)
library(cowplot)
library("ggsci")
library(data.table)
library(RColorBrewer)

sample_metadata = read.csv("data/metadata.csv")
flowcell_colors = colorRampPalette(brewer.pal(9,"Set1"))(length(unique(sample_metadata$flow_cell_id)))
names(flowcell_colors) = unique(sample_metadata$flow_cell_id)

unsucessful_runs = unique(sample_metadata[sample_metadata$run_fail == T,] %>% select(flow_cell_id,seq_date))

# throughput ####
throughput = Sys.glob(file.path(getwd(), "data/rawdata/vanRsoilactino_*/*/thro*.csv"))
throughput = lapply(throughput, function(x)
  read.csv(x) %>% mutate(flow_cell_id = x)) %>%  bind_rows()
throughput = throughput %>% mutate(
  flow_cell_id = str_extract(flow_cell_id, pattern = "throughput_\\w+_") %>%  str_remove(pattern =
                                                                                           "throughput_") %>% str_remove("_"),
  experiment_time_hr = Experiment.Time..minutes. / 60,
  estimated_bases_mb = Estimated.Bases / 1000000
)



p_throughput = ggplot(data = throughput,
                      mapping = aes(x = experiment_time_hr,
                                    y = estimated_bases_mb)) +
  geom_line(aes(color = flow_cell_id), size = 1.4) +
  theme_bw() +
  labs(x = "Time (hours)", y = "Estimated bases (Mb)", title = "Cumulative output over time") +
  scale_x_continuous(limits = c(0,25), expand = c(0, 0)) + scale_color_manual(values = flowcell_colors) +
  theme(legend.position = "none", #legend.background = element_rect(color = "black"),
        text = element_text(size=16),plot.title = element_text(hjust = 0.5,face="bold"),
        panel.grid = element_blank()) 
  #guides(color = guide_legend(nrow = 2, title = "Flow cell"))

# duty_time ####

duty_time = Sys.glob(file.path(getwd(), "data/rawdata/vanRsoilactino_*/*/duty_time*.csv"))
duty_time = lapply(duty_time, function(x)
  read.csv(x) %>% mutate(flow_cell_id = x)) %>%  bind_rows()

duty_time = duty_time %>% mutate(
  flow_cell_id = str_extract(flow_cell_id, pattern = "duty_time_\\w+_") %>%  str_remove(pattern = "duty_time_") %>% str_remove("_"),
  experiment_time_hr = Experiment.Time..minutes. / 60,
  channel_status = case_when(
    Channel.State == "adapter" | Channel.State == "strand" ~ "sequencing",
    Channel.State == "pore" ~ "pore",
    Channel.State == "unavailable" ~ "recovering",
    Channel.State == "saturated" | Channel.State == "multiple" | Channel.State == "zero" ~ "inactive",
    Channel.State == "unclassified" | Channel.State == "unclassified_following_reset" |Channel.State == "pending_manual_reset" | Channel.State == "pending_mux_change"  ~ "unclassified",
    Channel.State == "no_pore" ~ "no_pore",
    TRUE ~ "unknown"
  )
)


p_duty_time = ggplot(data = duty_time,
       mapping = aes(x = experiment_time_hr,
                     y=State.Time..samples.,
                     fill=channel_status)) +
  geom_bar(position = "stack", stat="identity", alpha = 0.65) +
  facet_wrap(~flow_cell_id,ncol = 1, strip.position = "right") +
  theme_bw() +
  labs(x= "Time (hours)", y = "Proportion of pores", title = "Pore activity over time") +
  theme(axis.text.y = element_blank(),
        legend.position = "bottom", legend.background = element_rect(color = "black"),
        text = element_text(size=16),plot.title = element_text(hjust = 0.5,face="bold"),
        panel.grid = element_blank(), strip.background = element_blank(),strip.text = element_text(size=7)) +
  scale_x_continuous(limits = c(0,17),expand = c(0, 0)) +
  scale_fill_rickandmorty() +
  guides(fill = guide_legend(nrow = 1, title = "Channel status", title.position = "top", title.hjust = 0.5))


facet_strips = ggplot(data = duty_time,
                      mapping = aes(x = experiment_time_hr,
                                    y=State.Time..samples.,
                                    fill=channel_status))  +
                    facet_wrap(~flow_cell_id,ncol = 1, strip.position = "right")  +
                    geom_rect(aes(fill = flow_cell_id),  xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
                    theme(strip.background = element_blank(), strip.text = element_text(size = 8)) + scale_fill_manual(values = flowcell_colors)


source("src/func_change_facet_strip.R")
p_duty_time = change_strip_color(p_duty_time,facet_strips)


# purification ####
#sample_metadata = sample_metadata[!duplicated(sample_metadata$sample), ]
sample_metadata = sample_metadata %>% mutate(recovery = (bc_conc*13)/(dna_conc*ep_dna_vol)*100)
ggplot(data = sample_metadata) +
  geom_histogram(
    aes(x = recovery, fill = bc_size_selection),
    binwidth = 10,
    alpha = 0.5,
    color = "black"
  ) +
  scale_x_continuous(breaks = seq(10, 150, 10)) +
  labs(x = "% DNA recovered after \n0.5x standard + (0.5x standard/0.10-0.15 modified)", y = "frequency", title = "DNA recovery") +
  theme_bw() +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank()) 

# dna_purity ####
nanodrop = read.csv("data/nanodrop_1.tsv",sep = '\t')
sample_metadata = merge(sample_metadata,nanodrop,by.x = "sample", by.y = "Sample.ID",all.x = T,all.y = F)
p_dna_purity = ggplot(data = sample_metadata,
                      aes(x = flow_cell_id, y = X260.280, fill = flow_cell_id)) +
  geom_boxplot() + 
  # geom_point(size = 5,
  #            shape = 21,
  #            position = position_jitterdodge()) +
  scale_fill_manual(values = flowcell_colors) +
  theme_bw() + labs(y = "OD 260/280", x = "Flow cell id", title = "Purity of starting DNA") +
  theme(legend.position = "none",
        text = element_text(size=16), plot.title = element_text(hjust = 0.5,face="bold"),
        panel.grid = element_blank()) +
  coord_flip()




# read length N50 ####

read_length = Sys.glob(file.path(getwd(), "data/rawdata/vanRsoilactino_*/basecalled_fastq/sequencing_summary.txt"))
read_length = lapply(read_length, function(x) {
  s_red_len = read.csv(x, sep = '\t') %>% mutate(flow_cell_id = str_extract(filename, pattern = "\\w+_") %>% str_remove("_\\w+"))
  s_red_len = s_red_len[order(s_red_len$sequence_length_template), ]
  s_red_len = s_red_len %>% mutate(cumulative_bases = cumsum(sequence_length_template))
  s_red_len = s_red_len %>% mutate(
    sequence_length_template_kb = sequence_length_template / 1000,
    pct_bases = sequence_length_template_kb / sum(sequence_length_template_kb)
  )
  return(s_red_len)
}) %>%  bind_rows()



read_n50 = read_length %>% group_by(flow_cell_id) %>% filter(cumulative_bases/sum(sequence_length_template) >= 0.5) %>% group_by(flow_cell_id) %>%  slice(1) %>% select(flow_cell_id,N50 = sequence_length_template_kb)

p_read_length = ggplot(data = read_length, aes(x = sequence_length_template_kb)) +
  geom_histogram(aes(fill = flow_cell_id, weight = pct_bases),
                 binwidth = 3,
                 color = "black"
                 ) +
  scale_x_continuous(
    limits = c(0, 80),
    expand = c(0, 0),
    breaks = seq(0, 80, by = 20)
  ) +
  scale_y_continuous(
    limits = c(0, 0.3),
    expand = c(0, 0)
  ) +
  geom_vline(
    data = read_n50,
    aes(xintercept = N50),
    linetype = "dotted",
    size = 1,
    color = "brown"
  ) +
  facet_wrap(
    ~ flow_cell_id,
    nrow = 2,
    scales = "fixed",
    strip.position = "top"
  ) +
  scale_fill_manual(values = flowcell_colors) +
  geom_text(data = read_n50, aes(y = 0.2, x = N50+30, label = paste("N50",N50,sep="\n")), size =
              6) +
  labs(x = "Read length (Kb)", y = "Proportion of reads", title = "Read length distribution") +
  theme_bw() + theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),strip.text = element_text(size=14),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5,face="bold"),
    text = element_text(size=16)
  )


# demultiplexed ####

demultiplexed = Sys.glob("data/rawdata/vanRsoilactino_*/basecalled_fastq/demultiplexed.tsv")
demultiplexed = lapply(demultiplexed, function(x) {
  read.csv(x, sep = '\t') %>% mutate(seq_date = as.numeric(str_extract(x, pattern = "\\d{8}")) )
}) %>% bind_rows()

#printing demultiplexing summary table
# demultiplexed = merge(demultiplexed,sample_metadata,by=c("seq_date","barcode"))
# demultiplexed_summarytable = demultiplexed %>% group_by(flow_cell_id,barcode) %>% 
#   summarize(sample = unique(sample),
#             total_bases = sum(length)/1000000, 
#             readlength_max = max(length)/1000, 
#             readlength_n50 = which(cumsum(sort(length))>(sum(length)/2))[1]/1000,
#             bc_size_selection = unique(bc_size_selection),
#             time_of_selection = unique(time_of_selection),
#             run_fail = unique(run_fail)
#             ) 
# 
# write_csv(demultiplexed_summarytable, file = "tables/demultiplexed_summary.csv")


demultiplexed = demultiplexed[!demultiplexed$seq_date %in% unsucessful_runs$seq_date,]
demultiplexed_summary = demultiplexed %>% group_by(seq_date,barcode) %>% summarize(total_reads = n(), total_bases = sum(length))
demultiplexed_seqindx = paste(demultiplexed_summary$seq_date,demultiplexed_summary$barcode) %in% paste(sample_metadata$seq_date,sample_metadata$barcode) 
demultiplexed_summary = demultiplexed_summary[demultiplexed_seqindx | demultiplexed_summary$barcode == "none",]
demultiplexed_summary = demultiplexed_summary %>% mutate(none = ifelse(barcode == "none",T,F))


ggplot(demultiplexed_summary, aes(y=total_bases,x=barcode, fill = seq_date)) + geom_bar(stat = "identity") + facet_wrap(~seq_date, scales = "free_x")

# cumulative ####
p_cumulative =  plot_grid(plot_grid(p_throughput, p_dna_purity, rel_heights = c(0.5, 0.5), ncol = 1, labels = c("A","B")),
                          plot_grid(p_read_length,p_duty_time, ncol = 1, rel_heights = c(0.35, 0.65), labels = c("C","D")),
                ncol = 2, rel_widths = c(0.3,0.7))
    

ggsave(p_cumulative,filename = "figures/cumulative_output.png",
       height = 12, width = 20, dpi = 300)