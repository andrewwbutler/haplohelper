library(cowplot)
library(dplyr)
library(stringr)

segment = "HA"
colorPalette <- c("A" = '#e41a1c',"C" = '#377eb8',"G" = '#4daf4a',"T" = '#984ea3'," " = '#ff7f00')


haplotype_file = paste("~/Projects/flu_project/pacbio_haplotypes/", segment, "_illumina_linkage_info.csv", sep="")
haplotypes = read.csv(haplotype_file)

samples_compared = unique(haplotypes$sample)

for (s in samples_compared) {
  haplotypes %>% filter(sample == s) -> sample_data 
  days = unique(sample_data$day)
  for (d in days){
    sample_data %>% filter(day == d) -> sample_data
    hap_length=nchar(toString(haplotypes$pb_haplotype[1]))
    sample_data$number.of.n <- -1 * str_count(sample_data$il_haplotype, "N") + -1000 * sample_data$pb_count
  
    hap_num = 0
    for (pb_hap in unique(sample_data$pb_haplotype)){
      hap_num = hap_num + 1
      sample_data %>% filter(pb_haplotype == pb_hap) -> sample_data_filtered 
      sample_data_filtered %>% mutate(hap_count = as.factor(paste(il_haplotype, "---", il_count))) -> sample_data_filtered
      sample_data_filtered$hap_count_ordered <- reorder(sample_data_filtered$hap_count, sample_data_filtered$number.of.n)
      
      p1 <-ggplot(sample_data_filtered, aes(x = as.factor(variant_position), y = hap_count_ordered, group=il_haplotype)) + geom_line(linetype=2) +geom_point(aes(color=nt), size=7) + 
        scale_color_manual(values=colorPalette, guide=F) + labs(y = "Haplotype --- Count") +
        theme(axis.line=element_blank(),
              #axis.text.x=element_blank(),
              axis.text.y = element_text(hjust=0),
              #axis.text.y=element_blank(),
              axis.ticks=element_blank(),
              axis.title.x=element_blank(),
              # axis.title.y=element_blank(),
              legend.position="none",
              panel.background=element_blank(),
              panel.border=element_blank(),
              panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              plot.background=element_blank(),
              text=element_text(family="Courier"))
      file_name = paste("~/Projects/flu_project/haplotype_figures/", "Sample ", s, "_", d, "_", "haplo", hap_num, "_illumina_linkages.pdf", sep="")
      save_plot(file_name, p1, base_height = 7, base_width = 12)
    }
  }
}