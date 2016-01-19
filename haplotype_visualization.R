library(cowplot)
library(dplyr)

segment = "HA"
haplotype_file = paste("~/Projects/flu_project/pacbio_haplotypes/", segment, "_tidy_haplotypes.csv", sep="")
haplotypes = read.csv(haplotype_file)
samples_compared = unique(haplotypes$sample)
main_title = paste("Sample", samples_compared[1], "vs", samples_compared[2])

haplotypes <- transform(haplotypes, haplotype.order  = factor(haplotype, levels=sort(unique(haplotypes$haplotype)),
                                                              ordered =TRUE))

cbPalette1 <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928')

p1 <- ggplot(haplotypes, aes(x = as.factor(day), y = freq, order = haplotype.order)) + geom_bar(stat="identity", position="stack", aes(fill=haplotype.order)) +
  labs(title=main_title, x = "Day", y = "Haplotype Frequency") + scale_fill_manual(values=cbPalette1, name ="Haplotypes") + facet_wrap(~sample) +
  theme(text=element_text(family="Courier"))
  
file_name = paste("~/Projects/flu_project/haplotype_figures/", "Sample_", samples_compared[1], "_vs_", samples_compared[2], ".pdf", sep="")
n_days = length(unique(haplotypes$day))
save_plot(file_name, p1, base_height =6 , base_width = 6.0+1.5*n_days)


### supplemental illumina snps figure
illumina_file = paste("~/Projects/flu_project/pacbio_haplotypes/", segment, "_illumina_minor_variants.csv", sep="")
illumina = read.csv(illumina_file)

illumina <- transform(illumina, nt.order  = factor(nt, levels=c('A','C', 'T', 'G', ""),
                      ordered =TRUE))

colorPalette <- c("A" = '#e41a1c',"C" = '#377eb8',"G" = '#4daf4a',"T" = '#984ea3'," " = '#ff7f00')

samples_compared = unique(illumina$sample)

for(i in 1:length(samples_compared)){
  illumina %>% filter(sample == samples_compared[i]) -> sample_data
  main_title = paste("Illumina frequencies", "for sample", samples_compared[i])
  p2 <- ggplot(sample_data, aes(x = as.factor(day), y = freq, order = nt.order)) + geom_bar(stat="identity", position="stack", aes(fill=nt.order)) +
      facet_wrap(~position, nrow =1) + labs(title=main_title, x = "Day", y = "Nucleotide Frequency") + scale_fill_manual(values=colorPalette, name ="nt") +
      theme(text=element_text(family="Courier"))

  # line plot option 
  #p2 <- ggplot(sample_data, aes(x = as.factor(day), y = freq)) + geom_line(aes(color=nt, group=nt)) + geom_point(aes(color=nt)) +
   # facet_wrap(~position, nrow =1) + labs(title=main_title, x = "Day", y = "Nucleotide Frequency") + scale_fill_manual(values=colorPalette) + scale_color_manual(values=colorPalette)

  file_name = paste("~/Projects/flu_project/haplotype_figures/", "Sample ", samples_compared[i], "_illumina_frequencies.pdf", sep="")
  save_plot(file_name, p2, base_height = 6, base_width = 9)
}


####
haplotype_file = paste("~/Projects/flu_project/pacbio_haplotypes/", segment, "_tidier_haplotypes.csv", sep="")
haplotypes = read.csv(haplotype_file)
samples_compared = unique(haplotypes$sample)

for(i in 1:length(samples_compared)){
  haplotypes %>% filter(sample == samples_compared[i]) -> sample_data
  hap_length=nchar(toString(haplotypes$haplotype[1]))
  sample_data %>% group_by(haplotype) %>% mutate(total_count = sum(count)/hap_length)%>% mutate(haplotype_count = paste(substr(toString(haplotype),1,hap_length), "---", total_count)) -> sample_data
  sample_data$haplotype_count <- reorder(sample_data$haplotype_count, sample_data$total_count)
  p3 <-ggplot(sample_data, aes(x = as.factor(variant_position), y = haplotype_count, group=haplotype_count)) + geom_line(linetype=2) +geom_point(aes(color=nt), size=7) + 
    scale_color_manual(values=colorPalette, guide=F) + labs(y = "Haplotype --- Count") +
    theme(axis.line=element_blank(),
                                 axis.text.x=element_blank(),
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
  
  illumina %>% filter(sample == samples_compared[i]) -> sample_data2
  main_title = paste("Illumina frequencies", "for sample", samples_compared[i])
  p2 <- ggplot(sample_data2, aes(x = as.factor(day), y = freq, order = nt.order)) + geom_bar(stat="identity", position="stack", aes(fill=nt.order)) +
    facet_wrap(~position, nrow =1) + labs(title=main_title, x = "Day", y = "Nucleotide Frequency") + scale_fill_manual(values=colorPalette, name ="nt") +
    theme(text=element_text(family="Courier"))
  
  
  file_name = paste("~/Projects/flu_project/haplotype_figures/", "Sample ", samples_compared[i], "_hap_freq_composite.pdf", sep="")
  save_plot(file_name, (ggdraw() + draw_plot(p3, 0, .5, 0.95, .5) + draw_plot(p2, 0, 0, 1, .5)), base_height=6, base_width = 9)
  
}