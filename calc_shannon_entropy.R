library(dplyr)
dir = "~/Projects/flu_project/illumina_snplists/"
files = list.files(path = dir, pattern = ".csv")
meta_data = read.csv("~/Projects/flu_project/metadata.csv")
meta_data$Ferret.ID = as.character(meta_data$Ferret.ID)
snp_cutoff = 0.03

shannon_entropy = data.frame(sample = NA, day = NA, segment = NA, generation = NA, exposure = NA, exposure_type = NA, entropy = NA)

trim <- function (x) gsub("^\\s+|\\s+$", "", x)
relevant_info= data.frame(segment = NA, generation = NA, t_mean_naive = NA, t_mean_pre = NA , t_pvalue= NA, tp_p1= NA, tp_p2= NA, tp_p3= NA, tp_p4= NA, stringsAsFactors=FALSE) 
names(relevant_info) <- c("segment", "generation", "mean naive entropy", "mean pre-immunized entropy", "p-value (naive vs pre-immunized)", "p-value (naive vs A/Brisbane/59/07)", "p-value (naive vs A/Denver/1/57)", "p-value (naive vs A/PR/8/34)", "p-value (naive vs A/Texas/36/91)")
row.names(relevant_info) <- NULL

for(file in files){
    data = read.csv(paste(dir, file, sep=""), na.strings = "" )
    if (grepl("Stock", file)){
      next
    }
    else{
      sample = toString(unlist(strsplit(toString(file), "_"))[1])
      day = as.integer(unlist(strsplit(toString(unlist(strsplit(toString(file), "\\."))[1]), "_"))[2])
      generation = toString(meta_data$Generation[which(meta_data$Ferret.ID == sample)[1]])
      exposure = toString(meta_data$Pre.Challenge.Exposure.1[which(meta_data$Ferret.ID == sample)[1]])
      exposure_type = trim(toString(meta_data$Pre.Challenge.Exposure[which(meta_data$Ferret.ID == sample)[1]]))
    }
    data %>% filter(majorfreq != 1 & minorfreq >= snp_cutoff) %>% group_by(segment) %>% summarize(entropy = sum(-(majorfreq * log2(majorfreq) + minorfreq * log2(minorfreq)))) -> entropy
    for (r in 1:nrow(entropy)){
      newRow = c(sample, day, entropy[r, 1], generation, exposure, exposure_type, entropy[r, 2])
      names(newRow) <- c("sample", "day", "segment", "generation", "exposure", "exposure_type","entropy")
      shannon_entropy = rbind(shannon_entropy, newRow)
    }
}

shannon_entropy <- shannon_entropy[-1, ]
generations = c("F0", "F1")
segments = c("HA", "MP",  "NA",  "NP",  "NS",  "PA",  "PB1", "PB2")
for (seg in segments){
  for (gen in generations){
    shannon_entropy %>% filter(generation == gen & segment == seg) -> shannon_entropy_filtered
    shannon_entropy_filtered %>% filter(exposure == "Naive") %>% select(entropy) -> naive
    shannon_entropy_filtered %>% filter(exposure == "Preimm") %>% select(entropy) -> preimm
    t <- t.test(naive, preimm)
    pt <- pairwise.t.test(shannon_entropy_filtered$entropy, shannon_entropy_filtered$exposure_type)
    t_pvalue = t$p.value
    t_mean_naive = t$estimate[1]
    t_mean_pre = t$estimate[2]
    tp_p1 = pt$p.value[4]
    tp_p2 = pt$p.value[8]
    tp_p3= pt$p.value[12]
    tp_p4 = pt$p.value[16]
    
    newRow = c(seg, gen, t_mean_naive, t_mean_pre, t_pvalue, tp_p1, tp_p2, tp_p3, tp_p4) 
    relevant_info = rbind(relevant_info, newRow)
  }
  ### both generations combined
  gen = "Both"
  shannon_entropy %>% filter(segment == seg) -> shannon_entropy_filtered
  shannon_entropy_filtered %>% filter(exposure == "Naive") %>% select(entropy) -> naive
  shannon_entropy_filtered %>% filter(exposure == "Preimm") %>% select(entropy) -> preimm
  t <- t.test(naive, preimm)
  pt <- pairwise.t.test(shannon_entropy_filtered$entropy, shannon_entropy_filtered$exposure_type)
  t_pvalue = t$p.value
  t_mean_naive = t$estimate[1]
  t_mean_pre = t$estimate[2]
  tp_p1 = pt$p.value[4]
  tp_p2 = pt$p.value[8]
  tp_p3= pt$p.value[12]
  tp_p4 = pt$p.value[16]
  
  newRow = c(seg, gen, t_mean_naive, t_mean_pre, t_pvalue, tp_p1, tp_p2, tp_p3, tp_p4) 
  relevant_info = rbind(relevant_info, newRow)
}

### All segments combined 
seg = "ALL"
for (gen in generations){
  shannon_entropy %>% filter(generation == gen) -> shannon_entropy_filtered
  shannon_entropy_filtered %>% filter(exposure == "Naive") %>% select(entropy) -> naive
  shannon_entropy_filtered %>% filter(exposure == "Preimm") %>% select(entropy) -> preimm
  t <- t.test(naive, preimm)
  pt <- pairwise.t.test(shannon_entropy_filtered$entropy, shannon_entropy_filtered$exposure_type)
  t_pvalue = t$p.value
  t_mean_naive = t$estimate[1]
  t_mean_pre = t$estimate[2]
  tp_p1 = pt$p.value[4]
  tp_p2 = pt$p.value[8]
  tp_p3= pt$p.value[12]
  tp_p4 = pt$p.value[16]
  
  newRow = c(seg, gen, t_mean_naive, t_mean_pre, t_pvalue, tp_p1, tp_p2, tp_p3, tp_p4) 
  relevant_info = rbind(relevant_info, newRow)
}
### both generations combined, all segments
gen = "Both"
shannon_entropy -> shannon_entropy_filtered
shannon_entropy_filtered %>% filter(exposure == "Naive") %>% select(entropy) -> naive
shannon_entropy_filtered %>% filter(exposure == "Preimm") %>% select(entropy) -> preimm
t <- t.test(naive, preimm)
pt <- pairwise.t.test(shannon_entropy_filtered$entropy, shannon_entropy_filtered$exposure_type)
t_pvalue = t$p.value
t_mean_naive = t$estimate[1]
t_mean_pre = t$estimate[2]
tp_p1 = pt$p.value[4]
tp_p2 = pt$p.value[8]
tp_p3= pt$p.value[12]
tp_p4 = pt$p.value[16]

newRow = c(seg, gen, t_mean_naive, t_mean_pre, t_pvalue, tp_p1, tp_p2, tp_p3, tp_p4) 
relevant_info = rbind(relevant_info, newRow)

relevant_info <- relevant_info[-1, ]
write.csv(relevant_info, file = "~/Projects/flu_project/entropy_analysis.csv", row.names = FALSE)

write.csv(shannon_entropy, file = "~/Projects/flu_project/entropy_analysis_0.03_full.csv", row.names = FALSE)
