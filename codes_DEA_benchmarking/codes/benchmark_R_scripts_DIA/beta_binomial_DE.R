library(countdata)
library(MSnbase)
library(dplyr)
source('D:/data/benchmark/codes/benchmark_R_scripts/preprocessing_pro_counts.R')
source('D:/data/benchmark/codes/benchmark_R_scripts/utils.R')
args <- commandArgs(trailingOnly = TRUE)
print(args)
maxtrix_folder<-args[1]
mstats_file<-args[2]
evidence_file<-args[3]
msstats_design_file<-args[4]
main_output_protein_file<-args[5]
main_output_peptide_file<-args[6]
platform<-args[7]
inten_type<-args[8]
imput<-args[9]
normal<-args[10]
dataset<-args[11]
save_fold<-args[12]
print_lab<-args[13]
true_organism<-args[14]
DE_organism<-args[15]
true_lgfc<-args[16]
logT<-args[17]
conss<-args[18]

# maxtrix_folder<-"D:/data/benchmark/data/DDA/FragPipe/"
# mstats_file<-"E:/MS_data/PXD002099/FragPipe/MSstats.csv"
# evidence_file<-"NULL"
# msstats_design_file<-"D:/data/benchmark/data/DDA/FragPipe/YUltq099_LFQ_FragPipe_design_msstats.tsv"
# main_output_protein_file<-"E:/MS_data/PXD002099/FragPipe/combined_protein.tsv"
# main_output_peptide_file<-"E:/MS_data/PXD002099/FragPipe/combined_peptide.tsv"
# platform<-"FragPipe"
# inten_type<-"LFQ"
# imput<-"Impseq"
# normal<-"lossf"
# dataset<-"YUltq099_LFQ"
# save_fold<-"D:/data/benchmark/benchmark_res/DDA/FragPipe/YUltq099_LFQ/"
# print_lab<-"T"
# true_organism<-"YEAST;UPS"
# DE_organism<-"UPS"
# true_lgfc<-"D:/data/benchmark/data/dataset_info/PXD002099_true_fc.txt"

if(imput=='blank'){
  imput=''
}

if(normal=='blank'){
  normal=''
}

file_name = paste0(maxtrix_folder, dataset, '_', platform, '_pro_count.tsv')


design_file = paste0(maxtrix_folder, '', dataset, '_', platform, '_design.tsv') 
#print(design_file)
designs = read.table(design_file, sep = '\t', header = TRUE)

set.seed(123)

prepro_res = preprocessing_raw_count(file_name, designs, platform, inten_type, imput = imput, normal = normal)

all_conts<-function(design_data){
  conds<-design_data$condition
  uni_con<-unique(conds)
  conts<-c()
  groups<-vector()
  for (i in 1:(length(uni_con)-1)) {
    for(j in (i+1):length(uni_con)){
      conts<-c(conts, paste0('condition', uni_con[j], '-', 'condition', uni_con[i]))
      groups<-rbind(groups, c(uni_con[i], uni_con[j]))
    }
  }
  return(list(conts=conts, groups=groups))
}

consts<-all_conts(designs)
conss<-strsplit(conss, ';', fixed = T)[[1]]

compute_logFC<-function(group1, group2){
  logFC<-c()
  for (i in 1:length(group1[,1])) {
    logfc<-mean(group2[i,][which(group2[i,]!=0 & !is.na(group2[i,]))])-mean(group1[i,][which(group1[i,]!=0 & !is.na(group1[i,]))])
    logFC<-c(logFC, logfc)
  }
  return(logFC)
}

res_all<-vector()
for (i in 1:length(consts$conts)) {
  if(length(grep(consts$conts[i], conss))==1){
  groups<-consts$groups
  g1<-groups[i,][1]
  g2<-groups[i,][2]
  
  sample_names1<-designs$sample_name[grep(g1,designs$condition)]
  sample_names2<-designs$sample_name[grep(g2,designs$condition)]
  idx_g1=match(sample_names1,colnames(prepro_res$normed))
  idx_g2=match(sample_names2,colnames(prepro_res$normed))
  
  counts<-prepro_res$normed[,c(idx_g1, idx_g2)]
  condition = as.factor(designs[c(grep(g1, designs$condition), grep(g2, designs$condition)),]$condition)
  
  counts<-as.data.frame(counts)
  
  counts$na_g1 = apply(counts[,grep(g1,condition)],1,function(x) sum(is.na(x)))
  counts$na_g2 = apply(counts[,grep(g2,condition)],1,function(x) sum(is.na(x)))
  
  # Filter protein table. DEqMS require minimum two values for each group.
  #filter_idx = which(counts$na_g1<2 & counts$na_g2<2)
  filter_idx = which((length(idx_g1)-counts$na_g1)>=2 & (length(idx_g2)-counts$na_g2)>=2)
  counts.filter = counts[filter_idx,][,1:(length(counts)-2)]
  ## edgeR can't process NA, use 0 instead
  counts.filter<-as.matrix(counts.filter)
  counts.filter[is.na(counts.filter)] = 0
  counts.filter[which(counts.filter<0)] = 0
  counts.filter<-counts.filter+1
  
  x <- counts.filter
  tx <- colSums(counts.filter)
  group <- as.character(condition)
  res.beta.binom <- bb.test(x, tx, group) %>% as.data.frame
  row.names(res.beta.binom)<-row.names(counts.filter)
  res.beta.binom$log2FC<-compute_logFC(log2(counts.filter[,which(condition==g1)]), log2(counts.filter[,which(condition==g2)]))
  res.beta.binom$qvalue <- p.adjust(res.beta.binom$p.value, method = "BH")
  res.beta.binom$contrast <- consts$conts[i]
  res.beta.binom <- res.beta.binom %>% cbind(protein = rownames(.),.)
  
  res_all<-rbind(res_all, res.beta.binom)
  }
}
colnames(res_all)[1:5]<-c('protein', 'pvalue', 'logFC', 'adj.pvalue', 'contrast')
res_all<-res_all[order(res_all$pvalue),]
if(print_lab=='T'){
  true_organisms=strsplit(true_organism,';',fixed = T)[[1]]
  DE_organisms=strsplit(DE_organism,';',fixed = T)[[1]]
  labs<-get_labels(res_all, true_organisms, DE_organisms, true_lgfc)
  res_all$organism<-labs$Organism
  res_all$label<-labs$DEP
  res_all$TlogFC<-labs$TlogFC
}
write.table(res_all, paste0(save_fold, dataset, '_beta_binomial_', platform, '_', inten_type, '_', imput, '_', normal, '.csv'), sep = ',', row.names = TRUE, col.names = TRUE)

#write.table(res_all, paste0(save_fold, 'beta_binomial_', in_type, '_', inten_type, '_', imput, '_', normal, '.csv'), sep = ',', row.names = TRUE, col.names = TRUE)