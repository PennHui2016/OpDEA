library(countdata)
library(MSnbase)
library(dplyr)
source('/home/hui/PycharmProjects/DL_proteomics/quant_funcs/preprocessing_counts.R')

args <- commandArgs(trailingOnly = TRUE)
#print(args)
protein_file<-args[1]
peptide_file<-args[2]
design_file<-args[3]
msstats_file<-args[4]
in_type<-args[5]
inten_type<-args[6]
imput<-args[7]
normal<-args[8]
save_fold<-args[9]

if(imput=='blank'){
  imput=''
}

if(normal=='blank'){
  normal=''
}

if(inten_type=='basic'){
  inten_type=''
}
set.seed(123)
#design_file = '/home/hui/Documents/VM_share/DIA_data/yeast/design_yeast.txt'
file_name = protein_file #'/home/hui/Documents/VM_share/DIA_data/yeast/maxquant/txt/proteinGroups.txt'
#file_name = '/home/hui/Documents/VM_share/DIA_data/yeast/frager/combined_protein.tsv'
#save_fold = '/home/hui/PycharmProjects/DL_proteomics/benchmark_res/yeast_'
designs = read.table(design_file, quote = "", sep = '\t', header = TRUE)

# in_type = 'maxquant'
# #in_type = 'fragpipe'
# #inten_type = ''
# imput = ''
# normal = ''
prepro_res = preprocessing_raw_count(file_name, design_file, in_type, imput = imput, normal = normal)

#condition = as.factor(designs$condition)
#design = model.matrix(~0+condition) # fitting without intercept

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
  groups<-consts$groups
  g1<-groups[i,][1]
  g2<-groups[i,][2]
  
  counts<-prepro_res$normed[,c(grep(g1, designs$condition), grep(g2, designs$condition))]
  condition = as.factor(designs[c(grep(g1, designs$condition), grep(g2, designs$condition)),]$condition)
  #design = model.matrix(~0+condition) # fitting without intercept
  counts<-as.data.frame(counts)
  
  counts$na_g1 = apply(counts[,c(grep(g1, colnames(counts)))],1,function(x) sum(is.na(x)))
  counts$na_g2 = apply(counts[,c(grep(g2, colnames(counts)))],1,function(x) sum(is.na(x)))
  # Filter protein table. DEqMS require minimum two values for each group.
  filter_idx = which(counts$na_g1<2 & counts$na_g2<2)
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
colnames(res_all)[1:5]<-c('protein', 'pvalue', 'logFC', 'adj.pvalue', 'contrast')
res_all<-res_all[order(res_all$pvalue),]
write.table(res_all, paste0(save_fold, 'beta_binomial_', in_type, '_', inten_type, '_', imput, '_', normal, '.csv'), sep = ',', row.names = TRUE, col.names = TRUE)