library(siggenes)
source('D:/data/benchmark/codes/benchmark_R_scripts_TMT/preprocessing_pro_intensity.R')
source('D:/data/benchmark/codes/benchmark_R_scripts_TMT/utils.R')
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

# maxtrix_folder<- "D:/data/benchmark/data/DDA/FragPipe/"
# mstats_file<- "E:/MS_data/PXD028735/FragPipe/5600/MSstats.csv"
# evidence_file<- "NULL"
# msstats_design_file<-"D:/data/benchmark/data/DDA/FragPipe/HYE5600735_LFQ_FragPipe_design_msstats.tsv"
# main_output_protein_file<-"E:/MS_data/PXD028735/FragPipe/5600/combined_protein.tsv"
# main_output_peptide_file<-"E:/MS_data/PXD028735/FragPipe/5600/combined_peptide.tsv"
# platform<-"FragPipe"
# inten_type<-"top0"
# imput<-"blank"
# normal<-"blank"
# dataset<-"HYE5600735_LFQ"
# save_fold<-"D:/data/benchmark/benchmark_res/DDA/FragPipe/HYE5600735_LFQ/"
# print_lab<-"T"
# true_organism<-"HUMAN;YEAST;ECOLI"
# DE_organism<-"YEAST;ECOLI"
# true_lgfc<-"D:/data/benchmark/data/dataset_info/PXD028735_true_fc.txt"
# logT<-"T"

if(logT=='T'){
  logT=T
}else if(logT=='F'){
  logT=F
}

if(imput=='blank'){
  imput=''
}

if(normal=='blank'){
  normal=''
}

if(inten_type == 'abd'){
  file_name = paste0(maxtrix_folder, dataset, '_', platform, '_tmt_abd.tsv')
}else if(inten_type == 'ratio'){
  file_name = paste0(maxtrix_folder, dataset, '_', platform, '_tmt_ratio.tsv')
}else if(inten_type == 'phi'){
  file_name = paste0(maxtrix_folder, dataset, '_', platform, '_tmt_phi.tsv')
}else if(inten_type == 'intensity'){
  file_name = paste0(maxtrix_folder, dataset, '_', platform, '_tmt_intensity.tsv')
}

design_file = paste0(maxtrix_folder, '', dataset, '_', platform, '_tmt_design.tsv') 
#print(design_file)
designs = read.table(design_file, sep = '\t', header = TRUE)

set.seed(123)

prepro_res = preprocessing_raw(file_name, designs, platform, inten_type, imput = imput, normal = normal, log2=logT)

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
  
  intens<-prepro_res$normed[,c(idx_g1, idx_g2)]
  condition = as.factor(designs[c(grep(g1, designs$condition), grep(g2, designs$condition)),]$condition)
  design = model.matrix(~0+condition) # fitting without intercept
  intens<-as.data.frame(intens)
  
  intens$na_g1 = apply(intens[,grep(g1,condition)],1,function(x) sum(is.na(x)))
  intens$na_g2 = apply(intens[,grep(g2,condition)],1,function(x) sum(is.na(x)))
  # Filter protein table. DEqMS require minimum two values for each group.
  #filter_idx = which(intens$na_g1<2 & intens$na_g2<2)
  filter_idx = which((length(idx_g1)-intens$na_g1)>=2 & (length(idx_g2)-intens$na_g2)>=2)
  intens.filter = intens[filter_idx,][,1:(length(intens)-2)]
  
  label<-c(rep(0, length(condition)))
  label[which(condition==g2)]=1
  
  intens.filter[is.na(intens.filter)] = 0
  #ebam.out<-ebam(intens.filter, label, method=wilc.ebam, rand =123, gene.names = row.names(intens.filter))
  sam.out <- sam(intens.filter, label, rand = 123, gene.names = row.names(intens.filter))
  sam.res<-cbind(row.names(intens.filter),sam.out@d,sam.out@d.bar,sam.out@vec.false,sam.out@p.value,sam.out@s,sam.out@q.value,log2(sam.out@fold),rep(consts$conts[i],length(sam.out@d)))
  colnames(sam.res)<-c('protein','d','d.bar','vec.fasle','p.value','s','q.value','log2FC','contrast')
  res_all<-as.data.frame(rbind(res_all, sam.res))
  }
}
#colnames(res_all)[length(res_all[1,])]<-'contrast'
colnames(res_all)[c(5,7,8)]<-c('pvalue','adj.pvalue','logFC')
res_all<-res_all[order(res_all$pvalue),]
#res_all$protein<-row.names(res_all)
if(print_lab=='T'){
  true_organisms=strsplit(true_organism,';',fixed = T)[[1]]
  DE_organisms=strsplit(DE_organism,';',fixed = T)[[1]]
  labs<-get_labels(res_all, true_organisms, DE_organisms, true_lgfc)
  res_all$organism<-labs$Organism
  res_all$label<-labs$DEP
  res_all$TlogFC<-labs$TlogFC
}
write.table(res_all, paste0(save_fold, dataset, '_siggenes_', platform, '_', inten_type, '_', imput, '_', normal, '.csv'), sep = ',', row.names = TRUE, col.names = TRUE)

#write.table(res_all, paste0(save_fold, 'siggenes_', in_type, '_', inten_type, '_', imput, '_', normal, '.csv'), sep = ',', row.names = TRUE, col.names = TRUE)
