library(DEqMS)
#source('/home/hui/PycharmProjects/DL_proteomics/quant_funcs/preprocessing.R')
library(matrixStats)

source('D:/data/benchmark/codes/benchmark_R_scripts/preprocessing_pro_intensity_deqms.R')
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

# maxtrix_folder<-'D:/data/benchmark/data/DDA/Maxquant/'
# mstats_file<-'NULL'
# evidence_file<-'E:/MS_data/PXD028735/maxquant/5600/new/combined/txt/evidence.txt'
# msstats_design_file<-'D:/data/benchmark/data/DDA/Maxquant/HYE5600735_LFQ_Maxquant_design_msstats.tsv'
# main_output_protein_file<-'E:/MS_data/PXD028735/maxquant/5600/new/combined/txt/proteinGroups.txt'
# main_output_peptide_file<-'E:/MS_data/PXD028735/maxquant/5600/new/combined/txt/peptides.txt'
# platform<-'Maxquant'
# inten_type<-'dlfq'
# imput<-'blank'
# normal<-'blank'
# dataset<-'HYE5600735_LFQ'
# save_fold<-'D:/data/benchmark/benchmark_res/DDA/Maxquant/HYE5600735_LFQ/'
# print_lab<-'T'
# true_organism<-'HUMAN;YEAST;ECOLI'
# DE_organism<-'YEAST;ECOLI'
# true_lgfc<-'D:/data/benchmark/data/dataset_info/PXD028735_true_fc.txt'
# logT<-'T'
# conss<-'conditionB-conditionA'


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

if(inten_type == 'top0'){
  file_name = paste0(maxtrix_folder, dataset, '_', platform, '_pro_intensity.tsv')
}else if(inten_type == 'top3'){
  file_name = paste0(maxtrix_folder, dataset, '_', platform, '_top3_pro_intensity.tsv')
}else if(inten_type == 'LFQ'){
  file_name = paste0(maxtrix_folder, dataset, '_', platform, '_pro_maxlfq.tsv')
}else if(inten_type == 'dlfq'){
  file_name = paste0(maxtrix_folder, dataset, '_', platform, '_dlfq_pro_intensity.tsv')
}
count_file = paste0(maxtrix_folder, dataset, '_', platform, '_pro_count.tsv')

design_file = paste0(maxtrix_folder, '', dataset, '_', platform, '_design.tsv') 
#print(design_file)
designs = read.table(design_file, sep = '\t', header = TRUE)

set.seed(123)

prepro_res = preprocessing_raw(file_name, count_file, designs, platform, inten_type, imput = imput, normal = normal, log2=logT)

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
  
  # use minimum peptide count among  samples
  # count unique+razor peptides used for quantification
  counts<-prepro_res$counts
  
  col_idx_pep_c<-c(match(sample_names1,colnames(counts)), match(sample_names2,colnames(counts)))
  count.table = data.frame(count = rowMins(as.matrix(counts[filter_idx,][,col_idx_pep_c])),
                                                               row.names = counts[,1][filter_idx])
  count.table$count = count.table$count+1
  # if(in_type=='maxquant'){
  #   col_idx_pep_c<-c(grep(paste0('Razor...unique.peptides.', groups[i,][1]), colnames(prepro_res$filtered)), grep(paste0('Razor...unique.peptides.', groups[i,][2]), colnames(prepro_res$filtered)))
  #   pep.count.table = data.frame(count = rowMins(as.matrix(prepro_res$filtered[filter_idx,][,col_idx_pep_c])),
  #                                row.names = prepro_res$filtered$Majority.protein.IDs[filter_idx])
  #   # Minimum peptide count of some proteins can be 0
  #   # add pseudocount 1 to all proteins
  #   pep.count.table$count = pep.count.table$count+1
  # }else if(in_type=='fragpipe'){
  #   spc_idx<-c()
  #   for (j in 1:length(designs$sample.name)) {
  #     idx<-which(colnames(prepro_res$filtered)==paste0(designs$sample.name[j],'.Total.Spectral.Count'))
  #     spc_idx<-c(spc_idx, idx)
  #   }
  #   counts<-prepro_res$filtered[filter_idx,][,spc_idx]
  #   col_idx_pep_c<-c(grep(groups[i,][1], colnames(counts)), grep(groups[i,][2], colnames(counts)))
  #   pep.count.table = data.frame(count = rowMins(as.matrix(counts[,col_idx_pep_c])),
  #                                row.names = prepro_res$filtered$Protein[filter_idx])
  #   # Minimum peptide count of some proteins can be 0
  #   # add pseudocount 1 to all proteins
  #   pep.count.table$count = pep.count.table$count+1
  # }
  fit1 = lmFit(intens.filter,design = design)
  cont <- makeContrasts(contrasts =  consts$conts[i], levels = design)
  fit2 = contrasts.fit(fit1,contrasts = cont)
  fit3 <- eBayes(fit2)
  fit3$count = count.table[rownames(fit3$coefficients),"count"]
  fit3$sigma[which(fit3$sigma==0)]=1e-14
  fit4 = spectraCounteBayes(fit3)
  DEqMS.results = outputResult(fit4, coef_col = 1)
  #limma.results = topTable(fit3, adjust="BH", sort.by = 'logFC', n=Inf)
  DEqMS.results<-cbind(DEqMS.results, rep(consts$conts[i], length(DEqMS.results[,1])))
  res_all<-rbind(res_all, DEqMS.results)
  }
}
colnames(res_all)[length(res_all[1,])]<-'contrast'
colnames(res_all)[c(7, 10,11)]<-c('protein', 'pvalue','adj.pvalue')
res_all<-res_all[order(res_all$pvalue),]
if(print_lab=='T'){
  true_organisms=strsplit(true_organism,';',fixed = T)[[1]]
  DE_organisms=strsplit(DE_organism,';',fixed = T)[[1]]
  labs<-get_labels(res_all, true_organisms, DE_organisms, true_lgfc)
  res_all$organism<-labs$Organism
  res_all$label<-labs$DEP
  res_all$TlogFC<-labs$TlogFC
}
write.table(res_all, paste0(save_fold, dataset, '_DEqMS_', platform, '_', inten_type, '_', imput, '_', normal, '.csv'), sep = ',', row.names = TRUE, col.names = TRUE)

#write.table(res_all, paste0(save_fold, 'DEqMS_', in_type, '_', inten_type, '_', imput, '_', normal, '.csv'), sep = ',', row.names = TRUE, col.names = TRUE)

