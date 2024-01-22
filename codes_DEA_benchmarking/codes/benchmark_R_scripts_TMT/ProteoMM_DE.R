library(ProteoMM)
source('D:/data/benchmark/codes/benchmark_R_scripts/preprocessing_peptides.R')
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

if(inten_type == 'top0'){
  file_name = paste0(maxtrix_folder, dataset, '_', platform, '_pep_intensity.tsv')
}else if(inten_type == 'LFQ'){
  file_name = paste0(maxtrix_folder, dataset, '_', platform, '_pep_maxlfq.tsv')
}
design_file = paste0(maxtrix_folder, '', dataset, '_', platform, '_design.tsv') 
#print(design_file)
designs = read.table(design_file, sep = '\t', header = TRUE)

set.seed(123)

prepro_res = preprocessing_raw_pep(file_name, designs, platform, inten_type, imput = imput, normal = normal)

all_conts<-function(design_data){
  conds<-design_data$condition
  uni_con<-unique(conds)
  conts<-c()
  groups<-vector()
  #levels<-c()
  levels<-c(paste0('condition', uni_con))
  for (i in 1:(length(uni_con)-1)) {
    for(j in (i+1):length(uni_con)){
      conts<-c(conts, paste0('condition', uni_con[j], '-', 'condition', uni_con[i]))
      groups<-rbind(groups, c(uni_con[i], uni_con[j]))
    }
  }
  return(list(conts=conts, groups=groups, levels=levels))
}

consts<-all_conts(designs)


res_all<-vector()
for (k in 1:length(consts$conts)) {
  g1 = consts$groups[k,][1]
  g2 = consts$groups[k,][2]
  sample_names1<-designs$sample_name[grep(g1,designs$condition)]
  sample_names2<-designs$sample_name[grep(g2,designs$condition)]
  idx_g1=match(sample_names1,colnames(prepro_res$normed))
  idx_g2=match(sample_names2,colnames(prepro_res$normed))
  test_col = c(idx_g1, idx_g2)
  data_processed<-cbind(2^prepro_res$normed,prepro_res$filtered[,1:3])
  intsCols = test_col
  metaCols = colnames(data_processed[c((length(colnames(data_processed))-2):length(colnames(data_processed)))])
  # if(in_type=='maxquant'){
  #   metaCols = c(grep('Sequence', colnames(prepro_res$processed)), which(colnames(prepro_res$processed)=='Proteins'))
  # }else if(in_type=='fragpipe'){
  #   metaCols = c(grep('Peptide.Sequence', colnames(prepro_res$processed)), which(colnames(prepro_res$processed)=='Protein'))
  # }
  
  m_logInts = make_intencities(data_processed, intsCols)
  m_prot.info = make_meta(data_processed, metaCols)
  m_logInts = convert_log2(m_logInts)
  
  m_logInts$na_g1 = apply(m_logInts[,c(grep(g1, colnames(m_logInts)))],1,function(x) sum(is.na(x)))
  m_logInts$na_g2 = apply(m_logInts[,c(grep(g2, colnames(m_logInts)))],1,function(x) sum(is.na(x)))
  # Filter protein table. DEqMS require minimum two values for each group.
  filter_idx = which(m_logInts$na_g1<length(idx_g1) & m_logInts$na_g2<length(idx_g2))
  m_logInts = m_logInts[filter_idx,][,1:(length(m_logInts[1,])-2)]
  m_prot.info = m_prot.info[filter_idx,]
  
  grps = as.factor(designs$condition[c(idx_g1, idx_g2)])
  
  set.seed(135) # results rarely vary due to the random seed for EigenMS
  

  mm_m_ints_norm = list(normalized=cbind(m_prot.info, m_logInts))
  
  #already imputed
  data_imp = mm_m_ints_norm$normalized[,4:length(mm_m_ints_norm$normalized[1,])]
  data_info = mm_m_ints_norm$normalized[,1:3]
  
  
  res.ProteoMM = peptideLevel_DE(data_imp,
                           grps, data_info,
                           pr_ppos=2)
  
  res.ProteoMM$contrast<-consts$conts[k]
  res_all<-rbind(res_all, res.ProteoMM)
}

# logFC<-c()
# for (i in 1:length(res_all$FC)) {
#   if (res_all$FC[i]<0){
#     logFC<-c(logFC, -abs(log2(abs(res_all$FC[i]))))
#   }else{
#     logFC<-c(logFC, abs(log2(abs(res_all$FC[i]))))
#   }
# }
# 
# res_all$logFC<-logFC
colnames(res_all)[c(1,2,3,4)]<-c('protein', 'logFC','pvalue','adj.pvalue')
res_all<-res_all[order(res_all$pvalue),]
if(print_lab=='T'){
  true_organisms=strsplit(true_organism,';',fixed = T)[[1]]
  DE_organisms=strsplit(DE_organism,';',fixed = T)[[1]]
  labs<-get_labels(res_all, true_organisms, DE_organisms, true_lgfc)
  res_all$organism<-labs$Organism
  res_all$label<-labs$DEP
  res_all$TlogFC<-labs$TlogFC
}
write.table(res_all, paste0(save_fold, dataset, '_ProteoMM_', platform, '_', inten_type, '_', imput, '_', normal, '.csv'), sep = ',', row.names = TRUE, col.names = TRUE)

#write.table(res_all, paste0(save_fold, 'ProteoMM_', in_type, '_', inten_type, '_', imput, '_', normal, '.csv'), sep = ',', row.names = TRUE, col.names = TRUE)

