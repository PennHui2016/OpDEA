library(tidyverse)
library(QFeatures)
library(msqrob2)

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
conss<-args[18]

# maxtrix_folder<- "D:/data/benchmark/data/DDA/FragPipe/"
# mstats_file<- "E:/MS_data/PXD028735/FragPipe/5600/MSstats.csv"
# evidence_file<- "NULL"
# msstats_design_file<-"D:/data/benchmark/data/DDA/FragPipe/HYE5600735_LFQ_FragPipe_design_msstats.tsv"
# main_output_protein_file<-"E:/MS_data/PXD028735/FragPipe/5600/combined_protein.tsv"
# main_output_peptide_file<-"E:/MS_data/PXD028735/FragPipe/5600/combined_peptide.tsv"
# platform<-"FragPipe"
# inten_type<-"top3"
# imput<-"blank"
# normal<-"blank"
# dataset<-"HYE5600735_LFQ"
# save_fold<-"D:/data/benchmark/benchmark_res/DDA/FragPipe/HYE5600735_LFQ/"
# print_lab<-"T"
# true_organism<-"HUMAN;YEAST;ECOLI"
# DE_organism<-"YEAST;ECOLI"
# true_lgfc<-"D:/data/benchmark/data/dataset_info/PXD028735_true_fc.txt"
# logT<-"T"

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
conss<-strsplit(conss, ';', fixed = T)[[1]]

res_all<-vector()
for (k in 1:length(consts$conts)) {
  if(length(grep(consts$conts[i], conss))==1){
  g1 = consts$groups[k,][1]
  g2 = consts$groups[k,][2]
  
  sample_names1<-designs$sample_name[grep(g1,designs$condition)]
  sample_names2<-designs$sample_name[grep(g2,designs$condition)]
  idx_g1=match(sample_names1,colnames(prepro_res$normed))
  idx_g2=match(sample_names2,colnames(prepro_res$normed))
  #idx_g1 = grep(g1, colnames(prepro_res$normed))
  #idx_g2 = grep(g2, colnames(prepro_res$normed))
  test_col = c(idx_g1, idx_g2)
  data_processed<-cbind(2^prepro_res$normed, prepro_res$filtered[,1:3])
  pe1 <- readQFeatures(
    table = data_processed, ecol = test_col,name = "peptideRaw"
  )
  colData(pe1)$condition <-as.factor(designs[c(grep(g1, designs$condition), grep(g2, designs$condition)),]$condition)
  rowData(pe1[["peptideRaw"]])$nNonZero <- rowSums(assay(pe1[["peptideRaw"]]) > 0)
  pe1 <- zeroIsNA(pe1, "peptideRaw") # convert 0 to NA
  pe1 <- logTransform(pe1, base = 2, i = "peptideRaw", name = "peptideLog")
  pe1 <- filterFeatures(pe1, ~ nNonZero >= 2)
  
  pe1 <- aggregateFeatures(pe1,
                                                        i = "peptideLog", fcol = "Protein", na.rm = TRUE,
                                                        name = "protein"
                               )
  
  # if (normal!=''){
  #   pe1 <- normalize(pe1,
  #                    i = "peptideLog",
  #                    name = "peptideNorm",
  #                    method = normal
  #   )
  #   if(in_type=='maxquant'){
  #     pe1 <- aggregateFeatures(pe1,
  #                              i = "peptideNorm", fcol = "Proteins", na.rm = TRUE,
  #                              name = "protein"
  #     )
  #   }else if(in_type=='fragpipe'){
  #     pe1 <- aggregateFeatures(pe1,
  #                              i = "peptideNorm", fcol = "Protein", na.rm = TRUE,
  #                              name = "protein"
  #     )
  #   }
  # }else{
  #   if(in_type=='maxquant'){
  #     pe1 <- aggregateFeatures(pe1,
  #                              i = "peptideLog", fcol = "Proteins", na.rm = TRUE,
  #                              name = "protein"
  #     )
  #   }else if(in_type=='fragpipe'){
  #     pe1 <- aggregateFeatures(pe1,
  #                              i = "peptideLog", fcol = "Protein", na.rm = TRUE,
  #                              name = "protein"
  #     )
  #   }
  # }
  
  
  pe <- msqrob(object = pe1, i = "protein", formula = ~condition)
  
  
  
  L <- makeContrast(contrasts=paste0('condition',g2,'=0'),
                     parameterNames = c(paste0('condition',g2)))
  pe <- hypothesisTest(object = pe, i = "protein", contrast = L)
  res.msqrob2 = pe@ExperimentList@listData[["protein"]]@elementMetadata@listData[[paste0('condition',g2)]]
  res.msqrob2$contrast<-consts$conts[k]
  res.msqrob2$protein<-row.names(res.msqrob2)
  res_all<-rbind(res_all, res.msqrob2)
  }
}
colnames(res_all)[5:6]<-c('pvalue','adj.pvalue')
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
write.table(res_all, paste0(save_fold, dataset, '_msqrob2_', platform, '_', inten_type, '_', imput, '_', normal, '.csv'), sep = ',', row.names = TRUE, col.names = TRUE)

#write.table(res_all, paste0(save_fold, 'msqrob2_', in_type, '_', inten_type, '_', imput, '_', normal, '.csv'), sep = ',', row.names = TRUE, col.names = TRUE)
