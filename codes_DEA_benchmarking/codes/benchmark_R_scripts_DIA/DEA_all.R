#library(limma)
source('D:/data/benchmark/codes/benchmark_R_scripts/preprocessing_pro_intensity.R')
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
dea<-args[19]

if(imput=='blank'){
  imput=''
}

if(normal=='blank'){
  normal=''
}

if(logT=='T'){
  logT=T
}else if(logT=='F'){
  logT=F
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

design_file = paste0(maxtrix_folder, '', dataset, '_', platform, '_design.tsv') 
#print(design_file)
designs = read.table(design_file, sep = '\t', header = TRUE)

set.seed(123)

rds_file = paste0('D:/data/benchmark/benchmark_res/RDS_file/', dataset, '_', platform, '_', inten_type, '_',
                  normal, '_', imput, '_', logT, '.rds')
if(file.exists(rds_file)){
  prepro_res<-readRDS(file=rds_file)
}else{
  prepro_res = preprocessing_raw(file_name, designs, platform, inten_type, imput = imput, normal = normal, log2=logT)
  
  saveRDS(prepro_res, file=rds_file)
}


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

# 'limma', 'ttest'
# 'DEP', 'ANOVA'
# 'SAM', 'ROTS'
# 'proDA'


if(dea=='limma'){
  source('D:/data/benchmark/codes/benchmark_R_scripts_DIA/limma_DEA.R')
  
  
}else if(dea == 'ttest'){
  
}else if(dea == 'DEP'){
  
}else if(dea == 'ANOVA'){
  
}else if(dea == 'SAM'){
  
}else if(dea == 'ROTS'){
  
}else if(dea == 'proDA'){
  
}
