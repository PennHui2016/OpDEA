library(MSstats)
library(MSnbase)
library(readr)
#library(MSqRob)
library(tidyverse)

#source('D:/data/benchmark/codes/benchmark_R_scripts/preprocessing_pro_intensity.R')
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

# maxtrix_folder<-'D:/data/benchmark/data/DDA/FragPipe/'
# mstats_file<-'E:/MS_data/PXD028735/FragPipe/5600/TOP0/MSstats.csv'
# evidence_file<-'NULL'
# msstats_design_file<-'D:/data/benchmark/data/DDA/FragPipe/HYE5600735_LFQ_FragPipe_design_msstats.tsv'
# main_output_protein_file<-'E:/MS_data/PXD028735/FragPipe/5600/combined_protein.tsv'
# main_output_peptide_file<-'E:/MS_data/PXD028735/FragPipe/5600/combined_peptide.tsv'
# platform<-'FragPipe'
# inten_type<-'all'
# imput<-'FALSE'
# normal<-'globalStandards'
# dataset<-'HYE5600735_LFQ'
# save_fold<-'D:/data/benchmark/benchmark_res/DDA/FragPipe/HYE5600735_LFQ/'
# print_lab<-'T'
# true_organism<-'HUMAN;YEAST;ECOLI'
# DE_organism<-'YEAST;ECOLI'
# true_lgfc<-'D:/data/benchmark/data/dataset_info/PXD028735_true_fc.txt'
# logT<-'T'

# if(imput=='blank'){
#   imput=''
# }
# 
# if(normal=='blank'){
#   normal=''
# }

# if(inten_type=='basic' | inten_type=='blank'){
#   inten_type=''
# }

set.seed(123)

if(platform=='FragPipe'){
  evidence_file_name = ''
  protein_file_name = mstats_file
}else if(platform=='Maxquant'){
  evidence_file_name = evidence_file
  protein_file_name = main_output_protein_file
}else if(platform=='DIANN' | platform=='Spectronaut'){
  evidence_file_name = evidence_file
}

designs = read.table(msstats_design_file, sep = '\t', header = TRUE)

if(normal=='FALSE'){
  normal=FALSE
}

if(imput=='FALSE'){
  imput=FALSE
}else if(imput=='TRUE'){
  imput=TRUE
}

all_conts<-function(design_data){
  conds<-design_data$Condition
  uni_con<-unique(conds)
  conts<-c()
  groups<-vector()
  #levels<-c()
  levels<-c(paste0('condition', uni_con))
  code_cont<-matrix(data=0, nrow = round(length(uni_con)*(length(uni_con)-1)/2), ncol = length(uni_con))
  flag=1
  for (i in 1:(length(uni_con)-1)) {
    for(j in (i+1):length(uni_con)){
      conts<-c(conts, paste0('condition', uni_con[j], '-', 'condition', uni_con[i]))
      groups<-rbind(groups, c(uni_con[i], uni_con[j]))
      code_cont[flag,i]=-1
      code_cont[flag,j]=1
      flag=flag+1
    }
  }
  colnames(code_cont)=uni_con
  row.names(code_cont) = conts
  return(list(conts=conts, groups=groups, levels=levels, code_cont=code_cont))
}

consts<-all_conts(designs)
conss<-strsplit(conss, ';', fixed = T)[[1]]
contss_r<-c()
for (cts in conss) {
  cds<-strsplit(cts, '-', fixed = T)[[1]]
  ctr_r<-paste0(cds[2], '-', cds[1])
  contss_r<-c(contss_r, ctr_r)
}

if(platform=='Maxquant'){
  #if(length(grep(consts$conts[i], conss))==1){
  # Read in MaxQuant files
  proteinGroups <- read.table(protein_file_name, quote = "", sep = '\t', header = TRUE)
  infile <- read.table(evidence_file_name, quote = "", sep = '\t', header = TRUE)
  annot<-designs
  
  # for (run in designs$Run) {
  #   infile$Experiment[grep(run, infile$Raw.file)]=paste0(designs$Condition[grep(run, designs$Run)], designs$BioReplicate[grep(run, designs$Run)])
  #   #raw$Condition[grep(run, raw$Run)]=designs$Condition[grep(run, designs$Run)]
  # }
  
  raw <- MaxQtoMSstatsFormat(evidence = infile,
                             annotation = annot,
                             proteinGroups = proteinGroups)
  
  raw$Run<-gsub('\\+', '', raw$Run)
  designs$Run<-gsub('\\+', '', designs$Run)
  for (run in designs$Run) {
    raw$BioReplicate[grep(run, raw$Run)]=designs$BioReplicate[grep(run, designs$Run)]
    raw$Condition[grep(run, raw$Run)]=designs$Condition[grep(run, designs$Run)]
  }
  
  QuantData <- dataProcess(raw, censoredInt="NA", normalization=normal,
                           MBimpute = imput,
                           featureSubset = inten_type)
}else if (platform=='FragPipe'){
  raw <- read_csv(protein_file_name, na = c("", "NA", "0"))
  raw$BioReplicate[is.na(raw$BioReplicate)]=0
  raw$ProteinName <- factor(raw$ProteinName)
  raw$PeptideSequence <- factor(raw$PeptideSequence)
  ##
  #raw$Condition[which(raw$Condition=="A1" | raw$Condition=="A2" | raw$Condition=="A2")]='A'
  # for (m in 1:length(raw$ProteinName)) {
  #   raw$BioReplicate[m] = substr(raw$Condition[m],2,2)
  #   raw$Condition[m] = substr(raw$Condition[m],1,1)
  # }
  ##
  #uni_run<-unique(designs$Run)
  raw$Run<-gsub('\\+', '', raw$Run)
  designs$Run<-gsub('\\+', '', designs$Run)
  for (run in designs$Run) {
    raw$BioReplicate[grep(run, raw$Run)]=designs$BioReplicate[grep(run, designs$Run)]
    raw$Condition[grep(run, raw$Run)]=designs$Condition[grep(run, designs$Run)]
  }
  
  QuantData <- dataProcess(raw, censoredInt="NA", normalization=normal,
                           MBimpute = imput,
                           featureSubset = inten_type)
}
levels(QuantData$ProteinLevelData$GROUP)

testResultOneComparison <- groupComparison(contrast.matrix = consts$code_cont, data = QuantData)
res.MSstats<-testResultOneComparison$ComparisonResult%>% arrange(pvalue)
colnames(res.MSstats)[1:3]<-c('protein','contrast', 'logFC')

res.MSstats_rem<-vector()

for (m in 1:length(conss)) {
  ctr_idx<-which(res.MSstats$contrast==conss[m])
  if(length(ctr_idx)>0){
    ctr_res<-res.MSstats[ctr_idx,]
    if(length(res.MSstats_rem)==0){
      res.MSstats_rem<-ctr_res
    }else{
      res.MSstats_rem<-rbind(res.MSstats_rem, ctr_res)}
  }
}

for (m in 1:length(contss_r)) {
  ctr_idx<-which(res.MSstats$contrast==contss_r[m])
  if(length(ctr_idx)>0){
    ctr_res<-res.MSstats[ctr_idx,]
    ctr_res$contrast=conss[m]
    if(length(res.MSstats_rem)==0){
      res.MSstats_rem<-ctr_res
    }else{
      res.MSstats_rem<-rbind(res.MSstats_rem, ctr_res)}
  }
}
res.MSstats<-as.data.frame(res.MSstats_rem)


if(print_lab=='T'){
  true_organisms=strsplit(true_organism,';',fixed = T)[[1]]
  DE_organisms=strsplit(DE_organism,';',fixed = T)[[1]]
  labs<-get_labels(res.MSstats, true_organisms, DE_organisms, true_lgfc)
  res.MSstats$organism<-labs$Organism
  res.MSstats$label<-labs$DEP
  res.MSstats$TlogFC<-labs$TlogFC
}
write.table(res.MSstats, paste0(save_fold, dataset, '_MSstats_', platform, '_', inten_type, '_', imput, '_', normal, '.csv'), sep = ',', row.names = TRUE, col.names = TRUE)

#write.table(res.MSstats, paste0(save_fold, 'MSstats_', in_type, '_', inten_type, '_', imput, '_', normal, '.csv'), sep = ',', row.names = TRUE, col.names = TRUE)
