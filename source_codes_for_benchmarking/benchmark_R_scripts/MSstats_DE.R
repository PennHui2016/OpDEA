library(MSstats)
library(MSnbase)
library(readr)
library(MSqRob)
library(tidyverse)

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

if(inten_type=='basic' | inten_type=='blank'){
  inten_type=''
}
set.seed(123)
evidence_file_name = msstats_file #'/home/hui/Documents/VM_share/DIA_data/yeast/maxquant/txt/evidence.txt'
protein_file_name = protein_file #'/home/hui/Documents/VM_share/DIA_data/yeast/maxquant/txt/proteinGroups.txt'
#mstats_file_name = '/home/hui/Documents/VM_share/DIA_data/yeast/frager/MSstats.csv'
#design_file = '/home/hui/Documents/VM_share/DIA_data/yeast/design_yeast_msstats.txt'
#save_fold = '/home/hui/PycharmProjects/DL_proteomics/benchmark_res/yeast_'
designs = read.table(design_file, sep = '\t', header = TRUE)

#in_type = 'maxquant'
# in_type = 'fragpipe'
# normal = FALSE # 'equalizeMedians', 'quantile', 'globalStandards', FALSE
if(normal=='FALSE'){
  normal=FALSE
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

if(in_type=='maxquant'){
  # Read in MaxQuant files
  proteinGroups <- read.table(protein_file_name, quote = "", sep = '\t', header = TRUE)
  infile <- read.table(evidence_file_name, quote = "", sep = '\t', header = TRUE)
  #proteinCPTAC<- import2MSnSet(protein_file_name, filetype="MaxQuant", remove_pattern=TRUE)
  annot<-designs
  # runs <- sampleNames(proteinCPTAC)
  # condition <- factor(substr(runs, 1,1))
  # #lab <- factor(rep(rep(1:3,each=3),3))
  # exp_annotation <- data.frame(run=runs, condition=condition)
  # annot <- exp_annotation
  # 
  # annot$Raw.file <- c("110714_yeast_ups1_2fmol_r1", "110714_yeast_ups1_2fmol_r2",
  #                     "110714_yeast_ups1_2fmol_r3", "110714_yeast_ups1_4fmol_r2",
  #                     "110714_yeast_ups1_4fmol_r3", "110714_yeast_ups1_4fmol_r4",
  #                     "110616_yeast_ups_10fmol", "110616_yeast_ups_10fmol_r2",
  #                     "110616_yeast_ups_10fmol_r3", "110618_yeast_ups_25fmol_r1",
  #                     "110618_yeast_ups_25fmol_r2", "110618_yeast_ups_25fmol_r3",
  #                     "110618_yeast_ups_50fmol_r1", "110618_yeast_ups_50fmol_r2",
  #                     "110618_yeast_ups_50fmol_r3")
  # 
  # #annot$run <- NULL
  # annot$IsotopeLabelType <- "L"
  # colnames(annot)[c(1,2)] <- c("BioReplicate", "Condition")
  
  raw <- MaxQtoMSstatsFormat(evidence = infile,
                             annotation = annot,
                             proteinGroups = proteinGroups)
  
  #unique(raw[grepl("ups", raw$ProteinName),]$ProteinName)
  
  QuantData <- dataProcess(raw, censoredInt="NA", normalization=normal)
}else if (in_type=='fragpipe'){
raw <- read_csv(evidence_file_name, na = c("", "NA", "0"))
raw$BioReplicate[is.na(raw$BioReplicate)]=0
raw$ProteinName <- factor(raw$ProteinName)
raw$PeptideSequence <- factor(raw$PeptideSequence)
##
#raw$Condition[which(raw$Condition=="A1" | raw$Condition=="A2" | raw$Condition=="A2")]='A'
for (m in 1:length(raw$ProteinName)) {
  raw$BioReplicate[m] = substr(raw$Condition[m],2,2)
  raw$Condition[m] = substr(raw$Condition[m],1,1)
}
##
QuantData <- dataProcess(raw, censoredInt="NA", normalization=normal)
}
levels(QuantData$ProteinLevelData$GROUP)

# comparison <- t(matrix(c(-1,1,0,0,0,
#                          -1,0,1,0,0,
#                          0,-1,1,0,0,
#                          -1,0,0,1,0,
#                          0,-1,0,1,0,
#                          0,0,-1,1,0,
#                          -1,0,0,0,1,
#                          0,-1,0,0,1,
#                          0,0,-1,0,1,
#                          0,0,0,-1,1
#                          
# ), ncol=10))
# row.names(comparison) <- c("B-A","C-A","C-B",'D-A', 'D-B','D-C', 'E-A', 'E-B', 'E-C', 'E-D')
# colnames(comparison)<-c('A', 'B', 'C', 'D', 'E')

testResultOneComparison <- groupComparison(contrast.matrix = consts$code_cont, data = QuantData)
res.MSstats<-testResultOneComparison$ComparisonResult%>% arrange(pvalue)
colnames(res.MSstats)[1:3]<-c('protein','contrast', 'logFC')
write.table(res.MSstats, paste0(save_fold, 'MSstats_', in_type, '_', inten_type, '_', imput, '_', normal, '.csv'), sep = ',', row.names = TRUE, col.names = TRUE)
