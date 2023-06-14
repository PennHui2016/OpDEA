library(tidyverse)
#library(limma)
library(QFeatures)
library(msqrob2)
source('/home/hui/PycharmProjects/DL_proteomics/quant_funcs/preprocessing_peptides.R')

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
#file_name = '/home/hui/Documents/VM_share/DIA_data/yeast/maxquant/txt/peptides.txt'
file_name = peptide_file #'/home/hui/Documents/VM_share/DIA_data/yeast/frager/combined_peptide.tsv'
#design_file = '/home/hui/Documents/VM_share/DIA_data/yeast/design_yeast.txt'
#save_fold = '/home/hui/PycharmProjects/DL_proteomics/benchmark_res/yeast_'
designs = read.table(design_file, quote = "", sep = '\t', header = TRUE)

# in_type = 'maxquant'
# # in_type = 'fragpipe'
# inten_type = ''
# imput = ''
# ## imput = 'bpca', 'knn', 'QRILC', 'MLE', 'MinDet', 'MinProb', 'min', 'zero', 'mixed', 'nbavg', 'with', 'none', 'mice', 'missForest', 'mi'
# normal = "center.median"
# ## normal = "sum", "max", "center.mean", "center.median", "div.mean", "div.median", "diff.meda", "quantiles", "quantiles.robust" or "vsn"

prepro_res = preprocessing_raw_pep(file_name, design_file, in_type, inten_type, imput = imput, normal = normal)

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
  idx_g1 = grep(g1, designs$condition)
  idx_g2 = grep(g2, designs$condition)
  test_col = prepro_res$col[c(idx_g1, idx_g2)]
  pe1 <- readQFeatures(
    table = prepro_res$processed, fnames = 1, ecol = test_col,
    name = "peptideRaw", sep = "\t"
  )
  colData(pe1)$condition <-as.factor(unlist(designs$condition[c(idx_g1, idx_g2)]))
  rowData(pe1[["peptideRaw"]])$nNonZero <- rowSums(assay(pe1[["peptideRaw"]]) > 0)
  pe1 <- zeroIsNA(pe1, "peptideRaw") # convert 0 to NA
  pe1 <- logTransform(pe1, base = 2, i = "peptideRaw", name = "peptideLog")
  pe1 <- filterFeatures(pe1, ~ nNonZero >= 2)
  
  if (normal!=''){
    pe1 <- normalize(pe1,
                     i = "peptideLog",
                     name = "peptideNorm",
                     method = normal
    )
    if(in_type=='maxquant'){
      pe1 <- aggregateFeatures(pe1,
                               i = "peptideNorm", fcol = "Proteins", na.rm = TRUE,
                               name = "protein"
      )
    }else if(in_type=='fragpipe'){
      pe1 <- aggregateFeatures(pe1,
                               i = "peptideNorm", fcol = "Protein", na.rm = TRUE,
                               name = "protein"
      )
    }
  }else{
    if(in_type=='maxquant'){
      pe1 <- aggregateFeatures(pe1,
                               i = "peptideLog", fcol = "Proteins", na.rm = TRUE,
                               name = "protein"
      )
    }else if(in_type=='fragpipe'){
      pe1 <- aggregateFeatures(pe1,
                               i = "peptideLog", fcol = "Protein", na.rm = TRUE,
                               name = "protein"
      )
    }
  }
  
  
  pe <- msqrob(object = pe1, i = "protein", formula = ~condition)
  
  
  
  L <- makeContrast(contrasts=paste0('condition',g2,'=0'),
                     parameterNames = c(paste0('condition',g2)))
  pe <- hypothesisTest(object = pe, i = "protein", contrast = L)
  res.msqrob2 = pe@ExperimentList@listData[["protein"]]@elementMetadata@listData[[paste0('condition',g2)]]
  res.msqrob2$contrast<-consts$conts[k]
  res_all<-rbind(res_all, res.msqrob2)
}
colnames(res_all)[5:6]<-c('pvalue','adj.pvalue')
res_all<-res_all[order(res_all$pvalue),]
res_all$protein<-row.names(res_all)
write.table(res_all, paste0(save_fold, 'msqrob2_', in_type, '_', inten_type, '_', imput, '_', normal, '.csv'), sep = ',', row.names = TRUE, col.names = TRUE)
