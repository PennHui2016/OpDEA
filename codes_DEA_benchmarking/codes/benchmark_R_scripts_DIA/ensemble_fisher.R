### ensemble differential expression analysis
# method1: chisq
# method2: intersection with smallest
# method3: union with smallest
# method4: vote with smallest, more is correct
# method5: intersection with highest
# method6: union with highest
# method7: vote with highest, more is correct

####
## test1: combine 3 best workflows with three different inputs
## test2: combine top K workflows

############ start test1 ##################
#DEA = c('_plgem', '_DEP', '_ProteoMM')
#norm = c('_div.mean', '_', '_max')
#imput = c('_QRILC', '_MinProb', '_MinDet')
#platform = fragpipe
#inputs=c('', 'MaxLFQ', 'MaxLFQ')

library(dplyr)
source('D:/data/benchmark/codes/benchmark_R_scripts_DIA/utils.R')
args <- commandArgs(trailingOnly = TRUE)
#Rscript /home/hui/PycharmProjects/DL_proteomics/benchmark_R_scripts/ensemble_fisher.R /data/res_DE_files/CPTAC_DDA_maxquant_fisher_.csv /data/res_DE_files/CPTAC_DDA/LTQ86/plgem/plgem_maxquant__MinDet_.csv /data/res_DE_files/CPTAC_DDA/LTQ86/DEqMS/DEqMS_maxquant__MLE_center.median.csv /data/res_DE_files/CPTAC_DDA/LTQ86/ProteoMM/ProteoMM_maxquant__bpca_vsn.csv 

# args<-c('/data/res_DE_files/CPTAC_DDA_maxquant_fisher_.csv',
#         '/data/res_DE_files/CPTAC_DDA/LTQ86/plgem/plgem_maxquant__MinDet_.csv',
#         '/data/res_DE_files/CPTAC_DDA/LTQ86/DEqMS/DEqMS_maxquant__MLE_center.median.csv',
#         '/data/res_DE_files/CPTAC_DDA/LTQ86/ProteoMM/ProteoMM_maxquant__bpca_vsn.csv')
N<-length(args)

save_fold<-args[1]
print_lab<-args[2]
true_organism<-args[3]
DE_organism<-args[4]
true_lgfc<-args[5]
#logT<-args[6]
conss<-args[6]
r1<-args[7]

dnum<-function(X){
  n<-nchar(X)
  #if(grepl('_HUMAN', X) | grepl('_YEAST', X) | grepl('_ECOLI', X) | grepl('_UPS', X)){
  if(grepl('_', X)){
    sps<-strsplit(X, '_', fixed = TRUE)
    sps[[1]][length(sps[[1]])] = gsub("\\d","",sps[[1]][length(sps[[1]])])
    outX = sps[[1]][1]
    for (i in 2:length(sps[[1]])) {
      outX = paste0(outX, '_', sps[[1]][i])
    }
  }else{
    outX = X
  }
  #outX = X
  #outX = gsub("\\d","",X)
  # if((substr(X, n, n)!="T") & (substr(X, n, n)!='N') & (substr(X, n, n)!='S') & (substr(X, n, n)!='I')){
  #   outX<-substr(X, 1, n-1)
  # }
  return(outX)
}

res1<-read.table(r1, sep = ',', header = TRUE)#, quote = "")
res1$protein<-apply(as.array(res1$protein), 1, function(x) dnum(x))
base.vars <- c("protein", "contrast")
# Join everything together
b1<-data.frame(cbind(res1$protein, res1$contrast))
colnames(b1)<-c('protein', 'contrast')
res.Base<-b1

ress<-list()

for (i in 8:N) {
  res<-read.table(args[i], sep = ',', header = TRUE)#, quote = "")
  res$protein<-apply(as.array(res$protein), 1, function(x) dnum(x))
  ress[[i]]<-res
  b<-data.frame(cbind(apply(as.array(res$protein), 1, function(x) dnum(x)), res$contrast))
  colnames(b)<-c('protein', 'contrast')
  #base.vars <- c("protein", "UPS", "contrast")
  res.Base<-full_join(res.Base,b,base.vars)
}

###
# plgem pvalue logFC
# ANOVA pvalue logFC
# ProteoMM pvalue logFC
# beta_binomial pvalue logFC
# DEP DEqMS edgeR limma msqrob2 MSStats proDA ROTS SAM siggenes ttest
res11 <- full_join(res.Base, res1)
#z11<- -qnorm(res11$pvalue/2)*sign(res11$logFC)

res.Hurdle <- res.Base
#res.Hurdle$logFC_r1<-res11$logFC
logFCs<-res11$logFC
pvalues<-res11$pvalue
#dfs<-(!is.na(z11))

#all_z<-z11^2
for (i in 8:N) {
  resi <- full_join(res.Base, ress[[i]])
  #zi<- -qnorm(resi$pvalue/2)*sign(resi$logFC)
  #all_z<-cbind(all_z,zi^2)
  logFCs<-cbind(logFCs, resi$logFC)
  pvalues<-cbind(pvalues,resi$pvalue)
  #dfs = dfs + (!is.na(zi))
}

res.Hurdle$pensemble <- apply(pvalues,1, function(x) aggregation::fisher(x))
#res.Hurdle$pchisq <- 1 - pchisq(res.Hurdle$chisq, df = dfs)
#res.Hurdle$logFC_avg<-apply(logFCs,1,function(x) mean(x[which(!is.na(x) & !is.infinite(x))]))
# biggest |logFC|
#res.Hurdle$logFC_avg<-apply(logFCs,1,function(x) sign(x[which(abs(x[which(!is.na(x) & !is.infinite(x))])==max(abs(x[which(!is.na(x) & !is.infinite(x))])))]) * max(abs(x[which(!is.na(x) & !is.infinite(x))])))
res.Hurdle$logFC_avg<-apply(logFCs,1,function(x) sign(x[which(abs(x[which(!is.na(x) & !is.infinite(x))])==max(abs(x[which(!is.na(x) & !is.infinite(x))])))[1]]) * max(abs(x[which(!is.na(x) & !is.infinite(x))]))[1])

res.Hurdle <- res.Hurdle %>% group_by(contrast) %>% mutate(qensemble = pensemble %>%  p.adjust(method = "BH")) %>% ungroup()
res.Hurdle <- res.Hurdle %>% arrange(pensemble)
#res.Hurdle$UPS<-grepl('_HUMAN',res.Hurdle$protein)  
#save_file<-gsub('cv_quant_c.csv', 'cv_quant_h.csv', r1)
#save_file<-paste0(save_fold, 'hurdle.csv')
colnames(res.Hurdle)[c(3,4,5)]<-c('pvalue', 'logFC', 'adj.pvalue')
res_all = res.Hurdle
if(print_lab=='T'){
  true_organisms=strsplit(true_organism,';',fixed = T)[[1]]
  DE_organisms=strsplit(DE_organism,';',fixed = T)[[1]]
  labs<-get_labels(res_all, true_organisms, DE_organisms, true_lgfc)
  res_all$organism<-labs$Organism
  res_all$label<-labs$DEP
  res_all$TlogFC<-labs$TlogFC
}
res_all$logFC[is.na(res_all$logFC)]=0
conss<-strsplit(conss, ';', fixed = T)[[1]]
retain_idx<-c()
for (cons in conss) {
  retain_idx<-c(retain_idx, grep(cons, res_all$contrast))
}
res_all<-res_all[retain_idx,]
write.table(res_all, save_fold, sep = ',', row.names = TRUE, col.names = TRUE)
