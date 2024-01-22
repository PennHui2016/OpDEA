############################################
###
### Figure 5 ensmeble inference and cross-setting, cross-instrument comparisons
###
########################################

library("readxl")
library(ggplot2)
library(reshape2)
library(ggpubr)
library(cowplot)

data_fold<-paste0(dirname(rstudioapi::getSourceEditorContext()$path),'/data/')
save_fold<-paste0(dirname(rstudioapi::getSourceEditorContext()$path),'/save_data/')

plot_ens<-function(i, metric, st, en, wb){
  platforms<-c('FragPipe', 'Maxquant', 'Maxquant', 'FragPipe','DIANN', 'spt')
  acqs<-c('DDA', 'DDA','TMT', 'TMT', 'DIA', 'DIA')
  settings<-c('FG_DDA','MQ_DDA','MQ_TMT','FG_TMT','DIANN_DIA','spt_DIA')

  platform = platforms[i]
  setting = settings[i]
  acq = acqs[i]

  values=read_excel(paste0(data_fold, 'metrics_', platform, '_',acq,'.xlsx'), sheet = metric)
  ranks_all = read_excel(paste0(data_fold, 'ranks_all_', platform, '_', acq,'.xlsx'), sheet = 'ranking_all')

  if (setting!='MQ_TMT'){
    values_ens_mv = read_excel(paste0(data_fold, 'metrics_', platform, '_',acq,'_all_ensemble_mv.xlsx'), sheet = metric)
    ranks_ens_mv = read_excel(paste0(data_fold, 'ranks_all_', platform, '_', acq,'_ensemble_mv.xlsx'), sheet = 'ranking_all')
    T1_mv=ranks_ens_mv$workflow[which(ranks_ens_mv$avg_rank_mean==min(ranks_ens_mv$avg_rank_mean))]

    if(length(T1_mv)>1){
      T1_mv=T1_mv[1]
    }
  }
  values_ens_topk = read_excel(paste0(data_fold, 'metrics_', platform, '_',acq,'_all_ensemble_topk.xlsx'), sheet = metric)

  ranks_ens_topk = read_excel(paste0(data_fold, 'ranks_all_', platform, '_', acq,'_ensemble_topk.xlsx'), sheet = 'ranking_all')

  T1_wf=ranks_all$workflow[which(ranks_all$avg_rank_mean==min(ranks_all$avg_rank_mean))]
  if(length(T1_wf)>1){
    T1_wf=T1_wf[1]
  }

  T1_topk=ranks_ens_topk$workflow[which(ranks_ens_topk$avg_rank_mean==min(ranks_ens_topk$avg_rank_mean))]
  if(length(T1_topk)>1){
    T1_topk=T1_topk[1]
  }
  T1_wf<-gsub('\\|\\|',paste0('\\|',platform,'\\|'),T1_wf)

  comb_all<-cbind(as.numeric(values[which(values$workflow==T1_wf),][st:en]), rep(metric, en-st+1),
                  rep('TOP1',en-st+1))

  if (setting!='MQ_TMT'){
    comb_all<-rbind(comb_all,cbind(
      as.numeric(values_ens_mv[which(values_ens_mv$workflow==T1_mv),(st+1):(en+1)]), rep(metric, en-st+1),
      rep('ens_multi-quant',en-st+1)
    ))}

  comb_all<-rbind(comb_all,cbind(
    as.numeric(values_ens_topk[which(values_ens_topk$workflow==T1_topk),(st+1):(en+1)]), rep(metric, en-st+1),
    rep('ens_topk',en-st+1)
  ))


  comb_all<-as.data.frame(comb_all)
  colnames(comb_all)<-c('value','metric','method')

  my_comparisons=list(c('TOP1','ens_multi-quant'),
                      c('TOP1','ens_topk'),
                      c('ens_multi-quant', 'ens_topk'))
  melted_data <- comb_all
  melted_data$value<-as.numeric(melted_data$value)

  cols=c("#1f78b4","#b2df8a","#fb9a99","#fdbf6f","#e31a1c","#a6cee3","#33a02c")
  if(metric=='pAUC0.01' | metric=='pAUC0.05' | metric=='pAUC0.1'){
    ylab='pAUC(0.01)'
    ymin=0.5
    ymax=1
  }else if(metric=='G-mean0.05' | metric=='nMCC0.05' | metric=='G-mean0.01' | metric=='nMCC0.01'){
    ylab='G-mean'
    ymin=0
    ymax=1
  }

  if(setting!='MQ_TMT'){
    col=cols[c(1:3)]
  }else{
    col=cols[c(1,3)]
  }

  p1<-ggboxplot(melted_data, "method", "value",width = 0.7, size=0.8,
                color="method", palette = col, title = setting,
                add = "jitter")+
    labs(x = '', y = ylab)+scale_x_discrete(limits=unique(melted_data$method))+
    scale_y_continuous(limits = c(ymin, ymax))+
    theme_classic() +
    #theme(axis.text.x = element_text(size = 10), axis.text.x=element_blank(),axis.text.y = element_text(size = 10))+
    theme(axis.title.y= element_text(size=15))+theme(axis.title.x = element_text(size = 15), axis.text.y = element_text(size = 15))+
    theme(legend.title=element_text(size=14),legend.text=element_text(size=14))+
    #stat_compare_means(method="t.test",hide.ns = F,comparisons = my_comparisons,label="p.signif")+
    theme(legend.position = "right")+stat_summary(fun=mean,
                                                  geom="point",
                                                  shape=17, size=3, color="red", fill="red")+
    theme(plot.title = element_text(size = 16,hjust = 0.5, face = "bold"))+theme(axis.text.x=element_blank())

  p1
  sheet=paste0(setting,'_', metric)

  addWorksheet(wb,sheet)
  writeData(wb, sheet, melted_data, startRow = 1, startCol = 1, colNames = TRUE, rowNames = TRUE)
  saveWorkbook(wb, file = paste0(save_fold, 'metrics_top_ens.xlsx'), overwrite = TRUE)
  return(list(p1=p1,dt=melted_data))
}

library(openxlsx)
wb <- createWorkbook()

p1<-plot_ens(1, 'pAUC0.01', 7, 28, wb)
p2<-plot_ens(2, 'pAUC0.01', 7, 28, wb)
p3<-plot_ens(3, 'pAUC0.01', 7, 21, wb)
p4<-plot_ens(4, 'pAUC0.01', 7, 21, wb)
p5<-plot_ens(5, 'pAUC0.01', 7, 23, wb)
p6<-plot_ens(6, 'pAUC0.01', 7, 23, wb)

p11<-plot_ens(1, 'G-mean0.05', 7, 28, wb)
p21<-plot_ens(2, 'G-mean0.05', 7, 28, wb)
p31<-plot_ens(3, 'G-mean0.05', 7, 21, wb)
p41<-plot_ens(4, 'G-mean0.05', 7, 21, wb)
p51<-plot_ens(5, 'G-mean0.05', 7, 23, wb)
p61<-plot_ens(6, 'G-mean0.05', 7, 23, wb)

p12<-plot_ens(1, 'pAUC0.05', 7, 28, wb)
p22<-plot_ens(2, 'pAUC0.05', 7, 28, wb)
p32<-plot_ens(3, 'pAUC0.05', 7, 21, wb)
p42<-plot_ens(4, 'pAUC0.05', 7, 21, wb)
p52<-plot_ens(5, 'pAUC0.05', 7, 23, wb)
p62<-plot_ens(6, 'pAUC0.05', 7, 23, wb)

p13<-plot_ens(1, 'pAUC0.1', 7, 28, wb)
p23<-plot_ens(2, 'pAUC0.1', 7, 28, wb)
p33<-plot_ens(3, 'pAUC0.1', 7, 21, wb)
p43<-plot_ens(4, 'pAUC0.1', 7, 21, wb)
p53<-plot_ens(5, 'pAUC0.1', 7, 23, wb)
p63<-plot_ens(6, 'pAUC0.1', 7, 23, wb)

p14<-plot_ens(1, 'nMCC0.05', 7, 28, wb)
p24<-plot_ens(2, 'nMCC0.05', 7, 28, wb)
p34<-plot_ens(3, 'nMCC0.05', 7, 21, wb)
p44<-plot_ens(4, 'nMCC0.05', 7, 21, wb)
p54<-plot_ens(5, 'nMCC0.05', 7, 23, wb)
p64<-plot_ens(6, 'nMCC0.05', 7, 23, wb)

extract_legend <- function(my_ggp) {
  step1 <- ggplot_gtable(ggplot_build(my_ggp))
  step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
  step3 <- step1$grobs[[step2]]
  return(step3)
}

shared_legend <- extract_legend(p1$p1)

## Figure 5A

library(ggpubr)
library("gridExtra")

pall<-grid.arrange(p1$p1+theme(legend.position = 'none'),
                   p2$p1+theme(legend.position = 'none', axis.title.y = element_blank(),axis.text.y = element_blank()),
                   p3$p1+theme(legend.position = 'none', axis.title.y = element_blank(),axis.text.y = element_blank()),
                   p4$p1+theme(legend.position = 'none', axis.title.y = element_blank(),axis.text.y = element_blank()),
                   p5$p1+theme(legend.position = 'none', axis.title.y = element_blank(),axis.text.y = element_blank()),
                   p6$p1+theme(legend.position = 'none', axis.title.y = element_blank(),axis.text.y = element_blank()),
                   p11$p1+theme(legend.position = 'none', plot.title = element_blank()),
                   p21$p1+theme(legend.position = 'none', axis.title.y = element_blank(), plot.title = element_blank(),axis.text.y = element_blank()),
                   p31$p1+theme(legend.position = 'none', axis.title.y = element_blank(), plot.title = element_blank(),axis.text.y = element_blank()),
                   p41$p1+theme(legend.position = 'none', axis.title.y = element_blank(), plot.title = element_blank(),axis.text.y = element_blank()),
                   p51$p1+theme(legend.position = 'none', axis.title.y = element_blank(), plot.title = element_blank(),axis.text.y = element_blank()),
                   p61$p1+theme(legend.position = 'none', axis.title.y = element_blank(), plot.title = element_blank(),axis.text.y = element_blank()),
                   nrow = 2, widths=c(1.5,1,1,1,1,1))
grid.arrange(pall,
             shared_legend, nrow = 1, widths = c(10, 1.5))
pall

write.table(p1$dt, paste0(save_fold, 'Figure5A_1.csv'), col.names = T, row.names = T)
write.table(p2$dt, paste0(save_fold, 'Figure5A_2.csv'), col.names = T, row.names = T)
write.table(p3$dt, paste0(save_fold, 'Figure5A_3.csv'), col.names = T, row.names = T)
write.table(p4$dt, paste0(save_fold, 'Figure5A_4.csv'), col.names = T, row.names = T)
write.table(p5$dt, paste0(save_fold, 'Figure5A_5.csv'), col.names = T, row.names = T)
write.table(p6$dt, paste0(save_fold, 'Figure5A_6.csv'), col.names = T, row.names = T)

write.table(p11$dt, paste0(save_fold, 'Figure5A_7.csv'), col.names = T, row.names = T)
write.table(p21$dt, paste0(save_fold, 'Figure5A_8.csv'), col.names = T, row.names = T)
write.table(p31$dt, paste0(save_fold, 'Figure5A_9.csv'), col.names = T, row.names = T)
write.table(p41$dt, paste0(save_fold, 'Figure5A_10.csv'), col.names = T, row.names = T)
write.table(p51$dt, paste0(save_fold, 'Figure5A_11.csv'), col.names = T, row.names = T)
write.table(p61$dt, paste0(save_fold, 'Figure5A_12.csv'), col.names = T, row.names = T)

#######################
## cross_setting compare
library(openxlsx)
wb <- createWorkbook()

data_fold<-paste0(dirname(rstudioapi::getSourceEditorContext()$path),'/data/')
save_fold<-paste0(dirname(rstudioapi::getSourceEditorContext()$path),'/save_data/')

top1wf<-function(dt, platform, acq){
  base_fold<-data_fold
  dea_res_fold<-paste0(base_fold, acq, '/',platform, '/',dt,'/')
  ranking_all_FG_DDA<-read_excel(paste0(data_fold, 'ranks_all_', platform, '_', acq, '.xlsx'), sheet = 'ranking_all')

  #ranking_all<-ranking_all[which(ranking_all$runTime!=0 & ranking_all$mean_pauc001!=0 & ranking_all$rank_median_pauc001!=0),]
  wftop1_FG_DDA<-ranking_all_FG_DDA[which(ranking_all_FG_DDA$avg_rank_mean==min(ranking_all_FG_DDA$avg_rank_mean)),]
  if(length((wftop1_FG_DDA$workflow))>1){
    idx=which(wftop1_FG_DDA$avg_rank_median==min(wftop1_FG_DDA$avg_rank_median))
    wftop1_FG_DDA = wftop1_FG_DDA[idx,]
  }
  dea_res_FG_DDA_wf1<-read.table(paste0(dea_res_fold,
                                        dt,'_',wftop1_FG_DDA$DEA,'_',
                                        wftop1_FG_DDA$Platform,'_',
                                        wftop1_FG_DDA$Matrix,'_',gsub('blank','',wftop1_FG_DDA$Imput),
                                        '_',gsub('blank','',wftop1_FG_DDA$normalization),'.csv'),
                                 sep = ',',header = T)
  return(list(wf=wftop1_FG_DDA, dea_res=dea_res_FG_DDA_wf1))
}

top1wf_ens<-function(dt, platform, acq, ens_ty){
  base_fold<-data_fold
  if(ens_ty=='ensemble_mv'){
    dea_res_fold<-paste0(base_fold, acq, '/',platform, '/','mv','/',dt,'/')
  }else if(ens_ty=='ensemble_topk'){
    dea_res_fold<-paste0(base_fold, acq, '/',platform, '/','topk','/',dt,'/')
  }
  #dea_res_fold<-paste0(base_fold, acq, '/',platform, '/',dt,'/')
  ranking_all_FG_DDA<-read_excel(paste0(data_fold, 'ranks_all_', platform, '_', acq, '_', ens_ty, '.xlsx'), sheet = 'ranking_all')
  #ranking_all<-ranking_all[which(ranking_all$runTime!=0 & ranking_all$mean_pauc001!=0 & ranking_all$rank_median_pauc001!=0),]
  wftop1_FG_DDA<-ranking_all_FG_DDA[which(ranking_all_FG_DDA$avg_rank_mean==min(ranking_all_FG_DDA$avg_rank_mean)),]
  if(length((wftop1_FG_DDA$workflow))>1){
    idx=which(wftop1_FG_DDA$avg_rank_median==min(wftop1_FG_DDA$avg_rank_median))
    wftop1_FG_DDA = wftop1_FG_DDA[idx,]
  }
  if(wftop1_FG_DDA$method=='set'){
    dea_res_FG_DDA_wf1<-read.table(paste0(dea_res_fold,
                                          dt,'_',
                                          wftop1_FG_DDA$Platform,'_',
                                          wftop1_FG_DDA$method, '_',
                                          wftop1_FG_DDA$operation,
                                          '_',gsub('\\|','-',wftop1_FG_DDA$cbn),'.csv'),
                                   sep = ',',header = T)
  }else{
    dea_res_FG_DDA_wf1<-read.table(paste0(dea_res_fold,
                                          dt,'_',
                                          wftop1_FG_DDA$Platform,'_',
                                          wftop1_FG_DDA$method,'_',gsub('\\|','-',wftop1_FG_DDA$cbn),'.csv'),
                                   sep = ',',header = T)
  }

  return(list(wf=wftop1_FG_DDA, dea_res=dea_res_FG_DDA_wf1))
}

dt='HYEtims735'

all_protein_FG_DDA<-read.table(paste0(data_fold,'HYEtims735_LFQ_FragPipe_all_proteins.tsv'), sep = '\t', header = T)
all_protein_MQ_DDA<-read.table(paste0(data_fold,'HYEtims735_LFQ_Maxquant_all_proteins.tsv'), sep = '\t', header = T)
all_protein_DIANN_DIA<-read.table(paste0(data_fold,'HYEtims735_DIA_DIANN_all_proteins.tsv'), sep = '\t', header = T)
all_protein_spt_DIA<-read.table(paste0(data_fold,'HYEtims735_DIA_spt_all_proteins.tsv'), sep = '\t', header = T)

top1_FG_DDA<-top1wf('HYEtims735_LFQ', 'FragPipe', 'DDA')
top1_MQ_DDA<-top1wf('HYEtims735_LFQ', 'Maxquant', 'DDA')
top1_DIANN_DIA<-top1wf('HYEtims735_DIA', 'DIANN', 'DIA')
top1_spt_DIA<-top1wf('HYEtims735_DIA', 'spt', 'DIA')
top1_FG_DDA_ens=top1wf_ens('HYEtims735_LFQ', 'FragPipe', 'DDA', 'ensemble_mv')
top1_MQ_DDA_ens=top1wf_ens('HYEtims735_LFQ', 'Maxquant', 'DDA', 'ensemble_mv')
top1_DIANN_DIA_ens=top1wf_ens('HYEtims735_DIA', 'DIANN', 'DIA', 'ensemble_mv')
top1_spt_DIA_ens=top1wf_ens('HYEtims735_DIA', 'spt', 'DIA', 'ensemble_mv')



# dt='HEqe408'
#
# all_protein_FG_DDA<-read.table('D:/data/benchmark/data/DDA/FragPipe/HEqe408_LFQ_FragPipe_all_proteins.tsv', sep = '\t', header = T)
# all_protein_MQ_DDA<-read.table('D:/data/benchmark/data/DDA/Maxquant/HEqe408_LFQ_Maxquant_all_proteins.tsv', sep = '\t', header = T)
# all_protein_DIANN_DIA<-read.table('D:/data/benchmark/data/DIA/DIANN/HEqe408_DIA_DIANN_all_proteins.tsv', sep = '\t', header = T)
# all_protein_spt_DIA<-read.table('D:/data/benchmark/data/DIA/Spectronaut/HEqe408_DIA_spt_all_proteins.tsv', sep = '\t', header = T)
#
# top1_FG_DDA<-top1wf('HEqe408_LFQ', 'FragPipe', 'DDA')
# top1_MQ_DDA<-top1wf('HEqe408_LFQ', 'Maxquant', 'DDA')
# top1_DIANN_DIA<-top1wf('HEqe408_DIA', 'DIANN', 'DIA')
# top1_spt_DIA<-top1wf('HEqe408_DIA', 'spt', 'DIA')
# top1_FG_DDA_ens=top1wf_ens('HEqe408_LFQ', 'FragPipe', 'DDA', 'ensemble_mv')
# top1_MQ_DDA_ens=top1wf_ens('HEqe408_LFQ', 'Maxquant', 'DDA', 'ensemble_mv')
# top1_DIANN_DIA_ens=top1wf_ens('HEqe408_DIA', 'DIANN', 'DIA', 'ensemble_mv')
# top1_spt_DIA_ens=top1wf_ens('HEqe408_DIA', 'spt', 'DIA', 'ensemble_mv')



FG_DDA<-top1_FG_DDA$dea_res$protein[which(abs(top1_FG_DDA$dea_res$logFC)>=log2(1.5) & top1_FG_DDA$dea_res$adj.pvalue<=0.05)]
MQ_DDA<-top1_MQ_DDA$dea_res$protein[which(abs(top1_MQ_DDA$dea_res$logFC)>=log2(1.5) & top1_MQ_DDA$dea_res$adj.pvalue<=0.05)]
DIANN_DIA<-top1_DIANN_DIA$dea_res$protein[which(abs(top1_DIANN_DIA$dea_res$logFC)>=log2(1.5) & top1_DIANN_DIA$dea_res$adj.pvalue<=0.05)]
spt_DIA<-top1_spt_DIA$dea_res$protein[which(abs(top1_spt_DIA$dea_res$logFC)>=log2(1.5) & top1_spt_DIA$dea_res$adj.pvalue<=0.05)]

FG_DDA_ens<-top1_FG_DDA_ens$dea_res$protein[which(abs(top1_FG_DDA_ens$dea_res$logFC)>=log2(1.5) & top1_FG_DDA_ens$dea_res$adj.pvalue<=0.05)]
MQ_DDA_ens<-top1_MQ_DDA_ens$dea_res$protein[which(abs(top1_MQ_DDA_ens$dea_res$logFC)>=log2(1.5) & top1_MQ_DDA_ens$dea_res$adj.pvalue<=0.05)]
DIANN_DIA_ens<-top1_DIANN_DIA_ens$dea_res$protein[which(abs(top1_DIANN_DIA_ens$dea_res$logFC)>=log2(1.5) & top1_DIANN_DIA_ens$dea_res$adj.pvalue<=0.05)]
spt_DIA_ens<-top1_spt_DIA_ens$dea_res$protein[which(abs(top1_spt_DIA_ens$dea_res$logFC)>=log2(1.5) & top1_spt_DIA_ens$dea_res$adj.pvalue<=0.05)]


T_YEAST<-unique(union(union(union(all_protein_FG_DDA$Protein[which(all_protein_FG_DDA$Organism=='YEAST')],
                                  all_protein_MQ_DDA$Protein[which(all_protein_MQ_DDA$Organism=='YEAST')]),
                            all_protein_DIANN_DIA$Protein[which(all_protein_DIANN_DIA$Organism=='YEAST')]),
                      all_protein_spt_DIA$Protein[which(all_protein_spt_DIA$Organism=='YEAST')]))
T_ECOLI<-unique(union(union(union(all_protein_FG_DDA$Protein[which(all_protein_FG_DDA$Organism=='ECOLI')],
                                  all_protein_MQ_DDA$Protein[which(all_protein_MQ_DDA$Organism=='ECOLI')]),
                            all_protein_DIANN_DIA$Protein[which(all_protein_DIANN_DIA$Organism=='ECOLI')]),
                      all_protein_spt_DIA$Protein[which(all_protein_spt_DIA$Organism=='ECOLI')]))
T_HUMAN<-unique(union(union(union(all_protein_FG_DDA$Protein[which(all_protein_FG_DDA$Organism=='HUMAN')],
                                  all_protein_MQ_DDA$Protein[which(all_protein_MQ_DDA$Organism=='HUMAN')]),
                            all_protein_DIANN_DIA$Protein[which(all_protein_DIANN_DIA$Organism=='HUMAN')]),
                      all_protein_spt_DIA$Protein[which(all_protein_spt_DIA$Organism=='HUMAN')]))

DEP_HUMAN<-unique(union(union(union(top1_FG_DDA$dea_res$protein[which(abs(top1_FG_DDA$dea_res$logFC)>=log2(1.5) &
                                                                        top1_FG_DDA$dea_res$adj.pvalue<=0.05 &
                                                                        top1_FG_DDA$dea_res$organism=='HUMAN')],
                                    top1_MQ_DDA$dea_res$protein[which(abs(top1_MQ_DDA$dea_res$logFC)>=log2(1.5) &
                                                                        top1_MQ_DDA$dea_res$adj.pvalue<=0.05 &
                                                                        top1_MQ_DDA$dea_res$organism=='HUMAN')]),
                              top1_DIANN_DIA$dea_res$protein[which(abs(top1_DIANN_DIA$dea_res$logFC)>=log2(1.5) &
                                                                     top1_DIANN_DIA$dea_res$adj.pvalue<=0.05 &
                                                                     top1_DIANN_DIA$dea_res$organism=='HUMAN')]),
                        top1_spt_DIA$dea_res$protein[which(abs(top1_spt_DIA$dea_res$logFC)>=log2(1.5) &
                                                             top1_spt_DIA$dea_res$adj.pvalue<=0.05 &
                                                             top1_spt_DIA$dea_res$organism=='HUMAN')]))


DEP_HUMAN_ens<-unique(union(union(union(top1_FG_DDA_ens$dea_res$protein[which(abs(top1_FG_DDA_ens$dea_res$logFC)>=log2(1.5) &
                                                                                top1_FG_DDA_ens$dea_res$adj.pvalue<=0.05 &
                                                                                top1_FG_DDA_ens$dea_res$organism=='HUMAN')],
                                        top1_MQ_DDA_ens$dea_res$protein[which(abs(top1_MQ_DDA_ens$dea_res$logFC)>=log2(1.5) &
                                                                                top1_MQ_DDA_ens$dea_res$adj.pvalue<=0.05 &
                                                                                top1_MQ_DDA_ens$dea_res$organism=='HUMAN')]),
                                  top1_DIANN_DIA_ens$dea_res$protein[which(abs(top1_DIANN_DIA_ens$dea_res$logFC)>=log2(1.5) &
                                                                             top1_DIANN_DIA_ens$dea_res$adj.pvalue<=0.05 &
                                                                             top1_DIANN_DIA_ens$dea_res$organism=='HUMAN')]),
                            top1_spt_DIA_ens$dea_res$protein[which(abs(top1_spt_DIA_ens$dea_res$logFC)>=log2(1.5) &
                                                                     top1_spt_DIA_ens$dea_res$adj.pvalue<=0.05 &
                                                                     top1_spt_DIA_ens$dea_res$organism=='HUMAN')]))

DEP_HUMAN<-unique(union(DEP_HUMAN, DEP_HUMAN_ens))

DEPs<-data.frame(Proteins=unique(union(union(union(all_protein_FG_DDA$Protein,
                                                   all_protein_MQ_DDA$Protein),
                                             all_protein_DIANN_DIA$Protein),
                                       all_protein_spt_DIA$Protein)))
DEPs$T_YEAST=c(rep(0, length(DEPs$Proteins)))
DEPs$T_YEAST[match(T_YEAST, DEPs$Proteins)]=1
DEPs$T_ECOLI=c(rep(0, length(DEPs$Proteins)))
DEPs$T_ECOLI[match(T_ECOLI, DEPs$Proteins)]=1
# DEPs$T_HUMAN=c(rep(0, length(DEPs$Proteins)))
# DEPs$T_HUMAN[match(T_HUMAN, DEPs$Proteins)]=1
DEPs$DEP_HUMAN=c(rep(0, length(DEPs$Proteins)))
DEPs$DEP_HUMAN[match(DEP_HUMAN, DEPs$Proteins)]=1
DEPs$FG_DDA=c(rep(0, length(DEPs$Proteins)))
DEPs$FG_DDA[match(FG_DDA, DEPs$Proteins)]=1
DEPs$MQ_DDA=c(rep(0, length(DEPs$Proteins)))
DEPs$MQ_DDA[match(MQ_DDA, DEPs$Proteins)]=1
DEPs$DIANN_DIA=c(rep(0, length(DEPs$Proteins)))
DEPs$DIANN_DIA[match(DIANN_DIA, DEPs$Proteins)]=1
DEPs$spt_DIA=c(rep(0, length(DEPs$Proteins)))
DEPs$spt_DIA[match(spt_DIA, DEPs$Proteins)]=1

DEPs$FG_DDA_ens=c(rep(0, length(DEPs$Proteins)))
DEPs$FG_DDA_ens[match(FG_DDA_ens, DEPs$Proteins)]=1
DEPs$MQ_DDA_ens=c(rep(0, length(DEPs$Proteins)))
DEPs$MQ_DDA_ens[match(MQ_DDA_ens, DEPs$Proteins)]=1
DEPs$DIANN_DIA_ens=c(rep(0, length(DEPs$Proteins)))
DEPs$DIANN_DIA_ens[match(DIANN_DIA_ens, DEPs$Proteins)]=1
DEPs$spt_DIA_ens=c(rep(0, length(DEPs$Proteins)))
DEPs$spt_DIA_ens[match(spt_DIA_ens, DEPs$Proteins)]=1

# simply remove protein groups
DEPs<-DEPs[setdiff(c(1:length(DEPs$Proteins)),grep(';', DEPs$Proteins)),]

sheet=paste(dt, '_cp_DEPs')

addWorksheet(wb,sheet)
writeData(wb, sheet, DEPs, startRow = 1, startCol = 1, colNames = TRUE, rowNames = TRUE)
saveWorkbook(wb, file = paste0(save_fold, 'cross_platform_compare.xlsx'), overwrite = TRUE)

library(UpSetR)
library(UpSetR)
library(grid)
library(gridExtra)
library(ggplot2)

NoAttBasePlot <- function (legend, size_plot_height, Main_bar_plot, Matrix_plot,
                           hratios, Size_plot, query_legend, set_metadata, set_metadata_plots,
                           newpage) {
  top <- 1
  bottom <- 100
  if ((!is.null(legend)) && (query_legend != tolower("none"))) {
    if (query_legend == tolower("top")) {
      top <- 3
      bottom <- 102
      legend_top <- 1
      legend_bottom <- 3
      size_plot_height <- (size_plot_height + 2)
    }
    else if (query_legend == tolower("bottom")) {
      legend_top <- 101
      legend_bottom <- 103
    }
  }
  if (is.null(set_metadata)) {
    matrix_and_mainbar_right <- 100
    matrix_and_mainbar_left <- 21
    size_bar_right <- 20
    size_bar_left <- 1
  }
  else if (!is.null(set_metadata)) {
    matrix_and_mainbar_right <- set_metadata$ncols + 100
    matrix_and_mainbar_left <- set_metadata$ncols + 21
    size_bar_right <- set_metadata$ncols + 20
    size_bar_left <- set_metadata$ncols + 1
    metadata_right <- set_metadata$ncols
    metadata_left <- 1
  }
  if (newpage) {
    grid::grid.newpage()
  }
  if ((!is.null(legend)) && (query_legend != tolower("none"))) {
    if (query_legend == tolower("top")) {
      pushViewport(viewport(layout = grid.layout(102, matrix_and_mainbar_right)))
    }
    else if (query_legend == tolower("bottom")) {
      pushViewport(viewport(layout = grid.layout(103, matrix_and_mainbar_right)))
    }
  }
  else if ((is.null(legend)) || (query_legend == tolower("none"))) {
    pushViewport(viewport(layout = grid.layout(100, matrix_and_mainbar_right)))
  }
  # Modified
  vp = UpSetR:::vplayout(top:bottom, 1:(matrix_and_mainbar_right-matrix_and_mainbar_left))
  pushViewport(vp)
  grid.draw(arrangeGrob(Main_bar_plot, Matrix_plot, heights = hratios))
  popViewport()
  # Modified
  vp = UpSetR:::vplayout(size_plot_height:bottom, (matrix_and_mainbar_right-matrix_and_mainbar_left-1):96)
  pushViewport(vp)
  grid.draw(arrangeGrob(Size_plot))
  popViewport()
  if (!is.null(set_metadata)) {
    for (i in 1:length(set_metadata_plots)) {
      if (i != 1) {
        metadata_left <- 1 + metadata_right
        metadata_right <- metadata_right + set_metadata$plots[[i]]$assign
      }
      else {
        metadata_left <- 1
        metadata_right <- set_metadata$plots[[i]]$assign
      }
      vp = UpSetR:::vplayout(size_plot_height:bottom, metadata_left:metadata_right)
      pushViewport(vp)
      grid.draw(arrangeGrob(set_metadata_plots[[i]]))
      popViewport()
    }
  }
  if ((!is.null(legend)) && (query_legend != tolower("none"))) {
    vp = UpSetR:::vplayout(legend_top:legend_bottom, matrix_and_mainbar_left:matrix_and_mainbar_right)
    pushViewport(vp)
    grid.draw(arrangeGrob(legend))
    popViewport()
  }
}

Make_size_plot <- function (Set_size_data, sbar_color, ratios, ylabel, scale_sets,
                            text_scale, set_size_angle, set_size.show, set_size.scale_max,
                            set_size.number_size) {
  if (length(text_scale) > 1 && length(text_scale) <= 6) {
    x_axis_title_scale <- text_scale[3]
    x_axis_tick_label_scale <- text_scale[4]
  }
  else {
    x_axis_title_scale <- text_scale
    x_axis_tick_label_scale <- text_scale
  }
  if (ylabel == "Set Size" && scale_sets != "identity") {
    ylabel <- paste("Set Size", paste0("( ",
                                       scale_sets, " )"))
    if (scale_sets == "log2") {
      Set_size_data$y <- log2(Set_size_data$y)
    }
    if (scale_sets == "log10") {
      Set_size_data$y <- log10(Set_size_data$y)
    }
  }
  if (!is.null(set_size.number_size)) {
    num.size <- (set_size.number_size/2.845276) * x_axis_tick_label_scale
  }
  else {
    num.size <- (7/2.845276) * x_axis_tick_label_scale
  }
  Size_plot <- (ggplot(data = Set_size_data, aes_string(x = "x",
                                                        y = "y")) + geom_bar(stat = "identity", colour = sbar_color,
                                                                             width = 0.4, fill = sbar_color, position = "identity") +
                  scale_x_continuous(limits = c(0.5, (nrow(Set_size_data) +
                                                        0.5)), breaks = c(0, max(Set_size_data)), expand = c(0,
                                                                                                             0)) + theme(panel.background = element_rect(fill = "white"),
                                                                                                                         plot.margin = unit(c(-0.11, -1.3, 0.5, 0.5), "lines"),
                                                                                                                         axis.title.x = element_text(size = 8.3 * x_axis_title_scale),
                                                                                                                         axis.text.x = element_text(size = 7 * x_axis_tick_label_scale,
                                                                                                                                                    vjust = 1, hjust = 0.5), axis.line = element_line(colour = "gray0"),
                                                                                                                         axis.line.y = element_blank(), axis.line.x = element_line(colour = "gray0",
                                                                                                                                                                                   size = 0.3), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                                                                                                                         panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
                  xlab(NULL) + ylab(ylabel) + coord_flip())
  if (set_size.show == TRUE) {
    Size_plot <- (Size_plot + geom_text(aes(label = y, vjust = 0.5,
                                            hjust = 1.2, angle = set_size_angle), size = num.size))
  }
  if (scale_sets == "log10") {
    if (!is.null(set_size.scale_max)) {
      Size_plot <- (Size_plot + scale_y_continuous(limits = c(set_size.scale_max,
                                                              0), trans = log10_reverse_trans()))
    }
    else {
      Size_plot <- (Size_plot + scale_y_continuous(trans = log10_reverse_trans()))
    }
  }
  else if (scale_sets == "log2") {
    if (!is.null(set_size.scale_max)) {
      Size_plot <- (Size_plot + scale_y_continuous(limits = c(set_size.scale_max,
                                                              0), trans = log2_reverse_trans()))
    }
    else {
      Size_plot <- (Size_plot + scale_y_continuous(trans = log2_reverse_trans()))
    }
  }
  else {
    if (!is.null(set_size.scale_max)) {
      Size_plot <- (Size_plot + scale_y_continuous(limits = c(set_size.scale_max,
                                                              0), trans = "reverse"))
    }
    else {
      # Modified
      #Size_plot <- (Size_plot + scale_y_continuous(trans = "reverse"))
    }
  }
  Size_plot <- ggplot_gtable(ggplot_build(Size_plot))
  return(Size_plot)
}

assignInNamespace(x="NoAttBasePlot", value=NoAttBasePlot, ns="UpSetR")
assignInNamespace(x="Make_size_plot", value=Make_size_plot, ns="UpSetR")

UpSetR::upset(DEPs,
              sets = c("T_YEAST", "T_ECOLI", "DEP_HUMAN", "FG_DDA", "MQ_DDA", "DIANN_DIA","spt_DIA","FG_DDA_ens","MQ_DDA_ens", "DIANN_DIA_ens","spt_DIA_ens"),
              order.by="freq", matrix.color="black", point.size=3,
              sets.bar.color=c("#15d08a","#f2ca01","#a45ee3","#1f11b4","#c3642c",
                               "#f078b4", '#5494cc',"gray","pink",'red','orange'),mb.ratio=c(0.5, 0.5), text.scale = 1.5,
              number.angles =0)
#,mainbar.y.max =900
# if (!require(devtools)) install.packages("devtools")
# devtools::install_github("yanlinlin82/ggvenn")
# library(ComplexHeatmap)
# lists<-list(T_YEAST=DEPs$Proteins[which(DEPs$T_YEAST==1)],
#             T_ECOLI=DEPs$Proteins[which(DEPs$T_ECOLI==1)],
#             DEP_HUMAN=DEPs$Proteins[which(DEPs$DEP_HUMAN==1)],
#             FG_DDA=DEPs$Proteins[which(DEPs$FG_DDA==1)],
#             MQ_DDA=DEPs$Proteins[which(DEPs$MQ_DDA==1)],
#             DIANN_DIA=DEPs$Proteins[which(DEPs$DIANN_DIA==1)],
#             spt_DIA=DEPs$Proteins[which(DEPs$spt_DIA==1)],
#             FG_DDA_ens=DEPs$Proteins[which(DEPs$FG_DDA_ens==1)],
#             MQ_DDA_ens=DEPs$Proteins[which(DEPs$MQ_DDA_ens==1)],
#             DIANN_DIA_ens=DEPs$Proteins[which(DEPs$DIANN_DIA_ens==1)],
#             spt_DIA_ens=DEPs$Proteins[which(DEPs$spt_DIA_ens==1)])
#
# m = make_comb_mat(lists, mode = "distinct")
# UpSet(m)

library(pROC)
library(ggpubr)

#all_pvals<-read.table('E:/proteomics/manu4_1/codes/reproduce_0_05/temp0_limma_min_FragPipe/HYEtims735_LFQ_top0_maxlfq_mv_outputs.csv',header = T, sep = ',')

cross_platform_proteins<-data.frame(Proteins=DEPs$Proteins, labels=0)
cross_platform_proteins$labels[grep('YEAST', cross_platform_proteins$Proteins)]=1
cross_platform_proteins$labels[grep('ECOLI', cross_platform_proteins$Proteins)]=1
cross_platform_proteins$qval_FG_DDA<-1
cross_platform_proteins$logFC_FG_DDA<-0
com=intersect(cross_platform_proteins$Proteins, top1_FG_DDA$dea_res$protein)
idx1=match(com, cross_platform_proteins$Proteins)
idx2=match(com, top1_FG_DDA$dea_res$protein)
cross_platform_proteins$qval_FG_DDA[idx1]=top1_FG_DDA$dea_res$adj.pvalue[idx2]
cross_platform_proteins$logFC_FG_DDA[idx1]=top1_FG_DDA$dea_res$logFC[idx2]

cross_platform_proteins$qval_MQ_DDA<-1
cross_platform_proteins$logFC_MQ_DDA<-0
com=intersect(cross_platform_proteins$Proteins, top1_MQ_DDA$dea_res$protein)
idx1=match(com, cross_platform_proteins$Proteins)
idx2=match(com, top1_MQ_DDA$dea_res$protein)
cross_platform_proteins$qval_MQ_DDA[idx1]=top1_MQ_DDA$dea_res$adj.pvalue[idx2]
cross_platform_proteins$logFC_MQ_DDA[idx1]=top1_MQ_DDA$dea_res$logFC[idx2]

cross_platform_proteins$qval_DIANN_DIA<-1
cross_platform_proteins$logFC_DIANN_DIA<-0
com=intersect(cross_platform_proteins$Proteins, top1_DIANN_DIA$dea_res$protein)
idx1=match(com, cross_platform_proteins$Proteins)
idx2=match(com, top1_DIANN_DIA$dea_res$protein)
cross_platform_proteins$qval_DIANN_DIA[idx1]=top1_DIANN_DIA$dea_res$adj.pvalue[idx2]
cross_platform_proteins$logFC_DIANN_DIA[idx1]=top1_DIANN_DIA$dea_res$logFC[idx2]

cross_platform_proteins$qval_spt_DIA<-1
cross_platform_proteins$logFC_spt_DIA<-0
com=intersect(cross_platform_proteins$Proteins, top1_spt_DIA$dea_res$protein)
idx1=match(com, cross_platform_proteins$Proteins)
idx2=match(com, top1_spt_DIA$dea_res$protein)
cross_platform_proteins$qval_spt_DIA[idx1]=top1_spt_DIA$dea_res$adj.pvalue[idx2]
cross_platform_proteins$logFC_spt_DIA[idx1]=top1_spt_DIA$dea_res$logFC[idx2]

cross_platform_proteins$qval_FG_DDA_ens<-1
cross_platform_proteins$logFC_FG_DDA_ens<-0
com=intersect(cross_platform_proteins$Proteins, top1_FG_DDA_ens$dea_res$protein)
idx1=match(com, cross_platform_proteins$Proteins)
idx2=match(com, top1_FG_DDA_ens$dea_res$protein)
cross_platform_proteins$qval_FG_DDA_ens[idx1]=top1_FG_DDA_ens$dea_res$adj.pvalue[idx2]
cross_platform_proteins$logFC_FG_DDA_ens[idx1]=top1_FG_DDA_ens$dea_res$logFC[idx2]

cross_platform_proteins$qval_MQ_DDA_ens<-1
cross_platform_proteins$logFC_MQ_DDA_ens<-0
com=intersect(cross_platform_proteins$Proteins, top1_MQ_DDA_ens$dea_res$protein)
idx1=match(com, cross_platform_proteins$Proteins)
idx2=match(com, top1_MQ_DDA_ens$dea_res$protein)
cross_platform_proteins$qval_MQ_DDA_ens[idx1]=top1_MQ_DDA_ens$dea_res$adj.pvalue[idx2]
cross_platform_proteins$logFC_MQ_DDA_ens[idx1]=top1_MQ_DDA_ens$dea_res$logFC[idx2]

cross_platform_proteins$qval_DIANN_DIA_ens<-1
cross_platform_proteins$logFC_DIANN_DIA_ens<-0
com=intersect(cross_platform_proteins$Proteins, top1_DIANN_DIA_ens$dea_res$protein)
idx1=match(com, cross_platform_proteins$Proteins)
idx2=match(com, top1_DIANN_DIA_ens$dea_res$protein)
cross_platform_proteins$qval_DIANN_DIA_ens[idx1]=top1_DIANN_DIA_ens$dea_res$adj.pvalue[idx2]
cross_platform_proteins$logFC_DIANN_DIA_ens[idx1]=top1_DIANN_DIA_ens$dea_res$logFC[idx2]

cross_platform_proteins$qval_spt_DIA_ens<-1
cross_platform_proteins$logFC_spt_DIA_ens<-0
com=intersect(cross_platform_proteins$Proteins, top1_spt_DIA_ens$dea_res$protein)
idx1=match(com, cross_platform_proteins$Proteins)
idx2=match(com, top1_spt_DIA_ens$dea_res$protein)
cross_platform_proteins$qval_spt_DIA_ens[idx1]=top1_spt_DIA_ens$dea_res$adj.pvalue[idx2]
cross_platform_proteins$logFC_spt_DIA_ens[idx1]=top1_spt_DIA_ens$dea_res$logFC[idx2]

sheet=paste(dt, '_cp_proteins')

addWorksheet(wb,sheet)
writeData(wb, sheet, cross_platform_proteins, startRow = 1, startCol = 1, colNames = TRUE, rowNames = TRUE)
saveWorkbook(wb, file = paste0(save_fold, 'cross_platform_compare.xlsx'), overwrite = TRUE)

roc_FG_DDA <- roc(cross_platform_proteins$labels, 1-cross_platform_proteins$qval_FG_DDA)
roc_MQ_DDA <- roc(cross_platform_proteins$labels, 1-cross_platform_proteins$qval_MQ_DDA)
roc_DIANN_DIA <- roc(cross_platform_proteins$labels, 1-cross_platform_proteins$qval_DIANN_DIA)
roc_spt_DIA <- roc(cross_platform_proteins$labels, 1-cross_platform_proteins$qval_spt_DIA)

roc_FG_DDA_ens <- roc(cross_platform_proteins$labels, 1-cross_platform_proteins$qval_FG_DDA_ens)
roc_MQ_DDA_ens <- roc(cross_platform_proteins$labels, 1-cross_platform_proteins$qval_MQ_DDA_ens)
roc_DIANN_DIA_ens <- roc(cross_platform_proteins$labels, 1-cross_platform_proteins$qval_DIANN_DIA_ens)
roc_spt_DIA_ens <- roc(cross_platform_proteins$labels, 1-cross_platform_proteins$qval_spt_DIA_ens)

# ms_bench_auc <- function(FPR, TPR, fpr_threshold = 1) {
#   # make sure that sorted.
#   oFPR <- order(FPR)
#   FPR <- FPR[oFPR]
#   TPR <- TPR[oFPR]
#
#   idx <- FPR < fpr_threshold
#   TPR <- TPR[idx]
#   FPR <- FPR[idx]
#   #integrate
#   res <- 1 / 2 * sum(diff(FPR) * (head(TPR,-1) + tail(TPR,-1)))
#   return(res / fpr_threshold * 100)
# }
# ms_bench_auc(1-roc_FG_DDA$specificities,roc_FG_DDA$sensitivities, fpr_threshold = 0.01)
# auc(roc_FG_DDA, partial.auc=c(1, 0.99), partial.auc.correct=TRUE)
spe_sens<-rbind(cbind(roc_FG_DDA$specificities, roc_FG_DDA$sensitivities,
                      c(rep(paste0('FG_DDA(',as.character(round(auc(roc_FG_DDA,
                                                                    partial.auc=c(1, 0.95),
                                                                    partial.auc.correct=TRUE),4)),')'),
                            length(roc_FG_DDA$specificities)))),
                cbind(roc_MQ_DDA$specificities, roc_MQ_DDA$sensitivities,
                      c(rep(paste0('MQ_DDA(',as.character(round(auc(roc_MQ_DDA,
                                                                    partial.auc=c(1, 0.95),
                                                                    partial.auc.correct=TRUE),4)),')'),
                            length(roc_MQ_DDA$specificities)))),
                cbind(roc_DIANN_DIA$specificities, roc_DIANN_DIA$sensitivities,
                      c(rep(paste0('DIANN_DIA(',as.character(round(auc(roc_DIANN_DIA,
                                                                       partial.auc=c(1, 0.95),
                                                                       partial.auc.correct=TRUE),4)),')'),
                            length(roc_DIANN_DIA$specificities)))),
                cbind(roc_spt_DIA$specificities, roc_spt_DIA$sensitivities,
                      c(rep(paste0('spt_DIA(',as.character(round(auc(roc_spt_DIA,
                                                                     partial.auc=c(1, 0.95),
                                                                     partial.auc.correct=TRUE),4)),')'),
                            length(roc_spt_DIA$specificities)))),
                cbind(roc_FG_DDA_ens$specificities, roc_FG_DDA_ens$sensitivities,
                      c(rep(paste0('FG_DDA_ens(',as.character(round(auc(roc_FG_DDA_ens,
                                                                        partial.auc=c(1, 0.95),
                                                                        partial.auc.correct=TRUE),4)),')'),
                            length(roc_FG_DDA_ens$specificities)))),
                cbind(roc_MQ_DDA_ens$specificities, roc_MQ_DDA_ens$sensitivities,
                      c(rep(paste0('MQ_DDA_ens(',as.character(round(auc(roc_MQ_DDA_ens,
                                                                        partial.auc=c(1, 0.95),
                                                                        partial.auc.correct=TRUE),4)),')'),
                            length(roc_MQ_DDA_ens$specificities)))),
                cbind(roc_DIANN_DIA_ens$specificities, roc_DIANN_DIA_ens$sensitivities,
                      c(rep(paste0('DIANN_DIA_ens(',as.character(round(auc(roc_DIANN_DIA_ens,
                                                                           partial.auc=c(1, 0.95),
                                                                           partial.auc.correct=TRUE),4)),')'),
                            length(roc_DIANN_DIA_ens$specificities)))),
                cbind(roc_spt_DIA_ens$specificities, roc_spt_DIA_ens$sensitivities,
                      c(rep(paste0('spt_DIA_ens(',as.character(round(auc(roc_spt_DIA_ens,
                                                                         partial.auc=c(1, 0.95),
                                                                         partial.auc.correct=TRUE),4)),')'),
                            length(roc_spt_DIA_ens$specificities))))
)

colnames(spe_sens)<-c('Specificity', 'Sensitivity', 'Method(pAUC(FPR=0.05))')
spe_sens<-as.data.frame(spe_sens)
spe_sens$Specificity<-1-as.numeric(spe_sens$Specificity)
spe_sens$Sensitivity<-as.numeric(spe_sens$Sensitivity)
#roc.test(roc_v1, roc_v2)
breaks = seq(0,1,0.1)
ggplot(spe_sens, aes(x = Specificity, y = Sensitivity, colour = `Method(pAUC(FPR=0.05))`))+
  scale_color_manual(values=c("#15d08a","#f2ca01","#a45ee3","#1f11b4","#c3642c",
                              "#f078b4", '#5494cc',"gray","pink",'red','orange')) +
  #geom_segment(aes(x = 0, y = 1, xend = 1,yend = 0),linewidth=10, alpha = 0.5, colour = "gray") +
  geom_step(linewidth=1,alpha = 0.9) +
  scale_x_continuous(name = "False Positive Rate (1 - Specificity)",limits = c(0,0.1), breaks = breaks) +
  scale_y_continuous(name = "True Positive Rate (Sensitivity)", limits = c(0,0.7), breaks = breaks) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 1))+scale_x_continuous(limits = c(0, 1))+
  labs(
    x = "False Positive Rate (1-Specificity)",
    y = "True Positive Rate (Sensitivity)")+ theme(legend.position = c(0.7, 0.35)) +
  geom_vline(xintercept = 0.05,linetype="dashed")
#annotate("text", x = 0.1, y = 0.1, vjust = 0, label = paste("pAUC(0.05) =",sprintf("%.3f",aucAvg))) +
#guides(colour = guide_legend(legendTitel)) +
#theme(axis.ticks = element_line(color = "grey80"))

cross_platform_metrics<-function(cross_platform_proteins, setting, roc){
  lab=cross_platform_proteins$labels
  if(setting=='FG_DDA'){
    logFCs=cross_platform_proteins$logFC_FG_DDA
    qvals=cross_platform_proteins$qval_FG_DDA
    roc_s<-roc_FG_DDA
  }else if(setting=='MQ_DDA'){
    logFCs=cross_platform_proteins$logFC_MQ_DDA
    qvals=cross_platform_proteins$qval_MQ_DDA
    roc_s<-roc_MQ_DDA
  }else if(setting=='DIANN_DIA'){
    logFCs=cross_platform_proteins$logFC_DIANN_DIA
    qvals=cross_platform_proteins$qval_DIANN_DIA
    roc_s<-roc_DIANN_DIA
  }else if(setting=='spt_DIA'){
    logFCs=cross_platform_proteins$logFC_spt_DIA
    qvals=cross_platform_proteins$qval_spt_DIA
    roc_s<-roc_spt_DIA
  }else if(setting=='FG_DDA_ens'){
    logFCs=cross_platform_proteins$logFC_FG_DDA_ens
    qvals=cross_platform_proteins$qval_FG_DDA_ens
    roc_s<-roc_FG_DDA_ens
  }else if(setting=='MQ_DDA_ens'){
    logFCs=cross_platform_proteins$logFC_MQ_DDA_ens
    qvals=cross_platform_proteins$qval_MQ_DDA_ens
    roc_s<-roc_MQ_DDA_ens
  }else if(setting=='DIANN_DIA_ens'){
    logFCs=cross_platform_proteins$logFC_DIANN_DIA_ens
    qvals=cross_platform_proteins$qval_DIANN_DIA_ens
    roc_s<-roc_DIANN_DIA_ens
  }else if(setting=='spt_DIA_ens'){
    logFCs=cross_platform_proteins$logFC_spt_DIA_ens
    qvals=cross_platform_proteins$qval_spt_DIA_ens
    roc_s<-roc_spt_DIA_ens
  }
  TP=length(which(lab==1 & abs(logFCs)>=log2(1.5) & qvals<=0.05))
  TN=length(which(lab==0 & (abs(logFCs)<=log2(1.5) | qvals>0.05)))
  FP=length(which(lab==0 & abs(logFCs)>=log2(1.5) & qvals<=0.05))
  FN=length(which(lab==1 & (abs(logFCs)<=log2(1.5) | qvals>0.05)))

  Recall = TP/(TP+FN)
  Precision = TP/(TP+FP)
  Specificity = TN/(TN+FP)
  F1=2*Recall*Precision/(Recall+Precision)
  MCC = (TP*TN-FP*FN)/(((TP+FP)^0.5)*((TP+FN)^0.5)*((TN+FP)^0.5)*((TN+FN))^0.5)
  nMCC = (1+MCC)/2
  Gmean = (Recall*Specificity)^0.5

  pauc001=as.numeric(auc(roc_s, partial.auc=c(1, 0.99), partial.auc.correct=TRUE))
  pauc005=as.numeric(auc(roc_s, partial.auc=c(1, 0.95), partial.auc.correct=TRUE))
  pauc01=as.numeric(auc(roc_s, partial.auc=c(1, 0.90), partial.auc.correct=TRUE))
  metrics=data.frame(metric=c('pAUC(0.01)','pAUC(0.05)','pAUC(0.1)','TP', 'TN', 'FP', 'FN', 'Recall', 'Precision', 'Specificity', 'F1', 'MCC', 'G-mean', 'nMCC'),
                     value=c(pauc001,pauc005,pauc01,TP, TN, FP, FN, Recall, Precision, Specificity, F1, MCC, Gmean,nMCC))
  return(metrics)
}

metrics_FG_DDA=cross_platform_metrics(cross_platform_proteins, 'FG_DDA', roc_FG_DDA)
metrics_MQ_DDA=cross_platform_metrics(cross_platform_proteins, 'MQ_DDA', roc_MQ_DDA)
metrics_DIANN_DIA=cross_platform_metrics(cross_platform_proteins, 'DIANN_DIA', roc_DIANN_DIA)
metrics_spt_DIA=cross_platform_metrics(cross_platform_proteins, 'spt_DIA', roc_spt_DIA)

metrics_FG_DDA_ens=cross_platform_metrics(cross_platform_proteins, 'FG_DDA_ens', roc_FG_DDA_ens)
metrics_MQ_DDA_ens=cross_platform_metrics(cross_platform_proteins, 'MQ_DDA_ens', roc_MQ_DDA_ens)
metrics_DIANN_DIA_ens=cross_platform_metrics(cross_platform_proteins, 'DIANN_DIA_ens', roc_DIANN_DIA_ens)
metrics_spt_DIA_ens=cross_platform_metrics(cross_platform_proteins, 'spt_DIA_ens', roc_spt_DIA_ens)

rdar_dt<-rbind(metrics_FG_DDA$value, metrics_MQ_DDA$value, metrics_DIANN_DIA$value, metrics_spt_DIA$value,
               metrics_FG_DDA_ens$value, metrics_MQ_DDA_ens$value, metrics_DIANN_DIA_ens$value, metrics_spt_DIA_ens$value)
colnames(rdar_dt)<-metrics_FG_DDA$metric
row.names(rdar_dt)<-c('FG_DDA','MQ_DDA', 'DIANN_DIA', 'spt_DIA', 'FG_DDA_ens','MQ_DDA_ens', 'DIANN_DIA_ens', 'spt_DIA_ens')
rdar_dt<-as.data.frame(as.matrix(rdar_dt))
pt_rdar_dt<-data.frame(group=row.names(rdar_dt))

sheet=paste(dt, '_cp_metrics')

addWorksheet(wb,sheet)
writeData(wb, sheet, rdar_dt, startRow = 1, startCol = 1, colNames = TRUE, rowNames = TRUE)
saveWorkbook(wb, file = paste0(save_fold, 'cross_platform_compare.xlsx'), overwrite = TRUE)
pt_rdar_dt<-cbind(pt_rdar_dt, round(rdar_dt[,c(1,2,3,13,14)], 2)) #,4,6
#pt_rdar_dt$group<-row.names(pt_rdar_dt)
library(ggradar)
ggradar(pt_rdar_dt,group.colours=c("#a45ee3","#1f11b4","#f078b4","#f2ca01","gray","pink",'red',"orange"),legend.position = 'right')
#a=c('#5494cc','#e18283','#0d898a','#f9cc52')
col=c("#15d08a","#f2ca01","#a45ee3","#1f11b4","#c3642c",
      "#f078b4", '#5494cc',"gray","pink",'red','orange')
melt_dt<-melt(pt_rdar_dt)
melt_dt$value<-round(melt_dt$value,2)
colnames(melt_dt)[1]<-'method'
# p3=ggplot(melt_dt, aes(x = variable, y = value, colour = method)) +
#   geom_bar(aes(fill = method),position=position_dodge(width = 0.5), stat = "identity", width = 0.1) +
#   geom_point(position=position_dodge(width = 0.1),size=3, shape=19) +
#   scale_fill_manual(values=col[])
#   mytheme = theme_classic() + theme(axis.text.x = element_text(size = 16, angle = -20),axis.text.y = element_text(size = 16))+
#   theme(axis.title.y= element_text(size=16))+theme(axis.title.x = element_text(size = 16))+
#   theme(legend.title=element_text(size=16),legend.text=element_text(size=16), legend.position = "top")
#   p3=p3+mytheme
#   p3

ggdotchart(melt_dt, "variable", "value", group = "method", color = "method",
           palette = col[c(6,1,8,3,10,7,5,4)],

           add = "segments", #label = "value",
           add.params = list(color = "lightgray", size = 1),
           dot.size = 3,
           #label = 'value',
           #font.label = list(color = "black", size = 9, vjust = 0.5),
           position=position_dodge(width = 0.9)
) +
  #geom_text(aes(label=value), position = 'indentity')+
  theme(axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
        axis.title.x = element_blank())#+#+scale_y_continuous(limits = c(0.2, 1))+
#coord_flip()


count_TPs<-data.frame(method=row.names(rdar_dt))
count_TPs<-cbind(count_TPs, rdar_dt[,c(4:7)]) #,4,6
melt_dt<-melt(count_TPs)
melt_dt$label<-as.character(melt_dt$value)
melt_dt$label[which(melt_dt$variable=='TN' | melt_dt$variable=='FN')]=''

nums<-vector()


for (i in 1:length(unique(melt_dt$method))){
  num_i<-melt_dt[which(melt_dt$method==melt_dt$method[i]),]
  uni_var<-levels(factor(melt_dt$variable))
  re_order<-match(num_i$variable, uni_var)
  num_i<-num_i[re_order,]
  l_y<-c(rep(0,length(uni_var)))
  for (j in length(uni_var):1) {
    if(j==length(uni_var)){
      l_y[j]=num_i$value[j]/2
    }else{
      l_y_j=num_i$value[j]/2 + sum(num_i$value[c((j+1):length(uni_var))])
      l_y[j]=l_y_j
    }
  }
  num_i$l_y<-l_y
  nums<-rbind(nums,num_i)
}

#Figure 5B

ggbarplot(nums, "method", "value", group = "variable", color = "variable",
          fill="variable", #label = T,lab.pos='in',
          platette = c('#5494cc','#e18283','#0d898a','#f9cc52'),
) + geom_text(aes(y=l_y, label=label))+
  scale_x_discrete(limits=c('FG_DDA_ens','MQ_DDA_ens','DIANN_DIA_ens', 'spt_DIA_ens',
                            'FG_DDA','MQ_DDA','DIANN_DIA','spt_DIA')) +
  geom_vline(xintercept = 4.5, linetype=2, show.legend=FALSE) + theme(axis.line.y = element_blank())+
  coord_flip()+theme(axis.title.y = element_blank())

write.table(nums, paste0(save_fold, 'Figure5B.csv'), col.names = T, row.names = F)
#############  cross_platform TMT_DDA

all_protein_FG_DDA<-read.table(paste0(data_fold, 'HYqfl683_LFQ_FragPipe_all_proteins.tsv'), sep = '\t', header = T)
all_protein_MQ_DDA<-read.table(paste0(data_fold, 'HYqfl683_LFQ_Maxquant_all_proteins.tsv'), sep = '\t', header = T)
all_protein_FG_TMT<-read.table(paste0(data_fold, 'HYqfl683_TMT11_FragPipe_tmt_all_proteins.tsv'), sep = '\t', header = T)
all_protein_MQ_TMT<-read.table(paste0(data_fold, 'HYqfl683_TMT11_Maxquant_tmt_all_proteins.tsv'), sep = '\t', header = T)


top1wf<-function(dt, platform, acq, contrast){
  base_fold<-data_fold
  dea_res_fold<-paste0(base_fold, acq, '/',platform, '/',dt,'/')
  ranking_all_FG_DDA<-read_excel(paste0(data_fold, 'ranks_all_', platform, '_', acq, '.xlsx'), sheet = 'ranking_all')
  #ranking_all<-ranking_all[which(ranking_all$runTime!=0 & ranking_all$mean_pauc001!=0 & ranking_all$rank_median_pauc001!=0),]
  wftop1_FG_DDA<-ranking_all_FG_DDA[which(ranking_all_FG_DDA$avg_rank_mean==min(ranking_all_FG_DDA$avg_rank_mean)),]
  if(length((wftop1_FG_DDA$workflow))>1){
    idx=which(wftop1_FG_DDA$avg_rank_median==min(wftop1_FG_DDA$avg_rank_median))
    wftop1_FG_DDA = wftop1_FG_DDA[idx,]
  }
  dea_res_FG_DDA_wf1<-read.table(paste0(dea_res_fold,
                                        dt,'_',wftop1_FG_DDA$DEA,'_',
                                        wftop1_FG_DDA$Platform,'_',
                                        wftop1_FG_DDA$Matrix,'_',gsub('blank','',wftop1_FG_DDA$Imput),
                                        '_',gsub('blank','',wftop1_FG_DDA$normalization),'.csv'),
                                 sep = ',',header = T)
  dea_res_FG_DDA_wf1<-dea_res_FG_DDA_wf1[which(dea_res_FG_DDA_wf1$contrast==contrast),]
  return(list(wf=wftop1_FG_DDA, dea_res=dea_res_FG_DDA_wf1))
}

top1wf_ens<-function(dt, platform, acq, ens_ty, contrast){
  base_fold<-data_fold
  if(ens_ty=='ensemble_mv'){
    dea_res_fold<-paste0(base_fold, acq, '/',platform, '/','mv','/',dt,'/')
  }else if(ens_ty=='ensemble_topk'){
    dea_res_fold<-paste0(base_fold, acq, '/',platform, '/','topk','/',dt,'/')
  }
  #dea_res_fold<-paste0(base_fold, acq, '/',platform, '/',dt,'/')
  ranking_all_FG_DDA<-read_excel(paste0(data_fold, 'ranks_all_', platform, '_', acq, '_', ens_ty, '.xlsx'), sheet = 'ranking_all')
  #ranking_all<-ranking_all[which(ranking_all$runTime!=0 & ranking_all$mean_pauc001!=0 & ranking_all$rank_median_pauc001!=0),]
  wftop1_FG_DDA<-ranking_all_FG_DDA[which(ranking_all_FG_DDA$avg_rank_mean==min(ranking_all_FG_DDA$avg_rank_mean)),]
  if(length((wftop1_FG_DDA$workflow))>1){
    idx=which(wftop1_FG_DDA$avg_rank_median==min(wftop1_FG_DDA$avg_rank_median))
    wftop1_FG_DDA = wftop1_FG_DDA[idx,]
  }
  if(wftop1_FG_DDA$method=='set'){
    dea_res_FG_DDA_wf1<-read.table(paste0(dea_res_fold,
                                          dt,'_',
                                          wftop1_FG_DDA$Platform,'_',
                                          wftop1_FG_DDA$method, '_',
                                          wftop1_FG_DDA$operation,
                                          '_',gsub('top','',gsub('\\|','-',wftop1_FG_DDA$cbn)),'.csv'),
                                   sep = ',',header = T)
  }else{
    if(acq=='TMT'){
      dea_res_FG_DDA_wf1<-read.table(paste0(dea_res_fold,
                                            dt,'_',
                                            wftop1_FG_DDA$Platform,'_',
                                            wftop1_FG_DDA$method,'_',gsub("top",'',gsub('\\|','-',wftop1_FG_DDA$cbn)),'.csv'),
                                     sep = ',',header = T)
    }else{
      dea_res_FG_DDA_wf1<-read.table(paste0(dea_res_fold,
                                            dt,'_',
                                            wftop1_FG_DDA$Platform,'_',
                                            wftop1_FG_DDA$method,'_',gsub('\\|','-',wftop1_FG_DDA$cbn),'.csv'),
                                     sep = ',',header = T)
    }

  }
  dea_res_FG_DDA_wf1<-dea_res_FG_DDA_wf1[which(dea_res_FG_DDA_wf1$contrast==contrast),]
  return(list(wf=wftop1_FG_DDA, dea_res=dea_res_FG_DDA_wf1))
}
dt = 'HYqfl683'
contrast='conditionB-conditionA'
top1_FG_DDA<-top1wf('HYqfl683_LFQ', 'FragPipe', 'DDA',contrast)
top1_MQ_DDA<-top1wf('HYqfl683_LFQ', 'Maxquant', 'DDA',contrast)
top1_FG_TMT<-top1wf('HYqfl683_TMT11', 'FragPipe', 'TMT',contrast)
top1_MQ_TMT<-top1wf('HYqfl683_TMT11', 'Maxquant', 'TMT',contrast)

top1_FG_DDA_ens<-top1wf_ens('HYqfl683_LFQ', 'FragPipe', 'DDA', 'ensemble_mv',contrast)
top1_MQ_DDA_ens<-top1wf_ens('HYqfl683_LFQ', 'Maxquant', 'DDA', 'ensemble_mv',contrast)
top1_FG_TMT_ens<-top1wf_ens('HYqfl683_TMT11', 'FragPipe', 'TMT', 'ensemble_mv',contrast)
top1_MQ_TMT_ens<-top1wf_ens('HYqfl683_TMT11', 'Maxquant', 'TMT', 'ensemble_topk',contrast)

FG_DDA<-top1_FG_DDA$dea_res$protein[which(abs(top1_FG_DDA$dea_res$logFC)>=log2(1.5) & top1_FG_DDA$dea_res$adj.pvalue<=0.05)]
MQ_DDA<-top1_MQ_DDA$dea_res$protein[which(abs(top1_MQ_DDA$dea_res$logFC)>=log2(1.5) & top1_MQ_DDA$dea_res$adj.pvalue<=0.05)]
FG_TMT<-top1_FG_TMT$dea_res$protein[which(abs(top1_FG_TMT$dea_res$logFC)>=log2(1.5) & top1_FG_TMT$dea_res$adj.pvalue<=0.05)]
MQ_TMT<-top1_MQ_TMT$dea_res$protein[which(abs(top1_MQ_TMT$dea_res$logFC)>=log2(1.5) & top1_MQ_TMT$dea_res$adj.pvalue<=0.05)]

FG_DDA_ens<-top1_FG_DDA_ens$dea_res$protein[which(abs(top1_FG_DDA_ens$dea_res$logFC)>=log2(1.5) & top1_FG_DDA_ens$dea_res$adj.pvalue<=0.05)]
MQ_DDA_ens<-top1_MQ_DDA_ens$dea_res$protein[which(abs(top1_MQ_DDA_ens$dea_res$logFC)>=log2(1.5) & top1_MQ_DDA_ens$dea_res$adj.pvalue<=0.05)]
FG_TMT_ens<-top1_FG_TMT_ens$dea_res$protein[which(abs(top1_FG_TMT_ens$dea_res$logFC)>=log2(1.5) & top1_FG_TMT_ens$dea_res$adj.pvalue<=0.05)]
MQ_TMT_ens<-top1_MQ_TMT_ens$dea_res$protein[which(abs(top1_MQ_TMT_ens$dea_res$logFC)>=log2(1.5) & top1_MQ_TMT_ens$dea_res$adj.pvalue<=0.05)]

T_YEAST<-unique(union(union(union(all_protein_FG_DDA$Protein[which(all_protein_FG_DDA$Organism=='YEAST')],
                                  all_protein_MQ_DDA$Protein[which(all_protein_MQ_DDA$Organism=='YEAST')]),
                            all_protein_FG_TMT$Protein[which(all_protein_FG_TMT$Organism=='YEAST')]),
                      all_protein_MQ_TMT$Protein[which(all_protein_MQ_TMT$Organism=='YEAST')]))

T_HUMAN<-unique(union(union(union(all_protein_FG_DDA$Protein[which(all_protein_FG_DDA$Organism=='HUMAN')],
                                  all_protein_MQ_DDA$Protein[which(all_protein_MQ_DDA$Organism=='HUMAN')]),
                            all_protein_FG_TMT$Protein[which(all_protein_FG_TMT$Organism=='HUMAN')]),
                      all_protein_MQ_TMT$Protein[which(all_protein_MQ_TMT$Organism=='HUMAN')]))

DEP_HUMAN<-unique(union(union(union(top1_FG_DDA$dea_res$protein[which(abs(top1_FG_DDA$dea_res$logFC)>=log2(1.5) &
                                                                        top1_FG_DDA$dea_res$adj.pvalue<=0.05 &
                                                                        top1_FG_DDA$dea_res$organism=='HUMAN')],
                                    top1_MQ_DDA$dea_res$protein[which(abs(top1_MQ_DDA$dea_res$logFC)>=log2(1.5) &
                                                                        top1_MQ_DDA$dea_res$adj.pvalue<=0.05 &
                                                                        top1_MQ_DDA$dea_res$organism=='HUMAN')]),
                              top1_FG_TMT$dea_res$protein[which(abs(top1_FG_TMT$dea_res$logFC)>=log2(1.5) &
                                                                  top1_FG_TMT$dea_res$adj.pvalue<=0.05 &
                                                                  top1_FG_TMT$dea_res$organism=='HUMAN')]),
                        top1_MQ_TMT$dea_res$protein[which(abs(top1_MQ_TMT$dea_res$logFC)>=log2(1.5) &
                                                            top1_MQ_TMT$dea_res$adj.pvalue<=0.05 &
                                                            top1_MQ_TMT$dea_res$organism=='HUMAN')]))

DEP_HUMAN_ens<-unique(union(union(union(top1_FG_DDA_ens$dea_res$protein[which(abs(top1_FG_DDA_ens$dea_res$logFC)>=log2(1.5) &
                                                                                top1_FG_DDA_ens$dea_res$adj.pvalue<=0.05 &
                                                                                top1_FG_DDA_ens$dea_res$organism=='HUMAN')],
                                        top1_MQ_DDA_ens$dea_res$protein[which(abs(top1_MQ_DDA_ens$dea_res$logFC)>=log2(1.5) &
                                                                                top1_MQ_DDA_ens$dea_res$adj.pvalue<=0.05 &
                                                                                top1_MQ_DDA_ens$dea_res$organism=='HUMAN')]),
                                  top1_FG_TMT_ens$dea_res$protein[which(abs(top1_FG_TMT_ens$dea_res$logFC)>=log2(1.5) &
                                                                          top1_FG_TMT_ens$dea_res$adj.pvalue<=0.05 &
                                                                          top1_FG_TMT_ens$dea_res$organism=='HUMAN')]),
                            top1_MQ_TMT_ens$dea_res$protein[which(abs(top1_MQ_TMT_ens$dea_res$logFC)>=log2(1.5) &
                                                                    top1_MQ_TMT_ens$dea_res$adj.pvalue<=0.05 &
                                                                    top1_MQ_TMT_ens$dea_res$organism=='HUMAN')]))

DEP_HUMAN<-unique(union(DEP_HUMAN, DEP_HUMAN_ens))

DEPs<-data.frame(Proteins=unique(union(union(union(all_protein_FG_DDA$Protein,
                                                   all_protein_MQ_DDA$Protein),
                                             all_protein_FG_TMT$Protein),
                                       all_protein_MQ_TMT$Protein)))
DEPs$T_YEAST=c(rep(0, length(DEPs$Proteins)))
DEPs$T_YEAST[match(T_YEAST, DEPs$Proteins)]=1
# DEPs$T_ECOLI=c(rep(0, length(DEPs$Proteins)))
# DEPs$T_ECOLI[match(T_ECOLI, DEPs$Proteins)]=1
# DEPs$T_HUMAN=c(rep(0, length(DEPs$Proteins)))
# DEPs$T_HUMAN[match(T_HUMAN, DEPs$Proteins)]=1
DEPs$DEP_HUMAN=c(rep(0, length(DEPs$Proteins)))
DEPs$DEP_HUMAN[match(DEP_HUMAN, DEPs$Proteins)]=1
DEPs$FG_DDA=c(rep(0, length(DEPs$Proteins)))
DEPs$FG_DDA[match(FG_DDA, DEPs$Proteins)]=1
DEPs$MQ_DDA=c(rep(0, length(DEPs$Proteins)))
DEPs$MQ_DDA[match(MQ_DDA, DEPs$Proteins)]=1
DEPs$FG_TMT=c(rep(0, length(DEPs$Proteins)))
DEPs$FG_TMT[match(FG_TMT, DEPs$Proteins)]=1
DEPs$MQ_TMT=c(rep(0, length(DEPs$Proteins)))
DEPs$MQ_TMT[match(MQ_TMT, DEPs$Proteins)]=1

DEPs$FG_DDA_ens=c(rep(0, length(DEPs$Proteins)))
DEPs$FG_DDA_ens[match(FG_DDA_ens, DEPs$Proteins)]=1
DEPs$MQ_DDA_ens=c(rep(0, length(DEPs$Proteins)))
DEPs$MQ_DDA_ens[match(MQ_DDA_ens, DEPs$Proteins)]=1
DEPs$FG_TMT_ens=c(rep(0, length(DEPs$Proteins)))
DEPs$FG_TMT_ens[match(FG_TMT_ens, DEPs$Proteins)]=1
DEPs$MQ_TMT_ens=c(rep(0, length(DEPs$Proteins)))
DEPs$MQ_TMT_ens[match(MQ_TMT_ens, DEPs$Proteins)]=1
# simply remove protein groups
DEPs<-DEPs[setdiff(c(1:length(DEPs$Proteins)),grep(';', DEPs$Proteins)),]

sheet=paste(dt, '_cp_DEPs')

addWorksheet(wb,sheet)
writeData(wb, sheet, DEPs, startRow = 1, startCol = 1, colNames = TRUE, rowNames = TRUE)
saveWorkbook(wb, file = paste0(save_fold, 'cross_platform_compare.xlsx'), overwrite = TRUE)
# library(UpSetR)
# UpSetR::upset(DEPs,
#               sets = c("T_YEAST", "DEP_HUMAN", "FG_DDA", "MQ_DDA", "FG_TMT","MQ_TMT",
#                        "FG_DDA_ens", "MQ_DDA_ens", "FG_TMT_ens","MQ_TMT_ens"),
#               order.by="freq", matrix.color="black", point.size=3,
#               sets.bar.color=c("#15d08a","#f2ca01","#a45ee3","#1f11b4","#c3642c",
#                                "#f078b4", '#5494cc',"gray",'red','orange'),mb.ratio=c(0.7, 0.3), text.scale = 1.5,
#               number.angles =-20)



library(pROC)
library(ggpubr)

#all_pvals<-read.table('E:/proteomics/manu4_1/codes/reproduce_0_05/temp0_limma_min_FragPipe/HYEtims735_LFQ_top0_maxlfq_mv_outputs.csv',header = T, sep = ',')

cross_platform_proteins<-data.frame(Proteins=DEPs$Proteins, labels=0)
cross_platform_proteins$labels[grep('YEAST', cross_platform_proteins$Proteins)]=1
cross_platform_proteins$qval_FG_DDA<-1
cross_platform_proteins$logFC_FG_DDA<-0
com=intersect(cross_platform_proteins$Proteins, top1_FG_DDA$dea_res$protein)
idx1=match(com, cross_platform_proteins$Proteins)
idx2=match(com, top1_FG_DDA$dea_res$protein)
cross_platform_proteins$qval_FG_DDA[idx1]=top1_FG_DDA$dea_res$adj.pvalue[idx2]
cross_platform_proteins$logFC_FG_DDA[idx1]=top1_FG_DDA$dea_res$logFC[idx2]

cross_platform_proteins$qval_MQ_DDA<-1
cross_platform_proteins$logFC_MQ_DDA<-0
com=intersect(cross_platform_proteins$Proteins, top1_MQ_DDA$dea_res$protein)
idx1=match(com, cross_platform_proteins$Proteins)
idx2=match(com, top1_MQ_DDA$dea_res$protein)
cross_platform_proteins$qval_MQ_DDA[idx1]=top1_MQ_DDA$dea_res$adj.pvalue[idx2]
cross_platform_proteins$logFC_MQ_DDA[idx1]=top1_MQ_DDA$dea_res$logFC[idx2]

cross_platform_proteins$qval_FG_TMT<-1
cross_platform_proteins$logFC_FG_TMT<-0
com=intersect(cross_platform_proteins$Proteins, top1_FG_TMT$dea_res$protein)
idx1=match(com, cross_platform_proteins$Proteins)
idx2=match(com, top1_FG_TMT$dea_res$protein)
cross_platform_proteins$qval_FG_TMT[idx1]=top1_FG_TMT$dea_res$adj.pvalue[idx2]
cross_platform_proteins$logFC_FG_TMT[idx1]=top1_FG_TMT$dea_res$logFC[idx2]

cross_platform_proteins$qval_MQ_TMT<-1
cross_platform_proteins$logFC_MQ_TMT<-0
com=intersect(cross_platform_proteins$Proteins, top1_MQ_TMT$dea_res$protein)
idx1=match(com, cross_platform_proteins$Proteins)
idx2=match(com, top1_MQ_TMT$dea_res$protein)
cross_platform_proteins$qval_MQ_TMT[idx1]=top1_MQ_TMT$dea_res$adj.pvalue[idx2]
cross_platform_proteins$logFC_MQ_TMT[idx1]=top1_MQ_TMT$dea_res$logFC[idx2]

cross_platform_proteins$qval_FG_DDA_ens<-1
cross_platform_proteins$logFC_FG_DDA_ens<-0
com=intersect(cross_platform_proteins$Proteins, top1_FG_DDA_ens$dea_res$protein)
idx1=match(com, cross_platform_proteins$Proteins)
idx2=match(com, top1_FG_DDA_ens$dea_res$protein)
cross_platform_proteins$qval_FG_DDA_ens[idx1]=top1_FG_DDA_ens$dea_res$adj.pvalue[idx2]
cross_platform_proteins$logFC_FG_DDA_ens[idx1]=top1_FG_DDA_ens$dea_res$logFC[idx2]

cross_platform_proteins$qval_MQ_DDA_ens<-1
cross_platform_proteins$logFC_MQ_DDA_ens<-0
com=intersect(cross_platform_proteins$Proteins, top1_MQ_DDA_ens$dea_res$protein)
idx1=match(com, cross_platform_proteins$Proteins)
idx2=match(com, top1_MQ_DDA_ens$dea_res$protein)
cross_platform_proteins$qval_MQ_DDA_ens[idx1]=top1_MQ_DDA_ens$dea_res$adj.pvalue[idx2]
cross_platform_proteins$logFC_MQ_DDA_ens[idx1]=top1_MQ_DDA_ens$dea_res$logFC[idx2]

cross_platform_proteins$qval_FG_TMT_ens<-1
cross_platform_proteins$logFC_FG_TMT_ens<-0
com=intersect(cross_platform_proteins$Proteins, top1_FG_TMT_ens$dea_res$protein)
idx1=match(com, cross_platform_proteins$Proteins)
idx2=match(com, top1_FG_TMT_ens$dea_res$protein)
cross_platform_proteins$qval_FG_TMT_ens[idx1]=top1_FG_TMT_ens$dea_res$adj.pvalue[idx2]
cross_platform_proteins$logFC_FG_TMT_ens[idx1]=top1_FG_TMT_ens$dea_res$logFC[idx2]

cross_platform_proteins$qval_MQ_TMT_ens<-1
cross_platform_proteins$logFC_MQ_TMT_ens<-0
com=intersect(cross_platform_proteins$Proteins, top1_MQ_TMT_ens$dea_res$protein)
idx1=match(com, cross_platform_proteins$Proteins)
idx2=match(com, top1_MQ_TMT_ens$dea_res$protein)
cross_platform_proteins$qval_MQ_TMT_ens[idx1]=top1_MQ_TMT_ens$dea_res$adj.pvalue[idx2]
cross_platform_proteins$logFC_MQ_TMT_ens[idx1]=top1_MQ_TMT_ens$dea_res$logFC[idx2]

sheet=paste(dt, '_cp_proteins')

addWorksheet(wb,sheet)
writeData(wb, sheet, cross_platform_proteins, startRow = 1, startCol = 1, colNames = TRUE, rowNames = TRUE)
saveWorkbook(wb, file = paste0(save_fold, 'cross_platform_compare.xlsx'), overwrite = TRUE)

roc_FG_DDA <- roc(cross_platform_proteins$labels, 1-cross_platform_proteins$qval_FG_DDA)
roc_MQ_DDA <- roc(cross_platform_proteins$labels, 1-cross_platform_proteins$qval_MQ_DDA)
roc_FG_TMT <- roc(cross_platform_proteins$labels, 1-cross_platform_proteins$qval_FG_TMT)
roc_MQ_TMT <- roc(cross_platform_proteins$labels, 1-cross_platform_proteins$qval_MQ_TMT)

roc_FG_DDA_ens <- roc(cross_platform_proteins$labels, 1-cross_platform_proteins$qval_FG_DDA_ens)
roc_MQ_DDA_ens <- roc(cross_platform_proteins$labels, 1-cross_platform_proteins$qval_MQ_DDA_ens)
roc_FG_TMT_ens <- roc(cross_platform_proteins$labels, 1-cross_platform_proteins$qval_FG_TMT_ens)
roc_MQ_TMT_ens <- roc(cross_platform_proteins$labels, 1-cross_platform_proteins$qval_MQ_TMT_ens)

spe_sens<-rbind(cbind(roc_FG_DDA$specificities, roc_FG_DDA$sensitivities,
                      c(rep(paste0('FG_DDA(',as.character(round(auc(roc_FG_DDA,
                                                                    partial.auc=c(1, 0.95),
                                                                    partial.auc.correct=TRUE),4)),')'),
                            length(roc_FG_DDA$specificities)))),
                cbind(roc_MQ_DDA$specificities, roc_MQ_DDA$sensitivities,
                      c(rep(paste0('MQ_DDA(',as.character(round(auc(roc_MQ_DDA,
                                                                    partial.auc=c(1, 0.95),
                                                                    partial.auc.correct=TRUE),4)),')'),
                            length(roc_MQ_DDA$specificities)))),
                cbind(roc_FG_TMT$specificities, roc_FG_TMT$sensitivities,
                      c(rep(paste0('FG_TMT(',as.character(round(auc(roc_FG_TMT,
                                                                    partial.auc=c(1, 0.95),
                                                                    partial.auc.correct=TRUE),4)),')'),
                            length(roc_FG_TMT$specificities)))),
                cbind(roc_MQ_TMT$specificities, roc_MQ_TMT$sensitivities,
                      c(rep(paste0('MQ_TMT(',as.character(round(auc(roc_MQ_TMT,
                                                                    partial.auc=c(1, 0.95),
                                                                    partial.auc.correct=TRUE),4)),')'),
                            length(roc_MQ_TMT$specificities)))),
                cbind(roc_FG_DDA_ens$specificities, roc_FG_DDA_ens$sensitivities,
                      c(rep(paste0('FG_DDA_ens(',as.character(round(auc(roc_FG_DDA_ens,
                                                                        partial.auc=c(1, 0.95),
                                                                        partial.auc.correct=TRUE),4)),')'),
                            length(roc_FG_DDA_ens$specificities)))),
                cbind(roc_MQ_DDA_ens$specificities, roc_MQ_DDA_ens$sensitivities,
                      c(rep(paste0('MQ_DDA_ens(',as.character(round(auc(roc_MQ_DDA_ens,
                                                                        partial.auc=c(1, 0.95),
                                                                        partial.auc.correct=TRUE),4)),')'),
                            length(roc_MQ_DDA_ens$specificities)))),
                cbind(roc_FG_TMT_ens$specificities, roc_FG_TMT_ens$sensitivities,
                      c(rep(paste0('FG_TMT_ens(',as.character(round(auc(roc_FG_TMT_ens,
                                                                        partial.auc=c(1, 0.95),
                                                                        partial.auc.correct=TRUE),4)),')'),
                            length(roc_FG_TMT_ens$specificities)))),
                cbind(roc_MQ_TMT_ens$specificities, roc_MQ_TMT_ens$sensitivities,
                      c(rep(paste0('MQ_TMT_ens(',as.character(round(auc(roc_MQ_TMT_ens,
                                                                        partial.auc=c(1, 0.95),
                                                                        partial.auc.correct=TRUE),4)),')'),
                            length(roc_MQ_TMT_ens$specificities))))
)

colnames(spe_sens)<-c('Specificity', 'Sensitivity', 'Method(pAUC(FPR=0.05))')
spe_sens<-as.data.frame(spe_sens)
spe_sens$Specificity<-1-as.numeric(spe_sens$Specificity)
spe_sens$Sensitivity<-as.numeric(spe_sens$Sensitivity)

breaks = seq(0,1,0.1)
ggplot(spe_sens, aes(x = Specificity, y = Sensitivity, colour = `Method(pAUC(FPR=0.05))`))+
  scale_color_manual(values=c("#15d08a","#f2ca01","#a45ee3","#1f11b4","#c3642c",
                              "#f078b4", '#5494cc',"gray","pink",'red','orange')) +
  #geom_segment(aes(x = 0, y = 1, xend = 1,yend = 0),linewidth=10, alpha = 0.5, colour = "gray") +
  geom_step(linewidth=1,alpha = 0.9) +
  scale_x_continuous(name = "False Positive Rate (1 - Specificity)",limits = c(0,0.1), breaks = breaks) +
  scale_y_continuous(name = "True Positive Rate (Sensitivity)", limits = c(0,0.7), breaks = breaks) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 1))+scale_x_continuous(limits = c(0, 1))+
  labs(
    x = "False Positive Rate (1-Specificity)",
    y = "True Positive Rate (Sensitivity)")+ theme(legend.position = c(0.7, 0.35)) +
  geom_vline(xintercept = 0.05,linetype="dashed")

cross_platform_metrics<-function(cross_platform_proteins, setting, roc){
  lab=cross_platform_proteins$labels
  if(setting=='FG_DDA'){
    logFCs=cross_platform_proteins$logFC_FG_DDA
    qvals=cross_platform_proteins$qval_FG_DDA
    roc_s<-roc_FG_DDA
  }else if(setting=='MQ_DDA'){
    logFCs=cross_platform_proteins$logFC_MQ_DDA
    qvals=cross_platform_proteins$qval_MQ_DDA
    roc_s<-roc_MQ_DDA
  }else if(setting=='FG_TMT'){
    logFCs=cross_platform_proteins$logFC_FG_TMT
    qvals=cross_platform_proteins$qval_FG_TMT
    roc_s<-roc_FG_TMT
  }else if(setting=='MQ_TMT'){
    logFCs=cross_platform_proteins$logFC_MQ_TMT
    qvals=cross_platform_proteins$qval_MQ_TMT
    roc_s<-roc_MQ_TMT
  }else if(setting=='FG_DDA_ens'){
    logFCs=cross_platform_proteins$logFC_FG_DDA_ens
    qvals=cross_platform_proteins$qval_FG_DDA_ens
    roc_s<-roc_FG_DDA_ens
  }else if(setting=='MQ_DDA_ens'){
    logFCs=cross_platform_proteins$logFC_MQ_DDA_ens
    qvals=cross_platform_proteins$qval_MQ_DDA_ens
    roc_s<-roc_MQ_DDA_ens
  }else if(setting=='FG_TMT_ens'){
    logFCs=cross_platform_proteins$logFC_FG_TMT_ens
    qvals=cross_platform_proteins$qval_FG_TMT_ens
    roc_s<-roc_FG_TMT_ens
  }else if(setting=='MQ_TMT_ens'){
    logFCs=cross_platform_proteins$logFC_MQ_TMT_ens
    qvals=cross_platform_proteins$qval_MQ_TMT_ens
    roc_s<-roc_MQ_TMT_ens
  }
  TP=length(which(lab==1 & abs(logFCs)>=log2(1.5) & qvals<=0.05))
  TN=length(which(lab==0 & (abs(logFCs)<=log2(1.5) | qvals>0.05)))
  FP=length(which(lab==0 & abs(logFCs)>=log2(1.5) & qvals<=0.05))
  FN=length(which(lab==1 & (abs(logFCs)<=log2(1.5) | qvals>0.05)))

  Recall = TP/(TP+FN)
  Precision = TP/(TP+FP)
  Specificity = TN/(TN+FP)
  F1=2*Recall*Precision/(Recall+Precision)
  MCC = (TP*TN-FP*FN)/(((TP+FP)^0.5)*((TP+FN)^0.5)*((TN+FP)^0.5)*((TN+FN))^0.5)
  Gmean = (Recall*Specificity)^0.5
  nMCC = (1+MCC)/2

  pauc001=as.numeric(auc(roc_s, partial.auc=c(1, 0.99), partial.auc.correct=TRUE))
  pauc005=as.numeric(auc(roc_s, partial.auc=c(1, 0.95), partial.auc.correct=TRUE))
  pauc01=as.numeric(auc(roc_s, partial.auc=c(1, 0.90), partial.auc.correct=TRUE))
  metrics=data.frame(metric=c('pAUC(0.01)','pAUC(0.05)','pAUC(0.1)','TP', 'TN', 'FP', 'FN', 'Recall', 'Precision', 'Specificity', 'F1', 'MCC', 'G-mean','nMCC'),
                     value=c(pauc001,pauc005,pauc01,TP, TN, FP, FN, Recall, Precision, Specificity, F1, MCC, Gmean, nMCC))
  return(metrics)
}

metrics_FG_DDA_ens=cross_platform_metrics(cross_platform_proteins, 'FG_DDA_ens', roc_FG_DDA_ens)
metrics_MQ_DDA_ens=cross_platform_metrics(cross_platform_proteins, 'MQ_DDA_ens', roc_MQ_DDA_ens)
metrics_FG_TMT_ens=cross_platform_metrics(cross_platform_proteins, 'FG_TMT_ens', roc_FG_TMT_ens)
metrics_MQ_TMT_ens=cross_platform_metrics(cross_platform_proteins, 'MQ_TMT_ens', roc_MQ_TMT_ens)

metrics_FG_DDA=cross_platform_metrics(cross_platform_proteins, 'FG_DDA', roc_FG_DDA)
metrics_MQ_DDA=cross_platform_metrics(cross_platform_proteins, 'MQ_DDA', roc_MQ_DDA)
metrics_FG_TMT=cross_platform_metrics(cross_platform_proteins, 'FG_TMT', roc_FG_TMT)
metrics_MQ_TMT=cross_platform_metrics(cross_platform_proteins, 'MQ_TMT', roc_MQ_TMT)

rdar_dt<-rbind(metrics_FG_DDA$value, metrics_MQ_DDA$value, metrics_FG_TMT$value,
               metrics_MQ_TMT$value,metrics_FG_DDA_ens$value, metrics_MQ_DDA_ens$value,
               metrics_FG_TMT_ens$value,
               metrics_MQ_TMT_ens$value)
colnames(rdar_dt)<-metrics_FG_DDA$metric
row.names(rdar_dt)<-c('FG_DDA','MQ_DDA', 'FG_TMT', 'MQ_TMT', 'FG_DDA_ens',
                      'MQ_DDA_ens', 'FG_TMT_ens', 'MQ_TMT_ens')
rdar_dt<-as.data.frame(as.matrix(rdar_dt))
sheet=paste(dt, '_cp_metrics')

addWorksheet(wb,sheet)
writeData(wb, sheet, rdar_dt, startRow = 1, startCol = 1, colNames = TRUE, rowNames = TRUE)
saveWorkbook(wb, file = paste0(save_fold, 'cross_platform_compare.xlsx'), overwrite = TRUE)
pt_rdar_dt<-data.frame(group=row.names(rdar_dt))
pt_rdar_dt<-cbind(pt_rdar_dt, rdar_dt[,c(1,2,3,13,14)]) #,4,6

col=c("#15d08a","#f2ca01","#a45ee3","#1f11b4","#c3642c",
      "#f078b4", '#5494cc',"gray","pink",'red','orange')
melt_dt<-melt(pt_rdar_dt)
colnames(melt_dt)[1]<-'method'


ggdotchart(melt_dt, "variable", "value", group = "method", color = "method",
           palette = col[c(6,1,8,3,10,7,5,4)],

           add = "segments", #label = "value",
           add.params = list(color = "lightgray", size = 1),
           dot.size = 3,
           #label = 'value',
           #font.label = list(color = "black", size = 9, vjust = 0.5),
           position=position_dodge(width = 0.9)
) +
  #geom_text(aes(label=value), position = 'indentity')+
  theme(axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
        axis.title.x = element_blank())#+#+scale_y_continuous(limits = c(0.2, 1))+



count_TPs<-data.frame(method=row.names(rdar_dt))
count_TPs<-cbind(count_TPs, rdar_dt[,c(4:7)]) #,4,6
melt_dt<-melt(count_TPs)
melt_dt$label<-as.character(melt_dt$value)
melt_dt$label[which(melt_dt$variable=='TN' | melt_dt$variable=='FN')]=''

nums<-vector()


for (i in 1:length(unique(melt_dt$method))){
  num_i<-melt_dt[which(melt_dt$method==melt_dt$method[i]),]
  uni_var<-levels(factor(melt_dt$variable))
  re_order<-match(num_i$variable, uni_var)
  num_i<-num_i[re_order,]
  l_y<-c(rep(0,length(uni_var)))
  for (j in length(uni_var):1) {
    if(j==length(uni_var)){
      l_y[j]=num_i$value[j]/2
    }else{
      l_y_j=num_i$value[j]/2 + sum(num_i$value[c((j+1):length(uni_var))])
      l_y[j]=l_y_j
    }
  }
  num_i$l_y<-l_y
  nums<-rbind(nums,num_i)
}


## Figure 5C
ggbarplot(nums, "method", "value", group = "variable", color = "variable",
          fill="variable", #label = T,lab.pos='in',
          platette = c('#5494cc','#e18283','#0d898a','#f9cc52'),
) + geom_text(aes(y=l_y, label=label))+
  scale_x_discrete(limits=c('FG_DDA_ens', 'MQ_DDA_ens',
                            'FG_TMT_ens','MQ_TMT_ens', 'FG_DDA','MQ_DDA','FG_TMT','MQ_TMT')) +
  geom_vline(xintercept = 4.5, linetype=2, show.legend=FALSE) + theme(axis.line.y = element_blank())+
  coord_flip()+theme(axis.title.y = element_blank())

write.table(nums, paste0(save_fold, 'Figure5C.csv'), col.names = T, row.names = F)
###################################################
################# HEqe408
##############
################

top1wf<-function(dt, platform, acq){
  base_fold<-data_fold
  dea_res_fold<-paste0(base_fold, acq, '/',platform, '/',dt,'/')
  ranking_all_FG_DDA<-read_excel(paste0(data_fold, 'ranks_all_', platform, '_', acq, '.xlsx'), sheet = 'ranking_all')
  #ranking_all<-ranking_all[which(ranking_all$runTime!=0 & ranking_all$mean_pauc001!=0 & ranking_all$rank_median_pauc001!=0),]
  wftop1_FG_DDA<-ranking_all_FG_DDA[which(ranking_all_FG_DDA$avg_rank_mean==min(ranking_all_FG_DDA$avg_rank_mean)),]
  if(length((wftop1_FG_DDA$workflow))>1){
    idx=which(wftop1_FG_DDA$avg_rank_median==min(wftop1_FG_DDA$avg_rank_median))
    wftop1_FG_DDA = wftop1_FG_DDA[idx,]
  }
  dea_res_FG_DDA_wf1<-read.table(paste0(dea_res_fold,
                                        dt,'_',wftop1_FG_DDA$DEA,'_',
                                        wftop1_FG_DDA$Platform,'_',
                                        wftop1_FG_DDA$Matrix,'_',gsub('blank','',wftop1_FG_DDA$Imput),
                                        '_',gsub('blank','',wftop1_FG_DDA$normalization),'.csv'),
                                 sep = ',',header = T)
  return(list(wf=wftop1_FG_DDA, dea_res=dea_res_FG_DDA_wf1))
}

top1wf_ens<-function(dt, platform, acq, ens_ty){
  base_fold<-data_fold
  if(ens_ty=='ensemble_mv'){
    dea_res_fold<-paste0(base_fold, acq, '/',platform, '/','mv','/',dt,'/')
  }else if(ens_ty=='ensemble_topk'){
    dea_res_fold<-paste0(base_fold, acq, '/',platform, '/','topk','/',dt,'/')
  }
  #dea_res_fold<-paste0(base_fold, acq, '/',platform, '/',dt,'/')
  ranking_all_FG_DDA<-read_excel(paste0(data_fold, 'ranks_all_', platform, '_', acq, '_', ens_ty, '.xlsx'), sheet = 'ranking_all')
  #ranking_all<-ranking_all[which(ranking_all$runTime!=0 & ranking_all$mean_pauc001!=0 & ranking_all$rank_median_pauc001!=0),]
  wftop1_FG_DDA<-ranking_all_FG_DDA[which(ranking_all_FG_DDA$avg_rank_mean==min(ranking_all_FG_DDA$avg_rank_mean)),]
  if(length((wftop1_FG_DDA$workflow))>1){
    idx=which(wftop1_FG_DDA$avg_rank_median==min(wftop1_FG_DDA$avg_rank_median))
    wftop1_FG_DDA = wftop1_FG_DDA[idx,]
  }
  if(wftop1_FG_DDA$method=='set'){
    dea_res_FG_DDA_wf1<-read.table(paste0(dea_res_fold,
                                          dt,'_',
                                          wftop1_FG_DDA$Platform,'_',
                                          wftop1_FG_DDA$method, '_',
                                          wftop1_FG_DDA$operation,
                                          '_',gsub('\\|','-',wftop1_FG_DDA$cbn),'.csv'),
                                   sep = ',',header = T)
  }else{
    dea_res_FG_DDA_wf1<-read.table(paste0(dea_res_fold,
                                          dt,'_',
                                          wftop1_FG_DDA$Platform,'_',
                                          wftop1_FG_DDA$method,'_',gsub('\\|','-',wftop1_FG_DDA$cbn),'.csv'),
                                   sep = ',',header = T)
  }

  return(list(wf=wftop1_FG_DDA, dea_res=dea_res_FG_DDA_wf1))
}

dt='HEqe408'

all_protein_FG_DDA<-read.table(paste0(data_fold, 'HEqe408_LFQ_FragPipe_all_proteins.tsv'), sep = '\t', header = T)
all_protein_MQ_DDA<-read.table(paste0(data_fold, 'HEqe408_LFQ_Maxquant_all_proteins.tsv'), sep = '\t', header = T)
all_protein_DIANN_DIA<-read.table(paste0(data_fold, 'HEqe408_DIA_DIANN_all_proteins.tsv'), sep = '\t', header = T)
all_protein_spt_DIA<-read.table(paste0(data_fold, 'HEqe408_DIA_spt_all_proteins.tsv'), sep = '\t', header = T)

top1_FG_DDA<-top1wf('HEqe408_LFQ', 'FragPipe', 'DDA')
top1_MQ_DDA<-top1wf('HEqe408_LFQ', 'Maxquant', 'DDA')
top1_DIANN_DIA<-top1wf('HEqe408_DIA', 'DIANN', 'DIA')
top1_spt_DIA<-top1wf('HEqe408_DIA', 'spt', 'DIA')
top1_FG_DDA_ens=top1wf_ens('HEqe408_LFQ', 'FragPipe', 'DDA', 'ensemble_mv')
top1_MQ_DDA_ens=top1wf_ens('HEqe408_LFQ', 'Maxquant', 'DDA', 'ensemble_mv')
top1_DIANN_DIA_ens=top1wf_ens('HEqe408_DIA', 'DIANN', 'DIA', 'ensemble_mv')
top1_spt_DIA_ens=top1wf_ens('HEqe408_DIA', 'spt', 'DIA', 'ensemble_mv')



FG_DDA<-top1_FG_DDA$dea_res$protein[which(abs(top1_FG_DDA$dea_res$logFC)>=log2(1.5) & top1_FG_DDA$dea_res$adj.pvalue<=0.05)]
MQ_DDA<-top1_MQ_DDA$dea_res$protein[which(abs(top1_MQ_DDA$dea_res$logFC)>=log2(1.5) & top1_MQ_DDA$dea_res$adj.pvalue<=0.05)]
DIANN_DIA<-top1_DIANN_DIA$dea_res$protein[which(abs(top1_DIANN_DIA$dea_res$logFC)>=log2(1.5) & top1_DIANN_DIA$dea_res$adj.pvalue<=0.05)]
spt_DIA<-top1_spt_DIA$dea_res$protein[which(abs(top1_spt_DIA$dea_res$logFC)>=log2(1.5) & top1_spt_DIA$dea_res$adj.pvalue<=0.05)]

FG_DDA_ens<-top1_FG_DDA_ens$dea_res$protein[which(abs(top1_FG_DDA_ens$dea_res$logFC)>=log2(1.5) & top1_FG_DDA_ens$dea_res$adj.pvalue<=0.05)]
MQ_DDA_ens<-top1_MQ_DDA_ens$dea_res$protein[which(abs(top1_MQ_DDA_ens$dea_res$logFC)>=log2(1.5) & top1_MQ_DDA_ens$dea_res$adj.pvalue<=0.05)]
DIANN_DIA_ens<-top1_DIANN_DIA_ens$dea_res$protein[which(abs(top1_DIANN_DIA_ens$dea_res$logFC)>=log2(1.5) & top1_DIANN_DIA_ens$dea_res$adj.pvalue<=0.05)]
spt_DIA_ens<-top1_spt_DIA_ens$dea_res$protein[which(abs(top1_spt_DIA_ens$dea_res$logFC)>=log2(1.5) & top1_spt_DIA_ens$dea_res$adj.pvalue<=0.05)]


# T_YEAST<-unique(union(union(union(all_protein_FG_DDA$Protein[which(all_protein_FG_DDA$Organism=='YEAST')],
#                                   all_protein_MQ_DDA$Protein[which(all_protein_MQ_DDA$Organism=='YEAST')]),
#                             all_protein_DIANN_DIA$Protein[which(all_protein_DIANN_DIA$Organism=='YEAST')]),
#                       all_protein_spt_DIA$Protein[which(all_protein_spt_DIA$Organism=='YEAST')]))
T_ECOLI<-unique(union(union(union(all_protein_FG_DDA$Protein[which(all_protein_FG_DDA$Organism=='ECOLI')],
                                  all_protein_MQ_DDA$Protein[which(all_protein_MQ_DDA$Organism=='ECOLI')]),
                            all_protein_DIANN_DIA$Protein[which(all_protein_DIANN_DIA$Organism=='ECOLI')]),
                      all_protein_spt_DIA$Protein[which(all_protein_spt_DIA$Organism=='ECOLI')]))
T_HUMAN<-unique(union(union(union(all_protein_FG_DDA$Protein[which(all_protein_FG_DDA$Organism=='HUMAN')],
                                  all_protein_MQ_DDA$Protein[which(all_protein_MQ_DDA$Organism=='HUMAN')]),
                            all_protein_DIANN_DIA$Protein[which(all_protein_DIANN_DIA$Organism=='HUMAN')]),
                      all_protein_spt_DIA$Protein[which(all_protein_spt_DIA$Organism=='HUMAN')]))

DEP_HUMAN<-unique(union(union(union(top1_FG_DDA$dea_res$protein[which(abs(top1_FG_DDA$dea_res$logFC)>=log2(1.5) &
                                                                        top1_FG_DDA$dea_res$adj.pvalue<=0.05 &
                                                                        top1_FG_DDA$dea_res$organism=='HUMAN')],
                                    top1_MQ_DDA$dea_res$protein[which(abs(top1_MQ_DDA$dea_res$logFC)>=log2(1.5) &
                                                                        top1_MQ_DDA$dea_res$adj.pvalue<=0.05 &
                                                                        top1_MQ_DDA$dea_res$organism=='HUMAN')]),
                              top1_DIANN_DIA$dea_res$protein[which(abs(top1_DIANN_DIA$dea_res$logFC)>=log2(1.5) &
                                                                     top1_DIANN_DIA$dea_res$adj.pvalue<=0.05 &
                                                                     top1_DIANN_DIA$dea_res$organism=='HUMAN')]),
                        top1_spt_DIA$dea_res$protein[which(abs(top1_spt_DIA$dea_res$logFC)>=log2(1.5) &
                                                             top1_spt_DIA$dea_res$adj.pvalue<=0.05 &
                                                             top1_spt_DIA$dea_res$organism=='HUMAN')]))


DEP_HUMAN_ens<-unique(union(union(union(top1_FG_DDA_ens$dea_res$protein[which(abs(top1_FG_DDA_ens$dea_res$logFC)>=log2(1.5) &
                                                                                top1_FG_DDA_ens$dea_res$adj.pvalue<=0.05 &
                                                                                top1_FG_DDA_ens$dea_res$organism=='HUMAN')],
                                        top1_MQ_DDA_ens$dea_res$protein[which(abs(top1_MQ_DDA_ens$dea_res$logFC)>=log2(1.5) &
                                                                                top1_MQ_DDA_ens$dea_res$adj.pvalue<=0.05 &
                                                                                top1_MQ_DDA_ens$dea_res$organism=='HUMAN')]),
                                  top1_DIANN_DIA_ens$dea_res$protein[which(abs(top1_DIANN_DIA_ens$dea_res$logFC)>=log2(1.5) &
                                                                             top1_DIANN_DIA_ens$dea_res$adj.pvalue<=0.05 &
                                                                             top1_DIANN_DIA_ens$dea_res$organism=='HUMAN')]),
                            top1_spt_DIA_ens$dea_res$protein[which(abs(top1_spt_DIA_ens$dea_res$logFC)>=log2(1.5) &
                                                                     top1_spt_DIA_ens$dea_res$adj.pvalue<=0.05 &
                                                                     top1_spt_DIA_ens$dea_res$organism=='HUMAN')]))

DEP_HUMAN<-unique(union(DEP_HUMAN, DEP_HUMAN_ens))

DEPs<-data.frame(Proteins=unique(union(union(union(all_protein_FG_DDA$Protein,
                                                   all_protein_MQ_DDA$Protein),
                                             all_protein_DIANN_DIA$Protein),
                                       all_protein_spt_DIA$Protein)))
# DEPs$T_YEAST=c(rep(0, length(DEPs$Proteins)))
# DEPs$T_YEAST[match(T_YEAST, DEPs$Proteins)]=1
DEPs$T_ECOLI=c(rep(0, length(DEPs$Proteins)))
DEPs$T_ECOLI[match(T_ECOLI, DEPs$Proteins)]=1
# DEPs$T_HUMAN=c(rep(0, length(DEPs$Proteins)))
# DEPs$T_HUMAN[match(T_HUMAN, DEPs$Proteins)]=1
DEPs$DEP_HUMAN=c(rep(0, length(DEPs$Proteins)))
DEPs$DEP_HUMAN[match(DEP_HUMAN, DEPs$Proteins)]=1
DEPs$FG_DDA=c(rep(0, length(DEPs$Proteins)))
DEPs$FG_DDA[match(FG_DDA, DEPs$Proteins)]=1
DEPs$MQ_DDA=c(rep(0, length(DEPs$Proteins)))
DEPs$MQ_DDA[match(MQ_DDA, DEPs$Proteins)]=1
DEPs$DIANN_DIA=c(rep(0, length(DEPs$Proteins)))
DEPs$DIANN_DIA[match(DIANN_DIA, DEPs$Proteins)]=1
DEPs$spt_DIA=c(rep(0, length(DEPs$Proteins)))
DEPs$spt_DIA[match(spt_DIA, DEPs$Proteins)]=1

DEPs$FG_DDA_ens=c(rep(0, length(DEPs$Proteins)))
DEPs$FG_DDA_ens[match(FG_DDA_ens, DEPs$Proteins)]=1
DEPs$MQ_DDA_ens=c(rep(0, length(DEPs$Proteins)))
DEPs$MQ_DDA_ens[match(MQ_DDA_ens, DEPs$Proteins)]=1
DEPs$DIANN_DIA_ens=c(rep(0, length(DEPs$Proteins)))
DEPs$DIANN_DIA_ens[match(DIANN_DIA_ens, DEPs$Proteins)]=1
DEPs$spt_DIA_ens=c(rep(0, length(DEPs$Proteins)))
DEPs$spt_DIA_ens[match(spt_DIA_ens, DEPs$Proteins)]=1

# simply remove protein groups
DEPs<-DEPs[setdiff(c(1:length(DEPs$Proteins)),grep(';', DEPs$Proteins)),]

sheet=paste(dt, '_cp_DEPs')

addWorksheet(wb,sheet)
writeData(wb, sheet, DEPs, startRow = 1, startCol = 1, colNames = TRUE, rowNames = TRUE)
saveWorkbook(wb, file = paste0(save_fold, 'cross_platform_compare.xlsx'), overwrite = TRUE)


# UpSetR::upset(DEPs,
#               sets = c("T_YEAST", "T_ECOLI", "DEP_HUMAN", "FG_DDA", "MQ_DDA", "DIANN_DIA","spt_DIA","FG_DDA_ens","MQ_DDA_ens", "DIANN_DIA_ens","spt_DIA_ens"),
#               order.by="freq", matrix.color="black", point.size=3,
#               sets.bar.color=c("#15d08a","#f2ca01","#a45ee3","#1f11b4","#c3642c",
#                                "#f078b4", '#5494cc',"gray","pink",'red','orange'),mb.ratio=c(0.5, 0.5), text.scale = 1.5,
#               number.angles =0)

library(pROC)
library(ggpubr)

#all_pvals<-read.table('E:/proteomics/manu4_1/codes/reproduce_0_05/temp0_limma_min_FragPipe/HYEtims735_LFQ_top0_maxlfq_mv_outputs.csv',header = T, sep = ',')

cross_platform_proteins<-data.frame(Proteins=DEPs$Proteins, labels=0)
#cross_platform_proteins$labels[grep('YEAST', cross_platform_proteins$Proteins)]=1
cross_platform_proteins$labels[grep('ECOLI', cross_platform_proteins$Proteins)]=1
cross_platform_proteins$qval_FG_DDA<-1
cross_platform_proteins$logFC_FG_DDA<-0
com=intersect(cross_platform_proteins$Proteins, top1_FG_DDA$dea_res$protein)
idx1=match(com, cross_platform_proteins$Proteins)
idx2=match(com, top1_FG_DDA$dea_res$protein)
cross_platform_proteins$qval_FG_DDA[idx1]=top1_FG_DDA$dea_res$adj.pvalue[idx2]
cross_platform_proteins$logFC_FG_DDA[idx1]=top1_FG_DDA$dea_res$logFC[idx2]

cross_platform_proteins$qval_MQ_DDA<-1
cross_platform_proteins$logFC_MQ_DDA<-0
com=intersect(cross_platform_proteins$Proteins, top1_MQ_DDA$dea_res$protein)
idx1=match(com, cross_platform_proteins$Proteins)
idx2=match(com, top1_MQ_DDA$dea_res$protein)
cross_platform_proteins$qval_MQ_DDA[idx1]=top1_MQ_DDA$dea_res$adj.pvalue[idx2]
cross_platform_proteins$logFC_MQ_DDA[idx1]=top1_MQ_DDA$dea_res$logFC[idx2]

cross_platform_proteins$qval_DIANN_DIA<-1
cross_platform_proteins$logFC_DIANN_DIA<-0
com=intersect(cross_platform_proteins$Proteins, top1_DIANN_DIA$dea_res$protein)
idx1=match(com, cross_platform_proteins$Proteins)
idx2=match(com, top1_DIANN_DIA$dea_res$protein)
cross_platform_proteins$qval_DIANN_DIA[idx1]=top1_DIANN_DIA$dea_res$adj.pvalue[idx2]
cross_platform_proteins$logFC_DIANN_DIA[idx1]=top1_DIANN_DIA$dea_res$logFC[idx2]

cross_platform_proteins$qval_spt_DIA<-1
cross_platform_proteins$logFC_spt_DIA<-0
com=intersect(cross_platform_proteins$Proteins, top1_spt_DIA$dea_res$protein)
idx1=match(com, cross_platform_proteins$Proteins)
idx2=match(com, top1_spt_DIA$dea_res$protein)
cross_platform_proteins$qval_spt_DIA[idx1]=top1_spt_DIA$dea_res$adj.pvalue[idx2]
cross_platform_proteins$logFC_spt_DIA[idx1]=top1_spt_DIA$dea_res$logFC[idx2]

cross_platform_proteins$qval_FG_DDA_ens<-1
cross_platform_proteins$logFC_FG_DDA_ens<-0
com=intersect(cross_platform_proteins$Proteins, top1_FG_DDA_ens$dea_res$protein)
idx1=match(com, cross_platform_proteins$Proteins)
idx2=match(com, top1_FG_DDA_ens$dea_res$protein)
cross_platform_proteins$qval_FG_DDA_ens[idx1]=top1_FG_DDA_ens$dea_res$adj.pvalue[idx2]
cross_platform_proteins$logFC_FG_DDA_ens[idx1]=top1_FG_DDA_ens$dea_res$logFC[idx2]

cross_platform_proteins$qval_MQ_DDA_ens<-1
cross_platform_proteins$logFC_MQ_DDA_ens<-0
com=intersect(cross_platform_proteins$Proteins, top1_MQ_DDA_ens$dea_res$protein)
idx1=match(com, cross_platform_proteins$Proteins)
idx2=match(com, top1_MQ_DDA_ens$dea_res$protein)
cross_platform_proteins$qval_MQ_DDA_ens[idx1]=top1_MQ_DDA_ens$dea_res$adj.pvalue[idx2]
cross_platform_proteins$logFC_MQ_DDA_ens[idx1]=top1_MQ_DDA_ens$dea_res$logFC[idx2]

cross_platform_proteins$qval_DIANN_DIA_ens<-1
cross_platform_proteins$logFC_DIANN_DIA_ens<-0
com=intersect(cross_platform_proteins$Proteins, top1_DIANN_DIA_ens$dea_res$protein)
idx1=match(com, cross_platform_proteins$Proteins)
idx2=match(com, top1_DIANN_DIA_ens$dea_res$protein)
cross_platform_proteins$qval_DIANN_DIA_ens[idx1]=top1_DIANN_DIA_ens$dea_res$adj.pvalue[idx2]
cross_platform_proteins$logFC_DIANN_DIA_ens[idx1]=top1_DIANN_DIA_ens$dea_res$logFC[idx2]

cross_platform_proteins$qval_spt_DIA_ens<-1
cross_platform_proteins$logFC_spt_DIA_ens<-0
com=intersect(cross_platform_proteins$Proteins, top1_spt_DIA_ens$dea_res$protein)
idx1=match(com, cross_platform_proteins$Proteins)
idx2=match(com, top1_spt_DIA_ens$dea_res$protein)
cross_platform_proteins$qval_spt_DIA_ens[idx1]=top1_spt_DIA_ens$dea_res$adj.pvalue[idx2]
cross_platform_proteins$logFC_spt_DIA_ens[idx1]=top1_spt_DIA_ens$dea_res$logFC[idx2]

sheet=paste(dt, '_cp_proteins')

addWorksheet(wb,sheet)
writeData(wb, sheet, cross_platform_proteins, startRow = 1, startCol = 1, colNames = TRUE, rowNames = TRUE)
saveWorkbook(wb, file = paste0(save_fold, 'cross_platform_compare.xlsx'), overwrite = TRUE)

roc_FG_DDA <- roc(cross_platform_proteins$labels, 1-cross_platform_proteins$qval_FG_DDA)
roc_MQ_DDA <- roc(cross_platform_proteins$labels, 1-cross_platform_proteins$qval_MQ_DDA)
roc_DIANN_DIA <- roc(cross_platform_proteins$labels, 1-cross_platform_proteins$qval_DIANN_DIA)
roc_spt_DIA <- roc(cross_platform_proteins$labels, 1-cross_platform_proteins$qval_spt_DIA)

roc_FG_DDA_ens <- roc(cross_platform_proteins$labels, 1-cross_platform_proteins$qval_FG_DDA_ens)
roc_MQ_DDA_ens <- roc(cross_platform_proteins$labels, 1-cross_platform_proteins$qval_MQ_DDA_ens)
roc_DIANN_DIA_ens <- roc(cross_platform_proteins$labels, 1-cross_platform_proteins$qval_DIANN_DIA_ens)
roc_spt_DIA_ens <- roc(cross_platform_proteins$labels, 1-cross_platform_proteins$qval_spt_DIA_ens)


cross_platform_metrics<-function(cross_platform_proteins, setting, roc){
  lab=cross_platform_proteins$labels
  if(setting=='FG_DDA'){
    logFCs=cross_platform_proteins$logFC_FG_DDA
    qvals=cross_platform_proteins$qval_FG_DDA
    roc_s<-roc_FG_DDA
  }else if(setting=='MQ_DDA'){
    logFCs=cross_platform_proteins$logFC_MQ_DDA
    qvals=cross_platform_proteins$qval_MQ_DDA
    roc_s<-roc_MQ_DDA
  }else if(setting=='DIANN_DIA'){
    logFCs=cross_platform_proteins$logFC_DIANN_DIA
    qvals=cross_platform_proteins$qval_DIANN_DIA
    roc_s<-roc_DIANN_DIA
  }else if(setting=='spt_DIA'){
    logFCs=cross_platform_proteins$logFC_spt_DIA
    qvals=cross_platform_proteins$qval_spt_DIA
    roc_s<-roc_spt_DIA
  }else if(setting=='FG_DDA_ens'){
    logFCs=cross_platform_proteins$logFC_FG_DDA_ens
    qvals=cross_platform_proteins$qval_FG_DDA_ens
    roc_s<-roc_FG_DDA_ens
  }else if(setting=='MQ_DDA_ens'){
    logFCs=cross_platform_proteins$logFC_MQ_DDA_ens
    qvals=cross_platform_proteins$qval_MQ_DDA_ens
    roc_s<-roc_MQ_DDA_ens
  }else if(setting=='DIANN_DIA_ens'){
    logFCs=cross_platform_proteins$logFC_DIANN_DIA_ens
    qvals=cross_platform_proteins$qval_DIANN_DIA_ens
    roc_s<-roc_DIANN_DIA_ens
  }else if(setting=='spt_DIA_ens'){
    logFCs=cross_platform_proteins$logFC_spt_DIA_ens
    qvals=cross_platform_proteins$qval_spt_DIA_ens
    roc_s<-roc_spt_DIA_ens
  }
  TP=length(which(lab==1 & abs(logFCs)>=log2(1.5) & qvals<=0.05))
  TN=length(which(lab==0 & (abs(logFCs)<=log2(1.5) | qvals>0.05)))
  FP=length(which(lab==0 & abs(logFCs)>=log2(1.5) & qvals<=0.05))
  FN=length(which(lab==1 & (abs(logFCs)<=log2(1.5) | qvals>0.05)))

  Recall = TP/(TP+FN)
  Precision = TP/(TP+FP)
  Specificity = TN/(TN+FP)
  F1=2*Recall*Precision/(Recall+Precision)
  MCC = (TP*TN-FP*FN)/(((TP+FP)^0.5)*((TP+FN)^0.5)*((TN+FP)^0.5)*((TN+FN))^0.5)
  Gmean = (Recall*Specificity)^0.5
  nMCC = (1+MCC)/2

  pauc001=as.numeric(auc(roc_s, partial.auc=c(1, 0.99), partial.auc.correct=TRUE))
  pauc005=as.numeric(auc(roc_s, partial.auc=c(1, 0.95), partial.auc.correct=TRUE))
  pauc01=as.numeric(auc(roc_s, partial.auc=c(1, 0.90), partial.auc.correct=TRUE))
  metrics=data.frame(metric=c('pAUC(0.01)','pAUC(0.05)','pAUC(0.1)','TP', 'TN', 'FP', 'FN', 'Recall', 'Precision', 'Specificity', 'F1', 'MCC', 'G-mean','nMCC'),
                     value=c(pauc001,pauc005,pauc01,TP, TN, FP, FN, Recall, Precision, Specificity, F1, MCC, Gmean,nMCC))
  return(metrics)
}

metrics_FG_DDA=cross_platform_metrics(cross_platform_proteins, 'FG_DDA', roc_FG_DDA)
metrics_MQ_DDA=cross_platform_metrics(cross_platform_proteins, 'MQ_DDA', roc_MQ_DDA)
metrics_DIANN_DIA=cross_platform_metrics(cross_platform_proteins, 'DIANN_DIA', roc_DIANN_DIA)
metrics_spt_DIA=cross_platform_metrics(cross_platform_proteins, 'spt_DIA', roc_spt_DIA)

metrics_FG_DDA_ens=cross_platform_metrics(cross_platform_proteins, 'FG_DDA_ens', roc_FG_DDA_ens)
metrics_MQ_DDA_ens=cross_platform_metrics(cross_platform_proteins, 'MQ_DDA_ens', roc_MQ_DDA_ens)
metrics_DIANN_DIA_ens=cross_platform_metrics(cross_platform_proteins, 'DIANN_DIA_ens', roc_DIANN_DIA_ens)
metrics_spt_DIA_ens=cross_platform_metrics(cross_platform_proteins, 'spt_DIA_ens', roc_spt_DIA_ens)

rdar_dt<-rbind(metrics_FG_DDA$value, metrics_MQ_DDA$value, metrics_DIANN_DIA$value, metrics_spt_DIA$value,
               metrics_FG_DDA_ens$value, metrics_MQ_DDA_ens$value, metrics_DIANN_DIA_ens$value, metrics_spt_DIA_ens$value)
colnames(rdar_dt)<-metrics_FG_DDA$metric
row.names(rdar_dt)<-c('FG_DDA','MQ_DDA', 'DIANN_DIA', 'spt_DIA', 'FG_DDA_ens','MQ_DDA_ens', 'DIANN_DIA_ens', 'spt_DIA_ens')
rdar_dt<-as.data.frame(as.matrix(rdar_dt))
pt_rdar_dt<-data.frame(group=row.names(rdar_dt))

sheet=paste(dt, '_cp_metrics')

addWorksheet(wb,sheet)
writeData(wb, sheet, rdar_dt, startRow = 1, startCol = 1, colNames = TRUE, rowNames = TRUE)
saveWorkbook(wb, file = paste0(save_fold, 'cross_platform_compare.xlsx'), overwrite = TRUE)
pt_rdar_dt<-cbind(pt_rdar_dt, rdar_dt[,c(1,2,3,13,14)]) #,4,6
#pt_rdar_dt$group<-row.names(pt_rdar_dt)
library(ggradar)
ggradar(pt_rdar_dt,group.colours=c("#a45ee3","#1f11b4","#f078b4","#f2ca01","gray","pink",'red',"orange"),legend.position = 'right')
#a=c('#5494cc','#e18283','#0d898a','#f9cc52')
col=c("#15d08a","#f2ca01","#a45ee3","#1f11b4","#c3642c",
      "#f078b4", '#5494cc',"gray","pink",'red','orange')
melt_dt<-melt(pt_rdar_dt)
colnames(melt_dt)[1]<-'method'

ggdotchart(melt_dt, "variable", "value", group = "method", color = "method",
           palette = col[c(6,1,8,3,10,7,5,4)],

           add = "segments", #label = "value",
           add.params = list(color = "lightgray", size = 1),
           dot.size = 3,
           #label = 'value',
           #font.label = list(color = "black", size = 9, vjust = 0.5),
           position=position_dodge(width = 0.9)
) +
  #geom_text(aes(label=value), position = 'indentity')+
  theme(axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
        axis.title.x = element_blank())#+#+scale_y_continuous(limits = c(0.2, 1))+
#coord_flip()


count_TPs<-data.frame(method=row.names(rdar_dt))
count_TPs<-cbind(count_TPs, rdar_dt[,c(4:7)]) #,4,6
melt_dt<-melt(count_TPs)
melt_dt$label<-as.character(melt_dt$value)
melt_dt$label[which(melt_dt$variable=='TN' | melt_dt$variable=='FN')]=''

nums<-vector()


for (i in 1:length(unique(melt_dt$method))){
  num_i<-melt_dt[which(melt_dt$method==melt_dt$method[i]),]
  uni_var<-levels(factor(melt_dt$variable))
  re_order<-match(num_i$variable, uni_var)
  num_i<-num_i[re_order,]
  l_y<-c(rep(0,length(uni_var)))
  for (j in length(uni_var):1) {
    if(j==length(uni_var)){
      l_y[j]=num_i$value[j]/2
    }else{
      l_y_j=num_i$value[j]/2 + sum(num_i$value[c((j+1):length(uni_var))])
      l_y[j]=l_y_j
    }
  }
  num_i$l_y<-l_y
  nums<-rbind(nums,num_i)
}

# supp. Fig. 5

ggbarplot(nums, "method", "value", group = "variable", color = "variable",
          fill="variable", #label = T,lab.pos='in',
          platette = c('#5494cc','#e18283','#0d898a','#f9cc52'),
) + geom_text(aes(y=l_y, label=label))+
  scale_x_discrete(limits=c('FG_DDA_ens','MQ_DDA_ens','DIANN_DIA_ens', 'spt_DIA_ens',
                            'FG_DDA','MQ_DDA','DIANN_DIA','spt_DIA')) +
  geom_vline(xintercept = 4.5, linetype=2, show.legend=FALSE) + theme(axis.line.y = element_blank())+
  coord_flip()+theme(axis.title.y = element_blank())
###################################################
################# HYtims134
##############
################

top1wf<-function(dt, platform, acq, contrast){
  base_fold<-data_fold
  dea_res_fold<-paste0(base_fold, acq, '/',platform, '/',dt,'/')
  ranking_all_FG_DDA<-read_excel(paste0(data_fold, 'ranks_all_', platform, '_', acq, '.xlsx'), sheet = 'ranking_all')
  #ranking_all<-ranking_all[which(ranking_all$runTime!=0 & ranking_all$mean_pauc001!=0 & ranking_all$rank_median_pauc001!=0),]
  wftop1_FG_DDA<-ranking_all_FG_DDA[which(ranking_all_FG_DDA$avg_rank_mean==min(ranking_all_FG_DDA$avg_rank_mean)),]
  if(length((wftop1_FG_DDA$workflow))>1){
    idx=which(wftop1_FG_DDA$avg_rank_median==min(wftop1_FG_DDA$avg_rank_median))
    wftop1_FG_DDA = wftop1_FG_DDA[idx,]
  }
  dea_res_FG_DDA_wf1<-read.table(paste0(dea_res_fold,
                                        dt,'_',wftop1_FG_DDA$DEA,'_',
                                        wftop1_FG_DDA$Platform,'_',
                                        wftop1_FG_DDA$Matrix,'_',gsub('blank','',wftop1_FG_DDA$Imput),
                                        '_',gsub('blank','',wftop1_FG_DDA$normalization),'.csv'),
                                 sep = ',',header = T)
  dea_res_FG_DDA_wf1<-dea_res_FG_DDA_wf1[which(dea_res_FG_DDA_wf1$contrast==contrast),]
  return(list(wf=wftop1_FG_DDA, dea_res=dea_res_FG_DDA_wf1))
}

top1wf_ens<-function(dt, platform, acq, ens_ty, contrast){
  base_fold<-data_fold
  if(ens_ty=='ensemble_mv'){
    dea_res_fold<-paste0(base_fold, acq, '/',platform, '/','mv','/',dt,'/')
  }else if(ens_ty=='ensemble_topk'){
    dea_res_fold<-paste0(base_fold, acq, '/',platform, '/','topk','/',dt,'/')
  }
  #dea_res_fold<-paste0(base_fold, acq, '/',platform, '/',dt,'/')
  ranking_all_FG_DDA<-read_excel(paste0(data_fold, 'ranks_all_', platform, '_', acq, '_', ens_ty, '.xlsx'), sheet = 'ranking_all')
  #ranking_all<-ranking_all[which(ranking_all$runTime!=0 & ranking_all$mean_pauc001!=0 & ranking_all$rank_median_pauc001!=0),]
  wftop1_FG_DDA<-ranking_all_FG_DDA[which(ranking_all_FG_DDA$avg_rank_mean==min(ranking_all_FG_DDA$avg_rank_mean)),]
  if(length((wftop1_FG_DDA$workflow))>1){
    idx=which(wftop1_FG_DDA$avg_rank_median==min(wftop1_FG_DDA$avg_rank_median))
    wftop1_FG_DDA = wftop1_FG_DDA[idx,]
  }
  if(wftop1_FG_DDA$method=='set'){
    dea_res_FG_DDA_wf1<-read.table(paste0(dea_res_fold,
                                          dt,'_',
                                          wftop1_FG_DDA$Platform,'_',
                                          wftop1_FG_DDA$method, '_',
                                          wftop1_FG_DDA$operation,
                                          '_',gsub('top','',gsub('\\|','-',wftop1_FG_DDA$cbn)),'.csv'),
                                   sep = ',',header = T)
  }else{
    if(acq=='TMT'){
      dea_res_FG_DDA_wf1<-read.table(paste0(dea_res_fold,
                                            dt,'_',
                                            wftop1_FG_DDA$Platform,'_',
                                            wftop1_FG_DDA$method,'_',gsub("top",'',gsub('\\|','-',wftop1_FG_DDA$cbn)),'.csv'),
                                     sep = ',',header = T)
    }else{
      dea_res_FG_DDA_wf1<-read.table(paste0(dea_res_fold,
                                            dt,'_',
                                            wftop1_FG_DDA$Platform,'_',
                                            wftop1_FG_DDA$method,'_',gsub('\\|','-',wftop1_FG_DDA$cbn),'.csv'),
                                     sep = ',',header = T)
    }

  }
  dea_res_FG_DDA_wf1<-dea_res_FG_DDA_wf1[which(dea_res_FG_DDA_wf1$contrast==contrast),]
  return(list(wf=wftop1_FG_DDA, dea_res=dea_res_FG_DDA_wf1))
}

dt='HYtims134'

contrast='conditionD-conditionB'

all_protein_FG_DDA<-read.table(paste0(data_fold, 'HYtims134_LFQ_FragPipe_all_proteins.tsv'), sep = '\t', header = T)
all_protein_MQ_DDA<-read.table(paste0(data_fold, 'HYtims134_LFQ_Maxquant_all_proteins.tsv'), sep = '\t', header = T)
all_protein_DIANN_DIA<-read.table(paste0(data_fold, 'HYtims134_DIA_DIANN_all_proteins.tsv'), sep = '\t', header = T)
all_protein_spt_DIA<-read.table(paste0(data_fold, 'HYtims134_DIA_spt_all_proteins.tsv'), sep = '\t', header = T)

top1_FG_DDA<-top1wf('HYtims134_LFQ', 'FragPipe', 'DDA',contrast)
top1_MQ_DDA<-top1wf('HYtims134_LFQ', 'Maxquant', 'DDA',contrast)
top1_DIANN_DIA<-top1wf('HYtims134_DIA', 'DIANN', 'DIA',contrast)
top1_spt_DIA<-top1wf('HYtims134_DIA', 'spt', 'DIA',contrast)

top1_FG_DDA_ens<-top1wf_ens('HYtims134_LFQ', 'FragPipe', 'DDA', 'ensemble_mv',contrast)
top1_MQ_DDA_ens<-top1wf_ens('HYtims134_LFQ', 'Maxquant', 'DDA', 'ensemble_mv',contrast)
top1_DIANN_DIA_ens<-top1wf_ens('HYtims134_DIA', 'DIANN', 'DIA', 'ensemble_mv',contrast)
top1_spt_DIA_ens<-top1wf_ens('HYtims134_DIA', 'spt', 'DIA', 'ensemble_mv',contrast)

FG_DDA<-top1_FG_DDA$dea_res$protein[which(abs(top1_FG_DDA$dea_res$logFC)>=log2(1.5) & top1_FG_DDA$dea_res$adj.pvalue<=0.05)]
MQ_DDA<-top1_MQ_DDA$dea_res$protein[which(abs(top1_MQ_DDA$dea_res$logFC)>=log2(1.5) & top1_MQ_DDA$dea_res$adj.pvalue<=0.05)]
DIANN_DIA<-top1_DIANN_DIA$dea_res$protein[which(abs(top1_DIANN_DIA$dea_res$logFC)>=log2(1.5) & top1_DIANN_DIA$dea_res$adj.pvalue<=0.05)]
spt_DIA<-top1_spt_DIA$dea_res$protein[which(abs(top1_spt_DIA$dea_res$logFC)>=log2(1.5) & top1_spt_DIA$dea_res$adj.pvalue<=0.05)]

FG_DDA_ens<-top1_FG_DDA_ens$dea_res$protein[which(abs(top1_FG_DDA_ens$dea_res$logFC)>=log2(1.5) & top1_FG_DDA_ens$dea_res$adj.pvalue<=0.05)]
MQ_DDA_ens<-top1_MQ_DDA_ens$dea_res$protein[which(abs(top1_MQ_DDA_ens$dea_res$logFC)>=log2(1.5) & top1_MQ_DDA_ens$dea_res$adj.pvalue<=0.05)]
DIANN_DIA_ens<-top1_DIANN_DIA_ens$dea_res$protein[which(abs(top1_DIANN_DIA_ens$dea_res$logFC)>=log2(1.5) & top1_DIANN_DIA_ens$dea_res$adj.pvalue<=0.05)]
spt_DIA_ens<-top1_spt_DIA_ens$dea_res$protein[which(abs(top1_spt_DIA_ens$dea_res$logFC)>=log2(1.5) & top1_spt_DIA_ens$dea_res$adj.pvalue<=0.05)]


T_YEAST<-unique(union(union(union(all_protein_FG_DDA$Protein[which(all_protein_FG_DDA$Organism=='YEAST')],
                                  all_protein_MQ_DDA$Protein[which(all_protein_MQ_DDA$Organism=='YEAST')]),
                            all_protein_DIANN_DIA$Protein[which(all_protein_DIANN_DIA$Organism=='YEAST')]),
                      all_protein_spt_DIA$Protein[which(all_protein_spt_DIA$Organism=='YEAST')]))
# T_ECOLI<-unique(union(union(union(all_protein_FG_DDA$Protein[which(all_protein_FG_DDA$Organism=='ECOLI')],
#                                   all_protein_MQ_DDA$Protein[which(all_protein_MQ_DDA$Organism=='ECOLI')]),
#                             all_protein_DIANN_DIA$Protein[which(all_protein_DIANN_DIA$Organism=='ECOLI')]),
#                       all_protein_spt_DIA$Protein[which(all_protein_spt_DIA$Organism=='ECOLI')]))
T_HUMAN<-unique(union(union(union(all_protein_FG_DDA$Protein[which(all_protein_FG_DDA$Organism=='HUMAN')],
                                  all_protein_MQ_DDA$Protein[which(all_protein_MQ_DDA$Organism=='HUMAN')]),
                            all_protein_DIANN_DIA$Protein[which(all_protein_DIANN_DIA$Organism=='HUMAN')]),
                      all_protein_spt_DIA$Protein[which(all_protein_spt_DIA$Organism=='HUMAN')]))

DEP_HUMAN<-unique(union(union(union(top1_FG_DDA$dea_res$protein[which(abs(top1_FG_DDA$dea_res$logFC)>=log2(1.5) &
                                                                        top1_FG_DDA$dea_res$adj.pvalue<=0.05 &
                                                                        top1_FG_DDA$dea_res$organism=='HUMAN')],
                                    top1_MQ_DDA$dea_res$protein[which(abs(top1_MQ_DDA$dea_res$logFC)>=log2(1.5) &
                                                                        top1_MQ_DDA$dea_res$adj.pvalue<=0.05 &
                                                                        top1_MQ_DDA$dea_res$organism=='HUMAN')]),
                              top1_DIANN_DIA$dea_res$protein[which(abs(top1_DIANN_DIA$dea_res$logFC)>=log2(1.5) &
                                                                     top1_DIANN_DIA$dea_res$adj.pvalue<=0.05 &
                                                                     top1_DIANN_DIA$dea_res$organism=='HUMAN')]),
                        top1_spt_DIA$dea_res$protein[which(abs(top1_spt_DIA$dea_res$logFC)>=log2(1.5) &
                                                             top1_spt_DIA$dea_res$adj.pvalue<=0.05 &
                                                             top1_spt_DIA$dea_res$organism=='HUMAN')]))


DEP_HUMAN_ens<-unique(union(union(union(top1_FG_DDA_ens$dea_res$protein[which(abs(top1_FG_DDA_ens$dea_res$logFC)>=log2(1.5) &
                                                                                top1_FG_DDA_ens$dea_res$adj.pvalue<=0.05 &
                                                                                top1_FG_DDA_ens$dea_res$organism=='HUMAN')],
                                        top1_MQ_DDA_ens$dea_res$protein[which(abs(top1_MQ_DDA_ens$dea_res$logFC)>=log2(1.5) &
                                                                                top1_MQ_DDA_ens$dea_res$adj.pvalue<=0.05 &
                                                                                top1_MQ_DDA_ens$dea_res$organism=='HUMAN')]),
                                  top1_DIANN_DIA_ens$dea_res$protein[which(abs(top1_DIANN_DIA_ens$dea_res$logFC)>=log2(1.5) &
                                                                             top1_DIANN_DIA_ens$dea_res$adj.pvalue<=0.05 &
                                                                             top1_DIANN_DIA_ens$dea_res$organism=='HUMAN')]),
                            top1_spt_DIA_ens$dea_res$protein[which(abs(top1_spt_DIA_ens$dea_res$logFC)>=log2(1.5) &
                                                                     top1_spt_DIA_ens$dea_res$adj.pvalue<=0.05 &
                                                                     top1_spt_DIA_ens$dea_res$organism=='HUMAN')]))

DEP_HUMAN<-unique(union(DEP_HUMAN, DEP_HUMAN_ens))

DEPs<-data.frame(Proteins=unique(union(union(union(all_protein_FG_DDA$Protein,
                                                   all_protein_MQ_DDA$Protein),
                                             all_protein_DIANN_DIA$Protein),
                                       all_protein_spt_DIA$Protein)))
DEPs$T_YEAST=c(rep(0, length(DEPs$Proteins)))
DEPs$T_YEAST[match(T_YEAST, DEPs$Proteins)]=1
# DEPs$T_ECOLI=c(rep(0, length(DEPs$Proteins)))
# DEPs$T_ECOLI[match(T_ECOLI, DEPs$Proteins)]=1
# DEPs$T_HUMAN=c(rep(0, length(DEPs$Proteins)))
# DEPs$T_HUMAN[match(T_HUMAN, DEPs$Proteins)]=1
DEPs$DEP_HUMAN=c(rep(0, length(DEPs$Proteins)))
DEPs$DEP_HUMAN[match(DEP_HUMAN, DEPs$Proteins)]=1
DEPs$FG_DDA=c(rep(0, length(DEPs$Proteins)))
DEPs$FG_DDA[match(FG_DDA, DEPs$Proteins)]=1
DEPs$MQ_DDA=c(rep(0, length(DEPs$Proteins)))
DEPs$MQ_DDA[match(MQ_DDA, DEPs$Proteins)]=1
DEPs$DIANN_DIA=c(rep(0, length(DEPs$Proteins)))
DEPs$DIANN_DIA[match(DIANN_DIA, DEPs$Proteins)]=1
DEPs$spt_DIA=c(rep(0, length(DEPs$Proteins)))
DEPs$spt_DIA[match(spt_DIA, DEPs$Proteins)]=1

DEPs$FG_DDA_ens=c(rep(0, length(DEPs$Proteins)))
DEPs$FG_DDA_ens[match(FG_DDA_ens, DEPs$Proteins)]=1
DEPs$MQ_DDA_ens=c(rep(0, length(DEPs$Proteins)))
DEPs$MQ_DDA_ens[match(MQ_DDA_ens, DEPs$Proteins)]=1
DEPs$DIANN_DIA_ens=c(rep(0, length(DEPs$Proteins)))
DEPs$DIANN_DIA_ens[match(DIANN_DIA_ens, DEPs$Proteins)]=1
DEPs$spt_DIA_ens=c(rep(0, length(DEPs$Proteins)))
DEPs$spt_DIA_ens[match(spt_DIA_ens, DEPs$Proteins)]=1

# simply remove protein groups
DEPs<-DEPs[setdiff(c(1:length(DEPs$Proteins)),grep(';', DEPs$Proteins)),]

sheet=paste(dt, '_cp_DEPs')

addWorksheet(wb,sheet)
writeData(wb, sheet, DEPs, startRow = 1, startCol = 1, colNames = TRUE, rowNames = TRUE)
saveWorkbook(wb, file = paste0(save_fold, 'cross_platform_compare.xlsx'), overwrite = TRUE)


# UpSetR::upset(DEPs,
#               sets = c("T_YEAST", "DEP_HUMAN", "FG_DDA", "MQ_DDA", "DIANN_DIA","spt_DIA","FG_DDA_ens","MQ_DDA_ens", "DIANN_DIA_ens","spt_DIA_ens"),
#               order.by="freq", matrix.color="black", point.size=3,
#               sets.bar.color=c("#15d08a","#f2ca01","#a45ee3","#1f11b4","#c3642c",
#                                "#f078b4", '#5494cc',"gray","pink",'red','orange'),mb.ratio=c(0.5, 0.5), text.scale = 1.5,
#               number.angles =0)

library(pROC)
library(ggpubr)

#all_pvals<-read.table('E:/proteomics/manu4_1/codes/reproduce_0_05/temp0_limma_min_FragPipe/HYEtims735_LFQ_top0_maxlfq_mv_outputs.csv',header = T, sep = ',')

cross_platform_proteins<-data.frame(Proteins=DEPs$Proteins, labels=0)
cross_platform_proteins$labels[grep('YEAST', cross_platform_proteins$Proteins)]=1
#cross_platform_proteins$labels[grep('ECOLI', cross_platform_proteins$Proteins)]=1
cross_platform_proteins$qval_FG_DDA<-1
cross_platform_proteins$logFC_FG_DDA<-0
com=intersect(cross_platform_proteins$Proteins, top1_FG_DDA$dea_res$protein)
idx1=match(com, cross_platform_proteins$Proteins)
idx2=match(com, top1_FG_DDA$dea_res$protein)
cross_platform_proteins$qval_FG_DDA[idx1]=top1_FG_DDA$dea_res$adj.pvalue[idx2]
cross_platform_proteins$logFC_FG_DDA[idx1]=top1_FG_DDA$dea_res$logFC[idx2]

cross_platform_proteins$qval_MQ_DDA<-1
cross_platform_proteins$logFC_MQ_DDA<-0
com=intersect(cross_platform_proteins$Proteins, top1_MQ_DDA$dea_res$protein)
idx1=match(com, cross_platform_proteins$Proteins)
idx2=match(com, top1_MQ_DDA$dea_res$protein)
cross_platform_proteins$qval_MQ_DDA[idx1]=top1_MQ_DDA$dea_res$adj.pvalue[idx2]
cross_platform_proteins$logFC_MQ_DDA[idx1]=top1_MQ_DDA$dea_res$logFC[idx2]

cross_platform_proteins$qval_DIANN_DIA<-1
cross_platform_proteins$logFC_DIANN_DIA<-0
com=intersect(cross_platform_proteins$Proteins, top1_DIANN_DIA$dea_res$protein)
idx1=match(com, cross_platform_proteins$Proteins)
idx2=match(com, top1_DIANN_DIA$dea_res$protein)
cross_platform_proteins$qval_DIANN_DIA[idx1]=top1_DIANN_DIA$dea_res$adj.pvalue[idx2]
cross_platform_proteins$logFC_DIANN_DIA[idx1]=top1_DIANN_DIA$dea_res$logFC[idx2]

cross_platform_proteins$qval_spt_DIA<-1
cross_platform_proteins$logFC_spt_DIA<-0
com=intersect(cross_platform_proteins$Proteins, top1_spt_DIA$dea_res$protein)
idx1=match(com, cross_platform_proteins$Proteins)
idx2=match(com, top1_spt_DIA$dea_res$protein)
cross_platform_proteins$qval_spt_DIA[idx1]=top1_spt_DIA$dea_res$adj.pvalue[idx2]
cross_platform_proteins$logFC_spt_DIA[idx1]=top1_spt_DIA$dea_res$logFC[idx2]

cross_platform_proteins$qval_FG_DDA_ens<-1
cross_platform_proteins$logFC_FG_DDA_ens<-0
com=intersect(cross_platform_proteins$Proteins, top1_FG_DDA_ens$dea_res$protein)
idx1=match(com, cross_platform_proteins$Proteins)
idx2=match(com, top1_FG_DDA_ens$dea_res$protein)
cross_platform_proteins$qval_FG_DDA_ens[idx1]=top1_FG_DDA_ens$dea_res$adj.pvalue[idx2]
cross_platform_proteins$logFC_FG_DDA_ens[idx1]=top1_FG_DDA_ens$dea_res$logFC[idx2]

cross_platform_proteins$qval_MQ_DDA_ens<-1
cross_platform_proteins$logFC_MQ_DDA_ens<-0
com=intersect(cross_platform_proteins$Proteins, top1_MQ_DDA_ens$dea_res$protein)
idx1=match(com, cross_platform_proteins$Proteins)
idx2=match(com, top1_MQ_DDA_ens$dea_res$protein)
cross_platform_proteins$qval_MQ_DDA_ens[idx1]=top1_MQ_DDA_ens$dea_res$adj.pvalue[idx2]
cross_platform_proteins$logFC_MQ_DDA_ens[idx1]=top1_MQ_DDA_ens$dea_res$logFC[idx2]

cross_platform_proteins$qval_DIANN_DIA_ens<-1
cross_platform_proteins$logFC_DIANN_DIA_ens<-0
com=intersect(cross_platform_proteins$Proteins, top1_DIANN_DIA_ens$dea_res$protein)
idx1=match(com, cross_platform_proteins$Proteins)
idx2=match(com, top1_DIANN_DIA_ens$dea_res$protein)
cross_platform_proteins$qval_DIANN_DIA_ens[idx1]=top1_DIANN_DIA_ens$dea_res$adj.pvalue[idx2]
cross_platform_proteins$logFC_DIANN_DIA_ens[idx1]=top1_DIANN_DIA_ens$dea_res$logFC[idx2]

cross_platform_proteins$qval_spt_DIA_ens<-1
cross_platform_proteins$logFC_spt_DIA_ens<-0
com=intersect(cross_platform_proteins$Proteins, top1_spt_DIA_ens$dea_res$protein)
idx1=match(com, cross_platform_proteins$Proteins)
idx2=match(com, top1_spt_DIA_ens$dea_res$protein)
cross_platform_proteins$qval_spt_DIA_ens[idx1]=top1_spt_DIA_ens$dea_res$adj.pvalue[idx2]
cross_platform_proteins$logFC_spt_DIA_ens[idx1]=top1_spt_DIA_ens$dea_res$logFC[idx2]

sheet=paste(dt, '_cp_proteins')

addWorksheet(wb,sheet)
writeData(wb, sheet, cross_platform_proteins, startRow = 1, startCol = 1, colNames = TRUE, rowNames = TRUE)
saveWorkbook(wb, file = paste0(save_fold, 'cross_platform_compare.xlsx'), overwrite = TRUE)

roc_FG_DDA <- roc(cross_platform_proteins$labels, 1-cross_platform_proteins$qval_FG_DDA)
roc_MQ_DDA <- roc(cross_platform_proteins$labels, 1-cross_platform_proteins$qval_MQ_DDA)
roc_DIANN_DIA <- roc(cross_platform_proteins$labels, 1-cross_platform_proteins$qval_DIANN_DIA)
roc_spt_DIA <- roc(cross_platform_proteins$labels, 1-cross_platform_proteins$qval_spt_DIA)

roc_FG_DDA_ens <- roc(cross_platform_proteins$labels, 1-cross_platform_proteins$qval_FG_DDA_ens)
roc_MQ_DDA_ens <- roc(cross_platform_proteins$labels, 1-cross_platform_proteins$qval_MQ_DDA_ens)
roc_DIANN_DIA_ens <- roc(cross_platform_proteins$labels, 1-cross_platform_proteins$qval_DIANN_DIA_ens)
roc_spt_DIA_ens <- roc(cross_platform_proteins$labels, 1-cross_platform_proteins$qval_spt_DIA_ens)


cross_platform_metrics<-function(cross_platform_proteins, setting, roc){
  lab=cross_platform_proteins$labels
  if(setting=='FG_DDA'){
    logFCs=cross_platform_proteins$logFC_FG_DDA
    qvals=cross_platform_proteins$qval_FG_DDA
    roc_s<-roc_FG_DDA
  }else if(setting=='MQ_DDA'){
    logFCs=cross_platform_proteins$logFC_MQ_DDA
    qvals=cross_platform_proteins$qval_MQ_DDA
    roc_s<-roc_MQ_DDA
  }else if(setting=='DIANN_DIA'){
    logFCs=cross_platform_proteins$logFC_DIANN_DIA
    qvals=cross_platform_proteins$qval_DIANN_DIA
    roc_s<-roc_DIANN_DIA
  }else if(setting=='spt_DIA'){
    logFCs=cross_platform_proteins$logFC_spt_DIA
    qvals=cross_platform_proteins$qval_spt_DIA
    roc_s<-roc_spt_DIA
  }else if(setting=='FG_DDA_ens'){
    logFCs=cross_platform_proteins$logFC_FG_DDA_ens
    qvals=cross_platform_proteins$qval_FG_DDA_ens
    roc_s<-roc_FG_DDA_ens
  }else if(setting=='MQ_DDA_ens'){
    logFCs=cross_platform_proteins$logFC_MQ_DDA_ens
    qvals=cross_platform_proteins$qval_MQ_DDA_ens
    roc_s<-roc_MQ_DDA_ens
  }else if(setting=='DIANN_DIA_ens'){
    logFCs=cross_platform_proteins$logFC_DIANN_DIA_ens
    qvals=cross_platform_proteins$qval_DIANN_DIA_ens
    roc_s<-roc_DIANN_DIA_ens
  }else if(setting=='spt_DIA_ens'){
    logFCs=cross_platform_proteins$logFC_spt_DIA_ens
    qvals=cross_platform_proteins$qval_spt_DIA_ens
    roc_s<-roc_spt_DIA_ens
  }
  TP=length(which(lab==1 & abs(logFCs)>=log2(1.5) & qvals<=0.05))
  TN=length(which(lab==0 & (abs(logFCs)<=log2(1.5) | qvals>0.05)))
  FP=length(which(lab==0 & abs(logFCs)>=log2(1.5) & qvals<=0.05))
  FN=length(which(lab==1 & (abs(logFCs)<=log2(1.5) | qvals>0.05)))

  Recall = TP/(TP+FN)
  Precision = TP/(TP+FP)
  Specificity = TN/(TN+FP)
  F1=2*Recall*Precision/(Recall+Precision)
  MCC = (TP*TN-FP*FN)/(((TP+FP)^0.5)*((TP+FN)^0.5)*((TN+FP)^0.5)*((TN+FN))^0.5)
  Gmean = (Recall*Specificity)^0.5
  nMCC = (1+MCC)/2

  pauc001=as.numeric(auc(roc_s, partial.auc=c(1, 0.99), partial.auc.correct=TRUE))
  pauc005=as.numeric(auc(roc_s, partial.auc=c(1, 0.95), partial.auc.correct=TRUE))
  pauc01=as.numeric(auc(roc_s, partial.auc=c(1, 0.90), partial.auc.correct=TRUE))
  metrics=data.frame(metric=c('pAUC(0.01)','pAUC(0.05)','pAUC(0.1)','TP', 'TN', 'FP', 'FN', 'Recall', 'Precision', 'Specificity', 'F1', 'MCC', 'G-mean','nMCC'),
                     value=c(pauc001,pauc005,pauc01,TP, TN, FP, FN, Recall, Precision, Specificity, F1, MCC, Gmean,nMCC))
  return(metrics)
}

metrics_FG_DDA=cross_platform_metrics(cross_platform_proteins, 'FG_DDA', roc_FG_DDA)
metrics_MQ_DDA=cross_platform_metrics(cross_platform_proteins, 'MQ_DDA', roc_MQ_DDA)
metrics_DIANN_DIA=cross_platform_metrics(cross_platform_proteins, 'DIANN_DIA', roc_DIANN_DIA)
metrics_spt_DIA=cross_platform_metrics(cross_platform_proteins, 'spt_DIA', roc_spt_DIA)

metrics_FG_DDA_ens=cross_platform_metrics(cross_platform_proteins, 'FG_DDA_ens', roc_FG_DDA_ens)
metrics_MQ_DDA_ens=cross_platform_metrics(cross_platform_proteins, 'MQ_DDA_ens', roc_MQ_DDA_ens)
metrics_DIANN_DIA_ens=cross_platform_metrics(cross_platform_proteins, 'DIANN_DIA_ens', roc_DIANN_DIA_ens)
metrics_spt_DIA_ens=cross_platform_metrics(cross_platform_proteins, 'spt_DIA_ens', roc_spt_DIA_ens)

rdar_dt<-rbind(metrics_FG_DDA$value, metrics_MQ_DDA$value, metrics_DIANN_DIA$value, metrics_spt_DIA$value,
               metrics_FG_DDA_ens$value, metrics_MQ_DDA_ens$value, metrics_DIANN_DIA_ens$value, metrics_spt_DIA_ens$value)
colnames(rdar_dt)<-metrics_FG_DDA$metric
row.names(rdar_dt)<-c('FG_DDA','MQ_DDA', 'DIANN_DIA', 'spt_DIA', 'FG_DDA_ens','MQ_DDA_ens', 'DIANN_DIA_ens', 'spt_DIA_ens')
rdar_dt<-as.data.frame(as.matrix(rdar_dt))
pt_rdar_dt<-data.frame(group=row.names(rdar_dt))

sheet=paste(dt, '_cp_metrics')

addWorksheet(wb,sheet)
writeData(wb, sheet, rdar_dt, startRow = 1, startCol = 1, colNames = TRUE, rowNames = TRUE)
saveWorkbook(wb, file = paste0(save_fold, 'cross_platform_compare.xlsx'), overwrite = TRUE)
pt_rdar_dt<-cbind(pt_rdar_dt, rdar_dt[,c(1,2,3,13,14)]) #,4,6
#pt_rdar_dt$group<-row.names(pt_rdar_dt)
library(ggradar)
ggradar(pt_rdar_dt,group.colours=c("#a45ee3","#1f11b4","#f078b4","#f2ca01","gray","pink",'red',"orange"),legend.position = 'right')
#a=c('#5494cc','#e18283','#0d898a','#f9cc52')
col=c("#15d08a","#f2ca01","#a45ee3","#1f11b4","#c3642c",
      "#f078b4", '#5494cc',"gray","pink",'red','orange')
melt_dt<-melt(pt_rdar_dt)
colnames(melt_dt)[1]<-'method'

ggdotchart(melt_dt, "variable", "value", group = "method", color = "method",
           palette = col[c(6,1,8,3,10,7,5,4)],

           add = "segments", #label = "value",
           add.params = list(color = "lightgray", size = 1),
           dot.size = 3,
           #label = 'value',
           #font.label = list(color = "black", size = 9, vjust = 0.5),
           position=position_dodge(width = 0.9)
) +
  #geom_text(aes(label=value), position = 'indentity')+
  theme(axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
        axis.title.x = element_blank())#+#+scale_y_continuous(limits = c(0.2, 1))+
#coord_flip()


count_TPs<-data.frame(method=row.names(rdar_dt))
count_TPs<-cbind(count_TPs, rdar_dt[,c(4:7)]) #,4,6
melt_dt<-melt(count_TPs)
melt_dt$label<-as.character(melt_dt$value)
melt_dt$label[which(melt_dt$variable=='TN' | melt_dt$variable=='FN')]=''

nums<-vector()


for (i in 1:length(unique(melt_dt$method))){
  num_i<-melt_dt[which(melt_dt$method==melt_dt$method[i]),]
  uni_var<-levels(factor(melt_dt$variable))
  re_order<-match(num_i$variable, uni_var)
  num_i<-num_i[re_order,]
  l_y<-c(rep(0,length(uni_var)))
  for (j in length(uni_var):1) {
    if(j==length(uni_var)){
      l_y[j]=num_i$value[j]/2
    }else{
      l_y_j=num_i$value[j]/2 + sum(num_i$value[c((j+1):length(uni_var))])
      l_y[j]=l_y_j
    }
  }
  num_i$l_y<-l_y
  nums<-rbind(nums,num_i)
}
## supp. Fig. 6
ggbarplot(nums, "method", "value", group = "variable", color = "variable",
          fill="variable", #label = T,lab.pos='in',
          platette = c('#5494cc','#e18283','#0d898a','#f9cc52'),
) + geom_text(aes(y=l_y, label=label))+
  scale_x_discrete(limits=c('FG_DDA_ens','MQ_DDA_ens','DIANN_DIA_ens', 'spt_DIA_ens',
                            'FG_DDA','MQ_DDA','DIANN_DIA','spt_DIA')) +
  geom_vline(xintercept = 4.5, linetype=2, show.legend=FALSE) + theme(axis.line.y = element_blank())+
  coord_flip()+theme(axis.title.y = element_blank())

############
## cross instrument
library("readxl")
library(ggplot2)
library(reshape2)
library(ggpubr)
library(cowplot)

library(openxlsx)
wb <- createWorkbook()

top1wf<-function(dt, platform, acq){
  base_fold<-data_fold
  dea_res_fold<-paste0(base_fold, acq, '/',platform, '/',dt,'/')
  ranking_all_FG_DDA<-read_excel(paste0(data_fold, 'ranks_all_', platform, '_', acq, '.xlsx'), sheet = 'ranking_all')

  #ranking_all<-ranking_all[which(ranking_all$runTime!=0 & ranking_all$mean_pauc001!=0 & ranking_all$rank_median_pauc001!=0),]
  wftop1_FG_DDA<-ranking_all_FG_DDA[which(ranking_all_FG_DDA$avg_rank_mean==min(ranking_all_FG_DDA$avg_rank_mean)),]
  if(length((wftop1_FG_DDA$workflow))>1){
    idx=which(wftop1_FG_DDA$avg_rank_median==min(wftop1_FG_DDA$avg_rank_median))
    wftop1_FG_DDA = wftop1_FG_DDA[idx,]
  }
  dea_res_FG_DDA_wf1<-read.table(paste0(dea_res_fold,
                                        dt,'_',wftop1_FG_DDA$DEA,'_',
                                        wftop1_FG_DDA$Platform,'_',
                                        wftop1_FG_DDA$Matrix,'_',gsub('blank','',wftop1_FG_DDA$Imput),
                                        '_',gsub('blank','',wftop1_FG_DDA$normalization),'.csv'),
                                 sep = ',',header = T)
  return(list(wf=wftop1_FG_DDA, dea_res=dea_res_FG_DDA_wf1))
}

top1wf_ens<-function(dt, platform, acq, ens_ty){
  base_fold<-data_fold
  if(ens_ty=='ensemble_mv'){
    dea_res_fold<-paste0(base_fold, acq, '/',platform, '/','mv','/',dt,'/')
  }else if(ens_ty=='ensemble_topk'){
    dea_res_fold<-paste0(base_fold, acq, '/',platform, '/','topk','/',dt,'/')
  }
  #dea_res_fold<-paste0(base_fold, acq, '/',platform, '/',dt,'/')
  ranking_all_FG_DDA<-read_excel(paste0(data_fold, 'ranks_all_', platform, '_', acq, '_', ens_ty, '.xlsx'), sheet = 'ranking_all')
  #ranking_all<-ranking_all[which(ranking_all$runTime!=0 & ranking_all$mean_pauc001!=0 & ranking_all$rank_median_pauc001!=0),]
  wftop1_FG_DDA<-ranking_all_FG_DDA[which(ranking_all_FG_DDA$avg_rank_mean==min(ranking_all_FG_DDA$avg_rank_mean)),]
  if(length((wftop1_FG_DDA$workflow))>1){
    idx=which(wftop1_FG_DDA$avg_rank_median==min(wftop1_FG_DDA$avg_rank_median))
    wftop1_FG_DDA = wftop1_FG_DDA[idx,]
  }
  if(wftop1_FG_DDA$method=='set'){
    dea_res_FG_DDA_wf1<-read.table(paste0(dea_res_fold,
                                          dt,'_',
                                          wftop1_FG_DDA$Platform,'_',
                                          wftop1_FG_DDA$method, '_',
                                          wftop1_FG_DDA$operation,
                                          '_',gsub('\\|','-',wftop1_FG_DDA$cbn),'.csv'),
                                   sep = ',',header = T)
  }else{
    dea_res_FG_DDA_wf1<-read.table(paste0(dea_res_fold,
                                          dt,'_',
                                          wftop1_FG_DDA$Platform,'_',
                                          wftop1_FG_DDA$method,'_',gsub('\\|','-',wftop1_FG_DDA$cbn),'.csv'),
                                   sep = ',',header = T)
  }

  return(list(wf=wftop1_FG_DDA, dea_res=dea_res_FG_DDA_wf1))
}

dt='HYE'

all_protein_FG_DDA_tims<-read.table(paste0(data_fold,'HYEtims735_LFQ_FragPipe_all_proteins.tsv'), sep = '\t', header = T)
all_protein_MQ_DDA_tims<-read.table(paste0(data_fold,'HYEtims735_LFQ_Maxquant_all_proteins.tsv'), sep = '\t', header = T)
all_protein_DIANN_DIA_tims<-read.table(paste0(data_fold,'HYEtims735_DIA_DIANN_all_proteins.tsv'), sep = '\t', header = T)
all_protein_spt_DIA_tims<-read.table(paste0(data_fold,'HYEtims735_DIA_spt_all_proteins.tsv'), sep = '\t', header = T)
all_protein_FG_DDA_st5600<-read.table(paste0(data_fold,'HYE5600735_LFQ_FragPipe_all_proteins.tsv'), sep = '\t', header = T)
all_protein_MQ_DDA_st5600<-read.table(paste0(data_fold,'HYE5600735_LFQ_Maxquant_all_proteins.tsv'), sep = '\t', header = T)
all_protein_FG_DDA_st6600<-read.table(paste0(data_fold,'HYE6600735_LFQ_FragPipe_all_proteins.tsv'), sep = '\t', header = T)
all_protein_MQ_DDA_st6600<-read.table(paste0(data_fold,'HYE6600735_LFQ_Maxquant_all_proteins.tsv'), sep = '\t', header = T)
all_protein_FG_DDA_qe<-read.table(paste0(data_fold,'HYEqe735_LFQ_FragPipe_all_proteins.tsv'), sep = '\t', header = T)
all_protein_MQ_DDA_qe<-read.table(paste0(data_fold,'HYEqe735_LFQ_Maxquant_all_proteins.tsv'), sep = '\t', header = T)



top1_FG_DDA_tims<-top1wf('HYEtims735_LFQ', 'FragPipe', 'DDA')
top1_MQ_DDA_tims<-top1wf('HYEtims735_LFQ', 'Maxquant', 'DDA')
top1_DIANN_DIA<-top1wf('HYEtims735_DIA', 'DIANN', 'DIA')
top1_spt_DIA<-top1wf('HYEtims735_DIA', 'spt', 'DIA')

top1_FG_DDA_st5600<-top1wf('HYE5600735_LFQ', 'FragPipe', 'DDA')
top1_MQ_DDA_st5600<-top1wf('HYE5600735_LFQ', 'Maxquant', 'DDA')
top1_FG_DDA_st6600<-top1wf('HYE6600735_LFQ', 'FragPipe', 'DDA')
top1_MQ_DDA_st6600<-top1wf('HYE6600735_LFQ', 'Maxquant', 'DDA')
top1_FG_DDA_qe<-top1wf('HYEqe735_LFQ', 'FragPipe', 'DDA')
top1_MQ_DDA_qe<-top1wf('HYEqe735_LFQ', 'Maxquant', 'DDA')

top1_FG_DDA_tims_ens=top1wf_ens('HYEtims735_LFQ', 'FragPipe', 'DDA', 'ensemble_mv')
top1_MQ_DDA_tims_ens=top1wf_ens('HYEtims735_LFQ', 'Maxquant', 'DDA', 'ensemble_mv')
top1_DIANN_DIA_tims_ens=top1wf_ens('HYEtims735_DIA', 'DIANN', 'DIA', 'ensemble_mv')
top1_spt_DIA_tims_ens=top1wf_ens('HYEtims735_DIA', 'spt', 'DIA', 'ensemble_mv')

top1_FG_DDA_st5600_ens<-top1wf_ens('HYE5600735_LFQ', 'FragPipe', 'DDA','ensemble_mv')
top1_MQ_DDA_st5600_ens<-top1wf_ens('HYE5600735_LFQ', 'Maxquant', 'DDA','ensemble_mv')
top1_FG_DDA_st6600_ens<-top1wf_ens('HYE6600735_LFQ', 'FragPipe', 'DDA','ensemble_mv')
top1_MQ_DDA_st6600_ens<-top1wf_ens('HYE6600735_LFQ', 'Maxquant', 'DDA','ensemble_mv')
top1_FG_DDA_qe_ens<-top1wf_ens('HYEqe735_LFQ', 'FragPipe', 'DDA','ensemble_mv')
top1_MQ_DDA_qe_ens<-top1wf_ens('HYEqe735_LFQ', 'Maxquant', 'DDA','ensemble_mv')


FG_DDA_tims<-top1_FG_DDA_tims$dea_res$protein[which(abs(top1_FG_DDA_tims$dea_res$logFC)>=log2(1.5) & top1_FG_DDA_tims$dea_res$adj.pvalue<=0.05)]
MQ_DDA_tims<-top1_MQ_DDA_tims$dea_res$protein[which(abs(top1_MQ_DDA_tims$dea_res$logFC)>=log2(1.5) & top1_MQ_DDA_tims$dea_res$adj.pvalue<=0.05)]
DIANN_DIA<-top1_DIANN_DIA$dea_res$protein[which(abs(top1_DIANN_DIA$dea_res$logFC)>=log2(1.5) & top1_DIANN_DIA$dea_res$adj.pvalue<=0.05)]
spt_DIA<-top1_spt_DIA$dea_res$protein[which(abs(top1_spt_DIA$dea_res$logFC)>=log2(1.5) & top1_spt_DIA$dea_res$adj.pvalue<=0.05)]
FG_DDA_tt5600<-top1_FG_DDA_st5600$dea_res$protein[which(abs(top1_FG_DDA_st5600$dea_res$logFC)>=log2(1.5) & top1_FG_DDA_st5600$dea_res$adj.pvalue<=0.05)]
MQ_DDA_tt5600<-top1_MQ_DDA_st5600$dea_res$protein[which(abs(top1_MQ_DDA_st5600$dea_res$logFC)>=log2(1.5) & top1_MQ_DDA_st5600$dea_res$adj.pvalue<=0.05)]
FG_DDA_tt6600<-top1_FG_DDA_st6600$dea_res$protein[which(abs(top1_FG_DDA_st6600$dea_res$logFC)>=log2(1.5) & top1_FG_DDA_st6600$dea_res$adj.pvalue<=0.05)]
MQ_DDA_tt6600<-top1_MQ_DDA_st6600$dea_res$protein[which(abs(top1_MQ_DDA_st6600$dea_res$logFC)>=log2(1.5) & top1_MQ_DDA_st6600$dea_res$adj.pvalue<=0.05)]
FG_DDA_qe<-top1_FG_DDA_qe$dea_res$protein[which(abs(top1_FG_DDA_qe$dea_res$logFC)>=log2(1.5) & top1_FG_DDA_qe$dea_res$adj.pvalue<=0.05)]
MQ_DDA_qe<-top1_MQ_DDA_qe$dea_res$protein[which(abs(top1_MQ_DDA_qe$dea_res$logFC)>=log2(1.5) & top1_MQ_DDA_qe$dea_res$adj.pvalue<=0.05)]

FG_DDA_tims_ens<-top1_FG_DDA_tims_ens$dea_res$protein[which(abs(top1_FG_DDA_tims_ens$dea_res$logFC)>=log2(1.5) & top1_FG_DDA_tims_ens$dea_res$adj.pvalue<=0.05)]
MQ_DDA_tims_ens<-top1_MQ_DDA_tims_ens$dea_res$protein[which(abs(top1_MQ_DDA_tims_ens$dea_res$logFC)>=log2(1.5) & top1_MQ_DDA_tims_ens$dea_res$adj.pvalue<=0.05)]
DIANN_DIA_ens<-top1_DIANN_DIA_tims_ens$dea_res$protein[which(abs(top1_DIANN_DIA_tims_ens$dea_res$logFC)>=log2(1.5) & top1_DIANN_DIA_tims_ens$dea_res$adj.pvalue<=0.05)]
spt_DIA_ens<-top1_spt_DIA_tims_ens$dea_res$protein[which(abs(top1_spt_DIA_tims_ens$dea_res$logFC)>=log2(1.5) & top1_spt_DIA_tims_ens$dea_res$adj.pvalue<=0.05)]
FG_DDA_tt5600_ens<-top1_FG_DDA_st5600_ens$dea_res$protein[which(abs(top1_FG_DDA_st5600_ens$dea_res$logFC)>=log2(1.5) & top1_FG_DDA_st5600_ens$dea_res$adj.pvalue<=0.05)]
MQ_DDA_tt5600_ens<-top1_MQ_DDA_st5600_ens$dea_res$protein[which(abs(top1_MQ_DDA_st5600_ens$dea_res$logFC)>=log2(1.5) & top1_MQ_DDA_st5600_ens$dea_res$adj.pvalue<=0.05)]
FG_DDA_tt6600_ens<-top1_FG_DDA_st6600_ens$dea_res$protein[which(abs(top1_FG_DDA_st6600_ens$dea_res$logFC)>=log2(1.5) & top1_FG_DDA_st6600_ens$dea_res$adj.pvalue<=0.05)]
MQ_DDA_tt6600_ens<-top1_MQ_DDA_st6600_ens$dea_res$protein[which(abs(top1_MQ_DDA_st6600_ens$dea_res$logFC)>=log2(1.5) & top1_MQ_DDA_st6600_ens$dea_res$adj.pvalue<=0.05)]
FG_DDA_qe_ens<-top1_FG_DDA_qe_ens$dea_res$protein[which(abs(top1_FG_DDA_qe_ens$dea_res$logFC)>=log2(1.5) & top1_FG_DDA_qe_ens$dea_res$adj.pvalue<=0.05)]
MQ_DDA_qe_ens<-top1_MQ_DDA_qe_ens$dea_res$protein[which(abs(top1_MQ_DDA_qe_ens$dea_res$logFC)>=log2(1.5) & top1_MQ_DDA_qe_ens$dea_res$adj.pvalue<=0.05)]



T_YEAST<-unique(union(union(union(union(union(union(union(union(union(
  all_protein_FG_DDA_tims$Protein[which(all_protein_FG_DDA_tims$Organism=='YEAST')],
                                  all_protein_MQ_DDA_tims$Protein[which(all_protein_MQ_DDA_tims$Organism=='YEAST')]),
                            all_protein_DIANN_DIA_tims$Protein[which(all_protein_DIANN_DIA_tims$Organism=='YEAST')]),
                      all_protein_spt_DIA_tims$Protein[which(all_protein_spt_DIA_tims$Organism=='YEAST')]),
                      all_protein_FG_DDA_st5600$Protein[which(all_protein_FG_DDA_st5600$Organism=='YEAST')]),
                      all_protein_FG_DDA_st6600$Protein[which(all_protein_FG_DDA_st6600$Organism=='YEAST')]),
                      all_protein_FG_DDA_qe$Protein[which(all_protein_FG_DDA_qe$Organism=='YEAST')]),
                all_protein_MQ_DDA_st5600$Protein[which(all_protein_MQ_DDA_st5600$Organism=='YEAST')]),
                all_protein_MQ_DDA_st6600$Protein[which(all_protein_MQ_DDA_st6600$Organism=='YEAST')]),
                all_protein_MQ_DDA_qe$Protein[which(all_protein_MQ_DDA_qe$Organism=='YEAST')]))



T_ECOLI<-unique(union(union(union(union(union(union(union(union(union(
  all_protein_FG_DDA_tims$Protein[which(all_protein_FG_DDA_tims$Organism=='ECOLI')],
                                  all_protein_MQ_DDA_tims$Protein[which(all_protein_MQ_DDA_tims$Organism=='ECOLI')]),
                            all_protein_DIANN_DIA_tims$Protein[which(all_protein_DIANN_DIA_tims$Organism=='ECOLI')]),
                      all_protein_spt_DIA_tims$Protein[which(all_protein_spt_DIA_tims$Organism=='ECOLI')]),
                all_protein_FG_DDA_st5600$Protein[which(all_protein_FG_DDA_st5600$Organism=='ECOLI')]),
all_protein_FG_DDA_st6600$Protein[which(all_protein_FG_DDA_st6600$Organism=='ECOLI')]),
all_protein_FG_DDA_qe$Protein[which(all_protein_FG_DDA_qe$Organism=='ECOLI')]),
all_protein_MQ_DDA_st5600$Protein[which(all_protein_MQ_DDA_st5600$Organism=='ECOLI')]),
all_protein_MQ_DDA_st6600$Protein[which(all_protein_MQ_DDA_st6600$Organism=='ECOLI')]),
all_protein_MQ_DDA_qe$Protein[which(all_protein_MQ_DDA_qe$Organism=='ECOLI')]))

T_HUMAN<-unique(union(union(union(union(union(union(union(union(union(
  all_protein_FG_DDA_tims$Protein[which(all_protein_FG_DDA_tims$Organism=='HUMAN')],
  all_protein_MQ_DDA_tims$Protein[which(all_protein_MQ_DDA_tims$Organism=='HUMAN')]),
  all_protein_DIANN_DIA_tims$Protein[which(all_protein_DIANN_DIA_tims$Organism=='HUMAN')]),
  all_protein_spt_DIA_tims$Protein[which(all_protein_spt_DIA_tims$Organism=='HUMAN')]),
  all_protein_FG_DDA_st5600$Protein[which(all_protein_FG_DDA_st5600$Organism=='HUMAN')]),
  all_protein_FG_DDA_st6600$Protein[which(all_protein_FG_DDA_st6600$Organism=='HUMAN')]),
  all_protein_FG_DDA_qe$Protein[which(all_protein_FG_DDA_qe$Organism=='HUMAN')]),
  all_protein_MQ_DDA_st5600$Protein[which(all_protein_MQ_DDA_st5600$Organism=='HUMAN')]),
  all_protein_MQ_DDA_st6600$Protein[which(all_protein_MQ_DDA_st6600$Organism=='HUMAN')]),
  all_protein_MQ_DDA_qe$Protein[which(all_protein_MQ_DDA_qe$Organism=='HUMAN')]))

DEPs<-data.frame(Proteins=unique(union(union(union(union(union(union(union(union(union(
  all_protein_FG_DDA_tims$Protein, all_protein_MQ_DDA_tims$Protein),
                                             all_protein_DIANN_DIA_tims$Protein),
                                       all_protein_spt_DIA_tims$Protein),
                                 all_protein_FG_DDA_st5600$Protein),
                 all_protein_FG_DDA_st6600$Protein),
                 all_protein_FG_DDA_qe$Protein),
  all_protein_MQ_DDA_st5600$Protein),
  all_protein_MQ_DDA_st6600$Protein),
  all_protein_MQ_DDA_qe$Protein)))

DEPs$T_YEAST=c(rep(0, length(DEPs$Proteins)))
DEPs$T_YEAST[match(T_YEAST, DEPs$Proteins)]=1
DEPs$T_ECOLI=c(rep(0, length(DEPs$Proteins)))
DEPs$T_ECOLI[match(T_ECOLI, DEPs$Proteins)]=1
DEPs$T_HUMAN=c(rep(0, length(DEPs$Proteins)))
DEPs$T_HUMAN[match(T_HUMAN, DEPs$Proteins)]=1

# simply remove protein groups
DEPs<-DEPs[setdiff(c(1:length(DEPs$Proteins)),grep(';', DEPs$Proteins)),]

sheet=paste(dt, '_ci_DEPs')

addWorksheet(wb,sheet)
writeData(wb, sheet, DEPs, startRow = 1, startCol = 1, colNames = TRUE, rowNames = TRUE)
saveWorkbook(wb, file = paste0(save_fold, 'cross_instrument_compare.xlsx'), overwrite = TRUE)


library(pROC)
library(ggpubr)

#all_pvals<-read.table('E:/proteomics/manu4_1/codes/reproduce_0_05/temp0_limma_min_FragPipe/HYEtims735_LFQ_top0_maxlfq_mv_outputs.csv',header = T, sep = ',')

cross_platform_proteins<-data.frame(Proteins=DEPs$Proteins, labels=0)
cross_platform_proteins$labels[grep('YEAST', cross_platform_proteins$Proteins)]=1
cross_platform_proteins$labels[grep('ECOLI', cross_platform_proteins$Proteins)]=1
cross_platform_proteins$qval_FG_DDA_tims<-1
cross_platform_proteins$logFC_FG_DDA_tims<-0
com=intersect(cross_platform_proteins$Proteins, top1_FG_DDA_tims$dea_res$protein)
idx1=match(com, cross_platform_proteins$Proteins)
idx2=match(com, top1_FG_DDA_tims$dea_res$protein)
cross_platform_proteins$qval_FG_DDA_tims[idx1]=top1_FG_DDA_tims$dea_res$adj.pvalue[idx2]
cross_platform_proteins$logFC_FG_DDA_tims[idx1]=top1_FG_DDA_tims$dea_res$logFC[idx2]

cross_platform_proteins$qval_MQ_DDA_tims<-1
cross_platform_proteins$logFC_MQ_DDA_tims<-0
com=intersect(cross_platform_proteins$Proteins, top1_MQ_DDA_tims$dea_res$protein)
idx1=match(com, cross_platform_proteins$Proteins)
idx2=match(com, top1_MQ_DDA_tims$dea_res$protein)
cross_platform_proteins$qval_MQ_DDA_tims[idx1]=top1_MQ_DDA_tims$dea_res$adj.pvalue[idx2]
cross_platform_proteins$logFC_MQ_DDA_tims[idx1]=top1_MQ_DDA_tims$dea_res$logFC[idx2]

cross_platform_proteins$qval_DIANN_DIA<-1
cross_platform_proteins$logFC_DIANN_DIA<-0
com=intersect(cross_platform_proteins$Proteins, top1_DIANN_DIA$dea_res$protein)
idx1=match(com, cross_platform_proteins$Proteins)
idx2=match(com, top1_DIANN_DIA$dea_res$protein)
cross_platform_proteins$qval_DIANN_DIA[idx1]=top1_DIANN_DIA$dea_res$adj.pvalue[idx2]
cross_platform_proteins$logFC_DIANN_DIA[idx1]=top1_DIANN_DIA$dea_res$logFC[idx2]

cross_platform_proteins$qval_spt_DIA<-1
cross_platform_proteins$logFC_spt_DIA<-0
com=intersect(cross_platform_proteins$Proteins, top1_spt_DIA$dea_res$protein)
idx1=match(com, cross_platform_proteins$Proteins)
idx2=match(com, top1_spt_DIA$dea_res$protein)
cross_platform_proteins$qval_spt_DIA[idx1]=top1_spt_DIA$dea_res$adj.pvalue[idx2]
cross_platform_proteins$logFC_spt_DIA[idx1]=top1_spt_DIA$dea_res$logFC[idx2]

cross_platform_proteins$qval_FG_DDA_tims_ens<-1
cross_platform_proteins$logFC_FG_DDA_tims_ens<-0
com=intersect(cross_platform_proteins$Proteins, top1_FG_DDA_tims_ens$dea_res$protein)
idx1=match(com, cross_platform_proteins$Proteins)
idx2=match(com, top1_FG_DDA_tims_ens$dea_res$protein)
cross_platform_proteins$qval_FG_DDA_tims_ens[idx1]=top1_FG_DDA_tims_ens$dea_res$adj.pvalue[idx2]
cross_platform_proteins$logFC_FG_DDA_tims_ens[idx1]=top1_FG_DDA_tims_ens$dea_res$logFC[idx2]

cross_platform_proteins$qval_MQ_DDA_tims_ens<-1
cross_platform_proteins$logFC_MQ_DDA_tims_ens<-0
com=intersect(cross_platform_proteins$Proteins, top1_MQ_DDA_tims_ens$dea_res$protein)
idx1=match(com, cross_platform_proteins$Proteins)
idx2=match(com, top1_MQ_DDA_tims_ens$dea_res$protein)
cross_platform_proteins$qval_MQ_DDA_tims_ens[idx1]=top1_MQ_DDA_tims_ens$dea_res$adj.pvalue[idx2]
cross_platform_proteins$logFC_MQ_DDA_tims_ens[idx1]=top1_MQ_DDA_tims_ens$dea_res$logFC[idx2]

cross_platform_proteins$qval_DIANN_DIA_ens<-1
cross_platform_proteins$logFC_DIANN_DIA_ens<-0
com=intersect(cross_platform_proteins$Proteins, top1_DIANN_DIA_tims_ens$dea_res$protein)
idx1=match(com, cross_platform_proteins$Proteins)
idx2=match(com, top1_DIANN_DIA_tims_ens$dea_res$protein)
cross_platform_proteins$qval_DIANN_DIA_ens[idx1]=top1_DIANN_DIA_tims_ens$dea_res$adj.pvalue[idx2]
cross_platform_proteins$logFC_DIANN_DIA_ens[idx1]=top1_DIANN_DIA_tims_ens$dea_res$logFC[idx2]

cross_platform_proteins$qval_spt_DIA_ens<-1
cross_platform_proteins$logFC_spt_DIA_ens<-0
com=intersect(cross_platform_proteins$Proteins, top1_spt_DIA_tims_ens$dea_res$protein)
idx1=match(com, cross_platform_proteins$Proteins)
idx2=match(com, top1_spt_DIA_tims_ens$dea_res$protein)
cross_platform_proteins$qval_spt_DIA_ens[idx1]=top1_spt_DIA_tims_ens$dea_res$adj.pvalue[idx2]
cross_platform_proteins$logFC_spt_DIA_ens[idx1]=top1_spt_DIA_tims_ens$dea_res$logFC[idx2]

cross_platform_proteins$qval_FG_DDA_st5600<-1
cross_platform_proteins$logFC_FG_DDA_st5600<-0
com=intersect(cross_platform_proteins$Proteins, top1_FG_DDA_st5600$dea_res$protein)
idx1=match(com, cross_platform_proteins$Proteins)
idx2=match(com, top1_FG_DDA_st5600$dea_res$protein)
cross_platform_proteins$qval_FG_DDA_st5600[idx1]=top1_FG_DDA_st5600$dea_res$adj.pvalue[idx2]
cross_platform_proteins$logFC_FG_DDA_st5600[idx1]=top1_FG_DDA_st5600$dea_res$logFC[idx2]

cross_platform_proteins$qval_MQ_DDA_st5600<-1
cross_platform_proteins$logFC_MQ_DDA_st5600<-0
com=intersect(cross_platform_proteins$Proteins, top1_MQ_DDA_st5600$dea_res$protein)
idx1=match(com, cross_platform_proteins$Proteins)
idx2=match(com, top1_MQ_DDA_st5600$dea_res$protein)
cross_platform_proteins$qval_MQ_DDA_st5600[idx1]=top1_MQ_DDA_st5600$dea_res$adj.pvalue[idx2]
cross_platform_proteins$logFC_MQ_DDA_st5600[idx1]=top1_MQ_DDA_st5600$dea_res$logFC[idx2]

cross_platform_proteins$qval_FG_DDA_st6600<-1
cross_platform_proteins$logFC_FG_DDA_st6600<-0
com=intersect(cross_platform_proteins$Proteins, top1_FG_DDA_st6600$dea_res$protein)
idx1=match(com, cross_platform_proteins$Proteins)
idx2=match(com, top1_FG_DDA_st6600$dea_res$protein)
cross_platform_proteins$qval_FG_DDA_st6600[idx1]=top1_FG_DDA_st6600$dea_res$adj.pvalue[idx2]
cross_platform_proteins$logFC_FG_DDA_st6600[idx1]=top1_FG_DDA_st6600$dea_res$logFC[idx2]

cross_platform_proteins$qval_MQ_DDA_st6600<-1
cross_platform_proteins$logFC_MQ_DDA_st6600<-0
com=intersect(cross_platform_proteins$Proteins, top1_MQ_DDA_st6600$dea_res$protein)
idx1=match(com, cross_platform_proteins$Proteins)
idx2=match(com, top1_MQ_DDA_st6600$dea_res$protein)
cross_platform_proteins$qval_MQ_DDA_st6600[idx1]=top1_MQ_DDA_st6600$dea_res$adj.pvalue[idx2]
cross_platform_proteins$logFC_MQ_DDA_st6600[idx1]=top1_MQ_DDA_st6600$dea_res$logFC[idx2]

cross_platform_proteins$qval_FG_DDA_qe<-1
cross_platform_proteins$logFC_FG_DDA_qe<-0
com=intersect(cross_platform_proteins$Proteins, top1_FG_DDA_qe$dea_res$protein)
idx1=match(com, cross_platform_proteins$Proteins)
idx2=match(com, top1_FG_DDA_qe$dea_res$protein)
cross_platform_proteins$qval_FG_DDA_qe[idx1]=top1_FG_DDA_qe$dea_res$adj.pvalue[idx2]
cross_platform_proteins$logFC_FG_DDA_qe[idx1]=top1_FG_DDA_qe$dea_res$logFC[idx2]

cross_platform_proteins$qval_MQ_DDA_qe<-1
cross_platform_proteins$logFC_MQ_DDA_qe<-0
com=intersect(cross_platform_proteins$Proteins, top1_MQ_DDA_qe$dea_res$protein)
idx1=match(com, cross_platform_proteins$Proteins)
idx2=match(com, top1_MQ_DDA_qe$dea_res$protein)
cross_platform_proteins$qval_MQ_DDA_qe[idx1]=top1_MQ_DDA_qe$dea_res$adj.pvalue[idx2]
cross_platform_proteins$logFC_MQ_DDA_qe[idx1]=top1_MQ_DDA_qe$dea_res$logFC[idx2]


cross_platform_proteins$qval_FG_DDA_st5600_ens<-1
cross_platform_proteins$logFC_FG_DDA_st5600_ens<-0
com=intersect(cross_platform_proteins$Proteins, top1_FG_DDA_st5600_ens$dea_res$protein)
idx1=match(com, cross_platform_proteins$Proteins)
idx2=match(com, top1_FG_DDA_st5600_ens$dea_res$protein)
cross_platform_proteins$qval_FG_DDA_st5600_ens[idx1]=top1_FG_DDA_st5600_ens$dea_res$adj.pvalue[idx2]
cross_platform_proteins$logFC_FG_DDA_st5600_ens[idx1]=top1_FG_DDA_st5600_ens$dea_res$logFC[idx2]

cross_platform_proteins$qval_MQ_DDA_st5600_ens<-1
cross_platform_proteins$logFC_MQ_DDA_st5600_ens<-0
com=intersect(cross_platform_proteins$Proteins, top1_MQ_DDA_st5600_ens$dea_res$protein)
idx1=match(com, cross_platform_proteins$Proteins)
idx2=match(com, top1_MQ_DDA_st5600_ens$dea_res$protein)
cross_platform_proteins$qval_MQ_DDA_st5600_ens[idx1]=top1_MQ_DDA_st5600_ens$dea_res$adj.pvalue[idx2]
cross_platform_proteins$logFC_MQ_DDA_st5600_ens[idx1]=top1_MQ_DDA_st5600_ens$dea_res$logFC[idx2]

cross_platform_proteins$qval_FG_DDA_st6600_ens<-1
cross_platform_proteins$logFC_FG_DDA_st6600_ens<-0
com=intersect(cross_platform_proteins$Proteins, top1_FG_DDA_st6600_ens$dea_res$protein)
idx1=match(com, cross_platform_proteins$Proteins)
idx2=match(com, top1_FG_DDA_st6600_ens$dea_res$protein)
cross_platform_proteins$qval_FG_DDA_st6600_ens[idx1]=top1_FG_DDA_st6600_ens$dea_res$adj.pvalue[idx2]
cross_platform_proteins$logFC_FG_DDA_st6600_ens[idx1]=top1_FG_DDA_st6600_ens$dea_res$logFC[idx2]

cross_platform_proteins$qval_MQ_DDA_st6600_ens<-1
cross_platform_proteins$logFC_MQ_DDA_st6600_ens<-0
com=intersect(cross_platform_proteins$Proteins, top1_MQ_DDA_st6600_ens$dea_res$protein)
idx1=match(com, cross_platform_proteins$Proteins)
idx2=match(com, top1_MQ_DDA_st6600_ens$dea_res$protein)
cross_platform_proteins$qval_MQ_DDA_st6600_ens[idx1]=top1_MQ_DDA_st6600_ens$dea_res$adj.pvalue[idx2]
cross_platform_proteins$logFC_MQ_DDA_st6600_ens[idx1]=top1_MQ_DDA_st6600_ens$dea_res$logFC[idx2]

cross_platform_proteins$qval_FG_DDA_qe_ens<-1
cross_platform_proteins$logFC_FG_DDA_qe_ens<-0
com=intersect(cross_platform_proteins$Proteins, top1_FG_DDA_qe_ens$dea_res$protein)
idx1=match(com, cross_platform_proteins$Proteins)
idx2=match(com, top1_FG_DDA_qe_ens$dea_res$protein)
cross_platform_proteins$qval_FG_DDA_qe_ens[idx1]=top1_FG_DDA_qe_ens$dea_res$adj.pvalue[idx2]
cross_platform_proteins$logFC_FG_DDA_qe_ens[idx1]=top1_FG_DDA_qe_ens$dea_res$logFC[idx2]

cross_platform_proteins$qval_MQ_DDA_qe_ens<-1
cross_platform_proteins$logFC_MQ_DDA_qe_ens<-0
com=intersect(cross_platform_proteins$Proteins, top1_MQ_DDA_qe_ens$dea_res$protein)
idx1=match(com, cross_platform_proteins$Proteins)
idx2=match(com, top1_MQ_DDA_qe_ens$dea_res$protein)
cross_platform_proteins$qval_MQ_DDA_qe_ens[idx1]=top1_MQ_DDA_qe_ens$dea_res$adj.pvalue[idx2]
cross_platform_proteins$logFC_MQ_DDA_qe_ens[idx1]=top1_MQ_DDA_qe_ens$dea_res$logFC[idx2]
sheet=paste(dt, '_ci_proteins')

addWorksheet(wb,sheet)
writeData(wb, sheet, cross_platform_proteins, startRow = 1, startCol = 1, colNames = TRUE, rowNames = TRUE)
saveWorkbook(wb, file = paste0(save_fold, 'cross_instrument_compare.xlsx'), overwrite = TRUE)

roc_FG_DDA_tims <- roc(cross_platform_proteins$labels, 1-cross_platform_proteins$qval_FG_DDA_tims)
roc_MQ_DDA_tims <- roc(cross_platform_proteins$labels, 1-cross_platform_proteins$qval_MQ_DDA_tims)
roc_DIANN_DIA <- roc(cross_platform_proteins$labels, 1-cross_platform_proteins$qval_DIANN_DIA)
roc_spt_DIA <- roc(cross_platform_proteins$labels, 1-cross_platform_proteins$qval_spt_DIA)

roc_FG_DDA_tims_ens <- roc(cross_platform_proteins$labels, 1-cross_platform_proteins$qval_FG_DDA_tims_ens)
roc_MQ_DDA_tims_ens <- roc(cross_platform_proteins$labels, 1-cross_platform_proteins$qval_MQ_DDA_tims_ens)
roc_DIANN_DIA_ens <- roc(cross_platform_proteins$labels, 1-cross_platform_proteins$qval_DIANN_DIA_ens)
roc_spt_DIA_ens <- roc(cross_platform_proteins$labels, 1-cross_platform_proteins$qval_spt_DIA_ens)

roc_FG_DDA_st5600 <- roc(cross_platform_proteins$labels, 1-cross_platform_proteins$qval_FG_DDA_st5600)
roc_MQ_DDA_st5600 <- roc(cross_platform_proteins$labels, 1-cross_platform_proteins$qval_MQ_DDA_st5600)
roc_FG_DDA_st6600 <- roc(cross_platform_proteins$labels, 1-cross_platform_proteins$qval_FG_DDA_st6600)
roc_MQ_DDA_st6600 <- roc(cross_platform_proteins$labels, 1-cross_platform_proteins$qval_MQ_DDA_st6600)
roc_FG_DDA_qe <- roc(cross_platform_proteins$labels, 1-cross_platform_proteins$qval_FG_DDA_qe)
roc_MQ_DDA_qe <- roc(cross_platform_proteins$labels, 1-cross_platform_proteins$qval_MQ_DDA_qe)

roc_FG_DDA_st5600_ens <- roc(cross_platform_proteins$labels, 1-cross_platform_proteins$qval_FG_DDA_st5600_ens)
roc_MQ_DDA_st5600_ens <- roc(cross_platform_proteins$labels, 1-cross_platform_proteins$qval_MQ_DDA_st5600_ens)
roc_FG_DDA_st6600_ens <- roc(cross_platform_proteins$labels, 1-cross_platform_proteins$qval_FG_DDA_st6600_ens)
roc_MQ_DDA_st6600_ens <- roc(cross_platform_proteins$labels, 1-cross_platform_proteins$qval_MQ_DDA_st6600_ens)
roc_FG_DDA_qe_ens <- roc(cross_platform_proteins$labels, 1-cross_platform_proteins$qval_FG_DDA_qe_ens)
roc_MQ_DDA_qe_ens <- roc(cross_platform_proteins$labels, 1-cross_platform_proteins$qval_MQ_DDA_qe_ens)

cross_platform_metrics<-function(cross_platform_proteins, setting, roc){
  lab=cross_platform_proteins$labels
  if(setting=='FG_DDA_tims'){
    logFCs=cross_platform_proteins$logFC_FG_DDA_tims
    qvals=cross_platform_proteins$qval_FG_DDA_tims
    roc_s<-roc_FG_DDA_tims
  }else if(setting=='MQ_DDA_tims'){
    logFCs=cross_platform_proteins$logFC_MQ_DDA_tims
    qvals=cross_platform_proteins$qval_MQ_DDA_tims
    roc_s<-roc_MQ_DDA_tims
  }else if(setting=='DIANN_DIA'){
    logFCs=cross_platform_proteins$logFC_DIANN_DIA
    qvals=cross_platform_proteins$qval_DIANN_DIA
    roc_s<-roc_DIANN_DIA
  }else if(setting=='spt_DIA'){
    logFCs=cross_platform_proteins$logFC_spt_DIA
    qvals=cross_platform_proteins$qval_spt_DIA
    roc_s<-roc_spt_DIA
  }else if(setting=='FG_DDA_tims_ens'){
    logFCs=cross_platform_proteins$logFC_FG_DDA_tims_ens
    qvals=cross_platform_proteins$qval_FG_DDA_tims_ens
    roc_s<-roc_FG_DDA_tims_ens
  }else if(setting=='MQ_DDA_tims_ens'){
    logFCs=cross_platform_proteins$logFC_MQ_DDA_tims_ens
    qvals=cross_platform_proteins$qval_MQ_DDA_tims_ens
    roc_s<-roc_MQ_DDA_tims_ens
  }else if(setting=='DIANN_DIA_ens'){
    logFCs=cross_platform_proteins$logFC_DIANN_DIA_ens
    qvals=cross_platform_proteins$qval_DIANN_DIA_ens
    roc_s<-roc_DIANN_DIA_ens
  }else if(setting=='spt_DIA_ens'){
    logFCs=cross_platform_proteins$logFC_spt_DIA_ens
    qvals=cross_platform_proteins$qval_spt_DIA_ens
    roc_s<-roc_spt_DIA_ens
  }else if(setting=='FG_DDA_st5600'){
    logFCs=cross_platform_proteins$logFC_FG_DDA_st5600
    qvals=cross_platform_proteins$qval_FG_DDA_st5600
    roc_s<-roc_FG_DDA_st5600
  }else if(setting=='MQ_DDA_st5600'){
    logFCs=cross_platform_proteins$logFC_MQ_DDA_st5600
    qvals=cross_platform_proteins$qval_MQ_DDA_st5600
    roc_s<-roc_MQ_DDA_st5600
  }else if(setting=='FG_DDA_st6600'){
    logFCs=cross_platform_proteins$logFC_FG_DDA_st6600
    qvals=cross_platform_proteins$qval_FG_DDA_st6600
    roc_s<-roc_FG_DDA_st6600
  }else if(setting=='MQ_DDA_st6600'){
    logFCs=cross_platform_proteins$logFC_MQ_DDA_st6600
    qvals=cross_platform_proteins$qval_MQ_DDA_st6600
    roc_s<-roc_MQ_DDA_st6600
  }else if(setting=='FG_DDA_qe'){
    logFCs=cross_platform_proteins$logFC_FG_DDA_qe
    qvals=cross_platform_proteins$qval_FG_DDA_qe
    roc_s<-roc_FG_DDA_qe
  }else if(setting=='MQ_DDA_qe'){
    logFCs=cross_platform_proteins$logFC_MQ_DDA_qe
    qvals=cross_platform_proteins$qval_MQ_DDA_qe
    roc_s<-roc_MQ_DDA_qe
  }else if(setting=='FG_DDA_st5600_ens'){
    logFCs=cross_platform_proteins$logFC_FG_DDA_st5600_ens
    qvals=cross_platform_proteins$qval_FG_DDA_st5600_ens
    roc_s<-roc_FG_DDA_st5600_ens
  }else if(setting=='MQ_DDA_st5600_ens'){
    logFCs=cross_platform_proteins$logFC_MQ_DDA_st5600_ens
    qvals=cross_platform_proteins$qval_MQ_DDA_st5600_ens
    roc_s<-roc_MQ_DDA_st5600_ens
  }else if(setting=='FG_DDA_st6600_ens'){
    logFCs=cross_platform_proteins$logFC_FG_DDA_st6600_ens
    qvals=cross_platform_proteins$qval_FG_DDA_st6600_ens
    roc_s<-roc_FG_DDA_st6600_ens
  }else if(setting=='MQ_DDA_st6600_ens'){
    logFCs=cross_platform_proteins$logFC_MQ_DDA_st6600_ens
    qvals=cross_platform_proteins$qval_MQ_DDA_st6600_ens
    roc_s<-roc_MQ_DDA_st6600_ens
  }else if(setting=='FG_DDA_qe_ens'){
    logFCs=cross_platform_proteins$logFC_FG_DDA_qe_ens
    qvals=cross_platform_proteins$qval_FG_DDA_qe_ens
    roc_s<-roc_FG_DDA_qe_ens
  }else if(setting=='MQ_DDA_qe_ens'){
    logFCs=cross_platform_proteins$logFC_MQ_DDA_qe_ens
    qvals=cross_platform_proteins$qval_MQ_DDA_qe_ens
    roc_s<-roc_MQ_DDA_qe_ens
  }
  TP=length(which(lab==1 & abs(logFCs)>=log2(1.5) & qvals<=0.05))
  TN=length(which(lab==0 & (abs(logFCs)<=log2(1.5) | qvals>0.05)))
  FP=length(which(lab==0 & abs(logFCs)>=log2(1.5) & qvals<=0.05))
  FN=length(which(lab==1 & (abs(logFCs)<=log2(1.5) | qvals>0.05)))

  Recall = TP/(TP+FN)
  Precision = TP/(TP+FP)
  Specificity = TN/(TN+FP)
  F1=2*Recall*Precision/(Recall+Precision)
  MCC = (TP*TN-FP*FN)/(((TP+FP)^0.5)*((TP+FN)^0.5)*((TN+FP)^0.5)*((TN+FN))^0.5)
  Gmean = (Recall*Specificity)^0.5
  nMCC = (1+MCC)/2

  pauc001=as.numeric(auc(roc_s, partial.auc=c(1, 0.99), partial.auc.correct=TRUE))
  pauc005=as.numeric(auc(roc_s, partial.auc=c(1, 0.95), partial.auc.correct=TRUE))
  pauc01=as.numeric(auc(roc_s, partial.auc=c(1, 0.90), partial.auc.correct=TRUE))
  metrics=data.frame(metric=c('pAUC(0.01)','pAUC(0.05)','pAUC(0.1)','TP', 'TN', 'FP', 'FN', 'Recall', 'Precision', 'Specificity', 'F1', 'MCC', 'G-mean','nMCC'),
                     value=c(pauc001,pauc005,pauc01,TP, TN, FP, FN, Recall, Precision, Specificity, F1, MCC, Gmean,nMCC))
  return(metrics)
}

metrics_FG_DDA_tims=cross_platform_metrics(cross_platform_proteins, 'FG_DDA_tims', roc_FG_DDA_tims)
metrics_MQ_DDA_tims=cross_platform_metrics(cross_platform_proteins, 'MQ_DDA_tims', roc_MQ_DDA_tims)
metrics_DIANN_DIA=cross_platform_metrics(cross_platform_proteins, 'DIANN_DIA', roc_DIANN_DIA)
metrics_spt_DIA=cross_platform_metrics(cross_platform_proteins, 'spt_DIA', roc_spt_DIA)

metrics_FG_DDA_tims_ens=cross_platform_metrics(cross_platform_proteins, 'FG_DDA_tims_ens', roc_FG_DDA_tims_ens)
metrics_MQ_DDA_tims_ens=cross_platform_metrics(cross_platform_proteins, 'MQ_DDA_tims_ens', roc_MQ_DDA_tims_ens)
metrics_DIANN_DIA_ens=cross_platform_metrics(cross_platform_proteins, 'DIANN_DIA_ens', roc_DIANN_DIA_ens)
metrics_spt_DIA_ens=cross_platform_metrics(cross_platform_proteins, 'spt_DIA_ens', roc_spt_DIA_ens)

metrics_FG_DDA_st5600=cross_platform_metrics(cross_platform_proteins, 'FG_DDA_st5600', roc_FG_DDA_st5600)
metrics_MQ_DDA_st5600=cross_platform_metrics(cross_platform_proteins, 'MQ_DDA_st5600', roc_MQ_DDA_st5600)
metrics_FG_DDA_st6600=cross_platform_metrics(cross_platform_proteins, 'FG_DDA_st6600', roc_FG_DDA_st6600)
metrics_MQ_DDA_st6600=cross_platform_metrics(cross_platform_proteins, 'MQ_DDA_st6600', roc_MQ_DDA_st6600)
metrics_FG_DDA_qe=cross_platform_metrics(cross_platform_proteins, 'FG_DDA_qe', roc_FG_DDA_qe)
metrics_MQ_DDA_qe=cross_platform_metrics(cross_platform_proteins, 'MQ_DDA_qe', roc_MQ_DDA_qe)

metrics_FG_DDA_st5600_ens=cross_platform_metrics(cross_platform_proteins, 'FG_DDA_st5600_ens', roc_FG_DDA_st5600_ens)
metrics_MQ_DDA_st5600_ens=cross_platform_metrics(cross_platform_proteins, 'MQ_DDA_st5600_ens', roc_MQ_DDA_st5600_ens)
metrics_FG_DDA_st6600_ens=cross_platform_metrics(cross_platform_proteins, 'FG_DDA_st6600_ens', roc_FG_DDA_st6600_ens)
metrics_MQ_DDA_st6600_ens=cross_platform_metrics(cross_platform_proteins, 'MQ_DDA_st6600_ens', roc_MQ_DDA_st6600_ens)
metrics_FG_DDA_qe_ens=cross_platform_metrics(cross_platform_proteins, 'FG_DDA_qe_ens', roc_FG_DDA_qe_ens)
metrics_MQ_DDA_qe_ens=cross_platform_metrics(cross_platform_proteins, 'MQ_DDA_qe_ens', roc_MQ_DDA_qe_ens)

rdar_dt<-rbind(metrics_FG_DDA_tims$value, metrics_MQ_DDA_tims$value, metrics_DIANN_DIA$value, metrics_spt_DIA$value,
               metrics_FG_DDA_tims_ens$value, metrics_MQ_DDA_tims_ens$value, metrics_DIANN_DIA_ens$value,
               metrics_spt_DIA_ens$value,
               metrics_FG_DDA_st5600$value, metrics_MQ_DDA_st5600$value, metrics_FG_DDA_st5600_ens$value, metrics_MQ_DDA_st5600_ens$value,
               metrics_FG_DDA_st6600$value, metrics_MQ_DDA_st6600$value, metrics_FG_DDA_st6600_ens$value, metrics_MQ_DDA_st6600_ens$value,
               metrics_FG_DDA_qe$value, metrics_MQ_DDA_qe$value, metrics_FG_DDA_qe_ens$value, metrics_MQ_DDA_qe_ens$value
               )
colnames(rdar_dt)<-metrics_FG_DDA_tims$metric
row.names(rdar_dt)<-c('FG_DDA_tims','MQ_DDA_tims', 'DIANN_DIA_tims', 'spt_DIA_tims', 'FG_DDA_tims_ens',
                      'MQ_DDA_tims_ens', 'DIANN_DIA_tims_ens', 'spt_DIA_tims_ens',
                      'FG_DDA_st5600','MQ_DDA_st5600','FG_DDA_st5600_ens', 'MQ_DDA_st5600_ens',
                      'FG_DDA_st6600','MQ_DDA_st6600','FG_DDA_st6600_ens', 'MQ_DDA_st6600_ens',
                      'FG_DDA_qe','MQ_DDA_qe','FG_DDA_qe_ens', 'MQ_DDA_qe_ens')
rdar_dt<-as.data.frame(as.matrix(rdar_dt))
pt_rdar_dt<-data.frame(group=row.names(rdar_dt))

sheet=paste(dt, '_cp_metrics')

addWorksheet(wb,sheet)
writeData(wb, sheet, rdar_dt, startRow = 1, startCol = 1, colNames = TRUE, rowNames = TRUE)
saveWorkbook(wb, file = paste0(save_fold, 'cross_instrument_compare.xlsx'), overwrite = TRUE)
pt_rdar_dt<-cbind(pt_rdar_dt, rdar_dt[,c(1,2,3,13,14)]) #,4,6
#pt_rdar_dt$group<-row.names(pt_rdar_dt)
library(ggradar)
#ggradar(pt_rdar_dt,group.colours=c("#a45ee3","#1f11b4","#f078b4","#f2ca01","gray","pink",'red',"orange"),legend.position = 'right')
#a=c('#5494cc','#e18283','#0d898a','#f9cc52')
col=c("#15d08a","#f2ca01","#a45ee3","#1f11b4","#c3642c",
      "#f078b4", '#5494cc',"gray","pink",'red','orange',"#21708a","#13ca01","#a45673","#a232b4","#d4a42c",
      "#0a48b4", '#c1256c',"black","#ccffc1",'#ddcd11','#342456')
melt_dt<-melt(pt_rdar_dt)
colnames(melt_dt)[1]<-'method'
# p3=ggplot(melt_dt, aes(x = variable, y = value, colour = method)) +
#   geom_bar(aes(fill = method),position=position_dodge(width = 0.5), stat = "identity", width = 0.1) +
#   geom_point(position=position_dodge(width = 0.1),size=3, shape=19) +
#   scale_fill_manual(values=col[])
#   mytheme = theme_classic() + theme(axis.text.x = element_text(size = 16, angle = -20),axis.text.y = element_text(size = 16))+
#   theme(axis.title.y= element_text(size=16))+theme(axis.title.x = element_text(size = 16))+
#   theme(legend.title=element_text(size=16),legend.text=element_text(size=16), legend.position = "top")
#   p3=p3+mytheme
#   p3

ggdotchart(melt_dt, "variable", "value", group = "method", color = "method",
           palette = col,
           position=position_dodge(width = 0.9),
           add = "segments",
           add.params = list(color = "lightgray", size = 1),
           dot.size = 3) + theme(axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
                                 axis.title.x = element_blank()) + coord_flip()


count_TPs<-data.frame(method=row.names(rdar_dt))
count_TPs<-cbind(count_TPs, rdar_dt[,c(4:7)]) #,4,6
melt_dt<-melt(count_TPs)
melt_dt$label<-as.character(melt_dt$value)
melt_dt$label[which(melt_dt$variable=='TN' | melt_dt$variable=='FN')]=''

nums<-vector()


for (i in 1:length(unique(melt_dt$method))){
  num_i<-melt_dt[which(melt_dt$method==melt_dt$method[i]),]
  uni_var<-levels(factor(melt_dt$variable))
  re_order<-match(num_i$variable, uni_var)
  num_i<-num_i[re_order,]
  l_y<-c(rep(0,length(uni_var)))
  for (j in length(uni_var):1) {
    if(j==length(uni_var)){
      l_y[j]=num_i$value[j]/2
    }else{
      l_y_j=num_i$value[j]/2 + sum(num_i$value[c((j+1):length(uni_var))])
      l_y[j]=l_y_j
    }
  }
  num_i$l_y<-l_y
  nums<-rbind(nums,num_i)
}

## Figure 5D
ggbarplot(nums, "method", "value", group = "variable", color = "variable",fill="variable",
          platette = c('#5494cc','#e18283','#0d898a','#f9cc52') #label=as.character(melt_dt$value)),
) +geom_vline(xintercept = 10.5, linetype=2, show.legend=FALSE)+ scale_x_discrete(limits=c('FG_DDA_st5600_ens', 'MQ_DDA_st5600_ens',
                              'FG_DDA_st6600_ens','MQ_DDA_st6600_ens',
                              'FG_DDA_qe_ens','MQ_DDA_qe_ens',
                              'FG_DDA_tims_ens','MQ_DDA_tims_ens',
                              'DIANN_DIA_tims_ens', 'spt_DIA_tims_ens',
                              'FG_DDA_st5600', 'MQ_DDA_st5600', 'FG_DDA_st6600', 'MQ_DDA_st6600',
                              'FG_DDA_qe','MQ_DDA_qe','FG_DDA_tims','MQ_DDA_tims',
                              'DIANN_DIA_tims','spt_DIA_tims'
                              )) + geom_text(aes(y=l_y, label=label)) +
theme(axis.title.y = element_blank())+ theme(axis.line.y = element_blank())+ coord_flip()

write.table(nums, paste0(save_fold, 'Figure5D.csv'), col.names = T, row.names = F)
