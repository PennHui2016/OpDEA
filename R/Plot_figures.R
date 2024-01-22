library(shinydashboard)
library(shiny)
library(threejs)
library(DT)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(ggsci)
library("readxl")
library(ggalluvial)
library(golem)
#source('./R/app_config.R')
#addResourcePath('root_fold', 'inst/app/root_fold')
root_fold<-'inst/app/www/'
wf_frag<-read_excel(paste0(root_fold, 'workflows.xlsx'), sheet='FG_DDA', col_names = TRUE)
wf_mq<-read_excel(paste0(root_fold, 'workflows.xlsx'), sheet='MQ_DDA', col_names = TRUE)
wf_dia<-read_excel(paste0(root_fold, 'workflows.xlsx'), sheet='DIANN_DIA', col_names = TRUE)
wf_spt<-read_excel(paste0(root_fold, 'workflows.xlsx'), sheet='spt_DIA', col_names = TRUE)
wf_fg_tmt<-read_excel(paste0(root_fold, 'workflows.xlsx'), sheet='FG_TMT', col_names = TRUE)
wf_mq_tmt<-read_excel(paste0(root_fold, 'workflows.xlsx'), sheet='MQ_TMT', col_names = TRUE)

wf_fg_dda_ens<-read_excel(paste0(root_fold, 'workflows_ens.xlsx'), sheet='FG_DDA', col_names = TRUE)
wf_mq_dda_ens<-read_excel(paste0(root_fold, 'workflows_ens.xlsx'), sheet='MQ_DDA', col_names = TRUE)
wf_diann_dia_ens<-read_excel(paste0(root_fold, 'workflows_ens.xlsx'), sheet='DIANN_DIA', col_names = TRUE)
wf_spt_dia_ens<-read_excel(paste0(root_fold, 'workflows_ens.xlsx'), sheet='spt_DIA', col_names = TRUE)
#wf_file<-golem::app_sys('workflows.xlsx')

#use_internal_file(path = './inst/workflows.xlsx', name='workflows.xlsx')
#wf_mq<-read_excel('root_fold/workflows.xlsx', sheet='maxquant', col_names = TRUE)
#wf_dia<-read_excel(wf_file, sheet='DIA-NN', col_names = TRUE)

col=c("#1f78b4","#b2df8a","#a6cee3","#33a02c","#fb9a99","#fdbf6f")
metric_frag<-list(pauc001=read_excel(paste0(root_fold, 'metrics_FG_DDA.xlsx'), sheet='pAUC0.01'),
                  pauc005=read_excel(paste0(root_fold, 'metrics_FG_DDA.xlsx'), sheet='pAUC0.05'),
                  pauc01=read_excel(paste0(root_fold, 'metrics_FG_DDA.xlsx'), sheet='pAUC0.1'),
                  Gmean=read_excel(paste0(root_fold, 'metrics_FG_DDA.xlsx'), sheet='G-mean'),
                  nMCC=read_excel(paste0(root_fold, 'metrics_FG_DDA.xlsx'), sheet='nMCC'))
metric_mq<-list(pauc001=read_excel(paste0(root_fold, 'metrics_MQ_DDA.xlsx'), sheet='pAUC0.01'),
                pauc005=read_excel(paste0(root_fold, 'metrics_MQ_DDA.xlsx'), sheet='pAUC0.05'),
                pauc01=read_excel(paste0(root_fold, 'metrics_MQ_DDA.xlsx'), sheet='pAUC0.1'),
                Gmean=read_excel(paste0(root_fold, 'metrics_MQ_DDA.xlsx'), sheet='G-mean'),
                nMCC=read_excel(paste0(root_fold, 'metrics_MQ_DDA.xlsx'), sheet='nMCC'))
metric_dia<-list(pauc001=read_excel(paste0(root_fold, 'metrics_DIANN_DIA.xlsx'), sheet='pAUC0.01'),
                 pauc005=read_excel(paste0(root_fold, 'metrics_DIANN_DIA.xlsx'), sheet='pAUC0.05'),
                 pauc01=read_excel(paste0(root_fold, 'metrics_DIANN_DIA.xlsx'), sheet='pAUC0.1'),
                 Gmean=read_excel(paste0(root_fold, 'metrics_DIANN_DIA.xlsx'), sheet='G-mean'),
                 nMCC=read_excel(paste0(root_fold, 'metrics_DIANN_DIA.xlsx'), sheet='nMCC'))
metric_spt<-list(pauc001=read_excel(paste0(root_fold, 'metrics_spt_DIA.xlsx'), sheet='pAUC0.01'),
                 pauc005=read_excel(paste0(root_fold, 'metrics_spt_DIA.xlsx'), sheet='pAUC0.05'),
                 pauc01=read_excel(paste0(root_fold, 'metrics_spt_DIA.xlsx'), sheet='pAUC0.1'),
                 Gmean=read_excel(paste0(root_fold, 'metrics_spt_DIA.xlsx'), sheet='G-mean'),
                 nMCC=read_excel(paste0(root_fold, 'metrics_spt_DIA.xlsx'), sheet='nMCC'))
metric_fg_tmt<-list(pauc001=read_excel(paste0(root_fold, 'metrics_FG_TMT.xlsx'), sheet='pAUC0.01'),
                 pauc005=read_excel(paste0(root_fold, 'metrics_FG_TMT.xlsx'), sheet='pAUC0.05'),
                 pauc01=read_excel(paste0(root_fold, 'metrics_FG_TMT.xlsx'), sheet='pAUC0.1'),
                 Gmean=read_excel(paste0(root_fold, 'metrics_FG_TMT.xlsx'), sheet='G-mean'),
                 nMCC=read_excel(paste0(root_fold, 'metrics_FG_TMT.xlsx'), sheet='nMCC'))
metric_mq_tmt<-list(pauc001=read_excel(paste0(root_fold, 'metrics_MQ_TMT.xlsx'), sheet='pAUC0.01'),
                 pauc005=read_excel(paste0(root_fold, 'metrics_MQ_TMT.xlsx'), sheet='pAUC0.05'),
                 pauc01=read_excel(paste0(root_fold, 'metrics_MQ_TMT.xlsx'), sheet='pAUC0.1'),
                 Gmean=read_excel(paste0(root_fold, 'metrics_MQ_TMT.xlsx'), sheet='G-mean'),
                 nMCC=read_excel(paste0(root_fold, 'metrics_MQ_TMT.xlsx'), sheet='nMCC'))




plot_metrics<-function(wf, platform){
  library(ggplot2)
  library(reshape2)
  library(ggpubr)
  library(ggsci)
  library("readxl")
  library(ggalluvial)
  if(platform == 'FG_DDA'){
    filename = 'metrics_FG_DDA.xlsx'
    dt_st = 7
    dt_en = 28
    met_dt=metric_frag
  }else if(platform == 'MQ_DDA'){
    filename = 'metrics_MQ_DDA.xlsx'
    dt_st = 7
    dt_en = 28
    met_dt=metric_mq
  }else if(platform == 'DIANN_DIA'){
    filename = 'metrics_DIANN_DIA.xlsx'
    dt_st = 7
    dt_en = 23
    met_dt=metric_dia
  }else if(platform == 'spt_DIA'){
    filename = 'metrics_spt_DIA.xlsx'
    dt_st = 7
    dt_en = 23
    met_dt=metric_spt
  }else if(platform == 'FG_TMT'){
    filename = 'metrics_FG_TMT.xlsx'
    dt_st = 7
    dt_en = 21
    met_dt=metric_fg_tmt
  }else if(platform == 'MQ_TMT'){
    filename = 'metrics_MQ_TMT.xlsx'
    dt_st = 7
    dt_en = 21
    met_dt=metric_mq_tmt
  }
  pauc001<-met_dt$pauc001
  pauc005<-met_dt$pauc005
  pauc01<-met_dt$pauc01
  nMCC<-met_dt$nMCC
  Gmean<-met_dt$Gmean
  #print(wf)
  idx1<-which(pauc001$workflow == wf)
  idx2<-which(pauc005$workflow == wf)
  idx3<-which(pauc01$workflow == wf)
  idx4<-which(nMCC$workflow == wf)
  idx5<-which(Gmean$workflow == wf)
  dt<-data.frame(rbind(as.numeric(pauc001[idx1,][dt_st:dt_en]),
                       as.numeric(pauc005[idx2,][dt_st:dt_en]),
                       as.numeric(pauc01[idx3,][dt_st:dt_en]),
                       as.numeric(nMCC[idx4,][dt_st:dt_en]),
                       as.numeric(Gmean[idx5,][dt_st:dt_en])))
  dt<-t(dt)
  dt[is.na(dt)]<-0
  dt[is.infinite(dt)]<-0
  colnames(dt)<-c('pAUC0.01','pAUC0.05','pAUC0.1','nMCC','G-mean')
  dt_pl<-reshape2::melt(dt)
  colnames(dt_pl)<-c('cl', 'metric', 'value')


  p1<-ggpubr::ggboxplot(dt_pl, "metric", "value",width = 0.5, size=0.8, #outlier.shape=NA,
                color = "metric", title = wf,
                add = "none")+labs( y = 'performance')+
    scale_x_discrete(limits=unique(dt_pl$metric))+
    theme(plot.title = element_text(size = 14), axis.title.x = element_blank())+theme_classic() +
    theme(axis.text.x = element_text(size = 12, angle = 0), axis.text.y = element_text(size = 12))+
    theme(legend.position = "None")+
    theme(plot.title = element_text(size = 14,hjust = 0.5, face = "bold"))+ggsci::scale_color_d3()
  return(p1)
}

getTabFil_fg<-function(input){
    data<-wf_frag
  if (!is.null(input$filter1)) {
    data <- data[data$expression_matrix %in% input$filter1,]
  }
  if (!is.null(input$filter2)) {
    data <- data[data$normalization %in% input$filter2,]
  }
  if (!is.null(input$filter3)) {
    data <- data[data$imputation %in% input$filter3,]
  }
  if (!is.null(input$filter4)) {
    data <- data[data$DEA_tool %in% input$filter4,]
  }
  return(data)
}

getTabFil_mq<-function(input){

  data<-wf_mq

  if (!is.null(input$filter5)) {
    data <- data[data$expression_matrix %in% input$filter5,]
  }
  if (!is.null(input$filter6)) {
    data <- data[data$normalization %in% input$filter6,]
  }
  if (!is.null(input$filter7)) {
    data <- data[data$imputation %in% input$filter7,]
  }
  if (!is.null(input$filter8)) {
    data <- data[data$DEA_tool %in% input$filter8,]
  }
  return(data)
}

getTabFil_dia<-function(input){
    data<-wf_dia
  if (!is.null(input$filter9)) {
    data <- data[data$expression_matrix %in% input$filter9,]
  }
  if (!is.null(input$filter10)) {
    data <- data[data$normalization %in% input$filter10,]
  }
  if (!is.null(input$filter11)) {
    data <- data[data$imputation %in% input$filter11,]
  }
  if (!is.null(input$filter12)) {
    data <- data[data$DEA_tool %in% input$filter12,]
  }
  return(data)
}

getTabFil_spt<-function(input){
  data<-wf_spt
  if (!is.null(input$filter13)) {
    data <- data[data$expression_matrix %in% input$filter13,]
  }
  if (!is.null(input$filter14)) {
    data <- data[data$normalization %in% input$filter14,]
  }
  if (!is.null(input$filter15)) {
    data <- data[data$imputation %in% input$filter15,]
  }
  if (!is.null(input$filter16)) {
    data <- data[data$DEA_tool %in% input$filter16,]
  }
  return(data)
}

getTabFil_fg_tmt<-function(input){
  data<-wf_fg_tmt
  if (!is.null(input$filter17)) {
    data <- data[data$expression_matrix %in% input$filter17,]
  }
  if (!is.null(input$filter18)) {
    data <- data[data$normalization %in% input$filter18,]
  }
  if (!is.null(input$filter19)) {
    data <- data[data$imputation %in% input$filter19,]
  }
  if (!is.null(input$filter20)) {
    data <- data[data$DEA_tool %in% input$filter20,]
  }
  return(data)
}

getTabFil_mq_tmt<-function(input){
  data<-wf_mq_tmt
  if (!is.null(input$filter21)) {
    data <- data[data$expression_matrix %in% input$filter21,]
  }
  if (!is.null(input$filter22)) {
    data <- data[data$normalization %in% input$filter22,]
  }
  if (!is.null(input$filter23)) {
    data <- data[data$imputation %in% input$filter23,]
  }
  if (!is.null(input$filter24)) {
    data <- data[data$DEA_tool %in% input$filter24,]
  }
  return(data)
}

sprs_fg_DDA<-read.table(paste0(root_fold, 'Fragpipe-DDA_spr_avg_rank_test.csv'), sep = ',', header = TRUE)
sprs_mq_DDA<-read.table(paste0(root_fold, 'Maxquant-DDA_spr_avg_rank_test.csv'), sep = ',', header = TRUE)
sprs_diann_DIA<-read.table(paste0(root_fold, 'DIANN-DIA_spr_avg_rank_test.csv'), sep = ',', header = TRUE)
sprs_spt_DIA<-read.table(paste0(root_fold, 'spt-DIA_spr_avg_rank_test.csv'), sep = ',', header = TRUE)
sprs_fg_TMT<-read.table(paste0(root_fold, 'Fragpipe-TMT_spr_avg_rank_test.csv'), sep = ',', header = TRUE)
sprs_mq_TMT<-read.table(paste0(root_fold, 'Maxquant-TMT_spr_avg_rank_test.csv'), sep = ',', header = TRUE)

sprs_fg_DDA$setting<-c(rep('FG_DDA',length(sprs_fg_DDA$dt_name)))
sprs_mq_DDA$setting<-c(rep('MQ_DDA',length(sprs_mq_DDA$dt_name)))
sprs_fg_TMT$setting<-c(rep('FG_TMT',length(sprs_fg_TMT$dt_name)))
sprs_mq_TMT$setting<-c(rep('MQ_TMT',length(sprs_mq_TMT$dt_name)))
sprs_diann_DIA$setting<-c(rep('DIANN_DIA',length(sprs_diann_DIA$dt_name)))
sprs_spt_DIA$setting<-c(rep('spt_DIA',length(sprs_spt_DIA$dt_name)))

all_sprs<-rbind(sprs_fg_DDA,sprs_mq_DDA, sprs_fg_TMT, sprs_mq_TMT,
                sprs_diann_DIA, sprs_spt_DIA)

sprs_dat_mean<-cbind(all_sprs$mean_spr,all_sprs$setting,
                     c(rep('mean',length(all_sprs$dt_name))))


colnames(sprs_dat_mean)<-c('spearman','setting','method')
sprs_dat_mean<-as.data.frame(sprs_dat_mean)
sprs_dat_mean$spearman<-as.numeric(sprs_dat_mean$spearman)

col=c('#425AB2','#34df8a','#F45896','#c21000','#452904','#e2df8a')
p31 = ggplot(sprs_dat_mean, aes(x=setting, y=spearman, color=setting)) +
  geom_boxplot(position = position_dodge(width = 0.3), width=0.3,)+ geom_jitter(width = 0.2)+
  #geom_boxplot(alpha=1,outlier.size=0, size=0.3, width=0.3,fill="white") +
  scale_fill_manual(values = col)+
  labs(x="", y="spearman correlation", color=factor)+
  scale_x_discrete(limits=unique(all_sprs$setting))+stat_summary(fun=mean,
                                                                 geom="point",
                                                                 shape=17,
                                                                 aes(group=method),
                                                                 position=position_dodge(.25), size=2, color="red", fill="red")
mytheme = theme_classic() + theme(axis.text.x = element_text(size = 16, angle = -20),axis.text.y = element_text(size = 10))+
  theme(axis.title.y= element_text(size=16))+theme(axis.title.x = element_text(size = 14))+
  theme(legend.title=element_text(size=16),legend.text=element_text(size=16), legend.position = "none")
p31=p31+mytheme
#p31

lopocv_fg_fig<-p31
#lopocv_mq_fig<-lopocv(sprs_mq, 'maxquant')
#lopocv_dia_fig<-lopocv(sprs_dia, 'DIA-NN')

platforms<-c('Maxquant', 'FragPipe', 'FragPipe', 'Maxquant', 'DIANN', 'spt')
acqs<-c('TMT', 'TMT', 'DDA', 'DDA', 'DIA', 'DIA')
metric_names<-c('Acc', 'F1', 'rec', 'pre', 'Mcc')
setting<-c('MQ_TMT','FG_TMT','FG_DDA','MQ_DDA','DIANN_DIA','spt_DIA')
cls_res<-vector()

for (i in 1:length(platforms)){
  cls_cv_res<-read.table(paste0(root_fold, platforms[i], '_', acqs[i], '_', 'cv_cbt_cls_res_folds_str'),sep = ',')
  cls_res_i<-data.frame(setting=setting[i],metrics=metric_names)
  cls_res_i$value<-colMeans(cls_cv_res)
  cls_res_i$sd<-apply(cls_cv_res, 2, function(x) sd(x))
  cls_res<-rbind(cls_res,cls_res_i)
}

cls_res_fig<-as.data.frame(cls_res)
cls_res_fig$value<-round(cls_res_fig$value, 2)
cls_res_fig$sd<-round(cls_res_fig$sd, 2)
cls_res_fig<-cls_res_fig[which(cls_res_fig$metrics=='F1' | cls_res_fig$metrics=='Mcc'),]

col=c('#427AB2','#F09148','#FF9896','#DBDB8D','#C59D94','#b2df8a')
p5<-ggplot(cls_res_fig, aes(fill=setting, y=value, x=metrics)) +
  geom_bar(position="dodge", stat="identity") + coord_cartesian(ylim=c(0.82,1))+
  geom_text(aes(label=value), position=position_dodge(width=0.9), vjust=-0.25, hjust=-0.15)+
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.1,
                position=position_dodge(.9))+ scale_fill_manual(values=col)+
  labs(x="", y="value")#+coord_flip()

mytheme = theme_classic() + theme(axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14))+
  theme(axis.title.y= element_text(size=14))+theme(axis.title.x = element_text(size = 14))+
  theme(legend.title=element_text(size=12),legend.text=element_text(size=12), legend.position = "top")#+
#guides(color=guide_legend(nrow=2))

p5=p5+mytheme

fig_catboost=p5

##### feature importance
platforms<-c('Maxquant', 'Fragpipe', 'Fragpipe', 'Maxquant', 'DIANN', 'spt')
acqs<-c('TMT', 'TMT', 'DDA', 'DDA', 'DIA', 'DIA')
feature_names<-c('DEA', 'exp_ty', 'MVI', 'norm')
setting<-c('MQ_TMT','FG_TMT','FG_DDA','MQ_DDA','DIANN_DIA','spt_DIA')
importances<-vector()

for (i in 1:length(platforms)){
  cls_impt<-read.table(paste0(root_fold, platforms[i], '_', acqs[i], '_', 'fea_importance.csv'),sep = ',', header = T)
  impt_i<-data.frame(setting=setting[i],feature=feature_names)
  fea_imp<-c()
  for (fea in feature_names) {
    idx=which(cls_impt[,1]==fea)
    if(length(idx)==1){
      fea_imp<-c(fea_imp,cls_impt[,2][idx])
    }else{
      fea_imp<-c(fea_imp,0)
    }
  }
  impt_i$FeatureImportance<-fea_imp
  uni_fea<-factor(impt_i$feature)
  re_order<-match(impt_i$feature, uni_fea)
  impt_i<-impt_i[re_order,]
  l_y<-c(rep(0,length(uni_fea)))
  for (j in length(uni_fea):1) {
    if(j==length(uni_fea)){
      l_y[j]=impt_i$FeatureImportance[j]/2
    }else{
      l_y_j=impt_i$FeatureImportance[j]/2 + sum(impt_i$FeatureImportance[c((j+1):length(uni_fea))])
      l_y[j]=l_y_j
    }
  }
  impt_i$l_y<-l_y
  importances<-rbind(importances,impt_i)
}

#importances<-read_excel(paste0(data_fold, 'sum_DT1_final.xlsx'), sheet = 'feature_importance')
importances_fig<-as.data.frame(importances)
importances_fig$importance<-round(as.numeric(importances_fig$FeatureImportance), 2)
importances_fig$l_y<-round(as.numeric(importances_fig$l_y), 2)
importances_fig$feature<-gsub('exp_ty','Matrix',importances_fig$feature)

p51<-ggplot(importances_fig, aes(fill=feature, y=importance, x=setting,
                                 stratum = feature, alluvium = feature)) +
  geom_col(width = 0.4)+
  geom_alluvium( width = 0.4,alpha = 0.2,knot.pos = 0)+
  geom_text(aes(y=l_y, label=importance)) +
  scale_fill_manual(values = pal_npg()(4))+
  scale_x_discrete(limits=c('FG_DDA', 'MQ_DDA', 'DIANN_DIA', 'spt_DIA', 'FG_TMT','MQ_TMT'))+
  #geom_flow(width = 0.4,alpha = 0.2,knot.pos = 0)

  labs(x="", y="feature importance")

mytheme = theme_classic() + theme(axis.text.x = element_text(size = 13, angle = -20),axis.text.y = element_text(size = 14))+
  theme(axis.title.y= element_text(size=14))+theme(axis.title.x = element_text(size = 14))+
  theme(legend.title=element_text(size=13),legend.text=element_text(size=13), legend.position = "top")

p51=p51+mytheme


imp_catboost<-p51


# runDEA<-function(input){
#   if(!is.null(input$upload_fg_raw)){
#     raw = input$upload_fg_raw$datapath
#   }else{
#     print("no combine_protein.tsv provide!!!")
#   }
#   if(!is.null(input$upload_fg_ion)){
#     evid = input$upload_fg_ion$datapath
#   }else{
#     print("no combine_ion.tsv provide!!!")
#   }
#   if(!is.null(input$upload_fg_dg)){
#     design=input$upload_fg_dg$datapath
#   }else{
#     print("no design file provide!!!")
#   }
#   platform = 'FragPipe'
#   acq = 'DDA'
#   if(!is.null(input$rb11)){
#     exp = input$rb11
#   }
#   if(!is.null(input$rb12)){
#     norm = input$rb12
#   }
#   if(!is.null(input$rb13)){
#     imp = input$rb13
#   }
#   if(!is.null(input$rb14)){
#     dea = input$rb14
#   }
#   logFC=input$logFC
#   qval = input$adjp
#   source('run_DEA_single.R')
#   res_zipfile = run_DEA_fg_DDA(input)
#   #print(paste0(platform, acq, exp, norm, imp, dea,logFC, qval))
#   return(res_zipfile)
# }
