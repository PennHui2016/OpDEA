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
wf_frag<-read_excel(paste0(root_fold, 'workflows.xlsx'), sheet='FragPipe', col_names = TRUE)
wf_mq<-read_excel(paste0(root_fold, 'workflows.xlsx'), sheet='maxquant', col_names = TRUE)
wf_dia<-read_excel(paste0(root_fold, 'workflows.xlsx'), sheet='DIA-NN', col_names = TRUE)
#wf_file<-golem::app_sys('workflows.xlsx')

#use_internal_file(path = './inst/workflows.xlsx', name='workflows.xlsx')
#wf_mq<-read_excel('root_fold/workflows.xlsx', sheet='maxquant', col_names = TRUE)
#wf_dia<-read_excel(wf_file, sheet='DIA-NN', col_names = TRUE)

col=c("#1f78b4","#b2df8a","#a6cee3","#33a02c","#fb9a99","#fdbf6f")
metric_frag<-list(pauc001=read_excel(paste0(root_fold, 'metrics_fragpipe.xlsx'), sheet='pAUC0.01'),
                  pauc005=read_excel(paste0(root_fold, 'metrics_fragpipe.xlsx'), sheet='pAUC0.05'),
                  pauc01=read_excel(paste0(root_fold, 'metrics_fragpipe.xlsx'), sheet='pAUC0.1'),
                  Gmean=read_excel(paste0(root_fold, 'metrics_fragpipe.xlsx'), sheet='G-mean'),
                  nMCC=read_excel(paste0(root_fold, 'metrics_fragpipe.xlsx'), sheet='nMCC'))
metric_mq<-list(pauc001=read_excel(paste0(root_fold, 'metrics_maxquant.xlsx'), sheet='pAUC0.01'),
                pauc005=read_excel(paste0(root_fold, 'metrics_maxquant.xlsx'), sheet='pAUC0.05'),
                pauc01=read_excel(paste0(root_fold, 'metrics_maxquant.xlsx'), sheet='pAUC0.1'),
                Gmean=read_excel(paste0(root_fold, 'metrics_maxquant.xlsx'), sheet='G-mean'),
                nMCC=read_excel(paste0(root_fold, 'metrics_maxquant.xlsx'), sheet='nMCC'))
metric_dia<-list(pauc001=read_excel(paste0(root_fold, 'metrics_diann.xlsx'), sheet='pAUC0.01'),
                 pauc005=read_excel(paste0(root_fold, 'metrics_diann.xlsx'), sheet='pAUC0.05'),
                 pauc01=read_excel(paste0(root_fold, 'metrics_diann.xlsx'), sheet='pAUC0.1'),
                 Gmean=read_excel(paste0(root_fold, 'metrics_diann.xlsx'), sheet='G-mean'),
                 nMCC=read_excel(paste0(root_fold, 'metrics_diann.xlsx'), sheet='nMCC'))


plot_metrics<-function(wf, platform){

  if(platform == 'FragPipe'){
    filename = 'metrics_fragpipe.xlsx'
    dt_st = 7
    dt_en = 126
    met_dt=metric_frag
  }else if(platform == 'maxquant'){
    filename = 'metrics_maxquant.xlsx'
    dt_st = 7
    dt_en = 116
    met_dt=metric_mq
  }else if(platform == 'DIANN'){
    filename = 'metrics_diann.xlsx'
    dt_st = 7
    dt_en = 152
    met_dt=metric_dia
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
  dt<-data.frame(rbind(pauc001[idx1,][dt_st:dt_en], pauc005[idx2,][dt_st:dt_en],pauc01[idx3,][dt_st:dt_en],
                       nMCC[idx4,][dt_st:dt_en],Gmean[idx5,][dt_st:dt_en]))
  dt<-t(dt)
  colnames(dt)<-c('pAUC0.01','pAUC0.05','pAUC0.1','nMCC','G-mean')
  dt_pl<-melt(dt)
  colnames(dt_pl)<-c('cl', 'metric', 'value')


  p1<-ggboxplot(dt_pl, "metric", "value",width = 0.5, size=0.8, #outlier.shape=NA,
                color = "metric", title = wf,
                add = "none")+labs( y = 'performance')+
    scale_x_discrete(limits=unique(dt_pl$metric))+
    theme(plot.title = element_text(size = 14), axis.title.x = element_blank())+theme_classic() +
    theme(axis.text.x = element_text(size = 12, angle = 0), axis.text.y = element_text(size = 12))+
    theme(legend.position = "None")+
    theme(plot.title = element_text(size = 14,hjust = 0.5, face = "bold"))+scale_color_d3()
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

sprs_fg<-read.table(paste0(root_fold, 'fragpipe_spr_avg_rank_test.csv'), sep = ',', header = TRUE)
sprs_mq<-read.table(paste0(root_fold, 'maxquant_spr_avg_rank_test.csv'), sep = ',', header = TRUE)
sprs_dia<-read.table(paste0(root_fold, 'diann_spr_avg_rank_test.csv'), sep = ',', header = TRUE)

lopocv<-function(sprs_frag, tl){
  dt_mean<-data.frame(sprs_frag$dt_name,sprs_frag$mean_spr, c(rep('mean',length(sprs_frag$dt_name))))
  colnames(dt_mean)<-c('project', 'spr', 'method')
  dt_median<-data.frame(sprs_frag$dt_name,sprs_frag$median_spr, c(rep('median',length(sprs_frag$dt_name))))
  colnames(dt_median)<-c('project', 'spr', 'method')

  dt<-rbind(dt_mean, dt_median)
  dt$project[which(dt$project=='human')] = 'human_ecoli'

  p<-ggboxplot(dt, "project", "spr",width = 0.5, size=0.8, #outlier.shape=NA,
              color = "method", title = tl,
              add = "none")+labs(x = 'project', y = 'Spearman correlation')+
  scale_x_discrete(limits=unique(dt$project))+
  theme(plot.title = element_text(size = 14))+theme_classic() +
  theme(axis.text.x = element_text(size = 12, angle = 0), axis.text.y = element_text(size = 12))+
  theme(legend.position = "top")+
  theme(plot.title = element_text(size = 14,hjust = 0.5, face = "bold"))+scale_color_d3()
  return(p)
}

lopocv_fg_fig<-lopocv(sprs_fg, 'FragPipe')
lopocv_mq_fig<-lopocv(sprs_mq, 'maxquant')
lopocv_dia_fig<-lopocv(sprs_dia, 'DIA-NN')

cls_res<-read_excel(paste0(root_fold, 'sum_CatBoost.xlsx'), sheet = 'cls_cv_err')
cls_res_fig<-as.data.frame(cls_res)
cls_res_fig$value<-round(cls_res_fig$value, 2)
cls_res_fig$sd<-round(cls_res_fig$sd, 2)

fig_catboost<-ggplot(cls_res_fig, aes(fill=Platform, y=value, x=metrics)) +
    geom_bar(position="dodge", stat="identity") + coord_cartesian(ylim=c(0,1))+
    geom_text(aes(label=value), position=position_dodge(width=0.9), vjust=-0.25, hjust=-0.15)+
    geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.1,
                  position=position_dodge(.9))+ scale_fill_manual(values=col[c(1,3,5)])+
    labs(x="", y="value")

  mytheme = theme_classic() + theme(axis.text.x = element_text(size = 16),axis.text.y = element_text(size = 16))+
    theme(axis.title.y= element_text(size=16))+theme(axis.title.x = element_text(size = 16))+
    theme(legend.position = 'top', legend.title=element_text(size=16),legend.text=element_text(size=16))

fig_catboost=fig_catboost+mytheme

importances<-read_excel(paste0(root_fold, 'sum_CatBoost.xlsx'), sheet = 'importance')
importances_fig<-as.data.frame(importances)
importances_fig$importance<-round(importances_fig$importance, 2)
importances_fig$l_y<-round(importances_fig$l_y, 2)


imp_catboost<-ggplot(importances_fig, aes(fill=feature, y=importance, x=platform,
                                  stratum = feature, alluvium = feature)) +
    geom_col(width = 0.4)+
    geom_alluvium( width = 0.4,alpha = 0.2,knot.pos = 0)+
    geom_text(aes(y=l_y, label=importance)) +scale_fill_manual(values = pal_npg()(4))+
    #geom_flow(width = 0.4,alpha = 0.2,knot.pos = 0)

    labs(x="", y="feature importance")

  mytheme = theme_classic() + theme(axis.text.x = element_text(size = 16),axis.text.y = element_text(size = 16))+
    theme(axis.title.y= element_text(size=16))+theme(axis.title.x = element_text(size = 16))+
    theme(legend.title=element_text(size=16),legend.text=element_text(size=16), legend.position = "top")

imp_catboost=imp_catboost+mytheme


