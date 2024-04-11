############################
####
#### Figure3 workflow classification, linear model and frequent pattern mining
####
#############################
####### A classification performance

##### !!!!
## The classification of workflow performance levels are implemented with python
## Catboost_workflows_classification.py classifies the FG_DDA, MQ_DDA,FG_TMT, DIANN_DIA and spt_DIA
## Catboost_workflows_classification_TMT_mq.py classifies workflows of MQ_TMT
## run with command: python Catboost_workflows_classification.py
##  python Catboost_workflows_classification.py Catboost_workflows_classification_TMT_mq.py
##


library("readxl")
library(ggplot2)
library(reshape2)
library(ggpubr)
library(cowplot)
library(openxlsx)

library(openxlsx)

wb2 <- createWorkbook()

data_fold<-paste0(dirname(rstudioapi::getSourceEditorContext()$path),'/data/')
save_fold<-paste0(dirname(rstudioapi::getSourceEditorContext()$path),'/save_data/')

platforms<-c('Maxquant', 'FragPipe', 'FragPipe', 'Maxquant', 'DIANN', 'spt')
acqs<-c('TMT', 'TMT', 'DDA', 'DDA', 'DIA', 'DIA')
metric_names<-c('Acc', 'F1', 'rec', 'pre', 'Mcc')
setting<-c('MQ_TMT','FG_TMT','FG_DDA','MQ_DDA','DIANN_DIA','spt_DIA')
cls_res<-vector()

cls_res_all<-c()

for (i in 1:length(platforms)){
  cls_cv_res<-read.table(paste0(data_fold, platforms[i], '_', acqs[i], '_', 'cv_cbt_cls_res_folds_str'),sep = ',')
  cls_res_i<-data.frame(setting=setting[i],metrics=metric_names)
  cls_res_i$value<-colMeans(cls_cv_res)
  cls_res_i$sd<-apply(cls_cv_res, 2, function(x) sd(x))
  cls_res<-rbind(cls_res,cls_res_i)
  cls_res_p <-data.frame(setting=setting[i], metric = 'F1', value = cls_cv_res[,which(metric_names=='F1')])
  cls_res_all<-rbind(cls_res_all, cls_res_p)
  cls_res_p1 <-data.frame(setting=setting[i], metric = 'Mcc', value = cls_cv_res[,which(metric_names=='Mcc')])
  cls_res_all<-rbind(cls_res_all, cls_res_p1)
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

p5=p5+mytheme
p5

cls_res_all<-as.data.frame(cls_res_all)
cls_res_all$value<-round(as.numeric(cls_res_all$value), 2)

sheet='10-fold_cv_res_fig3a'
addWorksheet(wb2,sheet)
writeData(wb2, sheet, cls_res_all, startRow = 1, startCol = 1, colNames = TRUE, rowNames = TRUE)
#saveWorkbook(wb2, file = paste0(save_fold,'10-fold_cv_res.xlsx'), overwrite = TRUE)

ebtop<-function(x){
  return(mean(x)+sd(x))
}
ebbottom<-function(x){
  return(mean(x)-sd(x))
}
set.seed(2024)
col=c('#427AB2','#F09148','#FF9896','#DBDB8D','#C59D94','#b2df8a')
p5_1<-ggplot(cls_res_all, aes(colour=setting, y=value, x=metric)) +
  stat_summary(geom = "bar",fun = "mean", fill='white',
               position = position_dodge(0.9)) + coord_cartesian(ylim=c(0.8,1))+
  stat_summary(geom = "errorbar",
               fun.min = ebbottom,
               fun.max = ebtop,
               position = position_dodge(0.9),
               width=0.2) + #scale_color_manual(values=col)+
  geom_point(aes(colour=setting, y=value, x=metric),
             position = position_jitterdodge(dodge.width = 0.9, 
                                             jitter.width = 0.1, 
                                             jitter.height = 0),
             show.legend = FALSE)+scale_color_manual(values=col)+
  labs(x="", y="value")#+coord_flip()

mytheme = theme_classic() + theme(axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14))+
  theme(axis.title.y= element_text(size=14))+theme(axis.title.x = element_text(size = 14))+
  theme(legend.title=element_text(size=12),legend.text=element_text(size=12), legend.position = "top")#+

p5_1=p5_1+mytheme
p5_1

write.table(cls_res_all, paste0(save_fold, 'Figure3A.csv'), col.names = T, row.names = F)

################# cls feature importance
data_fold<-paste0(dirname(rstudioapi::getSourceEditorContext()$path),'/data/')
save_fold<-paste0(dirname(rstudioapi::getSourceEditorContext()$path),'/save_data/')
acqs<-c('TMT', 'TMT', 'DDA', 'DDA', 'DIA', 'DIA')
feature_names<-c('DEA', 'exp_ty', 'MVI', 'norm')
setting<-c('MQ_TMT','FG_TMT','FG_DDA','MQ_DDA','DIANN_DIA','spt_DIA')
importances<-vector()

for (i in 1:length(platforms)){
  cls_impt<-read.table(paste0(data_fold, platforms[i], '_', acqs[i], '_', 'fea_importance.csv'),sep = ',', header = T)
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

importances_fig<-as.data.frame(importances)
importances_fig$importance<-round(as.numeric(importances_fig$FeatureImportance), 2)
importances_fig$l_y<-round(as.numeric(importances_fig$l_y), 2)
importances_fig$feature<-gsub('exp_ty','Matrix',importances_fig$feature)

sheet='feature_importance'
addWorksheet(wb2,sheet)
writeData(wb2, sheet, importances_fig, startRow = 1, startCol = 1, colNames = TRUE, rowNames = TRUE)
#saveWorkbook(wb2, file = paste0(save_fold,'10-fold_cv_res.xlsx'), overwrite = TRUE)

### Figure 3B
library(ggalluvial)
library(ggsci)

p52<-ggplot(importances_fig, aes(fill=feature, y=importance, x=setting,
                                 stratum = feature, alluvium = feature)) +
  geom_col(width = 0.4)+
  geom_alluvium( width = 0.4,alpha = 0.2,knot.pos = 0)+
  scale_fill_manual(values = pal_npg()(4))+
  scale_x_discrete(limits=c('FG_DDA', 'MQ_DDA', 'DIANN_DIA', 'spt_DIA', 'FG_TMT','MQ_TMT'))+
  labs(x="", y="feature importance")

mytheme = theme_classic() + theme(axis.text.x = element_text(size = 13, angle = -20),axis.text.y = element_text(size = 14))+
  theme(axis.title.y= element_text(size=14))+theme(axis.title.x = element_text(size = 14))+
  theme(legend.title=element_text(size=13),legend.text=element_text(size=13), legend.position = "top")

p52=p52+mytheme
p52

saveWorkbook(wb2, file = paste0(save_fold,'10-fold_cv_res.xlsx'), overwrite = TRUE)

write.table(importances_fig, paste0(save_fold, 'Figure3B.csv'), col.names = T, row.names = F)

#################### C linear model checking interactions between variables
library("readxl")
library(ggplot2)
library(reshape2)
library(ggpubr)
library(cowplot)
#save_fold<-'E:/proteomics/maus2/submission_NC/revision2/Figures/Figure3/'

volcano_lm<-function(rank_frag_DDA, platform, acq, inter_num, wb_lm,wb_anova){
  
  rank_frag_DDA$idx<-c(1:length(rank_frag_DDA$workflow))
  rank_frag_DDA<-rank_frag_DDA[order(rank_frag_DDA$avg_rank_mean),]
  rank_frag_DDA$rank_mean_exact<-c(1:length(rank_frag_DDA$workflow))
  rank_frag_DDA<-rank_frag_DDA[order(rank_frag_DDA$idx),]
  rank_frag_DDA$score<-length(rank_frag_DDA$workflow)-rank_frag_DDA$rank_mean_exact
  
  rank_frag_DDA$Imput<- gsub('blank','no_imp', rank_frag_DDA$Imput)
  rank_frag_DDA$normalization<- gsub('blank','no_norm', rank_frag_DDA$normalization)
  
  if(platform=='Maxquant' & acq=='TMT'){
    if(inter_num==1){
      model <- lm(score ~ DEA+Imput+normalization, data = rank_frag_DDA)
    }else if(inter_num==2){
      model <- lm(score ~ DEA+Imput+normalization+DEA*normalization+
                    DEA*Imput+Imput*normalization, data = rank_frag_DDA)
    }else if(inter_num==3){
      model <- lm(score ~ DEA+Imput+normalization+DEA*normalization+
                    DEA*Imput+Imput*normalization+
                    DEA*Imput*normalization, data = rank_frag_DDA)
    }
  }else{
    if(inter_num==1){
      model <- lm(score ~ DEA+Imput+normalization+Matrix, data = rank_frag_DDA)
    }else if(inter_num==2){
      model <- lm(score ~ DEA+Imput+normalization+Matrix+DEA*normalization+
                    DEA*Imput+Imput*normalization+DEA*Matrix+
                    Matrix*Imput+Matrix*normalization, data = rank_frag_DDA)
    }else if(inter_num==3){
      model <- lm(score ~ DEA+Imput+normalization+Matrix+DEA*normalization+
                    DEA*Imput+Imput*normalization+DEA*Matrix+
                    Matrix*Imput+Matrix*normalization+DEA*Imput*normalization+
                    DEA*Imput*Matrix+DEA*normalization*Matrix+Imput*normalization*Matrix,
                  data = rank_frag_DDA)
    }
  }
  lm_res<-summary(model)$coef
  anova_res<-anova(model)
  
  lm_res<-as.data.frame(lm_res)
  lm_res$variable<-row.names(lm_res)
  lm_res$log10Est<-sign(lm_res$Estimate)*log10(abs(lm_res$Estimate))
  lm_res$`Pr(>|t|)`[which(lm_res$`Pr(>|t|)`==0)]=min(lm_res$`Pr(>|t|)`[which(lm_res$`Pr(>|t|)`!=0)])
  lm_res$log10p<--log10(lm_res$`Pr(>|t|)`)
  
  lm_res$effection<-as.factor(ifelse(lm_res$`Pr(>|t|)` <= 0.01,
                                     ifelse(lm_res$Estimate>0,'positive','negative'),'NA'))
  lm_res$label<-lm_res$variable
  lm_res$label[which(lm_res$`Pr(>|t|)`>0.01)]=''
  
  sheet=paste0(platform,'_',acq,'_',as.character(inter_num))
  addWorksheet(wb_lm,sheet)
  writeData(wb_lm, sheet, lm_res, startRow = 1, startCol = 1, colNames = TRUE, rowNames = TRUE)
  
  addWorksheet(wb_anova,sheet)
  writeData(wb_anova, sheet, anova_res, startRow = 1, startCol = 1, colNames = TRUE, rowNames = TRUE)
  saveWorkbook(wb_lm, file = paste0(save_fold, 'linear_model_interaction/','lm_res.xlsx'), overwrite = TRUE)
  saveWorkbook(wb_anova, file = paste0(save_fold, 'linear_model_interaction/','anova_res.xlsx'), overwrite = TRUE)
  
  library(ggrepel)
  p6<-ggplot(lm_res, aes(x=log10Est, y=log10p, color=effection)) +
    geom_point(alpha=1, size=2) +
    theme_bw(base_size = 12) +
    xlab("Log10(Estimate)") +
    ylab("-Log10(pvalue) of linear model") +coord_flip()+
    scale_colour_manual(values = c('gray','steelblue','brown')) +
    geom_hline(yintercept = -log10(0.01), lty = 4) +
    geom_vline(xintercept = c(-log10(100), log10(100)), lty = 4)+
    geom_label_repel(data = lm_res, aes(label = label),
                     size = 3.5,#box.padding = unit(0.5, "lines"),
                     point.padding = unit(0.8, "lines"),
                     segment.color = "black",force=10,label.size = NA,
                     show.legend = FALSE, max.overlaps = 100)+
    theme(legend.title=element_text(size=14),legend.text=element_text(size=14), legend.position = "top")+
    guides(color=guide_legend(nrow=1), shape=guide_legend(nrow=3))
  
  return(list(p6=p6,lm_res=lm_res, anova_res=anova_res))
}

rank_frag_DDA<-read_excel(paste0(data_fold,'ranks_all_', 'FragPipe', '_', 'DDA','.xlsx'), sheet = 'ranking_all')
rank_mq_DDA<-read_excel(paste0(data_fold,'ranks_all_', 'Maxquant', '_', 'DDA','.xlsx'), sheet = 'ranking_all')
rank_frag_TMT<-read_excel(paste0(data_fold,'ranks_all_', 'FragPipe', '_', 'TMT','.xlsx'), sheet = 'ranking_all')
rank_mq_TMT<-read_excel(paste0(data_fold,'ranks_all_', 'Maxquant', '_', 'TMT','.xlsx'), sheet = 'ranking_all')
rank_diann_DIA<-read_excel(paste0(data_fold,'ranks_all_', 'DIANN', '_', 'DIA','.xlsx'), sheet = 'ranking_all')
rank_spt_DIA<-read_excel(paste0(data_fold,'ranks_all_', 'spt', '_', 'DIA','.xlsx'), sheet = 'ranking_all')

library(openxlsx)
wb_lm <- createWorkbook()
wb_anova<-createWorkbook()

for (inter_num in 1:1){
  p6_1<-volcano_lm(rank_frag_DDA,'FragPipe','DDA',inter_num,wb_lm, wb_anova)
  p6_2<-volcano_lm(rank_mq_DDA,'Maxquant','DDA',inter_num,wb_lm, wb_anova)
  p6_3<-volcano_lm(rank_frag_TMT,'FragPipe','TMT',inter_num,wb_lm, wb_anova)
  p6_4<-volcano_lm(rank_mq_TMT,'Maxquant','TMT',inter_num,wb_lm, wb_anova)
  p6_5<-volcano_lm(rank_diann_DIA,'DIANN','DIA',inter_num,wb_lm, wb_anova)
  p6_6<-volcano_lm(rank_spt_DIA,'spt','DIA',inter_num,wb_lm, wb_anova)
}
p6_1$p6
p6_2$p6
p6_3$p6
p6_4$p6
p6_5$p6
p6_6$p6



write.table(p6_1$lm_res, paste0(save_fold, 'Figure3C.csv'), col.names = T, row.names = F)


################ D frequent pattern mining from high-performing workflows

data_fold<-paste0(dirname(rstudioapi::getSourceEditorContext()$path),'/data/')
save_fold<-paste0(dirname(rstudioapi::getSourceEditorContext()$path),'/save_data/')

platforms<-c('Maxquant', 'Fragpipe', 'Fragpipe', 'Maxquant', 'DIANN', 'spt')
acqs<-c('TMT', 'TMT', 'DDA', 'DDA', 'DIA', 'DIA')
setting<-c('MQ_TMT','FG_TMT','FG_DDA','MQ_DDA','DIANN_DIA','spt_DIA')
FPs<-data.frame()
library(openxlsx)
wb4 <- createWorkbook()
for(i in 1:length(platforms)){
  hfp<-read.table(paste0(data_fold, platforms[i],'_',acqs[i],'FP_freitem_H.csv'),sep=',', header=F)
  colnames(hfp)<-c('support_ratio','pattern')
  hfp$pattern<-gsub('frozenset\\(\\{', '', hfp$pattern)
  hfp$pattern<-gsub('\\})', '', hfp$pattern)
  hfp$pattern<-gsub("\\'",'',gsub(', ', '&', hfp$pattern))
  hfp$setting<-setting[i]
  
  lfp<-read.table(paste0(data_fold, platforms[i],'_',acqs[i],'FP_freitem_L.csv'),sep=',', header=F)
  colnames(lfp)<-c('support_ratio','pattern')
  lfp$pattern<-gsub('frozenset\\(\\{', '', lfp$pattern)
  lfp$pattern<-gsub('\\})', '', lfp$pattern)
  lfp$pattern<-gsub("\\'",'',gsub(', ', '&', lfp$pattern))
  lfp$setting<-setting[i]
  
  FPs<-rbind(FPs, hfp)
  FPs$support_ratio<-as.numeric(FPs$support_ratio)
  
  sheet=paste0(setting[i], '_H')
  addWorksheet(wb4,sheet)
  writeData(wb4, sheet, hfp, startRow = 1, startCol = 1, colNames = TRUE, rowNames = TRUE)
  sheet=paste0(setting[i], '_L')
  addWorksheet(wb4,sheet)
  writeData(wb4, sheet, lfp, startRow = 1, startCol = 1, colNames = TRUE, rowNames = TRUE)
  saveWorkbook(wb4, file = paste0(save_fold, 'FPs_all.xlsx'), overwrite = TRUE)
}

FPs<-FPs[which(FPs$support_ratio>=0.3),]
FPs<-FPs[order(-FPs$support_ratio),]

FP_fig<-as.data.frame(FPs)
FP_fig$pattern<-FP_fig$pattern
FP_fig$support_ratio<-round(FP_fig$support_ratio, 2)
FP_fig$idx<-c(1:length(FP_fig$setting))

write.table(FP_fig, paste0(save_fold, 'FP_plot.csv'), sep = ',', col.names = T, row.names = F)

col=c('#427AB2','#F09148','#FF9896','#DBDB8D','#C59D94','#b2df8a')

library(ggthemes)
p7<-ggbarplot(FP_fig, x="idx", y="support_ratio",
              fill="setting", color = NA, rotate=TRUE,
              palette=col,xlab='',
              ylab='support_ratio', sort.val='asc',
              group="setting", sort.by.groups = TRUE,
              label = TRUE,lab.pos = "out",lab.vjust = 0.5, lab.hjust = -0.2
) + scale_y_continuous(expand = c(0,0),
                       limits = c(0,0.85))
p71<-p7+
  scale_x_discrete(labels=c(FP_fig$pattern[as.numeric(ggplot_build(p7)$layout$panel_params[[1]]$y$get_labels())])) +
  geom_vline(xintercept = c(4.5,7.5,11.5,12.5,16.5),
             linetype='dashed')
mytheme = theme_classic() + theme(axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14))+
  theme(axis.title.y= element_text(size=14))+theme(axis.title.x = element_text(size = 14))+
  theme(legend.title=element_text(size=12),legend.text=element_text(size=12), legend.position = "right")
p71=p71+mytheme

p71

write.table(FP_fig, paste0(save_fold, 'Figure3D.csv'), col.names = T, row.names = F)


