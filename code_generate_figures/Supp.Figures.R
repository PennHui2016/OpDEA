######################
####
#### source code for generating Supplementary Figures in Supplementary Information
#### authored by Hui Peng: hui.peng@ntu.edu.sg
####
####
#############################

###############################
## supp. Fig. 1 The spearman correlations of the leave-one-dataset-out 
## cross-validations using mean or median performance for workflow benchmarking 
library(openxlsx)
wb_all<-createWorkbook()

data_fold<-paste0(dirname(rstudioapi::getSourceEditorContext()$path),'/data/')
save_fold<-paste0(dirname(rstudioapi::getSourceEditorContext()$path),'/save_data/')

rank_values<-function(rk_values, all_values){
  s<-data.frame(all_values$workflow, rk_values, c(1:length(all_values$workflow)))
  colnames(s)<-c('workflow','value', 'idx')
  s<-s[order(-s$value),]
  s$rank=c(1:length(s$value))
  s<-s[order(s$idx),]
  return(s)
}

spr_test<-function(platform, acq){
  library("readxl")
  
  if(acq=='DIA'){
    st=7
    ed=23
  }else if(acq=='DDA'){
    st=7
    ed=28
  }else if(acq=='TMT'){
    st=7
    ed=21
  }
  pauc001<-read_excel(paste0(data_fold, 'metrics_',platform, '_',acq, '.xlsx'), sheet='pAUC0.01')
  pauc005<-read_excel(paste0(data_fold, 'metrics_',platform, '_',acq, '.xlsx'), sheet='pAUC0.05')
  pauc01<-read_excel(paste0(data_fold, 'metrics_',platform, '_',acq, '.xlsx'), sheet='pAUC0.1')
  nMCC<-read_excel(paste0(data_fold, 'metrics_',platform, '_',acq, '.xlsx'), sheet='nMCC0.05')
  Gmean<-read_excel(paste0(data_fold, 'metrics_',platform, '_',acq, '.xlsx'), sheet='G-mean0.05')
  
  header<-read_excel(paste0(data_fold, 'metrics_',platform, '_',acq, '.xlsx'), sheet='header', col_names = FALSE)
  header_dataset<-paste0(header[1,],'_', header[2,])
  dts<-header[1,][st:ed]
  dt_n<-c()
  
  for(i in dts[1,]){
    s = strsplit(i, '_', fixed = TRUE)[[1]]
    if(length(s)>2){
      dt=s[1]
      for (j in 2:(length(s)-1)) {
        dt = paste0(dt, '_', s[j])
      }
      dt_n<-c(dt_n, dt)
    }else if(length(s)==2){
      dt_n<-c(dt_n, s[1])
    }
  }
  
  sprs<-vector()
  tr_ranks_mean<-vector()
  tr_ranks_median<-vector()
  te_ranks<-vector()
  train_dt<-c()
  test_dt<-c()
  
  for (dt in unique(dt_n)) {
    train_dt<-c(train_dt, dt)
    dt_pauc001=pauc001[,st:ed]
    colnames(dt_pauc001)<-header_dataset[st:ed]
    
    dt_pauc005=pauc005[,st:ed]
    colnames(dt_pauc005)<-header_dataset[st:ed]
    
    dt_pauc01=pauc01[,st:ed]
    colnames(dt_pauc01)<-header_dataset[st:ed]
    
    dt_nMCC=nMCC[,st:ed]
    colnames(dt_nMCC)<-header_dataset[st:ed]
    
    dt_Gmean=Gmean[,st:ed]
    colnames(dt_Gmean)<-header_dataset[st:ed]
    
    mean_tr_pauc001<-rowMeans(dt_pauc001[,which(dt_n!=dt)])
    rank_mean_tr_pauc001<-rank_values(mean_tr_pauc001, pauc001)
    median_tr_pauc001<-apply(dt_pauc001[,which(dt_n!=dt)], 1, median)
    rank_median_tr_pauc001<-rank_values(median_tr_pauc001, pauc001)
    
    mean_tr_pauc005<-rowMeans(dt_pauc005[,which(dt_n!=dt)])
    rank_mean_tr_pauc005<-rank_values(mean_tr_pauc005, pauc005)
    median_tr_pauc005<-apply(dt_pauc005[,which(dt_n!=dt)], 1, median)
    rank_median_tr_pauc005<-rank_values(median_tr_pauc005, pauc005)
    
    mean_tr_pauc01<-rowMeans(dt_pauc01[,which(dt_n!=dt)])
    rank_mean_tr_pauc01<-rank_values(mean_tr_pauc01, pauc01)
    median_tr_pauc01<-apply(dt_pauc01[,which(dt_n!=dt)], 1, median)
    rank_median_tr_pauc01<-rank_values(median_tr_pauc01, pauc01)
    
    mean_tr_nMCC<-rowMeans(dt_nMCC[,which(dt_n!=dt)])
    rank_mean_tr_nMCC<-rank_values(mean_tr_nMCC, nMCC)
    median_tr_nMCC<-apply(dt_nMCC[,which(dt_n!=dt)], 1, median)
    rank_median_tr_nMCC<-rank_values(median_tr_nMCC, nMCC)
    
    mean_tr_Gmean<-rowMeans(dt_Gmean[,which(dt_n!=dt)])
    rank_mean_tr_Gmean<-rank_values(mean_tr_Gmean, Gmean)
    median_tr_Gmean<-apply(dt_Gmean[,which(dt_n!=dt)], 1, median)
    rank_median_tr_Gmean<-rank_values(median_tr_Gmean, Gmean)
    
    ranks_avg_mean<-rowMeans(cbind(rank_mean_tr_pauc001$rank, rank_mean_tr_pauc005$rank,
                                   rank_mean_tr_pauc01$rank, rank_mean_tr_nMCC$rank,
                                   rank_mean_tr_Gmean$rank))
    
    ranks_avg_mean1<-data.frame(avg_rank=ranks_avg_mean,
                                idx=c(1:length(ranks_avg_mean)))
    ranks_avg_mean1<-ranks_avg_mean1[order(ranks_avg_mean1$avg_rank),]
    ranks_avg_mean1$final_rank<-c(1:length(ranks_avg_mean1$avg_rank))
    ranks_avg_mean1<-ranks_avg_mean1[order(ranks_avg_mean1$idx),]
    
    ranks_avg_median<-rowMeans(cbind(rank_median_tr_pauc001$rank, rank_median_tr_pauc005$rank,
                                     rank_median_tr_pauc01$rank, rank_median_tr_nMCC$rank,
                                     rank_median_tr_Gmean$rank))
    ranks_avg_median1<-data.frame(avg_rank=ranks_avg_median,
                                  idx=c(1:length(ranks_avg_median)))
    
    ranks_avg_median1<-ranks_avg_median1[order(ranks_avg_median1$avg_rank),]
    ranks_avg_median1$final_rank<-c(1:length(ranks_avg_median1$avg_rank))
    ranks_avg_median1<-ranks_avg_median1[order(ranks_avg_median1$idx),]
    
    tr_ranks_mean<-cbind(tr_ranks_mean, ranks_avg_mean1$final_rank)
    tr_ranks_median<-cbind(tr_ranks_median, ranks_avg_median1$final_rank)
    
    tests_pauc001<-dt_pauc001[,which(dt_n==dt)]
    tests_pauc005<-dt_pauc005[,which(dt_n==dt)]
    tests_pauc01<-dt_pauc01[,which(dt_n==dt)]
    tests_nMCC<-dt_nMCC[,which(dt_n==dt)]
    tests_Gmean<-dt_Gmean[,which(dt_n==dt)]
    
    for (i in 1:length(tests_pauc001[1,])) {
      test_dt<-c(test_dt, colnames(tests_pauc001)[i])
      rk_te_pauc001<-rank_values(tests_pauc001[,i], pauc001)
      rk_te_pauc005<-rank_values(tests_pauc005[,i], pauc005)
      rk_te_pauc01<-rank_values(tests_pauc01[,i], pauc01)
      rk_te_Gmean<-rank_values(tests_Gmean[,i], Gmean)
      rk_te_nMCC<-rank_values(tests_nMCC[,i], nMCC)
      
      rk_avg_te<-rowMeans(cbind(rk_te_pauc001$rank, rk_te_pauc005$rank,
                                rk_te_pauc01$rank, rk_te_Gmean$rank,
                                rk_te_nMCC$rank))
      
      te_ranks<-cbind(te_ranks, rk_avg_te)
      
      spr1 = cor(ranks_avg_mean1$final_rank, rk_avg_te, method = 'spearman')
      pr1=cor(ranks_avg_mean1$final_rank, rk_avg_te, method = 'pearson')
      spr2 = cor(ranks_avg_median1$final_rank, rk_avg_te, method = 'spearman')
      pr2 = cor(ranks_avg_median1$final_rank, rk_avg_te, method = 'pearson')
      sprs<-rbind(sprs, c(dt, 'avg_rank', i, spr1, spr2, pr1, pr2))
      
    }
  }
  colnames(tr_ranks_mean)<-train_dt
  colnames(tr_ranks_median)<-train_dt
  colnames(te_ranks)<-test_dt
  colnames(sprs)<-c('dt_name', 'metric', 'dt_col', 'mean_spr', 'median_spr', 'mean_pr', 'median_pr')
  #write.table(sprs, paste0(save_fold, platform, '-', acq, '_spr_avg_rank_test.csv'), sep = ',', col.names = TRUE, row.names = FALSE)
  return(list(tr_rk_mean=tr_ranks_mean, tr_rk_median=tr_ranks_median, te_ranks=te_ranks, sprs=sprs))
}

lopocv_fg_dda<-spr_test('Fragpipe', 'DDA')
lopocv_mq_dda<-spr_test('Maxquant', 'DDA')
lopocv_fg_tmt<-spr_test('Fragpipe', 'TMT')
lopocv_mq_tmt<-spr_test('Maxquant', 'TMT')
lopocv_diann_dia<-spr_test('DIANN', 'DIA')
lopocv_spt_dia<-spr_test('spt', 'DIA')

sprs_fg_DDA<-as.data.frame(lopocv_fg_dda$sprs)
sprs_mq_DDA<-as.data.frame(lopocv_mq_dda$sprs)
sprs_fg_TMT<-as.data.frame(lopocv_fg_tmt$sprs)
sprs_mq_TMT<-as.data.frame(lopocv_mq_tmt$sprs)
sprs_diann_DIA<-as.data.frame(lopocv_diann_dia$sprs)
sprs_spt_DIA<-as.data.frame(lopocv_spt_dia$sprs)

sprs_fg_DDA$setting<-c(rep('FG_DDA',length(sprs_fg_DDA$dt_name)))
sprs_mq_DDA$setting<-c(rep('MQ_DDA',length(sprs_mq_DDA$dt_name)))
sprs_fg_TMT$setting<-c(rep('FG_TMT',length(sprs_fg_TMT$dt_name)))
sprs_mq_TMT$setting<-c(rep('MQ_TMT',length(sprs_mq_TMT$dt_name)))
sprs_diann_DIA$setting<-c(rep('DIANN_DIA',length(sprs_diann_DIA$dt_name)))
sprs_spt_DIA$setting<-c(rep('spt_DIA',length(sprs_spt_DIA$dt_name)))

all_sprs<-rbind(sprs_fg_DDA,sprs_mq_DDA, sprs_fg_TMT, sprs_mq_TMT,
                sprs_diann_DIA, sprs_spt_DIA)

sprs_dat<-cbind(all_sprs$mean_spr,all_sprs$setting,
                c(rep('mean',length(all_sprs$dt_name))))
sprs_dat<-as.data.frame(rbind(sprs_dat,
                              cbind(all_sprs$median_spr,all_sprs$setting,
                                    c(rep('median',length(all_sprs$dt_name))))))
colnames(sprs_dat)<-c('spearman','setting','method')
sprs_dat$spearman<-as.numeric(sprs_dat$spearman)

col=c('#427AB2','#b2df8a','#FF9896','#DBDB8D','#C59D94','#b2df8a')
p1 = ggplot(sprs_dat, aes(x=setting, y=spearman, fill=method)) +
  geom_boxplot(position = position_dodge(width = 0.75), width=0.5)+ geom_jitter(position = position_jitterdodge(jitter.width=0.1))+
  #geom_boxplot(alpha=1,outlier.size=0, size=0.3, width=0.3,fill="white") +
  scale_fill_manual(values = col)+
  labs(x="", y="spearman correlation", color=factor)+
  scale_x_discrete(limits=unique(all_sprs$setting))+stat_summary(fun=mean,
                                                                 geom="point",
                                                                 shape=17,
                                                                 aes(group=method),
                                                                 position=position_dodge(.75), size=2, color="red", fill="red")
mytheme = theme_classic() + theme(axis.text.x = element_text(size = 16, angle = -20),axis.text.y = element_text(size = 10))+
  theme(axis.title.y= element_text(size=16))+theme(axis.title.x = element_text(size = 14))+
  theme(legend.title=element_text(size=16),legend.text=element_text(size=16), legend.position = "top")
p1=p1+mytheme
p1

sheet = 'supp. Fig.1'
addWorksheet(wb_all, sheet)
writeData(wb_all, sheet, sprs_dat, startRow = 2, startCol = 2)
###############################
## supp. Fig. 2 Top 2 workflows for each matrix type of settings FG_TMT and MQ_TMT

data_fold<-paste0(dirname(rstudioapi::getSourceEditorContext()$path),'/data/')
save_fold<-paste0(dirname(rstudioapi::getSourceEditorContext()$path),'/save_data/')

library(readxl)

rank_frag_TMT<-read_excel(paste0(data_fold,'ranks_all_', 'FragPipe', '_', 'TMT','.xlsx'), sheet = 'ranking_all')
rank_mq_TMT<-read_excel(paste0(data_fold,'ranks_all_', 'Maxquant', '_', 'TMT','.xlsx'), sheet = 'ranking_all')

rank_frag_TMT<-rank_frag_TMT[order(rank_frag_TMT$avg_rank_mean),]
rank_mq_TMT<-rank_mq_TMT[order(rank_mq_TMT$avg_rank_mean),]

plot_topwf<-function(i, cols){
  
  rank_res<-read_excel(paste0(data_fold,'ranks_all_', platforms[i], '_', acqs[i],'.xlsx'), sheet = 'ranking_all')
  rank_res<-rank_res[order(rank_res$avg_rank_mean),]
  rank_res$exact_rank<-c(1:length(rank_res$workflow))
  rank_res<-rank_res[which(rank_res$DEA!='MSstats'),]
  uni_matrix<-unique(rank_res$Matrix)
  top_ranks<-data.frame()
  for (j in 1:length(uni_matrix)) {
    idx_mat<-which(rank_res$Matrix==uni_matrix[j])
    rank_idx_mat<-rank_res[idx_mat,]
    rank_idx_mat<-rank_idx_mat[order(rank_idx_mat$avg_rank_mean),]
    rank_idx_mat$setting<-setting[i]
    top_rank<-rank_idx_mat[1:N,]
    top_rank$rank<-paste0('top',c(1:N))
    top_rank$Matrix<-gsub('intensity','reporter',top_rank$Matrix)
    top_rank$Imput<-gsub('blank','none',top_rank$Imput)
    top_rank$normalization<-gsub('blank','none',top_rank$normalization)
    
    top_rank$label<-paste0(top_rank$DEA,'|',top_rank$normalization,'|',top_rank$Imput,
                           '(',as.character(paste0(top_rank$exact_rank)),')')
    
    top_ranks<-rbind(top_ranks, top_rank)
  }
  library(ggrepel)
  p8<-ggplot(top_ranks, aes(x=mean_pauc001, y=mean_gmean005, color=Matrix)) +
    geom_point(alpha=1, size=2) +
    theme_bw(base_size = 12) +
    xlab("pAUC(0.01)") +
    ylab("G-mean") +
    labs(title = setting[i])+scale_colour_manual(values = c(cols))+
    geom_label_repel(data = top_ranks, aes(label = label),
                     size = 4,#box.padding = unit(0.5, "lines"),
                     label.size = NA, 
                     point.padding = unit(0.8, "lines"),
                     segment.color = "black",force=1000,
                     show.legend = FALSE, max.overlaps = 1000)+
    theme(legend.title=element_text(size=12),legend.text=element_text(size=12), legend.position = "top")+
    guides(color=guide_legend(nrow=1), shape=guide_legend(nrow=3))
  p8
  return(list(fig=p8, dt=top_ranks))
}

platforms<-c('Maxquant', 'Fragpipe')
acqs<-c('TMT', 'TMT')
setting<-c('MQ_TMT','FG_TMT')

DEA_sym<-c(1:10)
norm_sym<-LETTERS
imp_sym<-letters
matrix_sym<-c('alpha', 'beta', 'gamma', 'delta', 'epsilon')

N=2
col=c("#1f78b4","#33a02c","#fb9a99","orange","black","purple","green","darkgray","#a6cee3","red")
p81<-plot_topwf(1, col[7])
p82<-plot_topwf(2, col[8:10])

p81$fig
p82$fig

sheet = 'supp. Fig.2'
addWorksheet(wb_all, sheet)
writeData(wb_all, sheet, rbind(p81$dt, p82$dt), startRow = 2, startCol = 2)
##############################
## supp. Fig. 3 Comparison of performance metric values of workflows 
## ranked at top 5%, 5%-25%, 25%-50% and bottom 50% under setting FG_DDA.

library("readxl")
library(ggplot2)
library(reshape2)
library(ggpubr)
library(cowplot)
data_fold<-paste0(dirname(rstudioapi::getSourceEditorContext()$path),'/data/')
save_fold<-paste0(dirname(rstudioapi::getSourceEditorContext()$path),'/save_data/')
platform='FragPipe'
acq='DDA'
ranking_all<-read_excel(paste0(data_fold, 'ranks_all_', platform, '_', acq, '.xlsx'), sheet = 'ranking_all')
ranking_all<-ranking_all[order(ranking_all$avg_rank_mean),]
ranking_all$label<-c(rep('top5%',as.integer(length(ranking_all$workflow))))
ranking_all$label[(round(length(ranking_all$workflow)*0.05)+1):round(length(ranking_all$workflow)*0.25)] = '5-25%'
ranking_all$label[(round(length(ranking_all$workflow)*0.25)+1):round(length(ranking_all$workflow)*0.5)] = '25-50%'
ranking_all$label[(round(length(ranking_all$workflow)*0.5)+1):as.integer(length(ranking_all$workflow))] = '50-100%'
ranking_all<-ranking_all[which(ranking_all$mean_pauc001!=0 & ranking_all$median_pauc001!=0),]

mean_performances<-cbind(ranking_all$mean_pauc001, ranking_all$mean_pauc005, ranking_all$mean_pauc01,
                         ranking_all$mean_nmcc005, ranking_all$mean_gmean005)
colnames(mean_performances)<-c('pAUC(0.01)', 'pAUC(0.05)', 'pAUC(0.1)', 'nMCC','G-mean')
dat_mean<-melt(mean_performances)
colnames(dat_mean)<-c('idx','metrics','performance')

dat_mean$type<-'mean'
dat<-dat_mean

####### top5% vs. top5-25% vs. top25-50% vs bottom
col=c("#1f78b4","#b2df8a","#fb9a99","#fdbf6f")
dat_l<-rbind(cbind(mean_performances[,1],rep('pAUC(0.01)', length(mean_performances[,1])),ranking_all$label),
             cbind(mean_performances[,2],rep('pAUC(0.05)', length(mean_performances[,1])),ranking_all$label),
             cbind(mean_performances[,3],rep('pAUC(0.1)', length(mean_performances[,1])),ranking_all$label),
             cbind(mean_performances[,4],rep('nMCC', length(mean_performances[,1])),ranking_all$label),
             cbind(mean_performances[,5],rep('G-mean', length(mean_performances[,1])),ranking_all$label)
)
colnames(dat_l)<-c('performance','metrics','label')
dat_l<-as.data.frame(dat_l)
dat_l$performance<-as.numeric(dat_l$performance)
p3 = ggplot(dat_l, aes(x=metrics, y=performance, fill=label)) +
  geom_violin(position = position_dodge(width = 0.75), scale = 'width')+
  scale_fill_manual(values = col)+
  labs(x="", y="performance", color=factor)+
  scale_x_discrete(limits=c("pAUC(0.01)","pAUC(0.05)","pAUC(0.1)",'nMCC', 'G-mean'))
mytheme = theme_classic() + theme(axis.text.x = element_text(size = 16, angle = -20),axis.text.y = element_text(size = 16))+
  theme(axis.title.y= element_text(size=16))+theme(axis.title.x = element_text(size = 16))+
  theme(legend.title=element_text(size=16),legend.text=element_text(size=16), legend.position = "top")
p3=p3+mytheme
p3

sheet = 'supp. Fig.3'
addWorksheet(wb_all, sheet)
writeData(wb_all, sheet, dat_l, startRow = 2, startCol = 2)
##############################
## supp. Fig. 4 Linear model fitting for checking affections of options in each
## workflow on the workflow ranking score

library("readxl")
library(ggplot2)
library(reshape2)
library(ggpubr)
library(cowplot)
data_fold<-paste0(dirname(rstudioapi::getSourceEditorContext()$path),'/data/')
save_fold<-paste0(dirname(rstudioapi::getSourceEditorContext()$path),'/save_data/')

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

data_fold<-paste0(dirname(rstudioapi::getSourceEditorContext()$path),'/data/')
save_fold<-paste0(dirname(rstudioapi::getSourceEditorContext()$path),'/save_data/')

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


sheet = 'supp. Fig.41'
addWorksheet(wb_all, sheet)
writeData(wb_all, sheet, p6_1$lm_res, startRow = 2, startCol = 2)

sheet = 'supp. Fig.42'
addWorksheet(wb_all, sheet)
writeData(wb_all, sheet, p6_2$lm_res, startRow = 2, startCol = 2)

sheet = 'supp. Fig.43'
addWorksheet(wb_all, sheet)
writeData(wb_all, sheet, p6_3$lm_res, startRow = 2, startCol = 2)

sheet = 'supp. Fig.44'
addWorksheet(wb_all, sheet)
writeData(wb_all, sheet, p6_4$lm_res, startRow = 2, startCol = 2)

sheet = 'supp. Fig.45'
addWorksheet(wb_all, sheet)
writeData(wb_all, sheet, p6_5$lm_res, startRow = 2, startCol = 2)

sheet = 'supp. Fig.46'
addWorksheet(wb_all, sheet)
writeData(wb_all, sheet, p6_6$lm_res, startRow = 2, startCol = 2)
##############################
## supp. Fig. 5 comparison of DEPs detected by top ranked workflows
library(readxl)
library(UpSetR)
data_fold<-paste0(dirname(rstudioapi::getSourceEditorContext()$path),'/data/')
save_fold<-paste0(dirname(rstudioapi::getSourceEditorContext()$path),'/save_data/')

top1wf<-function(dt, platform, acq, mat_ty, rk){
  base_fold<-data_fold
  dea_res_fold<-paste0(base_fold, acq, '/',platform, '/',dt,'/')
  ranking_all_FG_DDA<-read_excel(paste0(data_fold, 'ranks_all_', platform, '_', acq, '.xlsx'), sheet = 'ranking_all')
  
  if(mat_ty==''){
    wfs<-ranking_all_FG_DDA[order(ranking_all_FG_DDA$avg_rank_mean),]
    wftop1_FG_DDA<-wfs[rk,]
  }else{
    wfs<-ranking_all_FG_DDA[which(ranking_all_FG_DDA$Matrix==mat_ty),]
    wfs=wfs[order(wfs$avg_rank_mean),]
    wftop1_FG_DDA<-wfs[rk,]
  }
 
  print(paste0(dea_res_fold,
               dt,'_',wftop1_FG_DDA$DEA,'_',
               wftop1_FG_DDA$Platform,'_',
               wftop1_FG_DDA$Matrix,'_',gsub('blank','',wftop1_FG_DDA$Imput),
               '_',gsub('blank','',wftop1_FG_DDA$normalization),'.csv'))
  dea_res_FG_DDA_wf1<-read.table(paste0(dea_res_fold,
                                        dt,'_',wftop1_FG_DDA$DEA,'_',
                                        wftop1_FG_DDA$Platform,'_',
                                        wftop1_FG_DDA$Matrix,'_',gsub('blank','',wftop1_FG_DDA$Imput),
                                        '_',gsub('blank','',wftop1_FG_DDA$normalization),'.csv'),
                                 sep = ',',header = T)
  return(list(wf=wftop1_FG_DDA, dea_res=dea_res_FG_DDA_wf1))
}

all_protein_FG_DDA<-read.table(paste0(data_fold,'HYEtims735_LFQ_FragPipe_all_proteins.tsv'), sep = '\t', header = T)
top0<-top1wf('HYEtims735_LFQ', 'FragPipe', 'DDA','top0',1)
top3<-top1wf('HYEtims735_LFQ', 'FragPipe', 'DDA','top3',1)
count<-top1wf('HYEtims735_LFQ', 'FragPipe', 'DDA','count',1)
dlfq<-top1wf('HYEtims735_LFQ', 'FragPipe', 'DDA','dlfq',1)
lfq<-top1wf('HYEtims735_LFQ', 'FragPipe', 'DDA','LFQ',1)

all_top1<-top1wf('HYEtims735_LFQ', 'FragPipe', 'DDA','',1)
all_top2<-top1wf('HYEtims735_LFQ', 'FragPipe', 'DDA','',2)
all_top3<-top1wf('HYEtims735_LFQ', 'FragPipe', 'DDA','',3)
all_top4<-top1wf('HYEtims735_LFQ', 'FragPipe', 'DDA','',4)
all_top5<-top1wf('HYEtims735_LFQ', 'FragPipe', 'DDA','',5)

top0_DEP<-top0$dea_res$protein[which(abs(top0$dea_res$logFC)>=log2(1.5) & top0$dea_res$adj.pvalue<=0.05)]
top3_DEP<-top3$dea_res$protein[which(abs(top3$dea_res$logFC)>=log2(1.5) & top3$dea_res$adj.pvalue<=0.05)]
count_DEP<-count$dea_res$protein[which(abs(count$dea_res$logFC)>=log2(1.5) & count$dea_res$adj.pvalue<=0.05)]
dlfq_DEP<-dlfq$dea_res$protein[which(abs(dlfq$dea_res$logFC)>=log2(1.5) & dlfq$dea_res$adj.pvalue<=0.05)]
lfq_DEP<-lfq$dea_res$protein[which(abs(lfq$dea_res$logFC)>=log2(1.5) & lfq$dea_res$adj.pvalue<=0.05)]

all_top1_DEP<-all_top1$dea_res$protein[which(abs(all_top1$dea_res$logFC)>=log2(1.5) & all_top1$dea_res$adj.pvalue<=0.05)]
all_top2_DEP<-all_top2$dea_res$protein[which(abs(all_top2$dea_res$logFC)>=log2(1.5) & all_top2$dea_res$adj.pvalue<=0.05)]
all_top3_DEP<-all_top3$dea_res$protein[which(abs(all_top3$dea_res$logFC)>=log2(1.5) & all_top3$dea_res$adj.pvalue<=0.05)]
all_top4_DEP<-all_top4$dea_res$protein[which(abs(all_top4$dea_res$logFC)>=log2(1.5) & all_top4$dea_res$adj.pvalue<=0.05)]
all_top5_DEP<-all_top5$dea_res$protein[which(abs(all_top5$dea_res$logFC)>=log2(1.5) & all_top5$dea_res$adj.pvalue<=0.05)]



T_DEP<-all_protein_FG_DDA$Protein[which(all_protein_FG_DDA$Organism=='YEAST' | 
                                          all_protein_FG_DDA$Organism=='ECOLI')]

Recall_top0 = length(intersect(top0_DEP,T_DEP))/length(T_DEP)
Recall_top3 = length(intersect(top3_DEP,T_DEP))/length(T_DEP)
Recall_count = length(intersect(count_DEP,T_DEP))/length(T_DEP)
Recall_dlfq = length(intersect(dlfq_DEP,T_DEP))/length(T_DEP)
Recall_lfq = length(intersect(lfq_DEP,T_DEP))/length(T_DEP)

Recall_Top1 = length(intersect(all_top1_DEP,T_DEP))/length(T_DEP)
Recall_Top2 = length(intersect(all_top2_DEP,T_DEP))/length(T_DEP)
Recall_Top3 = length(intersect(all_top3_DEP,T_DEP))/length(T_DEP)
Recall_Top4 = length(intersect(all_top4_DEP,T_DEP))/length(T_DEP)
Recall_Top5 = length(intersect(all_top5_DEP,T_DEP))/length(T_DEP)

DEPs_mat<- list(TOP1 = all_top1_DEP,
                TOP2 = all_top2_DEP,
                TOP3 = all_top3_DEP,
                TOP4 = all_top4_DEP,
                TOP5 = all_top5_DEP,
                T_DEP = T_DEP)
m=fromList(DEPs_mat)
all = unique(c(all_top1_DEP, all_top2_DEP, all_top3_DEP, all_top4_DEP,
               all_top4_DEP, all_top5_DEP, T_DEP))

out_DEPs<-data.frame(protein=all)
out_DEPs$Top1=c(rep(0, length(all)))
out_DEPs$Top1[match(all_top1_DEP, all)]=1

out_DEPs$Top2=c(rep(0, length(all)))
out_DEPs$Top2[match(all_top2_DEP, all)]=1

out_DEPs$Top3=c(rep(0, length(all)))
out_DEPs$Top3[match(all_top3_DEP, all)]=1

out_DEPs$Top4=c(rep(0, length(all)))
out_DEPs$Top4[match(all_top4_DEP, all)]=1

out_DEPs$Top5=c(rep(0, length(all)))
out_DEPs$Top5[match(all_top5_DEP, all)]=1

out_DEPs$T_DEP=c(rep(0, length(all)))
out_DEPs$T_DEP[match(T_DEP, all)]=1

sheet = 'supp. Fig.51'
addWorksheet(wb_all, sheet)
writeData(wb_all, sheet, out_DEPs, startRow = 2, startCol = 2)
  
UpSetR::upset(m,
              sets = c("TOP1", "TOP2", "TOP3", "TOP4", "TOP5", "T_DEP"),
              order.by="freq", matrix.color="black", point.size=3,
              nintersects = 15,
              sets.bar.color=c("#15d08a","#f2ca01","#a45ee3","#1f11b4","#c3642c",
                               "#f078b4"),mb.ratio=c(0.5, 0.5), text.scale = 1.5,
              number.angles =0)

DEPs_mat<- list(top0 = top0_DEP,
                top3 = top3_DEP,
                count = count_DEP,
                dlfq = dlfq_DEP,
                MaxLFQ = lfq_DEP,
                T_DEP = T_DEP)
m = fromList(DEPs_mat)

all = unique(c(top0_DEP, top3_DEP, count_DEP, dlfq_DEP,
               lfq_DEP, T_DEP))

out_DEPs<-data.frame(protein=all)
out_DEPs$top0=c(rep(0, length(all)))
out_DEPs$top0[match(top0_DEP, all)]=1

out_DEPs$top3=c(rep(0, length(all)))
out_DEPs$top3[match(top3_DEP, all)]=1

out_DEPs$count=c(rep(0, length(all)))
out_DEPs$count[match(count_DEP, all)]=1

out_DEPs$dlfq=c(rep(0, length(all)))
out_DEPs$dlfq[match(dlfq_DEP, all)]=1

out_DEPs$lfq=c(rep(0, length(all)))
out_DEPs$lfq[match(lfq_DEP, all)]=1

out_DEPs$T_DEP=c(rep(0, length(all)))
out_DEPs$T_DEP[match(T_DEP, all)]=1

UpSetR::upset(m,
              sets = c("top0", "top3", "count", "dlfq", "MaxLFQ", "T_DEP"),
              order.by="freq", matrix.color="black", point.size=3,
              nintersects = 20,
              sets.bar.color=c("#15d08a","#f2ca01","#a45ee3","#1f11b4","#c3642c",
                               "#f078b4"),mb.ratio=c(0.5, 0.5), text.scale = 1.5,
              number.angles =0)



sheet = 'supp. Fig.52'
addWorksheet(wb_all, sheet)
writeData(wb_all, sheet, out_DEPs, startRow = 2, startCol = 2)
##############################
## supp. Fig. 6 comparison of DEPs detected by top ranked workflows
################# HEqe408


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

dt='HEqe408'

all_protein_FG_DDA<-read.table(paste0(data_fold, 'HEqe408_LFQ_FragPipe_all_proteins.tsv'), sep = '\t', header = T)
all_protein_MQ_DDA<-read.table(paste0(data_fold, 'HEqe408_LFQ_Maxquant_all_proteins.tsv'), sep = '\t', header = T)
all_protein_DIANN_DIA<-read.table(paste0(data_fold,'HEqe408_DIA_DIANN_all_proteins.tsv'), sep = '\t', header = T)
all_protein_spt_DIA<-read.table(paste0(data_fold,'HEqe408_DIA_spt_all_proteins.tsv'), sep = '\t', header = T)

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

DEPs$T_ECOLI=c(rep(0, length(DEPs$Proteins)))
DEPs$T_ECOLI[match(T_ECOLI, DEPs$Proteins)]=1

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
wb=createWorkbook()
sheet=paste(dt, '_cp_DEPs')

addWorksheet(wb,sheet)
writeData(wb, sheet, DEPs, startRow = 1, startCol = 1, colNames = TRUE, rowNames = TRUE)
#saveWorkbook(wb, file = paste0(save_fold, 'cross_platform_compare.xlsx'), overwrite = TRUE)

library(pROC)
library(ggpubr)

cross_platform_proteins<-data.frame(Proteins=DEPs$Proteins, labels=0)
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
#saveWorkbook(wb, file = paste0(save_fold, 'cross_platform_compare.xlsx'), overwrite = TRUE)

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
#saveWorkbook(wb, file = paste0(save_fold, 'cross_platform_compare.xlsx'), overwrite = TRUE)
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


ggbarplot(nums, "method", "value", group = "variable", color = "variable",
          fill="variable", #label = T,lab.pos='in',
          platette = c('#5494cc','#e18283','#0d898a','#f9cc52'),
) + geom_text(aes(y=l_y, label=label))+
  scale_x_discrete(limits=c('FG_DDA_ens','MQ_DDA_ens','DIANN_DIA_ens', 'spt_DIA_ens',
                            'FG_DDA','MQ_DDA','DIANN_DIA','spt_DIA')) +
  geom_vline(xintercept = 4.5, linetype=2, show.legend=FALSE) + theme(axis.line.y = element_blank())+
  coord_flip()+theme(axis.title.y = element_blank())

sheet = 'supp. Fig.61'
addWorksheet(wb_all, sheet)
writeData(wb_all, sheet, nums, startRow = 2, startCol = 2)
###################################################
################# HYtims134


data_fold<-paste0(dirname(rstudioapi::getSourceEditorContext()$path),'/data/')
save_fold<-paste0(dirname(rstudioapi::getSourceEditorContext()$path),'/save_data/')


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

all_protein_FG_DDA<-read.table(paste0(data_fold,'HYtims134_LFQ_FragPipe_all_proteins.tsv'), sep = '\t', header = T)
all_protein_MQ_DDA<-read.table(paste0(data_fold,'HYtims134_LFQ_Maxquant_all_proteins.tsv'), sep = '\t', header = T)
all_protein_DIANN_DIA<-read.table(paste0(data_fold,'HYtims134_DIA_DIANN_all_proteins.tsv'), sep = '\t', header = T)
all_protein_spt_DIA<-read.table(paste0(data_fold,'HYtims134_DIA_spt_all_proteins.tsv'), sep = '\t', header = T)

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
#saveWorkbook(wb, file = paste0(save_fold, 'cross_platform_compare.xlsx'), overwrite = TRUE)



library(pROC)
library(ggpubr)

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
saveWorkbook(wb, file = paste0(save_fold, 'cross_platform_compare_supp.xlsx'), overwrite = TRUE)
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


ggbarplot(nums, "method", "value", group = "variable", color = "variable",
          fill="variable", #label = T,lab.pos='in',
          platette = c('#5494cc','#e18283','#0d898a','#f9cc52'),
) + geom_text(aes(y=l_y, label=label))+
  scale_x_discrete(limits=c('FG_DDA_ens','MQ_DDA_ens','DIANN_DIA_ens', 'spt_DIA_ens',
                            'FG_DDA','MQ_DDA','DIANN_DIA','spt_DIA')) +
  geom_vline(xintercept = 4.5, linetype=2, show.legend=FALSE) + theme(axis.line.y = element_blank())+
  coord_flip()+theme(axis.title.y = element_blank())


sheet = 'supp. Fig.62'
addWorksheet(wb_all, sheet)
writeData(wb_all, sheet, nums, startRow = 2, startCol = 2)
###################################################
################# 
## supp.fig.7 The Receiver Operating Characteristic (ROC) curves showing the 
## performances of the best single workflows and best ens_multi-quant (for MQ_TMT, the ens_topk was tested instead) workflows under different settings.
## supp.fig.8 pooling DEPs from different quantification platforms

data_fold<-paste0(dirname(rstudioapi::getSourceEditorContext()$path),'/data/')
save_fold<-paste0(dirname(rstudioapi::getSourceEditorContext()$path),'/save_data/')

top1wf<-function(dt, platform, acq){
  base_fold<-data_fold
  dea_res_fold<-paste0(base_fold, acq, '/',platform, '/',dt,'/')
  ranking_all_FG_DDA<-read_excel(paste0(data_fold, 'ranks_all_', platform, '_', acq, '.xlsx'), sheet = 'ranking_all')
  
  
  wftop1_FG_DDA<-ranking_all_FG_DDA[which(ranking_all_FG_DDA$avg_rank_mean==min(ranking_all_FG_DDA$avg_rank_mean)),]
  if(length((wftop1_FG_DDA$workflow))>1){
    idx=which(wftop1_FG_DDA$avg_rank_median==min(wftop1_FG_DDA$avg_rank_median))
    wftop1_FG_DDA = wftop1_FG_DDA[idx,]
  }
  print(paste0(dea_res_fold,
               dt,'_',wftop1_FG_DDA$DEA,'_',
               wftop1_FG_DDA$Platform,'_',
               wftop1_FG_DDA$Matrix,'_',gsub('blank','',wftop1_FG_DDA$Imput),
               '_',gsub('blank','',wftop1_FG_DDA$normalization),'.csv'))
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
  
  ranking_all_FG_DDA<-read_excel(paste0(data_fold, 'ranks_all_', platform, '_', acq, '_', ens_ty, '.xlsx'), sheet = 'ranking_all')
  
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


FG_DDA<-top1_FG_DDA$dea_res$protein[which(abs(top1_FG_DDA$dea_res$logFC)>=log2(1.5) & top1_FG_DDA$dea_res$adj.pvalue<=0.05)]
MQ_DDA<-top1_MQ_DDA$dea_res$protein[which(abs(top1_MQ_DDA$dea_res$logFC)>=log2(1.5) & top1_MQ_DDA$dea_res$adj.pvalue<=0.05)]
DIANN_DIA<-top1_DIANN_DIA$dea_res$protein[which(abs(top1_DIANN_DIA$dea_res$logFC)>=log2(1.5) & top1_DIANN_DIA$dea_res$adj.pvalue<=0.05)]
spt_DIA<-top1_spt_DIA$dea_res$protein[which(abs(top1_spt_DIA$dea_res$logFC)>=log2(1.5) & top1_spt_DIA$dea_res$adj.pvalue<=0.05)]

FG_DDA_ens<-top1_FG_DDA_ens$dea_res$protein[which(abs(top1_FG_DDA_ens$dea_res$logFC)>=log2(1.5) & top1_FG_DDA_ens$dea_res$adj.pvalue<=0.05)]
MQ_DDA_ens<-top1_MQ_DDA_ens$dea_res$protein[which(abs(top1_MQ_DDA_ens$dea_res$logFC)>=log2(1.5) & top1_MQ_DDA_ens$dea_res$adj.pvalue<=0.05)]
DIANN_DIA_ens<-top1_DIANN_DIA_ens$dea_res$protein[which(abs(top1_DIANN_DIA_ens$dea_res$logFC)>=log2(1.5) & top1_DIANN_DIA_ens$dea_res$adj.pvalue<=0.05)]
spt_DIA_ens<-top1_spt_DIA_ens$dea_res$protein[which(abs(top1_spt_DIA_ens$dea_res$logFC)>=log2(1.5) & top1_spt_DIA_ens$dea_res$adj.pvalue<=0.05)]

pool_FG_MQ_DDA<-unique(c(FG_DDA, MQ_DDA))
pool_DIANN_spt_DIA<-unique(c(DIANN_DIA, spt_DIA))
pool_DDA_DIA<-unique(c(pool_FG_MQ_DDA, pool_DIANN_spt_DIA))

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

DEPs$pool_DDA=c(rep(0, length(DEPs$Proteins)))
DEPs$pool_DDA[match(pool_FG_MQ_DDA, DEPs$Proteins)]=1
DEPs$pool_DIA=c(rep(0, length(DEPs$Proteins)))
DEPs$pool_DIA[match(pool_DIANN_spt_DIA, DEPs$Proteins)]=1
DEPs$pool_DDA_DIA=c(rep(0, length(DEPs$Proteins)))
DEPs$pool_DDA_DIA[match(pool_DDA_DIA, DEPs$Proteins)]=1
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

library(pROC)
library(ggpubr)

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

cross_platform_proteins$qval_pool_DDA<-apply(cbind(cross_platform_proteins$qval_FG_DDA,
                                                   cross_platform_proteins$qval_MQ_DDA),
                                             1, function(x) min(x))
cross_platform_proteins$logFC_pool_DDA<-unlist(apply(cbind(cross_platform_proteins$logFC_FG_DDA,
                                                           cross_platform_proteins$logFC_MQ_DDA),
                                                     1, function(x) x[which(abs(x)==max(abs(x)))[1]]))

cross_platform_proteins$qval_pool_DIA<-apply(cbind(cross_platform_proteins$qval_DIANN_DIA,
                                                   cross_platform_proteins$qval_spt_DIA),
                                             1, function(x) min(x))
cross_platform_proteins$logFC_pool_DIA<-unlist(apply(cbind(cross_platform_proteins$logFC_DIANN_DIA,
                                                           cross_platform_proteins$logFC_spt_DIA),
                                                     1, function(x) x[which(abs(x)==max(abs(x)))[1]]))

cross_platform_proteins$qval_pool_DDA_DIA<-apply(cbind(cross_platform_proteins$qval_pool_DDA,
                                                       cross_platform_proteins$qval_pool_DIA),
                                                 1, function(x) min(x))
cross_platform_proteins$logFC_pool_DDA_DIA<-unlist(apply(cbind(cross_platform_proteins$logFC_pool_DDA,
                                                               cross_platform_proteins$logFC_pool_DIA),
                                                         1, function(x) x[which(abs(x)==max(abs(x)))[1]]))
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

roc_pool_DDA<-roc(cross_platform_proteins$labels, 1-cross_platform_proteins$qval_pool_DDA)
roc_pool_DIA<-roc(cross_platform_proteins$labels, 1-cross_platform_proteins$qval_pool_DIA)
roc_pool_DDA_DIA<-roc(cross_platform_proteins$labels, 1-cross_platform_proteins$qval_pool_DDA_DIA)

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
#############
## supp.fig.7 (a)
breaks = seq(0,1,0.1)
ggplot(spe_sens, aes(x = Specificity, y = Sensitivity, colour = `Method(pAUC(FPR=0.05))`))+
  scale_color_manual(values=c("#15d08a","#f2ca01","#a45ee3","#1f11b4","#c3642c",
                              "#f078b4", '#5494cc',"gray","pink",'red','orange')) +
  
  geom_step(linewidth=1,alpha = 0.9) +
  scale_x_continuous(name = "False Positive Rate (1 - Specificity)",limits = c(0,0.1), breaks = breaks) +
  scale_y_continuous(name = "True Positive Rate (Sensitivity)", limits = c(0,0.7), breaks = breaks) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 1))+scale_x_continuous(limits = c(0, 1))+
  labs(
    x = "False Positive Rate (1-Specificity)",
    y = "True Positive Rate (Sensitivity)")+ theme(legend.position = c(0.7, 0.35)) +
  geom_vline(xintercept = 0.05,linetype="dashed")

sheet = 'supp. Fig.71'
addWorksheet(wb_all, sheet)
writeData(wb_all, sheet, spe_sens, startRow = 2, startCol = 2)

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
  }else if(setting == 'pool_DDA'){
    logFCs=cross_platform_proteins$logFC_pool_DDA
    qvals=cross_platform_proteins$qval_pool_DDA
    roc_s<-roc_pool_DDA
  }else if(setting == 'pool_DIA'){
    logFCs=cross_platform_proteins$logFC_pool_DIA
    qvals=cross_platform_proteins$qval_pool_DIA
    roc_s<-roc_pool_DIA
  }else if(setting == 'pool_DDA_DIA'){
    logFCs=cross_platform_proteins$logFC_pool_DDA_DIA
    qvals=cross_platform_proteins$qval_pool_DDA_DIA
    roc_s<-roc_pool_DDA_DIA
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

metrics_pool_DDA=cross_platform_metrics(cross_platform_proteins, 'pool_DDA', roc_pool_DDA)
metrics_pool_DIA=cross_platform_metrics(cross_platform_proteins, 'pool_DIA', roc_pool_DIA)
metrics_pool_DDA_DIA=cross_platform_metrics(cross_platform_proteins, 'pool_DDA_DIA', roc_pool_DDA_DIA)

rdar_dt<-rbind(metrics_FG_DDA$value, metrics_MQ_DDA$value, metrics_DIANN_DIA$value, metrics_spt_DIA$value,
               metrics_FG_DDA_ens$value, metrics_MQ_DDA_ens$value, 
               metrics_DIANN_DIA_ens$value, metrics_spt_DIA_ens$value,
metrics_pool_DDA$value, metrics_pool_DIA$value, metrics_pool_DDA_DIA$value)
colnames(rdar_dt)<-metrics_FG_DDA$metric
row.names(rdar_dt)<-c('FG_DDA','MQ_DDA', 'DIANN_DIA', 'spt_DIA', 'FG_DDA_ens',
                      'MQ_DDA_ens', 'DIANN_DIA_ens', 'spt_DIA_ens',
'pool_DDA','pool_DIA', 'pool_DDA_DIA')
rdar_dt<-as.data.frame(as.matrix(rdar_dt))
pt_rdar_dt<-data.frame(group=row.names(rdar_dt))

sheet=paste(dt, '_cp_metrics')

addWorksheet(wb,sheet)
writeData(wb, sheet, rdar_dt, startRow = 1, startCol = 1, colNames = TRUE, rowNames = TRUE)
saveWorkbook(wb, file = paste0(save_fold, 'cross_platform_compare.xlsx'), overwrite = TRUE)
pt_rdar_dt<-cbind(pt_rdar_dt, round(rdar_dt[,c(1,2,3,13,14)], 2)) #,4,6

library(ggradar)
ggradar(pt_rdar_dt,group.colours=c("#a45ee3","#1f11b4","#f078b4","#f2ca01","gray","pink",'red',"orange"),legend.position = 'right')

col=c("#15d08a","#f2ca01","#a45ee3","#1f11b4","#c3642c",
      "#f078b4", '#5494cc',"gray","pink",'red','orange')
melt_dt<-melt(pt_rdar_dt)
melt_dt$value<-round(melt_dt$value,2)
colnames(melt_dt)[1]<-'method'


ggdotchart(melt_dt, "variable", "value", group = "method", color = "method",
           palette = col[c(6,1,8,3,10,7,5,4, 2, 9, 11)],
           
           add = "segments", #label = "value",
           add.params = list(color = "lightgray", size = 1),
           dot.size = 3,
           position=position_dodge(width = 0.9)
) +
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

#############
## supp.fig.8
ggbarplot(nums, "method", "value", group = "variable", color = "variable",
          fill="variable", #label = T,lab.pos='in',
          platette = c('#5494cc','#e18283','#0d898a','#f9cc52'),
) + geom_text(aes(y=l_y, label=label))+
  scale_x_discrete(limits=c('FG_DDA_ens','MQ_DDA_ens','DIANN_DIA_ens', 'spt_DIA_ens',
                            'FG_DDA','MQ_DDA','DIANN_DIA','spt_DIA', 'pool_DDA', 'pool_DIA', 'pool_DDA_DIA')) +
  geom_vline(xintercept = 4.5, linetype=2, show.legend=FALSE) + theme(axis.line.y = element_blank())+
  coord_flip()+theme(axis.title.y = element_blank())


sheet = 'supp. Fig.8'
addWorksheet(wb_all, sheet)
writeData(wb_all, sheet, nums, startRow = 2, startCol = 2)
#############
## supp.fig.7 (b)

all_protein_FG_DDA<-read.table(paste0(data_fold,'HYqfl683_LFQ_FragPipe_all_proteins.tsv'), sep = '\t', header = T)
all_protein_MQ_DDA<-read.table(paste0(data_fold,'HYqfl683_LFQ_Maxquant_all_proteins.tsv'), sep = '\t', header = T)
all_protein_FG_TMT<-read.table(paste0(data_fold,'HYqfl683_TMT11_FragPipe_tmt_all_proteins.tsv'), sep = '\t', header = T)
all_protein_MQ_TMT<-read.table(paste0(data_fold,'Maxquant/HYqfl683_TMT11_Maxquant_tmt_all_proteins.tsv'), sep = '\t', header = T)



data_fold<-'E:/proteomics/maus2/submission_NC/revision2/Figures/data/combined_res/'
save_fold<-'E:/proteomics/maus2/submission_NC/revision2/Figures/Figure5/'



top1wf<-function(dt, platform, acq, contrast){
  base_fold<-data_fold
  dea_res_fold<-paste0(base_fold, acq, '/',platform, '/',dt,'/')
  ranking_all_FG_DDA<-read_excel(paste0(data_fold, 'ranks_all_', platform, '_', acq, '.xlsx'), sheet = 'ranking_all')
  
  wftop1_FG_DDA<-ranking_all_FG_DDA[which(ranking_all_FG_DDA$avg_rank_mean==min(ranking_all_FG_DDA$avg_rank_mean)),]
  if(length((wftop1_FG_DDA$workflow))>1){
    idx=which(wftop1_FG_DDA$avg_rank_median==min(wftop1_FG_DDA$avg_rank_median))
    wftop1_FG_DDA = wftop1_FG_DDA[idx,]
  }
  print(paste0(dea_res_fold,
               dt,'_',wftop1_FG_DDA$DEA,'_',
               wftop1_FG_DDA$Platform,'_',
               wftop1_FG_DDA$Matrix,'_',gsub('blank','',wftop1_FG_DDA$Imput),
               '_',gsub('blank','',wftop1_FG_DDA$normalization),'.csv'))
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
  
  ranking_all_FG_DDA<-read_excel(paste0(data_fold, 'ranks_all_', platform, '_', acq, '_', ens_ty, '.xlsx'), sheet = 'ranking_all')
  
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


library(pROC)
library(ggpubr)

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

## supp.fig.7(b)

breaks = seq(0,1,0.1)
ggplot(spe_sens, aes(x = Specificity, y = Sensitivity, colour = `Method(pAUC(FPR=0.05))`))+
  scale_color_manual(values=c("#15d08a","#f2ca01","#a45ee3","#1f11b4","#c3642c",
                              "#f078b4", '#5494cc',"gray","pink",'red','orange')) +
  geom_step(linewidth=1,alpha = 0.9) +
  scale_x_continuous(name = "False Positive Rate (1 - Specificity)",limits = c(0,0.1), breaks = breaks) +
  scale_y_continuous(name = "True Positive Rate (Sensitivity)", limits = c(0,0.7), breaks = breaks) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 1))+scale_x_continuous(limits = c(0, 1))+
  labs(
    x = "False Positive Rate (1-Specificity)",
    y = "True Positive Rate (Sensitivity)")+ theme(legend.position = c(0.7, 0.25)) +
  geom_vline(xintercept = 0.05,linetype="dashed")


sheet = 'supp. Fig.72'
addWorksheet(wb_all, sheet)
writeData(wb_all, sheet, spe_sens, startRow = 2, startCol = 2)

saveWorkbook(wb_all, file = paste0(save_fold,'supp.Figs_source_data.xlsx'), overwrite = TRUE)
