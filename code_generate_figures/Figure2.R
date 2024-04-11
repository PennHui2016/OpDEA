####################################
#
# Figure2 for benchmarking results, generalizability, predicatable, interaction, rule, top ranking
#
#####################

### Figure 2A. performance distribution (FragPipe - DDA)

library("readxl")
library(ggplot2)
library(reshape2)
library(ggpubr)
library(cowplot)
data_fold<-paste0(dirname(rstudioapi::getSourceEditorContext()$path),'/data/')
save_fold<-paste0(dirname(rstudioapi::getSourceEditorContext()$path),'/save_data/')
platform='FragPipe'
acq='DDA'
setting = 'FG_DDA'
ranking_all<-read_excel(paste0(data_fold, 'ranks_all_', platform, '_', acq, '.xlsx'), sheet = 'ranking_all')
ranking_all<-ranking_all[order(ranking_all$avg_rank_mean),]
ranking_all$label<-c(rep('top5%',round(length(ranking_all$workflow))))
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

col=c("#1f78b4","#b2df8a","#a6cee3","#33a02c","#fb9a99","#fdbf6f")
p1 = ggplot(dat, aes(x=metrics, y=performance, fill=metrics)) +
  geom_violin(position = position_dodge(width = 0.75), scale = 'width')+
  scale_fill_manual(values = col)+
  labs(x="", y="performance", color=factor)+
  scale_x_discrete(limits=c("pAUC(0.01)","pAUC(0.05)","pAUC(0.1)",'nMCC', 'G-mean'))
mytheme = theme_classic() + theme(axis.text.x = element_text(size = 16, angle = -20),axis.text.y = element_text(size = 16))+
  theme(axis.title.y= element_text(size=16))+theme(axis.title.x = element_text(size = 16))+
  theme(legend.title=element_text(size=16),legend.text=element_text(size=16), legend.position = "none")
p1=p1+mytheme
p1

write.table(dat_mean, paste0(save_fold, 'Figure2A.csv'), col.names = T, row.names = F)

### Figure 2B. example LODOCV with HYEtims735_LFQ dataset
### scatter + regression showing comparsion of ranking with single dataset and cross-dataset

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
  write.table(sprs, paste0(save_fold, platform, '-', acq, '_spr_avg_rank_test.csv'), sep = ',', col.names = TRUE, row.names = FALSE)
  return(list(tr_rk_mean=tr_ranks_mean, tr_rk_median=tr_ranks_median, te_ranks=te_ranks, sprs=sprs))
}

lopocv_fg_dda<-spr_test('Fragpipe', 'DDA')
lopocv_mq_dda<-spr_test('Maxquant', 'DDA')
lopocv_fg_tmt<-spr_test('Fragpipe', 'TMT')
lopocv_mq_tmt<-spr_test('Maxquant', 'TMT')
lopocv_diann_dia<-spr_test('DIANN', 'DIA')
lopocv_spt_dia<-spr_test('spt', 'DIA')

tr_rank_fg_dda_mean<-as.data.frame(lopocv_fg_dda$tr_rk_mean)
te_rank_fg_dda<-as.data.frame(lopocv_fg_dda$te_ranks)

scatter_dat<-cbind(te_rank_fg_dda$`HYEtims735_LFQ_conditionB-conditionA`,
                   tr_rank_fg_dda_mean$HYEtims735,
                   rep('mean',length(te_rank_fg_dda$`HYEtims735_LFQ_conditionB-conditionA`)))

scatter_dat<-data.frame(new=as.numeric(scatter_dat[,1]),
                        bechmarking=as.numeric(scatter_dat[,2]),
                        method=scatter_dat[,3])

colnames(scatter_dat)<-c('new','bechmarking','method')


library(ggplot2)
library(ggpubr)

p2<-ggplot(scatter_dat, aes(new, bechmarking, color = method)) +
  geom_point(alpha=0.1) +
  geom_smooth(method = "lm") +
  stat_cor(method = "spearman",aes(color = method, label = after_stat(r.label)),label.x = c(5000,5000), label.y = c(800,2000)) +
  scale_color_manual(values = c("black"), aesthetics = c("color", "fill")) +
  theme_pubr() +theme(legend.title=element_text(size=16),legend.text=element_text(size=16), legend.position = "none")+
  labs(x = "ranking by performance of new dataset", y = "ranking by benchmarking")
mytheme = theme_classic() + theme(axis.text.x = element_text(size = 14, angle = 0),axis.text.y = element_text(size = 14))+
  theme(axis.title.y= element_text(size=15))+theme(axis.title.x = element_text(size = 15))+
  theme(legend.title=element_text(size=16),legend.text=element_text(size=16), legend.position = "none")
p2<-p2+mytheme
p2

write.table(scatter_dat, paste0(save_fold, 'Figure2B.csv'), col.names = T, row.names = F)

######## C. distributions of LOPOCV sprs
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

## only show mean

sprs_dat_mean<-cbind(all_sprs$mean_spr,all_sprs$setting,
                     c(rep('mean',length(all_sprs$dt_name))))


colnames(sprs_dat_mean)<-c('spearman','setting','method')
sprs_dat_mean<-as.data.frame(sprs_dat_mean)
sprs_dat_mean$spearman<-as.numeric(sprs_dat_mean$spearman)

col=c('#425AB2','#34df8a','#F45896','#c21000','#452904','#e2df8a')
p31 = ggplot(sprs_dat_mean, aes(x=setting, y=spearman, color=setting)) +
  geom_boxplot(position = position_dodge(width = 0.3), width=0.3,)+ geom_jitter(width = 0.2)+
  scale_fill_manual(values = col)+
  labs(x="", y="spearman correlation", color=factor)+
  scale_x_discrete(limits=unique(all_sprs$setting))+stat_summary(fun=mean,
                                                                 geom="point",
                                                                 shape=17,
                                                                 aes(group=method),
                                                                 position=position_dodge(.25), size=2, color="red", fill="red")#+


mytheme = theme_classic() + theme(axis.text.x = element_text(size = 16, angle = -20),axis.text.y = element_text(size = 10))+
  theme(axis.title.y= element_text(size=16))+theme(axis.title.x = element_text(size = 14))+
  theme(legend.title=element_text(size=16),legend.text=element_text(size=16), legend.position = "none")

p31=p31+mytheme
p31

write.table(sprs_dat_mean, paste0(save_fold, 'Figure2C.csv'), col.names = T, row.names = F)

############## D. ANOVA analysis of dataset-workflow ranking dependency

library("readxl")
library(ggplot2)
library(reshape2)
library(ggpubr)
library(cowplot)

data_fold<-paste0(dirname(rstudioapi::getSourceEditorContext()$path),'/data/')
save_fold<-paste0(dirname(rstudioapi::getSourceEditorContext()$path),'/save_data/')
single_rank<-function(data_single){
  colnames(data_single)<-'V1'
  single_rank<-as.data.frame(cbind(data_single,c(1:length(as.list(data_single$V1)))))
  colnames(single_rank)<-c('value','idx')
  single_rank<-single_rank[order(-single_rank$value),]
  single_rank$rank<-c(1:length(data_single$V1))
  single_rank<-single_rank[order(single_rank$idx),]
  return(single_rank)
}

avg_rank<-function(platform, acq, st, en){
  pauc001<-read_excel(paste0(data_fold, 'metrics_',platform, '_',acq, '.xlsx'), sheet='pAUC0.01')
  pauc005<-read_excel(paste0(data_fold, 'metrics_',platform, '_',acq, '.xlsx'), sheet='pAUC0.05')
  pauc01<-read_excel(paste0(data_fold, 'metrics_',platform, '_',acq, '.xlsx'), sheet='pAUC0.1')
  nMCC<-read_excel(paste0(data_fold, 'metrics_',platform, '_',acq, '.xlsx'), sheet='nMCC0.05')
  Gmean<-read_excel(paste0(data_fold, 'metrics_',platform, '_',acq, '.xlsx'), sheet='G-mean0.05')
  
  avg_ranks<-vector()
  for (i in c(st:en)) {
    rank_pauc001<-single_rank(pauc001[,i])
    rank_pauc005<-single_rank(pauc005[,i])
    rank_pauc01<-single_rank(pauc01[,i])
    rank_nMCC<-single_rank(nMCC[,i])
    rank_Gmean<-single_rank(Gmean[,i])
    avg_rank<-(rank_pauc001$rank+rank_pauc005$rank+rank_pauc01$rank+rank_nMCC$rank+rank_Gmean$rank)/5
    avg_final<-as.data.frame(cbind(avg_rank,c(1:length(avg_rank))))
    colnames(avg_final)<-c('avg_rank_c','idx')
    avg_final<-avg_final[order(avg_final$avg_rank_c),]
    avg_final$final_rank<-c(1:length(avg_final$avg_rank_c))
    avg_final<-avg_final[order(avg_final$idx),]
    avg_ranks<-cbind(avg_ranks, avg_final$final_rank)
  }
  colnames(avg_ranks)<-colnames(pauc001)[st:en]
  row.names(avg_ranks)<-pauc001$workflow
  return(avg_ranks)
}

ANOVA_lm<-function(platform, acq, st, en, wb){
  
  ranks_pauc001<-avg_rank(platform, acq, st, en)
  write.table(ranks_pauc001, paste0(save_fold,
                                    platform, '_', acq, '_', '_all_ranks.csv'),
              sep = ',', col.names = T)
  header<-read_excel(paste0(data_fold, 'metrics_',platform, '_',acq, '.xlsx'), sheet='header', col_names = FALSE)
  pauc001<-read_excel(paste0(data_fold, 'metrics_',platform, '_',acq, '.xlsx'), sheet='pAUC0.01')
  
  p_kw<-vector()
  sta_kw<-vector()
  
  
  for (i in 1:length(ranks_pauc001[,1])) {
    melt_data<-as.data.frame(cbind(colnames(ranks_pauc001), as.numeric(ranks_pauc001[i,]), as.character(header[3,][st:en])))
    
    colnames(melt_data)<-c('dataset','performance','instrument')
    remain<-c()
    for (i in 1:length(unique(melt_data$instrument))) {
      
      idx<-which(melt_data$instrument==unique(melt_data$instrument)[i])
      if(length(idx)>=5){
        remain<-c(remain, idx)
      }
    }
    melt_data<-melt_data[remain,]
    res_kw<-kruskal.test(performance ~ instrument, data = melt_data)
    
    p_kw<-c(p_kw, res_kw$p.value)
    sta_kw<-c(sta_kw, res_kw$statistic)
    
  }
  
  all_tk<-as.data.frame(cbind(pauc001$workflow, sta_kw, p_kw))
  colnames(all_tk)<-c('workflow', 'statistic', 'pvalue')
  sheet=paste0(platform, '_', acq)
  addWorksheet(wb,sheet)
  writeData(wb, sheet, all_tk, startRow = 1, startCol = 1, colNames = TRUE, rowNames = TRUE)
  saveWorkbook(wb, file = paste0(save_fold, 'cross_instrument_avg_rank_Kruskal-Wallis_test.xlsx'), overwrite = TRUE)
  return(list(anova=all_tk))
}


library(openxlsx)
wb1 <- createWorkbook()

fg_avg_anova_DDA<-ANOVA_lm('FragPipe', 'DDA', 7, 28, wb1)
mq_avg_anova_DDA<-ANOVA_lm('Maxquant', 'DDA', 7, 28, wb1)
diann_avg_anova_DIA<-ANOVA_lm('DIANN', 'DIA', 7, 23, wb1)
spt_avg_anova_DIA<-ANOVA_lm('spt', 'DIA', 7, 23, wb1)

rank_fg_DDA<-read_excel(paste0(data_fold,'ranks_all_', 'FragPipe', '_', 'DDA', '.xlsx'), sheet = 'ranking_all')
rank_mq_DDA<-read_excel(paste0(data_fold,'ranks_all_', 'Maxquant', '_', 'DDA', '.xlsx'), sheet = 'ranking_all')
rank_fg_TMT<-read_excel(paste0(data_fold,'ranks_all_', 'FragPipe', '_', 'TMT', '.xlsx'), sheet = 'ranking_all')
rank_mq_TMT<-read_excel(paste0(data_fold,'ranks_all_', 'Maxquant', '_', 'TMT', '.xlsx'), sheet = 'ranking_all')
rank_diann_DIA<-read_excel(paste0(data_fold,'ranks_all_', 'DIANN', '_', 'DIA', '.xlsx'), sheet = 'ranking_all')
rank_spt_DIA<-read_excel(paste0(data_fold,'ranks_all_', 'spt', '_', 'DIA', '.xlsx'), sheet = 'ranking_all')

topwf_anova_p<-function(fg_avg_anova_DDA, rank_fg_DDA, setting){
  
  fg_DDA_ps<-data.frame(workflow=fg_avg_anova_DDA$anova$workflow,p=fg_avg_anova_DDA$anova$pvalue)
  
  
  idx_wf<-match(fg_DDA_ps$workflow, rank_fg_DDA$workflow)
  fg_DDA_ps$rank_mean<-rank_fg_DDA$avg_rank_mean[idx_wf]
  fg_DDA_ps$rank_median<-rank_fg_DDA$avg_rank_median[idx_wf]
  fg_DDA_ps$setting<-setting
  fg_DDA_ps$idx<-c(1:length(fg_DDA_ps$workflow))
  fg_DDA_ps<-fg_DDA_ps[order(fg_DDA_ps$rank_mean),]
  fg_DDA_ps$rank_mean_exact<-c(1:length(fg_DDA_ps$workflow))
  fg_DDA_ps<-fg_DDA_ps[order(fg_DDA_ps$rank_median),]
  fg_DDA_ps$rank_median_exact<-c(1:length(fg_DDA_ps$workflow))
  fg_DDA_ps<-fg_DDA_ps[order(fg_DDA_ps$idx),]
  return(fg_DDA_ps)
}


fg_DDA_ps<-topwf_anova_p(fg_avg_anova_DDA, rank_fg_DDA, "FG_DDA")
mq_DDA_ps<-topwf_anova_p(mq_avg_anova_DDA, rank_mq_DDA, "MQ_DDA")
diann_DIA_ps<-topwf_anova_p(diann_avg_anova_DIA, rank_diann_DIA, "DIANN_DIA")
spt_DIA_ps<-topwf_anova_p(spt_avg_anova_DIA, rank_spt_DIA, "spt_DIA")

top=30
idx_fg_dda<-which(fg_DDA_ps$rank_mean_exact<=top | fg_DDA_ps$rank_median_exact<=top)
idx_mq_dda<-which(mq_DDA_ps$rank_mean_exact<=top | mq_DDA_ps$rank_median_exact<=top)
idx_diann_dia<-which(diann_DIA_ps$rank_mean_exact<=top | diann_DIA_ps$rank_median_exact<=top)
idx_spt_dia<-which(spt_DIA_ps$rank_mean_exact<=top | spt_DIA_ps$rank_median_exact<=top)

all_anova_ps<-rbind(fg_DDA_ps[idx_fg_dda,],
                    mq_DDA_ps[idx_mq_dda,],
                    diann_DIA_ps[idx_diann_dia,],
                    spt_DIA_ps[idx_spt_dia,])


melt_anova_ps_dat<-cbind(all_anova_ps$workflow[which(all_anova_ps$rank_mean_exact<=top)],
                         all_anova_ps$p[which(all_anova_ps$rank_mean_exact<=top)],
                         all_anova_ps$rank_mean_exact[which(all_anova_ps$rank_mean_exact<=top)],
                         all_anova_ps$setting[which(all_anova_ps$rank_mean_exact<=top)],
                         c(rep('mean',top)))


melt_anova_ps_dat<-rbind(melt_anova_ps_dat,
                         cbind(all_anova_ps$workflow[which(all_anova_ps$rank_median_exact<=top)],
                               all_anova_ps$p[which(all_anova_ps$rank_median_exact<=top)],
                               all_anova_ps$rank_median_exact[which(all_anova_ps$rank_median_exact<=top)],
                               all_anova_ps$setting[which(all_anova_ps$rank_median_exact<=top)],
                               c(rep('median',top))))

melt_anova_ps_dat<-as.data.frame(melt_anova_ps_dat)
colnames(melt_anova_ps_dat)<-c('workflow','pvalue','rank','setting','method')
melt_anova_ps_dat$log10p<--log10(as.numeric(melt_anova_ps_dat$pvalue))
melt_anova_ps_dat$rank<-as.integer(melt_anova_ps_dat$rank)
melt_anova_ps_dat$stable<-as.factor(ifelse(as.numeric(melt_anova_ps_dat$pvalue) <= 0.05,
                                           ifelse(melt_anova_ps_dat$method=='mean',
                                                  'Non-STABLE & mean',
                                                  'Non-STABLE & median'),
                                           'STABLE'))
melt_anova_ps_dat$label<-paste0(melt_anova_ps_dat$setting,'|',melt_anova_ps_dat$method,'|',melt_anova_ps_dat$rank)
melt_anova_ps_dat$label[which(as.numeric(melt_anova_ps_dat$pvalue)>0.05)]=''

### only present mean
melt_anova_ps_dat<-cbind(all_anova_ps$workflow[which(all_anova_ps$rank_mean_exact<=top)],
                         all_anova_ps$p[which(all_anova_ps$rank_mean_exact<=top)],
                         all_anova_ps$rank_mean_exact[which(all_anova_ps$rank_mean_exact<=top)],
                         all_anova_ps$setting[which(all_anova_ps$rank_mean_exact<=top)],
                         c(rep('mean',top)))
melt_anova_ps_dat<-as.data.frame(melt_anova_ps_dat)
colnames(melt_anova_ps_dat)<-c('workflow', 'pvalue','rank','setting','method')
melt_anova_ps_dat$log10p<--log10(as.numeric(melt_anova_ps_dat$pvalue))
melt_anova_ps_dat$rank<-as.integer(melt_anova_ps_dat$rank)
melt_anova_ps_dat$stable<-as.factor(ifelse(as.numeric(melt_anova_ps_dat$pvalue) <= 0.05,
                                           'Non-STABLE',
                                           'STABLE'))
melt_anova_ps_dat$label<-paste0(melt_anova_ps_dat$setting,'|',melt_anova_ps_dat$rank)
melt_anova_ps_dat$label[which(as.numeric(melt_anova_ps_dat$pvalue)>0.05)]=''
library(ggrepel)
p4_1<-ggplot(melt_anova_ps_dat, aes(x=rank, y=log10p,color=stable, shape=setting)) +
  geom_point(alpha=1, size=3) +
  theme_bw(base_size = 12) +
  scale_shape_manual(values=c(1, 3,4, 8, 13, 24))+
  xlab("ranking of workflows") +
  ylab("-Log10(pvalue) of Kruskal-Wallis test") +
  scale_colour_manual(values = c("#ff6600",'gray')) +
  geom_hline(yintercept = -log10(0.05), lty = 4) +
  geom_label_repel(data = melt_anova_ps_dat, aes(label = label),
                   size = 4,box.padding = unit(0.5, "lines"),
                   point.padding = unit(0.8, "lines"),
                   segment.color = "black",force = 2,
                   show.legend = FALSE, max.overlaps = 10000)+
  theme(legend.title=element_text(size=13),legend.text=element_text(size=12), legend.position = "top")+
  guides(color=guide_legend(nrow=2), shape=guide_legend(nrow=2))

p4_1

write.table(melt_anova_ps_dat, paste0(save_fold, 'Figure2D.csv'), col.names = T, row.names = F)

######################## top ranking workflows
## Figure 2E
data_fold<-paste0(dirname(rstudioapi::getSourceEditorContext()$path),'/data/')
save_fold<-paste0(dirname(rstudioapi::getSourceEditorContext()$path),'/save_data/')
library(readxl)

rank_frag_DDA<-read_excel(paste0(data_fold,'ranks_all_', 'FragPipe', '_', 'DDA','.xlsx'), sheet = 'ranking_all')
rank_mq_DDA<-read_excel(paste0(data_fold,'ranks_all_', 'Maxquant', '_', 'DDA','.xlsx'), sheet = 'ranking_all')
rank_frag_TMT<-read_excel(paste0(data_fold,'ranks_all_', 'FragPipe', '_', 'TMT','.xlsx'), sheet = 'ranking_all')
rank_mq_TMT<-read_excel(paste0(data_fold,'ranks_all_', 'Maxquant', '_', 'TMT','.xlsx'), sheet = 'ranking_all')
rank_diann_DIA<-read_excel(paste0(data_fold,'ranks_all_', 'DIANN', '_', 'DIA','.xlsx'), sheet = 'ranking_all')
rank_spt_DIA<-read_excel(paste0(data_fold,'ranks_all_', 'spt', '_', 'DIA','.xlsx'), sheet = 'ranking_all')

rank_frag_DDA<-rank_frag_DDA[order(rank_frag_DDA$avg_rank_mean),]
rank_mq_DDA<-rank_mq_DDA[order(rank_mq_DDA$avg_rank_mean),]
rank_frag_TMT<-rank_frag_TMT[order(rank_frag_TMT$avg_rank_mean),]
rank_mq_TMT<-rank_mq_TMT[order(rank_mq_TMT$avg_rank_mean),]
rank_diann_DIA<-rank_diann_DIA[order(rank_diann_DIA$avg_rank_mean),]
rank_spt_DIA<-rank_spt_DIA[order(rank_spt_DIA$avg_rank_mean),]

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

platforms<-c('Maxquant', 'Fragpipe', 'Fragpipe', 'Maxquant', 'DIANN', 'spt')
acqs<-c('TMT', 'TMT', 'DDA', 'DDA', 'DIA', 'DIA')
setting<-c('MQ_TMT','FG_TMT','FG_DDA','MQ_DDA','DIANN_DIA','spt_DIA')

DEA_sym<-c(1:10)
norm_sym<-LETTERS
imp_sym<-letters
matrix_sym<-c('alpha', 'beta', 'gamma', 'delta', 'epsilon')

N=2
col=c("#1f78b4","#33a02c","#fb9a99","orange","black","purple","green","darkgray","#a6cee3","red")
p81<-plot_topwf(1, col[7])
p82<-plot_topwf(2, col[8:10])
p83<-plot_topwf(3, col[c(1:5)])
p84<-plot_topwf(4, col[c(1:5)])
p85<-plot_topwf(5, col[c(2,3,6,5)])
p86<-plot_topwf(6, col[c(2,3,6,5)])
# save and combine in PPT (500*350)

all_tops<-rbind(p81$dt,p82$dt,p83$dt,p84$dt,p85$dt,p86$dt)
write.table(all_tops, paste0(save_fold, 'all_top_workflows.csv'), sep = ',', col.names = T, row.names = F)

# save and combine in PPT (500*350)
write.table(p83$dt, paste0(save_fold, 'Figure2E_1.csv'), col.names = T, row.names = F)
write.table(p84$dt, paste0(save_fold, 'Figure2E_2.csv'), col.names = T, row.names = F)
write.table(p85$dt, paste0(save_fold, 'Figure2E_3.csv'), col.names = T, row.names = F)
write.table(p86$dt, paste0(save_fold, 'Figure2E_4.csv'), col.names = T, row.names = F)
