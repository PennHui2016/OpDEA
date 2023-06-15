library("readxl")
library(ComplexHeatmap)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(cowplot)

# DIANN
data_fold<-'F:/NTU/quantification/data/2022.10.18/source_data_for_figures/'
header<-read_excel(paste0(data_fold, 'metrics_diann.xlsx'), sheet='header', col_names = FALSE)
diann_pauc001<-read_excel(paste0(data_fold, 'metrics_diann.xlsx'), sheet='pAUC0.01')
datasets<-unique(t(header[1,][7:length(header[1,])]))

header_dataset<-paste0(header[1,],'_', header[2,])
datas<-diann_pauc001[,7:152]
colnames(datas)<-header_dataset[7:152]
library(corrplot)

corrmatrix = cor(datas, method = 'spearman')

rownames(corrmatrix) = paste(1:ncol(datas), names(datas),sep = ' ') # ?????????????????????????????? <?????? ?????????> ?????????
colnames(corrmatrix) = as.character(1:ncol(datas)) # ????????????????????????????????????
col3 <- colorRampPalette(c('DodgerBlue3','white', "OrangeRed")) # ??????????????????col3, ??????????????????????????????colorbar?????????
corrplot(corrmatrix, method = 'circle', diag = F, type = 'full', outline = F,
         col = col3(10), cl.lim = c(-0.2,1),addgrid.col = NA, 
         tl.pos = 'lb',tl.cex = 0.5, tl.col = 'black', tl.srt = 0, tl.offset = 0.5)
axis(1,at = 1:12, labels = NA, pos = 12.5, tck = -0.01) # ??????????????????????????????,????????????????????????,??????,????????????,???????????????
axis(4,at = 1:12, labels = NA, pos = 0.5, tck = -0.01) # ??????????????????????????????,??????????????????

#########################################
#checking average performance for prediction

metrics<-c('pAUC0.01', 'pAUC0.05', 'pAUC0.1', 'G-mean', 'nMCC')
data_fold<-'F:/NTU/quantification/data/2022.10.18/source_data_for_figures/'
header<-read_excel(paste0(data_fold, 'metrics_diann.xlsx'), sheet='header', col_names = FALSE)
header_dataset<-paste0(header[1,],'_', header[2,])
dts<-header[1,][7:152]
dt_n<-c()

for(i in dts[1,]){
  s = strsplit(i, '_', fixed = TRUE)[[1]]
  dt_n<-c(dt_n, s[1])
}

sprs<-vector()
for (metric in metrics) {
  #metric = 'pAUC0.01'
  diann_pauc001<-read_excel(paste0(data_fold, 'metrics_diann.xlsx'), sheet=metric)
  datas<-diann_pauc001[,7:152]
  colnames(datas)<-header_dataset[7:152]
  
  for (dt in unique(dt_n)) {
    mean_tr<-rowMeans(datas[,which(dt_n!=dt)])
    median_tr<-apply(datas[,which(dt_n!=dt)], 1, median)
    
    tests<-datas[,which(dt_n==dt)]
    for (i in 1:length(tests[1,])) {
      spr1 = cor(mean_tr, tests[, i], method = 'spearman')
      spr2 = cor(median_tr, tests[, i], method = 'spearman')
      sprs<-rbind(sprs, c(dt, metric, i, spr1, spr2))
    }
  }
}
colnames(sprs)<-c('dt_name', 'metric', 'dt_col', 'mean_spr', 'median_spr')
write.table(sprs, paste0(data_fold, 'DIANN_spr_test.csv'), sep = ',', col.names = TRUE, row.names = FALSE)

# fragpipe
header<-read_excel(paste0(data_fold, 'metrics_fragpipe.xlsx'), sheet='header', col_names = FALSE)
diann_pauc001<-read_excel(paste0(data_fold, 'metrics_fragpipe.xlsx'), sheet='G-mean')
datasets<-unique(t(header[1,][7:length(header[1,])]))

header_dataset<-paste0(header[1,],'_', header[2,])
datas<-diann_pauc001[,7:126]
colnames(datas)<-header_dataset[7:126]
library(corrplot)

corrmatrix = cor(datas, method = 'spearman')

rownames(corrmatrix) = paste(1:ncol(datas), names(datas),sep = ' ') # ?????????????????????????????? <?????? ?????????> ?????????
colnames(corrmatrix) = as.character(1:ncol(datas)) # ????????????????????????????????????
col3 <- colorRampPalette(c('DodgerBlue3','white', "OrangeRed")) # ??????????????????col3, ??????????????????????????????colorbar?????????
corrplot(corrmatrix, method = 'circle', diag = F, type = 'full', outline = F,
         col = col3(10), cl.lim = c(-0.2,1),addgrid.col = NA, 
         tl.pos = 'lb',tl.cex = 0.5, tl.col = 'black', tl.srt = 0, tl.offset = 0.5)
axis(1,at = 1:12, labels = NA, pos = 12.5, tck = -0.01) # ??????????????????????????????,????????????????????????,??????,????????????,???????????????
axis(4,at = 1:12, labels = NA, pos = 0.5, tck = -0.01) # ??????????????????????????????,??????????????????

#########################################
#checking average performance for prediction

metrics<-c('pAUC0.01', 'pAUC0.05', 'pAUC0.1', 'G-mean', 'nMCC')
data_fold<-'F:/NTU/quantification/data/2022.10.18/source_data_for_figures/'
header<-read_excel(paste0(data_fold, 'metrics_fragpipe.xlsx'), sheet='header', col_names = FALSE)
header_dataset<-paste0(header[1,],'_', header[2,])
dts<-header[1,][7:126]
dt_n<-c()

for(i in dts[1,]){
  s = strsplit(i, '_', fixed = TRUE)[[1]]
  dt_n<-c(dt_n, s[1])
}

sprs<-vector()
for (metric in metrics) {
  #metric = 'pAUC0.01'
  diann_pauc001<-read_excel(paste0(data_fold, 'metrics_fragpipe.xlsx'), sheet=metric)
  datas<-diann_pauc001[,7:126]
  colnames(datas)<-header_dataset[7:126]
  
  for (dt in unique(dt_n)) {
    mean_tr<-rowMeans(datas[,which(dt_n!=dt)])
    median_tr<-apply(datas[,which(dt_n!=dt)], 1, median)
    
    tests<-datas[,which(dt_n==dt)]
    for (i in 1:length(tests[1,])) {
      spr1 = cor(mean_tr, tests[, i], method = 'spearman')
      spr2 = cor(median_tr, tests[, i], method = 'spearman')
      sprs<-rbind(sprs, c(dt, metric, i, spr1, spr2))
    }
  }
}
colnames(sprs)<-c('dt_name', 'metric', 'dt_col', 'mean_spr', 'median_spr')
write.table(sprs, paste0(data_fold, 'fragpipe_spr_test.csv'), sep = ',', col.names = TRUE, row.names = FALSE)



# maxquant
header<-read_excel(paste0(data_fold, 'metrics_maxquant.xlsx'), sheet='header', col_names = FALSE)
diann_pauc001<-read_excel(paste0(data_fold, 'metrics_maxquant.xlsx'), sheet='G-mean')
datasets<-unique(t(header[1,][7:length(header[1,])]))

header_dataset<-paste0(header[1,],'_', header[2,])
datas<-diann_pauc001[,7:116]
datas[is.na(datas)]<-0
colnames(datas)<-header_dataset[7:116]
library(corrplot)

corrmatrix = cor(datas, method = 'spearman')

rownames(corrmatrix) = paste(1:ncol(datas), names(datas),sep = ' ') # ?????????????????????????????? <?????? ?????????> ?????????
colnames(corrmatrix) = as.character(1:ncol(datas)) # ????????????????????????????????????
col3 <- colorRampPalette(c('DodgerBlue3','white', "OrangeRed")) # ??????????????????col3, ??????????????????????????????colorbar?????????
corrplot(corrmatrix, method = 'circle', diag = F, type = 'full', outline = F,
         col = col3(10), cl.lim = c(-0.2,1),addgrid.col = NA, 
         tl.pos = 'lb',tl.cex = 0.5, tl.col = 'black', tl.srt = 0, tl.offset = 0.5)
axis(1,at = 1:12, labels = NA, pos = 12.5, tck = -0.01) # ??????????????????????????????,????????????????????????,??????,????????????,???????????????
axis(4,at = 1:12, labels = NA, pos = 0.5, tck = -0.01) # ??????????????????????????????,??????????????????

#########################################
#checking average performance for prediction

metrics<-c('pAUC0.01', 'pAUC0.05', 'pAUC0.1', 'G-mean', 'nMCC')
data_fold<-'F:/NTU/quantification/data/2022.10.18/source_data_for_figures/'
header<-read_excel(paste0(data_fold, 'metrics_maxquant.xlsx'), sheet='header', col_names = FALSE)
header_dataset<-paste0(header[1,],'_', header[2,])
dts<-header[1,][7:116]
dt_n<-c()

for(i in dts[1,]){
  s = strsplit(i, '_', fixed = TRUE)[[1]]
  dt_n<-c(dt_n, s[1])
}

sprs<-vector()
for (metric in metrics) {
  #metric = 'pAUC0.01'
  diann_pauc001<-read_excel(paste0(data_fold, 'metrics_maxquant.xlsx'), sheet=metric)
  datas<-diann_pauc001[,7:116]
  datas[is.na(datas)]<-0
  colnames(datas)<-header_dataset[7:116]
  
  for (dt in unique(dt_n)) {
    mean_tr<-rowMeans(datas[,which(dt_n!=dt)])
    median_tr<-apply(datas[,which(dt_n!=dt)], 1, median)
    
    tests<-datas[,which(dt_n==dt)]
    for (i in 1:length(tests[1,])) {
      spr1 = cor(mean_tr, tests[, i], method = 'spearman')
      spr2 = cor(median_tr, tests[, i], method = 'spearman')
      sprs<-rbind(sprs, c(dt, metric, i, spr1, spr2))
    }
  }
}
colnames(sprs)<-c('dt_name', 'metric', 'dt_col', 'mean_spr', 'median_spr')
write.table(sprs, paste0(data_fold, 'maxquant_spr_test.csv'), sep = ',', col.names = TRUE, row.names = FALSE)


#################################################################
#################

#fragpipe_feature<-read.table(paste0(data_fold, 'fragpipe_dataset_feature.csv'), sep = ',', header = FALSE)
fragpipe_feature<-read.table(paste0(data_fold, 'maxquant_dataset_feature.csv'), sep = ',', header = FALSE)
#fragpipe_feature<-read.table(paste0(data_fold, 'DIANN_dataset_feature.csv'), sep = ',', header = FALSE)

labels<-fragpipe_feature$V42
dt_n<-c()
dt<-c()
for (i in labels) {
  lbs<-strsplit(i, '_', fixed = TRUE)[[1]]
  if(lbs[2]=='human'){
    dt_n<-c(dt_n, paste0(lbs[3],'_', lbs[4]))
    dt<-c(dt, lbs[3])
  }else{
    dt<-c(dt, lbs[2])
  dt_n<-c(dt_n, paste0(lbs[2],'_', lbs[3]))}
}

fragpipe_feature<-fragpipe_feature[,1:length(fragpipe_feature[1,])-1]

library(pls)
library("FactoMineR")
library(factoextra)
library("corrplot")
library(gridExtra)

dtt<-fragpipe_feature
dtt[which(is.na(dtt))]=0
res.pca <- PCA(dtt, graph = FALSE,scale.unit = FALSE)
pcs<-res.pca$ind$coord
colnames(pcs)<-c('PC1','PC2','PC3','PC4','PC5')
eigs<-as.vector(res.pca$eig[1:5,2])
p1<-ggplot(data = as.data.frame(pcs),mapping = aes(x = PC1, y =PC2,colour=dt, shape=dt))+geom_point(size=2)+ theme(legend.position="none")+labs(title='PCA',x=paste0('PC1(',substr(as.character(eigs[1]),1,5),'%)'),y=paste0('PC2(',substr(as.character(eigs[2]),1,5),'%)'))#+geom_text_repel(aes(label=Name))
p2<-ggplot(data = as.data.frame(pcs),mapping = aes(x = PC2, y =PC3,colour=dt, shape=dt))+geom_point(size=2)+ theme(legend.position="none")+labs(x=paste0('PC2(',substr(as.character(eigs[2]),1,5),'%)'),y=paste0('PC3(',substr(as.character(eigs[3]),1,5),'%)'))#+geom_text(aes(label=Name))
p3<-ggplot(data = as.data.frame(pcs),mapping = aes(x = PC3, y =PC4,colour=dt, shape=dt))+geom_point(size=2)+ theme(legend.position="right")+labs(x=paste0('PC3(',substr(as.character(eigs[3]),1,5),'%)'),y=paste0('PC4(',substr(as.character(eigs[4]),1,5),'%)'))#+geom_text(aes(label=Name))
p<-grid.arrange(p1,p2,p3,nrow=1)


fea_name<-c('mean_mr_g1', 'std_mr_g1', 'mean_mr_g2', 'std_mr_g2', 'mean_pr_g1', 'std_pr_g1', 'mean_pr_g2',
            'std_pr_g2', 'mean_pr_g1g2', 'std_pr_g1g2', 'mean_mr_g1_blgfc', 'std_mr_g1_blgfc', 'mean_mr_g2_blgfc',
            'std_mr_g2_blgfc', 'mean_pr_g1_blgfc', 'std_pr_g1_blgfc', 'mean_pr_g2_blgfc', 'std_pr_g2_blgfc',
            'mean_pr_g1g2_blgfc', 'std_pr_g1g2_blgfc', 'cov_g1', 'cov_g2', 'mean_cov_r_g1', 'std_cov_r_g1',
            'mean_cov_r_g2', 'std_cov_r_g2', 'cov_g1_blgfc', 'cov_g2_blgfc', 'mean_cov_r_g1_blgfc',
            'std_cov_r_g1_blgfc', 'mean_cov_r_g2_blgfc', 'std_cov_r_g2_blgfc',
            'mean_logFC', 'std_logFC', 'min_logFC', 'max_logFC', 'per_logFC_b1', 'per_p_s005', 'per_p_s001',
            'per_q_s005', 'per_q_s001')


features<-fragpipe_feature
colnames(features)<-fea_name

features1<-features[,1:10]
features2<-features[,11:20]
features3<-features[,21:30]
features4<-features[,31:41]
features1$dt<-dt
features2$dt<-dt
features3$dt<-dt
features4$dt<-dt

dtt1<-melt(features1)
dtt2<-melt(features2)
dtt3<-melt(features3)
dtt4<-melt(features4)

colnames(dtt1)<-c('dataset', 'feature', 'value')
colnames(dtt2)<-c('dataset', 'feature', 'value')
colnames(dtt3)<-c('dataset', 'feature', 'value')
colnames(dtt4)<-c('dataset', 'feature', 'value')

p1<-ggboxplot(dtt1, "feature", "value",width = 0.5, size=0.8, #outlier.shape=NA,
              color = "dataset",
              add = "none")+labs(x = 'feature', y = '')+
  scale_x_discrete(limits=unique(dtt1$feature))+
  theme(plot.title = element_text(size = 14))+theme_classic() + 
  theme(axis.text.x = element_text(size = 12, angle = 0), axis.text.y = element_text(size = 12))+
  theme(legend.position = "none")+
  theme(plot.title = element_text(size = 14,hjust = 0.5, face = "bold"))+coord_flip()
p2<-ggboxplot(dtt2, "feature", "value",width = 0.5, size=0.8, #outlier.shape=NA,
              color = "dataset",
              add = "none")+labs(x = 'feature', y = '')+
  scale_x_discrete(limits=unique(dtt2$feature))+
  theme(plot.title = element_text(size = 14))+theme_classic() + 
  theme(axis.text.x = element_text(size = 12, angle = 0), axis.text.y = element_text(size = 12))+
  theme(legend.position = "none")+
  theme(plot.title = element_text(size = 14,hjust = 0.5, face = "bold"))+coord_flip()
p3<-ggboxplot(dtt3, "feature", "value",width = 0.5, size=0.8, #outlier.shape=NA,
              color = "dataset",
              add = "none")+labs(x = 'feature', y = '')+
  scale_x_discrete(limits=unique(dtt3$feature))+
  theme(plot.title = element_text(size = 14))+theme_classic() + 
  theme(axis.text.x = element_text(size = 12, angle = 0), axis.text.y = element_text(size = 12))+
  theme(legend.position = "none")+
  theme(plot.title = element_text(size = 14,hjust = 0.5, face = "bold"))+coord_flip()
p4<-ggboxplot(dtt4, "feature", "value",width = 0.5, size=0.8, #outlier.shape=NA,
              color = "dataset",
              add = "none")+labs(x = 'feature', y = '')+
  scale_x_discrete(limits=unique(dtt4$feature))+
  theme(plot.title = element_text(size = 14))+theme_classic() + 
  theme(axis.text.x = element_text(size = 12, angle = 0), axis.text.y = element_text(size = 12))+
  theme(legend.position = "none")+
  theme(plot.title = element_text(size = 14,hjust = 0.5, face = "bold"))+coord_flip()
  #geom_hline(aes(yintercept = 0, colour = "red"))
p<-grid.arrange(p1,p2,p3, p4,nrow=2)



###

#test_frag<-read.table(paste0(data_fold, 'fragpipe_spr_test.csv'), header = TRUE, sep = ',')
#test_frag<-read.table(paste0(data_fold, 'maxquant_spr_test.csv'), header = TRUE, sep = ',')
test_frag<-read.table(paste0(data_fold, 'DIANN_spr_test.csv'), header = TRUE, sep = ',')
p1<-ggboxplot(test_frag, "dt_name", "mean_spr",width = 0.5, size=0.8, #outlier.shape=NA,
          color = "metric", title = 'mean_performance',
          add = "none")+labs(x = 'dt_name', y = '')+
  scale_x_discrete(limits=unique(test_frag$dt_name))+
  theme(plot.title = element_text(size = 14))+theme_classic() + 
  theme(axis.text.x = element_text(size = 12, angle = 0), axis.text.y = element_text(size = 12))+
  theme(legend.position = "none")+
  theme(plot.title = element_text(size = 14,hjust = 0.5, face = "bold"))

p2<-ggboxplot(test_frag, "dt_name", "median_spr",width = 0.5, size=0.8, #outlier.shape=NA,
              color = "metric", title = 'median_performance',
              add = "none")+labs(x = 'dt_name', y = '')+
  scale_x_discrete(limits=unique(test_frag$dt_name))+
  theme(plot.title = element_text(size = 14))+theme_classic() + 
  theme(axis.text.x = element_text(size = 12, angle = 0), axis.text.y = element_text(size = 12))+
  theme(legend.position = "right")+
  theme(plot.title = element_text(size = 14,hjust = 0.5, face = "bold"))
p<-grid.arrange(p1,p2,nrow=1)

#########################
#################################################
rank_values<-function(rk_values, all_values){
  s<-data.frame(all_values$workflow, rk_values, c(1:length(all_values$workflow)))
  colnames(s)<-c('workflow','value', 'idx')
  s<-s[order(-s$value),]
  s$rank=c(1:length(s$value))
  s<-s[order(s$workflow),]
}

spr_test<-function(platform){
  library("readxl")
  data_fold<-'F:/NTU/quantification/data/2022.10.18/source_data_for_figures/'
  if(platform=='diann'){
    st=7
    ed=152
  }else if(platform=='fragpipe'){
    st=7
    ed=126
  }else if(platform=='maxquant'){
    st=7
    ed=116
  }
  pauc001<-read_excel(paste0(data_fold, 'metrics_',platform,'.xlsx'), sheet='pAUC0.01')
  pauc005<-read_excel(paste0(data_fold, 'metrics_',platform,'.xlsx'), sheet='pAUC0.05')
  pauc01<-read_excel(paste0(data_fold, 'metrics_',platform,'.xlsx'), sheet='pAUC0.1')
  nMCC<-read_excel(paste0(data_fold, 'metrics_',platform,'.xlsx'), sheet='nMCC')
  Gmean<-read_excel(paste0(data_fold, 'metrics_',platform,'.xlsx'), sheet='G-mean')
  
  header<-read_excel(paste0(data_fold, 'metrics_',platform,'.xlsx'), sheet='header', col_names = FALSE)
  header_dataset<-paste0(header[1,],'_', header[2,])
  dts<-header[1,][st:ed]
  dt_n<-c()
  
  for(i in dts[1,]){
    s = strsplit(i, '_', fixed = TRUE)[[1]]
    dt_n<-c(dt_n, s[1])
  }
  
  sprs<-vector()
  #leave-one-dataset-out-cv
  for (dt in unique(dt_n)) {
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
    
    ranks_avg_median<-rowMeans(cbind(rank_median_tr_pauc001$rank, rank_median_tr_pauc005$rank,
                                     rank_median_tr_pauc01$rank, rank_median_tr_nMCC$rank,
                                     rank_median_tr_Gmean$rank))
    
    tests_pauc001<-dt_pauc001[,which(dt_n==dt)]
    tests_pauc005<-dt_pauc005[,which(dt_n==dt)]
    tests_pauc01<-dt_pauc01[,which(dt_n==dt)]
    tests_nMCC<-dt_nMCC[,which(dt_n==dt)]
    tests_Gmean<-dt_Gmean[,which(dt_n==dt)]
    
    for (i in 1:length(tests_pauc001[1,])) {
      rk_te_pauc001<-rank_values(tests_pauc001[,i], pauc001)
      rk_te_pauc005<-rank_values(tests_pauc005[,i], pauc005)
      rk_te_pauc01<-rank_values(tests_pauc01[,i], pauc01)
      rk_te_Gmean<-rank_values(tests_Gmean[,i], Gmean)
      rk_te_nMCC<-rank_values(tests_nMCC[,i], nMCC)
      
      rk_avg_te<-rowMeans(cbind(rk_te_pauc001$rank, rk_te_pauc005$rank,
                                rk_te_pauc01$rank, rk_te_Gmean$rank,
                                rk_te_nMCC$rank))
      
      spr1 = cor(ranks_avg_mean, rk_avg_te, method = 'spearman')
      spr2 = cor(ranks_avg_median, rk_avg_te, method = 'spearman')
      sprs<-rbind(sprs, c(dt, 'avg_rank', i, spr1, spr2))
    }
  }
  colnames(sprs)<-c('dt_name', 'metric', 'dt_col', 'mean_spr', 'median_spr')
  write.table(sprs, paste0(data_fold, platform, '_spr_avg_rank_test.csv'), sep = ',', col.names = TRUE, row.names = FALSE)
  
}

spr_test('diann')
spr_test('fragpipe')
spr_test('maxquant')

library(ComplexHeatmap)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(cowplot)
library("ggsci")
data_fold<-'F:/NTU/quantification/data/2022.10.18/source_data_for_figures/'
sprs_frag<-read.table(paste0(data_fold, 'fragpipe_spr_avg_rank_test.csv'), sep = ',', header = TRUE)

dt_mean<-data.frame(sprs_frag$dt_name,sprs_frag$mean_spr, c(rep('mean',length(sprs_frag$dt_name))))
colnames(dt_mean)<-c('project', 'spr', 'method')
dt_median<-data.frame(sprs_frag$dt_name,sprs_frag$median_spr, c(rep('median',length(sprs_frag$dt_name))))
colnames(dt_median)<-c('project', 'spr', 'method')

dt<-rbind(dt_mean, dt_median)
dt$project[which(dt$project=='human')] = 'human_ecoli'

p1<-ggboxplot(dt, "project", "spr",width = 0.5, size=0.8, #outlier.shape=NA,
              color = "method", title = 'FragPipe',
              add = "none")+labs(x = 'project', y = 'Spearman correlation')+
  scale_x_discrete(limits=unique(dt$project))+
  theme(plot.title = element_text(size = 14))+theme_classic() + 
  theme(axis.text.x = element_text(size = 12, angle = 0), axis.text.y = element_text(size = 12))+
  theme(legend.position = "right")+
  theme(plot.title = element_text(size = 14,hjust = 0.5, face = "bold"))+scale_color_d3()
p1


sprs_frag<-read.table(paste0(data_fold, 'maxquant_spr_avg_rank_test.csv'), sep = ',', header = TRUE)

dt_mean<-data.frame(sprs_frag$dt_name,sprs_frag$mean_spr, c(rep('mean',length(sprs_frag$dt_name))))
colnames(dt_mean)<-c('project', 'spr', 'method')
dt_median<-data.frame(sprs_frag$dt_name,sprs_frag$median_spr, c(rep('median',length(sprs_frag$dt_name))))
colnames(dt_median)<-c('project', 'spr', 'method')

dt<-rbind(dt_mean, dt_median)
dt$project[which(dt$project=='human')] = 'human_ecoli'

p2<-ggboxplot(dt, "project", "spr",width = 0.5, size=0.8, #outlier.shape=NA,
              color = "method", title = 'maxquant',
              add = "none")+labs(x = 'project', y = 'Spearman correlation')+
  scale_x_discrete(limits=unique(dt$project))+
  theme(plot.title = element_text(size = 14))+theme_classic() + 
  theme(axis.text.x = element_text(size = 12, angle = 0), axis.text.y = element_text(size = 12))+
  theme(legend.position = "right")+
  theme(plot.title = element_text(size = 14,hjust = 0.5, face = "bold"))+scale_color_d3()
p2

sprs_frag<-read.table(paste0(data_fold, 'diann_spr_avg_rank_test.csv'), sep = ',', header = TRUE)

dt_mean<-data.frame(sprs_frag$dt_name,sprs_frag$mean_spr, c(rep('mean',length(sprs_frag$dt_name))))
colnames(dt_mean)<-c('project', 'spr', 'method')
dt_median<-data.frame(sprs_frag$dt_name,sprs_frag$median_spr, c(rep('median',length(sprs_frag$dt_name))))
colnames(dt_median)<-c('project', 'spr', 'method')

dt<-rbind(dt_mean, dt_median)
dt$project[which(dt$project=='human')] = 'human_ecoli'

p3<-ggboxplot(dt, "project", "spr",width = 0.5, size=0.8, #outlier.shape=NA,
              color = "method", title = 'DIANN',
              add = "none")+labs(x = 'project', y = 'Spearman correlation')+
  scale_x_discrete(limits=unique(dt$project))+
  theme(plot.title = element_text(size = 14))+theme_classic() + 
  theme(axis.text.x = element_text(size = 12, angle = 0), axis.text.y = element_text(size = 12))+
  theme(legend.position = "right")+
  theme(plot.title = element_text(size = 14,hjust = 0.5, face = "bold"))+scale_color_d3()
p3