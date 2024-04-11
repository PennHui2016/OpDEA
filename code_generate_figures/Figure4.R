#####################################
##
## Figure 4 compares options in each step of a workflow
##
###################

####################################################################
data_fold<-paste0(dirname(rstudioapi::getSourceEditorContext()$path),'/data/')
save_fold<-paste0(dirname(rstudioapi::getSourceEditorContext()$path),'/save_data/')
## compare mean differences
library("readxl")
library(ggplot2)
library(reshape2)
library(ggpubr)
library(cowplot)

rank_with_diff<-function(out_data){
  tools<-unique(colnames(out_data))
  
  win_pro<-c()
  for (i in 1:length(tools)) {
    num=0
    for (j in setdiff(c(1:length(tools)),i)) {
      med_diff=median(out_data[,i]-out_data[,j])
      mea_diff=mean(out_data[,i]-out_data[,j])
      if(mea_diff>0){
        num=num+1
      }
      if(mea_diff==0 & med_diff>0){
        num=num+1
      }
    }
    win_pro<-c(win_pro,num/(length(tools)-1))
  }
  all_tool=tools
  rank<-vector()
  r=1
  while (length(all_tool)>1) {
    ind<-which(win_pro==max(win_pro))
    if(length(ind)==1){
      remain_ind<-setdiff(c(1:length(all_tool)),ind)
      rank<-rbind(rank,c(all_tool[ind],win_pro[ind],r))
      r=r+1
      all_tool<-all_tool[remain_ind]
      win_pro<-win_pro[remain_ind]
    }else if(length(ind)>1){
      degree<-c()
      st<-all_tool[ind]
      for (m in 1:length(st)) {
        num=0
        for (n in setdiff(c(1:length(st)),m)) {
          med=median(out_data[,which(colnames(out_data)==st[m])]-out_data[,which(colnames(out_data)==st[n])])
          mea=mean(out_data[,which(colnames(out_data)==st[m])]-out_data[,which(colnames(out_data)==st[n])])
          if(mea>0){
            num=num+1
          }
          if(mea==0 & med>0){
            num=num+1
          }
        }
        degree<-c(degree,num)
      }
      strk<-data.frame(cbind(st,degree))
      strk<-strk[order(strk$degree, decreasing=TRUE),]
      for (s in 1:length(strk$st)) {
        idx<-which(all_tool==strk$st[s])
        remain<-setdiff(c(1:length(all_tool)),idx)
        rank<-rbind(rank,c(strk$st[s],paste0(as.character(win_pro[ind[1]]),'|',as.character(strk$degree[s])),r))
        r=r+1
        all_tool<-all_tool[remain]
        win_pro<-win_pro[remain]
      }
    }
  }
  rank<-rbind(rank,c(all_tool[1],win_pro[1],r))
  return(rank)
}

get_paired_data<-function(save_dir, data_fold, metric, procedure, st, en, platform, acq){
  if (metric == 'pauc001'){
    sheet_name = 'pAUC0.01'
  }else if(metric == 'pauc005'){
    sheet_name = 'pAUC0.05'
  }else if(metric == 'pauc01'){
    sheet_name = 'pAUC0.1'
  }else if(metric == 'nMCC'){
    sheet_name = 'nMCC0.05'
  }else if(metric == 'G-mean'){
    sheet_name = 'G-mean0.05'
  }
  res_data<-read_excel(paste0(data_fold,'metrics_', platform, '_', acq,'.xlsx'), sheet = sheet_name)
  res_data<-res_data[,1:en]
  res_data[is.na(res_data)]=0
  res_data<-res_data[which(res_data$DEA!='MSstats'),]
  
  if(procedure=='Matrix'){
    res_data<-res_data[which(res_data$Matrix!='count'),]
  }
  
  if(procedure=='DEA_count'){
    res_data<-res_data[which(res_data$Matrix=='count'),]
  }else if(procedure=='DEA'){
    res_data<-res_data[which(res_data$Matrix!='count'),]
  }
  
  if (!dir.exists(save_dir)){
    dir.create(save_dir)
  }
  if(procedure=='DEA'){
    uni_matrix<-sort(unique(res_data$DEA))
    mat0<-uni_matrix[1]
    ind_mat0<-which(res_data$DEA==mat0)
  }else if(procedure=='Matrix'){
    uni_matrix<-sort(unique(res_data$Matrix))
    mat0<-uni_matrix[1]
    ind_mat0<-which(res_data$Matrix==mat0)
  }else if(procedure=='Imput'){
    uni_matrix<-sort(unique(res_data$Imput))
    mat0<-uni_matrix[1]
    ind_mat0<-which(res_data$Imput==mat0)
  }else if(procedure=='normalization'){
    uni_matrix<-sort(unique(res_data$normalization))
    mat0<-uni_matrix[1]
    ind_mat0<-which(res_data$normalization==mat0)
  }else if(procedure=='DEA_count'){
    uni_matrix<-sort(unique(res_data$DEA))
    mat0<-uni_matrix[1]
    ind_mat0<-which(res_data$DEA==mat0)
  }
  
  new_wf_all<-paste0(res_data$DEA, res_data$Matrix, res_data$Imput,
                     res_data$normalization)
  sl_idxs<-vector()
  
  for (i in 1:length(ind_mat0)) {
    ext=1
    idxs<-c()
    idxs<-c(idxs,ind_mat0[i])
    wf = new_wf_all[ind_mat0][i]
    for (j in 2:length(uni_matrix)) {
      wf_2<-gsub(mat0, uni_matrix[j], wf)
      if(procedure=='DEA'){
        ind_matj<-which(res_data$DEA==uni_matrix[j])
      }else if(procedure=='Matrix'){
        ind_matj<-which(res_data$Matrix==uni_matrix[j])
      }else if(procedure=='Imput'){
        ind_matj<-which(res_data$Imput==uni_matrix[j])
      }else if(procedure=='normalization'){
        ind_matj<-which(res_data$normalization==uni_matrix[j])
      }else if(procedure=='DEA_count'){
        ind_matj<-which(res_data$DEA==uni_matrix[j])
      }
      
      idx<-which(new_wf_all[ind_matj]==wf_2)
      if(length(idx)==1){
        idxs<-c(idxs,ind_matj[idx])
        ext=ext+1
      }
    }
    if (ext==length(uni_matrix)){
      sl_idxs<-rbind(sl_idxs,idxs)
    }
  }
  
  
  out_data<-vector()
  rowname_all<-c()
  for (i in 1:length(uni_matrix)) {
    vl<-c()
    rowname<-c()
    for (j in 1:length(sl_idxs[,1])) {
      vl<-c(vl,as.numeric(res_data[sl_idxs[j,i],][,st:en]))
      rowname<-c(rowname,paste0(gsub(uni_matrix[i], '', res_data$workflow[sl_idxs[j,i]]),'_',colnames(res_data)[st:en]))
    }
    rowname_all<-cbind(rowname_all, rowname)
    out_data<-cbind(out_data,vl)
  }
  colnames(out_data)<-uni_matrix
  row.names(out_data)<-rowname_all[,1]
  
  return(out_data)
}

calculte_diffs<-function(save_dir, out_data, metric, procedure, platform, acq){
  uni_matrix=colnames(out_data)
  cormat <- round(cor(out_data), 2)
  abbrs_all<-c('A', 'B','C','D','E','F','G','H','I','J','K','L','M','N','O', 'P', 'Q', 'R', 'S', 'T')
  
  abbrs<-abbrs_all[1:length(uni_matrix)]
  diffs_ht<-matrix(nrow=length(uni_matrix),ncol=length(uni_matrix))
  ps<-matrix(nrow=length(uni_matrix),ncol=length(uni_matrix))
  for(i in 1:(length(uni_matrix)-1)){
    for (j in (i+1):length(uni_matrix)) {
      diffs_ht[j,i]=mean(out_data[,i]-out_data[,j])
      ps[j,i]<-t.test(out_data[,i],out_data[,j], paired = TRUE)$p.value
    }
  }
  colnames(diffs_ht)<-uni_matrix
  row.names(diffs_ht)<-uni_matrix
  colnames(ps)<-uni_matrix
  row.names(ps)<-uni_matrix
  diffs_ht<-round(diffs_ht,2)
  
  
  diffs<-vector()
  
  for(i in 1:(length(uni_matrix)-1)){
    for (j in (i+1):length(uni_matrix)) {
      p<-t.test(out_data[,j],out_data[,i], paired = TRUE)$p.value
      if(is.nan(p)){
        p<-NA
      }
      if(p>0.05 | is.na(p)){
        sigf='na'
      }else if(p<0.05 & p>=0.01){
        sigf='*'
      }else if(p<0.01 & p>=0.001){
        sigf='**'
      }else if(p<0.001 & p>=0.0001){
        sigf='***'
      }else if(p<=0.0001){
        sigf='****'
      }
      diffs<-rbind(diffs, c(uni_matrix[i],uni_matrix[j],abbrs[i],abbrs[j],
                            mean(out_data[,i]-out_data[,j]),
                            median(out_data[,i]-out_data[,j]),
                            p,sigf))
    }
  }
  colnames(diffs)<-c('mat1','mat2','abbr1','abbr2','mean','median','t_pval','significance')
  diffs<-as.data.frame(diffs)
  diffs$contrast<-paste0(diffs$mat1,'-',diffs$mat2)
  
  write.table(diffs, paste0(save_dir, metric, '_comp_', procedure,'_',platform,'_', acq ,'mean_diffs.csv'), sep = ',')
  diffs$mean<-round(as.numeric(diffs$mean),3)
  
  return(list(diffs=diffs, diff_ht=diffs_ht))
}

avg_rank<-function(platform, acq, procedure, wb){
  data_folder<-paste0(save_fold, platform,"_", acq,'/')
  ranks<-vector()
  for (metric in c('pauc001', 'pauc005', 'pauc01', 'nMCC','G-mean')) {
    rank_metric<-read.table(paste0(data_folder, metric, '_comp_', procedure,'_',platform,'_', acq ,'_rank.csv'), sep = ',')
    colnames(rank_metric)<-paste0(metric,c('feature','sr','rank'))
    if(length(ranks)==0){
      ranks<-rank_metric
    }else{
      idx<-match(ranks[,1], rank_metric[,1])
      ranks<-cbind(ranks,rank_metric[idx,])
    }
  }
  ranks_only<-ranks[,grep('rank',colnames(ranks))]
  avg_rank<-rowMeans(ranks_only)
  ranks<-cbind(ranks,avg_rank)
  colnames(ranks)[length(ranks[1,])]<-'avg_rank'
  sheet=procedure
  addWorksheet(wb,sheet)
  writeData(wb, sheet, ranks, startRow = 1, startCol = 1, colNames = TRUE, rowNames = TRUE)
  saveWorkbook(wb, file = paste0(data_folder, platform,'_', acq,'avg_rankings.xlsx'), overwrite = TRUE)
  write.table(ranks, paste0(data_folder, 'avg_rank_',procedure, '_',platform,'_',acq,'.csv'),sep = ',',col.names = T)
  return(ranks)
}


procedure_compare<-function(data_fold, procedure, st, en, platform, acq, wb){ #ranks, 
  save_dir<-paste0(save_fold, platform,"_", acq,'/')
  if (!dir.exists(save_dir)){
    dir.create(save_dir)
  }
  
  out_data_pauc001<-get_paired_data(save_dir, data_fold, 'pauc001', procedure, st, en, platform, acq)
  out_data_pauc005<-get_paired_data(save_dir, data_fold, 'pauc005', procedure, st, en, platform, acq)
  out_data_pauc01<-get_paired_data(save_dir, data_fold, 'pauc01', procedure, st, en, platform, acq)
  out_data_nMCC<-get_paired_data(save_dir, data_fold, 'nMCC', procedure, st, en, platform, acq)
  out_data_gmean<-get_paired_data(save_dir, data_fold, 'G-mean', procedure, st, en, platform, acq)
  
  ref_out<-out_data_pauc001
  
  abbrs_all<-c('A', 'B','C','D','E','F','G','H','I','J','K','L','M','N','O', 'P', 'Q', 'R', 'S', 'T')
  uni_matrix = colnames(out_data_pauc001)
  abbrs<-abbrs_all[1:length(uni_matrix)]
  out_data_pauc001<-out_data_pauc001[which(apply(ref_out, 1, function(x) min(x))!=0),]
  write.table(out_data_pauc001, paste0(save_dir, 'pauc001', '_comp_', procedure,'_',platform,'_', acq ,'_mat.csv'), sep = ',')
  rank_metric<-rank_with_diff(out_data_pauc001)
  write.table(rank_metric, paste0(save_dir, 'pauc001', '_comp_', procedure,'_',platform,'_', acq ,'_rank.csv'), sep = ',')
  calculte_diffs(save_dir, out_data_pauc001, 'pauc001', procedure, platform, acq)
  
  out_data_pauc005<-out_data_pauc005[which(apply(ref_out, 1, function(x) min(x))!=0),]
  write.table(out_data_pauc005, paste0(save_dir, 'pauc005', '_comp_', procedure,'_',platform,'_', acq ,'_mat.csv'), sep = ',')
  rank_metric<-rank_with_diff(out_data_pauc005)
  write.table(rank_metric, paste0(save_dir, 'pauc005', '_comp_', procedure,'_',platform,'_', acq ,'_rank.csv'), sep = ',')
  calculte_diffs(save_dir, out_data_pauc005, 'pauc005', procedure, platform, acq)
  
  out_data_pauc01<-out_data_pauc01[which(apply(ref_out, 1, function(x) min(x))!=0),]
  write.table(out_data_pauc01, paste0(save_dir, 'pauc01', '_comp_', procedure,'_',platform,'_', acq ,'_mat.csv'), sep = ',')
  rank_metric<-rank_with_diff(out_data_pauc01)
  write.table(rank_metric, paste0(save_dir, 'pauc01', '_comp_', procedure,'_',platform,'_', acq ,'_rank.csv'), sep = ',')
  calculte_diffs(save_dir, out_data_pauc01, 'pauc01', procedure, platform, acq)
  
  out_data_nMCC<-out_data_nMCC[which(apply(ref_out, 1, function(x) min(x))!=0),]
  write.table(out_data_nMCC, paste0(save_dir, 'nMCC', '_comp_', procedure,'_',platform,'_', acq ,'_mat.csv'), sep = ',')
  rank_metric<-rank_with_diff(out_data_nMCC)
  write.table(rank_metric, paste0(save_dir, 'nMCC', '_comp_', procedure,'_',platform,'_', acq ,'_rank.csv'), sep = ',')
  calculte_diffs(save_dir, out_data_nMCC, 'nMCC', procedure, platform, acq)
  
  out_data_gmean<-out_data_gmean[which(apply(ref_out, 1, function(x) min(x))!=0),]
  write.table(out_data_gmean, paste0(save_dir, 'G-mean', '_comp_', procedure,'_',platform,'_', acq ,'_mat.csv'), sep = ',')
  rank_metric<-rank_with_diff(out_data_gmean)
  write.table(rank_metric, paste0(save_dir, 'G-mean', '_comp_', procedure,'_',platform,'_', acq ,'_rank.csv'), sep = ',')
  calculte_diffs(save_dir, out_data_gmean, 'G-mean', procedure, platform, acq)
  
  
  avg_rank(platform, acq, procedure, wb)
  
}

library(openxlsx)
wb5 <- createWorkbook()
procedure_compare(data_fold, 'Matrix', 7, 28, 'FragPipe','DDA',wb5)
procedure_compare(data_fold, 'DEA', 7, 28, 'FragPipe','DDA',wb5)
procedure_compare(data_fold, 'DEA_count', 7, 28, 'FragPipe','DDA',wb5)
procedure_compare(data_fold, 'Imput', 7, 28, 'FragPipe','DDA',wb5)
procedure_compare(data_fold, 'normalization', 7, 28, 'FragPipe','DDA',wb5)

wb5 <- createWorkbook()
procedure_compare(data_fold, 'Matrix', 7, 28, 'Maxquant','DDA',wb5)
procedure_compare(data_fold, 'DEA', 7, 28, 'Maxquant','DDA',wb5)
procedure_compare(data_fold, 'DEA_count', 7, 28, 'Maxquant','DDA',wb5)
procedure_compare(data_fold, 'Imput', 7, 28, 'Maxquant','DDA',wb5)
procedure_compare(data_fold, 'normalization', 7, 28, 'Maxquant','DDA',wb5)

wb5 <- createWorkbook()
procedure_compare(data_fold, 'Matrix', 7, 21, 'FragPipe','TMT',wb5)
procedure_compare(data_fold, 'DEA', 7, 21, 'FragPipe','TMT',wb5)
procedure_compare(data_fold, 'Imput', 7, 21, 'FragPipe','TMT',wb5)
procedure_compare(data_fold, 'normalization', 7, 21, 'FragPipe','TMT',wb5)

wb5 <- createWorkbook()
procedure_compare(data_fold, 'DEA', 7, 21, 'Maxquant','TMT',wb5)
procedure_compare(data_fold, 'Imput', 7, 21, 'Maxquant','TMT',wb5)
procedure_compare(data_fold, 'normalization', 7, 21, 'Maxquant','TMT',wb5)

wb5 <- createWorkbook()
procedure_compare(data_fold, 'Matrix', 7, 23, 'DIANN','DIA',wb5)
procedure_compare(data_fold, 'DEA', 7, 23, 'DIANN','DIA',wb5)
procedure_compare(data_fold, 'Imput', 7, 23, 'DIANN','DIA',wb5)
procedure_compare(data_fold, 'normalization',7, 23, 'DIANN','DIA',wb5)

wb5 <- createWorkbook()
procedure_compare(data_fold, 'Matrix', 7, 23, 'spt','DIA',wb5)
procedure_compare(data_fold, 'DEA', 7, 23, 'spt','DIA',wb5)
procedure_compare(data_fold, 'Imput', 7, 23, 'spt','DIA',wb5)
procedure_compare(data_fold, 'normalization',7, 23, 'spt','DIA',wb5)


paired_data_DIANN_DIA<-read.table(paste0(save_fold, 'DIANN_DIA/','pauc001_comp_Matrix_DIANN_DIA_mat.csv'),sep = ',', header = T)
paired_data_spt_DIA<-read.table(paste0(save_fold, 'spt_DIA/','pauc001_comp_Matrix_spt_DIA_mat.csv'),sep = ',', header = T)

uni_sc<-unique(colnames(paired_data_DIANN_DIA))
cb_diff_sc_fg<-vector()
for (i in 1:(length(uni_sc)-1)) {
  for (j in (i+1):length(uni_sc)) {
    cb_diff_sc_fg<-rbind(cb_diff_sc_fg,data.frame(contrast=paste0(uni_sc[i],'-',uni_sc[j]),
                                                  diff=paired_data_DIANN_DIA[,which(colnames(paired_data_DIANN_DIA)==uni_sc[i])]-
                                                    paired_data_DIANN_DIA[,which(colnames(paired_data_DIANN_DIA)==uni_sc[j])],
                                                  setting='DIANN_DIA'))
  }
}

cb_diff_sc_mq<-vector()
for (i in 1:(length(uni_sc)-1)) {
  for (j in (i+1):length(uni_sc)) {
    cb_diff_sc_mq<-rbind(cb_diff_sc_mq,data.frame(contrast=paste0(uni_sc[i],'-',uni_sc[j]),
                                                  diff=paired_data_spt_DIA[,which(colnames(paired_data_spt_DIA)==uni_sc[i])]-
                                                    paired_data_spt_DIA[,which(colnames(paired_data_spt_DIA)==uni_sc[j])],
                                                  setting='spt_DIA'))
  }
}



cb_diff<-rbind(cb_diff_sc_fg,cb_diff_sc_mq)

rank_diann<-read.table(paste0(save_fold, 'DIANN_DIA/','pauc001_comp_Matrix_DIANN_DIA_rank.csv'),sep = ',', header = T)
rank_spt<-read.table(paste0(save_fold, 'spt_DIA/','pauc001_comp_Matrix_spt_DIA_rank.csv'),sep = ',', header = T)
diff_diann<-read.table(paste0(save_fold, 'DIANN_DIA/','pauc001_comp_Matrix_DIANN_DIAmean_diffs.csv'),sep = ',', header = T)
diff_spt<-read.table(paste0(save_fold, 'spt_DIA/','pauc001_comp_Matrix_spt_DIAmean_diffs.csv'),sep = ',', header = T)


col=c("#a6cee3","#33a02c","#fb9a99","#fdbf6f","#1f78b4","#b2df8a")
p8<-ggplot(cb_diff,aes(x=contrast, y=diff, color=setting))+ 
  geom_boxplot(outlier.alpha = 0.1)+ #outlier.shape = NA,varwidth = FALSE
  
  stat_summary(fun=mean, 
               geom="point",
               shape=20, 
               aes(group=setting, shape="mean"),
               position=position_dodge(.75), size=2, color="red", fill="red", show.legend = FALSE)+
  
  geom_hline(aes(yintercept = 0), linetype = 2,colour='gray',show.legend = FALSE)+
  guides(color=guide_legend(title=NULL, nrow=1, byrow=TRUE))+ylim(-0.25,0.5)+scale_color_manual(values = col)+#
  theme_pubr()+geom_point(aes(shape = "mean"), colour = 'red',alpha = 0)+  # <-- added this line of code and next
  guides(shape=guide_legend(title=NULL, override.aes = list(alpha = 1)))+
  xlab('')+ylab('pairwised difference of pAUC(0.01)')
p8

write.table(cb_diff, paste0(save_fold, 'Figure4A.csv'), col.names = T, row.names = F)
################ heatmap for comparing all procedures
# obtain running time

library(openxlsx)
library("readxl")

data_fold<-paste0(dirname(rstudioapi::getSourceEditorContext()$path),'/data/')
save_fold<-paste0(dirname(rstudioapi::getSourceEditorContext()$path),'/save_data/')

DEA_run_times<-read_excel(paste0(data_fold,'/times.xlsx'), sheet = 'DEA', col_names = F)
colnames(DEA_run_times)<-c('Method','time')
Imput_run_times<-read_excel(paste0(data_fold,'/times.xlsx'), sheet = 'Imput', col_names = F)
colnames(Imput_run_times)<-c('Method','time')
norm_run_times<-read_excel(paste0(data_fold,'/times.xlsx'), sheet = 'normalization', col_names = F)
colnames(norm_run_times)<-c('Method','time')

platforms<-c('Fragpipe', 'Maxquant', 'Maxquant', 'Fragpipe','DIANN', 'spt')
acqs<-c('DDA', 'DDA','TMT', 'TMT', 'DIA', 'DIA')
setting<-c('FG_DDA','MQ_DDA','MQ_TMT','FG_TMT','DIANN_DIA','spt_DIA')
procedures<-c('DEA','DEA_count','Imput','normalization','Matrix')

scale_values <- function(x){(x-min(x))/(max(x)-min(x))}

ranks<-list()
heatmap_ranks<-list()
label_heatmap_ranks<-list()
for (i in 1:length(procedures)) {
  procedure = procedures[i]
  rank<-vector()
  for (j in 1:length(setting)) {
    rk_file<-paste0(save_fold, platforms[j],'_',acqs[j],'/',
                    'avg_rank_',procedure,'_',platforms[j],'_',acqs[j],'.csv')
    if (file.exists(rk_file)){
      rank_procedure<-read.table(rk_file,sep = ',',header = T)
      rank_procedure$setting<-setting[j]
      rank<-rbind(rank,cbind(rank_procedure$pauc001feature, rank_procedure$avg_rank, rank_procedure$setting))}
  }
  rank<-as.data.frame(rank)
  colnames(rank)<-c('Method', 'avg_rank', 'setting')
  rank$Method<-gsub('blank','none',rank$Method)
  uni_method<-unique(rank$Method)
  uni_setting<-unique(rank$setting)
  heatmap_rank<-data.frame(Method=uni_method)
  label_heatmap_rank<-data.frame(Method=uni_method)
  for (k in uni_setting) {
    rank_setting<-c(rep(NA, length(heatmap_rank$Method)))
    mat_setting<-rank[which(rank$setting==k),]
    idx<-match(mat_setting$Method, uni_method)
    
    rank_setting[idx]<-mat_setting$avg_rank
    heatmap_rank<-cbind(heatmap_rank, as.numeric(rank_setting))
  }
  colnames(heatmap_rank)[2:length(heatmap_rank[1,])]<-uni_setting
  
  #heatmap_rank$average<-apply(as.matrix(heatmap_rank[, 2:length(heatmap_rank[1, ])]), 1, function(x) mean(x[!is.na(x)]))
  heatmap_rank$average<-apply(as.matrix(heatmap_rank[, 2:length(heatmap_rank[1, ])]), 1, function(x) mean(x))
  for (l in 2:length(heatmap_rank)) {
    rank_setting=heatmap_rank[,l]
    rank_exact<-rank_setting
    rank_ava<-as.numeric(rank_setting)
    rank_ava[is.na(rank_ava)]=100
    flag=1
    ranked<-c()
    
    while (flag<=length(rank_ava)) {
      rst<-which(rank_ava==min(setdiff(rank_ava,ranked)))
      if(length(which(is.na(rank_exact[rst])))==0){
        rank_exact[rst]=flag
        ranked<-c(ranked, rank_ava[rst])
        flag=flag+length(rst)
      }else{
        flag=flag+length(rst)
      }
      
    }
    label_heatmap_rank<-cbind(label_heatmap_rank, as.numeric(rank_exact))
    #heatmap_rank<-cbind(heatmap_rank, rank_setting)
    
    
  }
  colnames(label_heatmap_rank)<-colnames(heatmap_rank)
  
  #label_heatmap_rank$average<-apply(as.matrix(label_heatmap_rank[, 2:length(label_heatmap_rank[1, ])]), 1, function(x) mean(x[!is.na(x)]))
  if(procedure!='Matrix'){
    if(procedure=='DEA' | procedure=='DEA_count'){
      runtime=DEA_run_times
    }else if(procedure=='Imput'){
      runtime=Imput_run_times
    }else if(procedure=='normalization'){
      runtime=norm_run_times
    }
    idx<-match(heatmap_rank$Method, runtime$Method)
    heatmap_rank$time<-runtime$time[idx]
    label_heatmap_rank$time<-NA
  }
  
  ranks[[i]]=rank
  heatmap_ranks[[i]]<-heatmap_rank
  label_heatmap_ranks[[i]]<-label_heatmap_rank
}


library(openxlsx)
wb6 <- createWorkbook()
for (i in 1:length(procedures)) {
  ranks_procedure<-heatmap_ranks[[i]]
  ranks_procedure_label<-label_heatmap_ranks[[i]]
  
  sheet=procedures[i]
  
  addWorksheet(wb6,sheet)
  writeData(wb6, sheet, ranks_procedure, startRow = 1, startCol = 1, colNames = TRUE, rowNames = TRUE)
  saveWorkbook(wb6, file = paste0(save_fold, 'Heatmap.xlsx'), overwrite = TRUE)
  sheet_label=paste0(procedures[i], '_lab')
  addWorksheet(wb6,sheet_label)
  writeData(wb6, sheet_label, ranks_procedure_label, startRow = 1, startCol = 1, colNames = TRUE, rowNames = TRUE)
  saveWorkbook(wb6, file = paste0(save_fold, 'Heatmap.xlsx'), overwrite = TRUE)
  
}

library(ComplexHeatmap)
library(cowplot)
library(ggpubr)

heatmap_dt<-heatmap_ranks[[5]]
heatmap_dt<-rbind(heatmap_dt, c('reporter', rep(NA, length(heatmap_dt[1,])-1)))
heatmap_dt$MQ_TMT<-c(rep(NA, length(heatmap_dt$Method)-1), 1)

heatmap_lab<-label_heatmap_ranks[[5]]
heatmap_lab<-rbind(heatmap_lab, c('reporter intensity', rep(NA, length(heatmap_lab[1,])-1)))
heatmap_lab$MQ_TMT<-c(rep(NA, length(heatmap_lab$Method)-1), 1)

heatmap_dt<-heatmap_dt[c(1,2,4,3,8,5,6,7,9),]
heatmap_lab<-heatmap_lab[c(1,2,4,3,8,5,6,7,9),]
heatmap_dt<-heatmap_dt[,c(1:6,8,7)]
heatmap_lab<-heatmap_lab[,c(1:6,8,7)]

row.names(heatmap_dt)<-heatmap_dt$Method
row.names(heatmap_lab)<-heatmap_lab$Method
heatmap_dt<-heatmap_dt[,2:(length(heatmap_dt[1,])-1)]
heatmap_lab<-heatmap_lab[,2:(length(heatmap_lab[1,])-1)]

heatmap_dt$FG_DDA<-as.numeric(heatmap_dt$FG_DDA)
heatmap_dt$MQ_DDA<-as.numeric(heatmap_dt$MQ_DDA)
heatmap_dt$FG_TMT<-as.numeric(heatmap_dt$FG_TMT)
heatmap_dt$DIANN_DIA<-as.numeric(heatmap_dt$DIANN_DIA)
heatmap_dt$spt_DIA<-as.numeric(heatmap_dt$spt_DIA)
heatmap_dt$MQ_TMT<-as.numeric(heatmap_dt$MQ_TMT)

p9<-Heatmap(as.matrix(heatmap_dt), name = "rank", km = 0, cluster_rows = FALSE, 
            cluster_columns = FALSE,
            show_row_names = TRUE, row_names_side = 'right', show_column_names = TRUE, 
            column_names_side = 'bottom', heatmap_width=unit(130, "mm"),heatmap_height=unit(80, "mm"), row_names_centered = TRUE,
            column_names_rot=-20, column_names_centered = TRUE,
            col=circlize::colorRamp2(c(1,as.integer(max(as.numeric(heatmap_lab$FG_DDA[!is.na(heatmap_lab$FG_DDA)]))/2), max(as.numeric(heatmap_lab$FG_DDA[!is.na(heatmap_lab$FG_DDA)]))), c("#299D8F", 'white',"#D87659")),
            cell_fun = function(j, i, x, y, w, h, col) {
              grid.text(heatmap_lab[i, j], x, y) }, border = 'black',row_title = "Expression Matrix type")
p9

write.table(heatmap_dt, paste0(save_fold, 'Figure4B_1.csv'), col.names = T, row.names = T)
write.table(heatmap_lab, paste0(save_fold, 'Figure4B_2.csv'), col.names = T, row.names = T)


###################################################################################
### Figure 4c: comparing normalization

## normalization
heatmap_dt<-heatmap_ranks[[4]]
heatmap_lab<-label_heatmap_ranks[[4]]
ha = data.frame(Speed=as.numeric(heatmap_dt$time))
row.names(ha)<-heatmap_dt$Method
row.names(heatmap_dt)<-heatmap_dt$Method
row.names(heatmap_lab)<-heatmap_lab$Method
heatmap_dt<-heatmap_dt[,2:(length(heatmap_dt[1,])-1)]
heatmap_lab<-heatmap_lab[,2:(length(heatmap_lab[1,])-1)]

row.names(heatmap_dt)[which(row.names(heatmap_dt)=='blank')]='none'

p10<-Heatmap(as.matrix(heatmap_dt), name = "rank", km = 0, cluster_rows = FALSE, 
             cluster_columns = FALSE,
             show_row_names = TRUE, row_names_side = 'right', show_column_names = TRUE, 
             column_names_side = 'bottom', heatmap_width=unit(130, "mm"),heatmap_height=unit(80, "mm"), row_names_centered = TRUE,
             column_names_rot=-20, column_names_centered = TRUE,
             col=circlize::colorRamp2(c(1,as.integer(max(as.numeric(heatmap_lab$FG_DDA))/2), max(as.numeric(heatmap_lab$FG_DDA))), c("#299D8F", 'white',"#D87659")),
             cell_fun = function(j, i, x, y, w, h, col) {
               grid.text(heatmap_lab[i, j], x, y) }, border = 'black',heatmap_legend_param =list(legend_direction="horizontal"))
p10

p10_1<-Heatmap(as.matrix(ha), name = "time(s)", km = 0, cluster_rows = FALSE, 
               cluster_columns = FALSE,
               show_row_names = FALSE, show_column_names = TRUE, 
               column_names_side = 'bottom', heatmap_width=unit(16, "mm"),heatmap_height=unit(70, "mm"), row_names_centered = TRUE,
               column_names_rot=-20, column_names_centered = TRUE,
               col=circlize::colorRamp2(c(1,as.integer(max(as.numeric(ha$Speed))/2), max(as.numeric(ha$Speed))), c("#299D8F", 'white',"#D87659")),
               border = 'black',row_title = "Normalization methods",heatmap_legend_param =list(legend_direction="horizontal"))
p10_1

write.table(heatmap_dt, paste0(save_fold, 'Figure4C_1.csv'), col.names = T, row.names = T)
write.table(ha, paste0(save_fold, 'Figure4C_2.csv'), col.names = T, row.names = T)
write.table(heatmap_lab, paste0(save_fold, 'Figure4C_3.csv'), col.names = T, row.names = T)


###################################################################################
### Figure 4d: comparing imputation

# Imputation
heatmap_dt<-heatmap_ranks[[3]]
heatmap_lab<-label_heatmap_ranks[[3]]
ha = data.frame(Speed=as.numeric(heatmap_dt$time))
row.names(ha)<-heatmap_dt$Method
row.names(heatmap_dt)<-heatmap_dt$Method
row.names(heatmap_lab)<-heatmap_lab$Method
heatmap_dt<-heatmap_dt[,2:(length(heatmap_dt[1,])-1)]
heatmap_lab<-heatmap_lab[,2:(length(heatmap_lab[1,])-1)]


p11<-Heatmap(as.matrix(heatmap_dt), name = "rank", km = 0, cluster_rows = FALSE, 
             cluster_columns = FALSE,
             show_row_names = TRUE, row_names_side = 'right', show_column_names = TRUE, 
             column_names_side = 'bottom', heatmap_width=unit(130, "mm"),heatmap_height=unit(80, "mm"), row_names_centered = TRUE,
             column_names_rot=-20, column_names_centered = TRUE,
             col=circlize::colorRamp2(c(1,as.integer(max(as.numeric(heatmap_lab$FG_DDA))/2), max(as.numeric(heatmap_lab$FG_DDA))), c("#299D8F", 'white',"#D87659")),
             cell_fun = function(j, i, x, y, w, h, col) {
               grid.text(heatmap_lab[i, j], x, y) }, border = 'black', heatmap_legend_param =list(legend_direction="horizontal"))

p11

p11_1<-Heatmap(as.matrix(ha), name = "time(s)", km = 0, cluster_rows = FALSE, 
               cluster_columns = FALSE,
               show_row_names = FALSE, show_column_names = TRUE, 
               column_names_side = 'bottom', heatmap_width=unit(18, "mm"),heatmap_height=unit(75, "mm"), row_names_centered = TRUE,
               column_names_rot=-20, column_names_centered = TRUE,
               col=circlize::colorRamp2(c(1,as.integer(max(as.numeric(ha$Speed))/2), max(as.numeric(ha$Speed))), c("#299D8F", 'white',"#D87659")),
               border = 'black',row_title = "Imputation algorithms",heatmap_legend_param =list(legend_direction="horizontal"))
p11_1


write.table(heatmap_dt, paste0(save_fold, 'Figure4D_1.csv'), col.names = T, row.names = T)
write.table(ha, paste0(save_fold, 'Figure4D_2.csv'), col.names = T, row.names = T)
write.table(heatmap_lab, paste0(save_fold, 'Figure4D_3.csv'), col.names = T, row.names = T)
###################################################################################
### Figure 4e: comparing DEA

# DEA
heatmap_dt<-heatmap_ranks[[1]]
heatmap_lab<-label_heatmap_ranks[[1]]
ha = data.frame(Speed=as.numeric(heatmap_dt$time))
row.names(ha)<-heatmap_dt$Method
row.names(heatmap_dt)<-heatmap_dt$Method
row.names(heatmap_lab)<-heatmap_lab$Method
heatmap_dt<-heatmap_dt[,2:(length(heatmap_dt[1,])-1)]
heatmap_lab<-heatmap_lab[,2:(length(heatmap_lab[1,])-1)]

p12<-Heatmap(as.matrix(heatmap_dt), name = "rank", km = 0, cluster_rows = FALSE, 
             cluster_columns = FALSE,#row_order =as.matrix(data_fragpipe$workflow), #ranking_frag,
             show_row_names = TRUE, row_names_side = 'right', show_column_names = TRUE, 
             column_names_side = 'bottom', heatmap_width=unit(130, "mm"),heatmap_height=unit(50, "mm"), row_names_centered = TRUE,
             column_names_rot=-20, column_names_centered = TRUE,
             col=circlize::colorRamp2(c(1,as.integer(max(as.numeric(heatmap_lab$FG_DDA))/2), max(as.numeric(heatmap_lab$FG_DDA))), c("#299D8F", 'white',"#D87659")),
             cell_fun = function(j, i, x, y, w, h, col) {
               grid.text(heatmap_lab[i, j], x, y) }, border = 'black',heatmap_legend_param =list(legend_direction="horizontal"))
p12

p12_1<-Heatmap(as.matrix(ha), name = "time(s)", km = 0, cluster_rows = FALSE, 
               cluster_columns = FALSE,#row_order =as.matrix(data_fragpipe$workflow), #ranking_frag,
               show_row_names = FALSE, show_column_names = TRUE, 
               column_names_side = 'bottom', heatmap_width=unit(18, "mm"),heatmap_height=unit(45, "mm"), row_names_centered = TRUE,
               column_names_rot=-20, column_names_centered = TRUE,
               col=circlize::colorRamp2(c(1,as.integer(max(as.numeric(DEA_run_times$time))/2), max(as.numeric(DEA_run_times$time))), c("#299D8F", 'white',"#D87659")),
               border = 'black',row_title = "DEA with intensity",heatmap_legend_param =list(legend_direction="horizontal"))
p12_1

write.table(heatmap_dt, paste0(save_fold, 'Figure4E_1.csv'), col.names = T, row.names = T)
write.table(ha, paste0(save_fold, 'Figure4E_2.csv'), col.names = T, row.names = T)
write.table(heatmap_lab, paste0(save_fold, 'Figure4E_3.csv'), col.names = T, row.names = T)

# DEA count
heatmap_dt<-heatmap_ranks[[2]]
heatmap_lab<-label_heatmap_ranks[[2]]
ha = data.frame(Speed=as.numeric(heatmap_dt$time))
row.names(ha)<-heatmap_dt$Method
row.names(heatmap_dt)<-heatmap_dt$Method
row.names(heatmap_lab)<-heatmap_lab$Method
heatmap_dt<-heatmap_dt[,2:(length(heatmap_dt[1,])-2)]
heatmap_lab<-heatmap_lab[,2:(length(heatmap_lab[1,])-2)]

p12_2<-Heatmap(as.matrix(heatmap_dt), name = "rank", km = 0, cluster_rows = FALSE, 
               cluster_columns = FALSE,
               show_row_names = TRUE, row_names_side = 'right', show_column_names = TRUE, 
               column_names_side = 'top', heatmap_width=unit(60, "mm"),heatmap_height=unit(25, "mm"), row_names_centered = TRUE,
               column_names_rot=-20, column_names_centered = TRUE,
               col=circlize::colorRamp2(c(1,as.integer(max(as.numeric(heatmap_ranks[[1]]$FG_DDA))/2), max(as.numeric(heatmap_ranks[[1]]$FG_DDA))), c("#299D8F", 'white',"#D87659")),
               cell_fun = function(j, i, x, y, w, h, col) {
                 grid.text(heatmap_lab[i, j], x, y) }, border = 'black',heatmap_legend_param =list(legend_direction="horizontal"))
p12_2

p12_3<-Heatmap(as.matrix(ha), name = "time(s)", km = 0, cluster_rows = FALSE, 
               cluster_columns = FALSE,
               show_row_names = FALSE, show_column_names = FALSE, 
               column_names_side = 'bottom', heatmap_width=unit(18, "mm"),heatmap_height=unit(15, "mm"), row_names_centered = TRUE,
               column_names_rot=-20, column_names_centered = TRUE,
               col=circlize::colorRamp2(c(1,as.integer(max(as.numeric(DEA_run_times$time))/2), max(as.numeric(DEA_run_times$time))), c("#299D8F", 'white',"#D87659")),
               border = 'black',row_title = "count",heatmap_legend_param =list(legend_direction="horizontal"))
p12_3


write.table(heatmap_dt, paste0(save_fold, 'Figure4E_4.csv'), col.names = T, row.names = T)
write.table(ha, paste0(save_fold, 'Figure4E_5.csv'), col.names = T, row.names = T)
write.table(heatmap_lab, paste0(save_fold, 'Figure4E_6.csv'), col.names = T, row.names = T)

