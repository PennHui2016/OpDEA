res_root<-'benchmark_res/'

rank_metrics<-function(metrics, st){
  means<-apply(metrics, 1, function(x) mean(as.numeric(x[c(st:length(x))])))
  medians<-apply(metrics, 1, function(x) median(as.numeric(x[c(st:length(x))])))
  metrics$mean<-means
  metrics$median<-medians
  metrics$rank_mean[order(-as.numeric(metrics$mean))]<-c(1:length(metrics$mean))
  metrics$rank_median[order(-as.numeric(metrics$median))]<-c(1:length(metrics$median))
  metrics$workflow<-paste0(metrics$DEA,'|',metrics$Platform,'|',metrics$Matrix,'|',metrics$Imput,'|',metrics$normalization)
  metrics<-metrics[,c(length(metrics[1,]), 1:(length(metrics[1,])-1))]
  ranks<-cbind(metrics$mean,metrics$median, metrics$rank_mean, metrics$rank_median)
  return(list(all=metrics, rank=ranks))
}

for (platform in c('DIANN', 'spt')) {
    acq<-'DIA'
  if (acq=='DIA' & platform=='DIANN'){
    dataset_info_file<-'data/dataset_info/DIA_DIANN.txt'
  }else if(acq=='DIA' & platform=='spt'){
    dataset_info_file<-'data/dataset_info/DIA_spt.txt'
  }
  dataset_info<-read.table(dataset_info_file, header = T, sep = '\t')
  datasets<-dataset_info$dataset
  contrss<-dataset_info$test_contrs
  
  pAUC001s<-vector()
  pAUC005s<-vector()
  pAUC01s<-vector()
  nMCCs001<-vector()
  Gmeans001<-vector()
  nMCCs005<-vector()
  Gmeans005<-vector()
  
  runtimes<-vector()
  flag=0
  
  avg_ranks<-vector()
  coln_pauc001<-c()
  coln_pauc005<-c()
  coln_pauc01<-c()
  coln_nmcc001<-c()
  coln_gmean001<-c()
  coln_nmcc005<-c()
  coln_gmean005<-c()
  for (j in 1:length(datasets)) {

    dt<-datasets[j]
    
    res_file_cls001<-paste0(res_root, 'clear_res/', dt,'_', platform, '_cls_001.csv')
    res_file_cls005<-paste0(res_root, 'clear_res/', dt,'_', platform, '_cls_005.csv')
    res_file_pAUC<-paste0(res_root, 'clear_res/', dt,'_', platform,'_pAUCs.csv')
    
    res_cls001<-read.table(res_file_cls001, sep = ',', header = T)
    res_cls005<-read.table(res_file_cls005, sep = ',', header = T)
    res_pAUC<-read.table(res_file_pAUC, sep = ',', header = T)
    
    colnames(res_cls001)[7:length(res_cls001[1,])]<-paste0(colnames(res_cls001)[7:length(res_cls001[1,])],'_', dt)
    colnames(res_cls005)[6:length(res_cls005[1,])]<-paste0(colnames(res_cls005)[6:length(res_cls005[1,])],'_', dt)
    colnames(res_pAUC)[6:length(res_pAUC[1,])]<-paste0(colnames(res_pAUC)[6:length(res_pAUC[1,])],'_', dt)
    
    nmcc001_idx<-grep('nMcc', colnames(res_cls001))
    nmcc001_all_idx<-grep('all', colnames(res_cls001))
    nmcc001_idx<-setdiff(nmcc001_idx, nmcc001_all_idx)
    
    nmcc005_idx<-grep('nMcc', colnames(res_cls005))
    nmcc005_all_idx<-grep('all', colnames(res_cls005))
    nmcc005_idx<-setdiff(nmcc005_idx, nmcc005_all_idx)
    
    gmean001_idx<-grep('geomean', colnames(res_cls001))
    gmean001_idx<-setdiff(gmean001_idx, nmcc001_all_idx)
    
    gmean005_idx<-grep('geomean', colnames(res_cls005))
    gmean005_idx<-setdiff(gmean005_idx, nmcc005_all_idx)
    
    pauc001_idx<-grep('pauc001', colnames(res_pAUC))
    pauc_all<-grep('all', colnames(res_pAUC))
    pauc_lfq<-grep('lfq', colnames(res_pAUC))
    pauc001_idx<-setdiff(pauc001_idx, c(pauc_all,pauc_lfq))
    
    pauc005_idx<-grep('pauc005', colnames(res_pAUC))
    pauc005_idx<-setdiff(pauc005_idx, c(pauc_all,pauc_lfq))
    
    pauc01_idx<-grep('pauc01', colnames(res_pAUC))
    pauc01_idx<-setdiff(pauc01_idx, c(pauc_all,pauc_lfq))
    
    
    
    if (flag==0){
      pAUC001s<-res_pAUC[,pauc001_idx]
      pAUC005s<-res_pAUC[,pauc005_idx]
      pAUC01s<-res_pAUC[,pauc01_idx]
      nMCCs001<-res_cls001[,nmcc001_idx]
      Gmeans001<-res_cls001[,gmean001_idx]
      nMCCs005<-res_cls005[,nmcc005_idx]
      Gmeans005<-res_cls005[,gmean005_idx]
      coln_pauc001<-c(coln_pauc001, colnames(res_pAUC)[pauc001_idx])
      coln_pauc005<-c(coln_pauc005, colnames(res_pAUC)[pauc005_idx])
      coln_pauc01<-c(coln_pauc01, colnames(res_pAUC)[pauc01_idx])
      coln_nmcc001<-c(coln_nmcc001, colnames(res_cls001)[nmcc001_idx])
      coln_gmean001<-c(coln_gmean001, colnames(res_cls001)[gmean001_idx])
      coln_nmcc005<-c(coln_nmcc005, colnames(res_cls005)[nmcc005_idx])
      coln_gmean005<-c(coln_gmean005, colnames(res_cls005)[gmean005_idx])
      runtimes<-res_cls001[, 6]
      flag = flag+1
    }else{
    pAUC001s<-cbind(pAUC001s, res_pAUC[,pauc001_idx])
    pAUC005s<-cbind(pAUC005s, res_pAUC[,pauc005_idx])
    pAUC01s<-cbind(pAUC01s, res_pAUC[,pauc01_idx])
    nMCCs001<-cbind(nMCCs001, res_cls001[,nmcc001_idx])
    Gmeans001<-cbind(Gmeans001, res_cls001[,gmean001_idx])
    nMCCs005<-cbind(nMCCs005, res_cls005[,nmcc005_idx])
    Gmeans005<-cbind(Gmeans005, res_cls005[,gmean005_idx])
    runtimes<-cbind(runtimes, res_cls001[, 6])
    coln_pauc001<-c(coln_pauc001, colnames(res_pAUC)[pauc001_idx])
    coln_pauc005<-c(coln_pauc005, colnames(res_pAUC)[pauc005_idx])
    coln_pauc01<-c(coln_pauc01, colnames(res_pAUC)[pauc01_idx])
    coln_nmcc001<-c(coln_nmcc001, colnames(res_cls001)[nmcc001_idx])
    coln_gmean001<-c(coln_gmean001, colnames(res_cls001)[gmean001_idx])
    coln_nmcc005<-c(coln_nmcc005, colnames(res_cls005)[nmcc005_idx])
    coln_gmean005<-c(coln_gmean005, colnames(res_cls005)[gmean005_idx])
    }
    
  }
  
  colnames(pAUC001s)<-coln_pauc001
  colnames(pAUC005s)<-coln_pauc005
  colnames(pAUC01s)<-coln_pauc01
  colnames(nMCCs001)<-coln_nmcc001
  colnames(Gmeans001)<-coln_gmean001
  colnames(nMCCs005)<-coln_nmcc005
  colnames(Gmeans005)<-coln_gmean005

  pAUC001s<-as.data.frame(cbind(res_cls001[,c(1:5)], pAUC001s))
  pAUC005s<-as.data.frame(cbind(res_cls001[,c(1:5)], pAUC005s))
  pAUC01s<-as.data.frame(cbind(res_cls001[,c(1:5)], pAUC01s))
  nMCCs001<-as.data.frame(cbind(res_cls001[,c(1:5)], nMCCs001))
  Gmeans001<-as.data.frame(cbind(res_cls001[,c(1:5)], Gmeans001))
  nMCCs005<-as.data.frame(cbind(res_cls001[,c(1:5)], nMCCs005))
  Gmeans005<-as.data.frame(cbind(res_cls001[,c(1:5)], Gmeans005))
  
  pAUC001s<-rank_metrics(pAUC001s, 6)
  pAUC005s<-rank_metrics(pAUC005s, 6)
  pAUC01s<-rank_metrics(pAUC01s, 6)
  nMCCs001<-rank_metrics(nMCCs001, 6)
  Gmeans001<-rank_metrics(Gmeans001, 6)
  nMCCs005<-rank_metrics(nMCCs005, 6)
  Gmeans005<-rank_metrics(Gmeans005, 6)
  
  ranks_all<-cbind(pAUC001s$rank, pAUC005s$rank, pAUC01s$rank, nMCCs005$rank, Gmeans005$rank)
  colnames(ranks_all)<-c('mean_pauc001', 'median_pauc001', 'rank_mean_pauc001', 'rank_median_pauc001',
                         'mean_pauc005', 'median_pauc005', 'rank_mean_pauc005', 'rank_median_pauc005',
                         'mean_pauc01', 'median_pauc01', 'rank_mean_pauc01', 'rank_median_pauc01',
                         'mean_nmcc005', 'median_nmcc005', 'rank_mean_nmcc005', 'rank_median_nmcc005',
                         'mean_gmean005', 'median_gmean005', 'rank_mean_gmean005', 'rank_median_gmean005')
  ranks_all<-as.data.frame(ranks_all)
  ranks_all$avg_rank_mean<-rowMeans(cbind(ranks_all$rank_mean_pauc001,
                                          ranks_all$rank_mean_pauc005,
                                          ranks_all$rank_mean_pauc01,
                                          ranks_all$rank_mean_nmcc005,
                                          ranks_all$rank_mean_gmean005))
  ranks_all$avg_rank_median<-rowMeans(cbind(ranks_all$rank_median_pauc001,
                                          ranks_all$rank_median_pauc005,
                                          ranks_all$rank_median_pauc01,
                                          ranks_all$rank_median_nmcc005,
                                          ranks_all$rank_median_gmean005))
  
  ranks_all<-cbind(res_cls001[,c(1:6)], ranks_all)
  ranks_all$workflow<-paste0(ranks_all$DEA,'|',ranks_all$Platform,'|',ranks_all$Matrix,'|',ranks_all$Imput,'|',ranks_all$normalization)
  ranks_all<-ranks_all[,c(length(ranks_all[1,]), 1:(length(ranks_all[1,])-1))]

  library(openxlsx)
  
  sheets<-c('pAUC0.01','pAUC0.05','pAUC0.1','nMCC0.01','G-mean0.01','nMCC0.05','G-mean0.05')
  out_Data<-list()
  out_Data[[1]]<-as.data.frame(pAUC001s$all)
  out_Data[[2]]<-as.data.frame(pAUC005s$all)
  out_Data[[3]]<-as.data.frame(pAUC01s$all)
  out_Data[[4]]<-as.data.frame(nMCCs001$all)
  out_Data[[5]]<-as.data.frame(Gmeans001$all)
  out_Data[[6]]<-as.data.frame(nMCCs005$all)
  out_Data[[7]]<-as.data.frame(Gmeans005$all)
  wb <- createWorkbook()
  for (i in 1:length(sheets) ) {
    addWorksheet(wb, sheets[i])
  }
  for (i in 1:length(sheets) ) {
    writeData(wb, sheets[i], out_Data[[i]], startRow = 1, startCol = 1)
  }
  
  headers<-rbind(colnames(pAUC001s$all),colnames(pAUC001s$all))
  for (s in 1:length(headers[1,])) {
    strs<-strsplit(headers[1,s],'_',fixed = T)[[1]]
    if(length(strs)>1 & length(grep('rank',strs))==0){
      cons<-strsplit(strs[2],'.',fixed = T)[[1]]
      constr_name<-paste0('condition',cons[1],'-','condition',cons[2])
      dt_name<-strs[3]
      for (x in 4:length(strs)) {
        dt_name<-paste0(dt_name,'_',strs[x])
      }
      headers[1,s]<-dt_name
      headers[2,s]<-constr_name
    }
  }
  addWorksheet(wb, 'header')
  writeData(wb, 'header', as.data.frame(headers), startRow = 1, startCol = 1, colNames = FALSE)
  
  saveWorkbook(wb, file = paste0(res_root, 'combined_res/metrics_', platform, '_', acq,'.xlsx'), overwrite = TRUE)
  
  openxlsx::write.xlsx(as.data.frame(ranks_all), file=paste0(res_root, 'combined_res/ranks_all_', platform, '_', acq, '.xlsx'), sheetName = paste0('ranking_all'),append = T, rowNames = F)
  
  colnames(runtimes)<-dataset_info$dataset
  runtimes<-cbind(res_cls001[,c(1:5)], runtimes)
  write.table(as.data.frame(runtimes), file=paste0(res_root, 'combined_res/', platform, '_', acq, '_runtimes_all.csv'),sep = ',', row.names=FALSE, col.names = T)
}

