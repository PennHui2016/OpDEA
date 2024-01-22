res_root<-'benchmark_res/'



rank_metrics<-function(metrics){
  metrics$mean<-apply(metrics, 1, function(x) mean(as.numeric(x[c(7:length(x))])))
  metrics$median<-apply(metrics, 1, function(x) median(as.numeric(x[c(7:length(x))])))
  metrics$rank_mean[order(-as.numeric(metrics$mean))]<-c(1:length(metrics$mean))
  metrics$rank_median[order(-as.numeric(metrics$median))]<-c(1:length(metrics$median))
  metrics$workflow<-paste0(metrics$ensemble_method,'||',metrics$cbn)
  metrics<-metrics[,c(length(metrics[1,]), 1:(length(metrics[1,])-1))]
  ranks<-cbind(metrics$mean,metrics$median, metrics$rank_mean, metrics$rank_median)
  return(list(all=metrics, rank=ranks))
}

data_fold<-'benchmark_res/'

for (i in c(c(1:6))) {
  if(i==1){
    platform='FragPipe'
    acq='DDA'
  }else if(i==2){
    platform='Maxquant'
    acq='DDA'
  }else if(i==3){
    platform='FragPipe'
    acq='TMT'
  }else if(i==4){
    platform='Maxquant'
    acq='TMT'
  }else if(i==5){
    platform='DIANN'
    acq='DIA'
  }else if(i==6){
    platform='spt'
    acq='DIA'
  }
  if(i==4){
    ensem_types<-c('topk')
  }else{
    ensem_types<-c('mv', 'topk')
  }
  for (ensem_type in ensem_types) {
    if (acq=='TMT' & platform=='FragPipe'){
      dataset_info_file<-'data/dataset_info/TMT_frag.txt'
      intens_ty = c('abd', 'ratio', 'phi')
    }else if(acq=='TMT' & platform=='Maxquant'){
      dataset_info_file<-'data/dataset_info/TMT_mq.txt'
      intens_ty = c('intensity')
    }else if (acq=='DDA' & platform=='FragPipe'){
      dataset_info_file<-'data/dataset_info/DDA_frag.txt'
      intens_ty = c('top0', 'top3', 'LFQ', 'dlfq')
    }else if(acq=='DDA' & platform=='Maxquant'){
      dataset_info_file<-'data/dataset_info/DDA_mq.txt'
      intens_ty = c('top0', 'top3', 'LFQ', 'dlfq')
    }else if (acq=='DIA' & platform=='DIANN'){
      dataset_info_file<-'data/dataset_info/DIA_DIANN.txt'
      intens_ty = c('top1', 'top3', 'LFQ', 'dlfq')
    }else if(acq=='DIA' & platform=='spt'){
      dataset_info_file<-'data/dataset_info/DIA_spt.txt'
      intens_ty = c('top1', 'top3', 'LFQ', 'dlfq')
    }
    
    dataset_info<-read.table(dataset_info_file, header = T, sep = '\t')
    if(acq=='TMT'){
      dataset_info<-dataset_info[c(1,2,5,6,7),]
    }
    
    datasets<-dataset_info$dataset
    contrss<-dataset_info$test_contrs
    res_root<-data_fold
    
    pAUC001s<-vector()
    pAUC005s<-vector()
    pAUC01s<-vector()
    nMCCs001<-vector()
    Gmeans001<-vector()
    nMCCs005<-vector()
    Gmeans005<-vector()
    flag=0
    
    coln_pauc001<-c()
    coln_pauc005<-c()
    coln_pauc01<-c()
    coln_nmcc001<-c()
    coln_gmean001<-c()
    coln_nmcc005<-c()
    coln_gmean005<-c()
    avg_ranks<-vector()
    
    for (j in c(1:length(datasets))){
      dt<-datasets[j]
      
      res_file_cls001<-paste0(res_root, 'clear_res/', dt,'_', platform, '_cls_001_ensemble_',ensem_type,'.csv')
      res_file_cls005<-paste0(res_root, 'clear_res/', dt,'_', platform, '_cls_005_ensemble_',ensem_type,'.csv')
      res_file_pAUC<-paste0(res_root, 'clear_res/', dt,'_', platform,'_pAUCs_ensemble_',ensem_type,'.csv')
      
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
    
    pAUC001s<-as.data.frame(cbind(res_cls001[,c(1:6)], pAUC001s))
    pAUC005s<-as.data.frame(cbind(res_cls001[,c(1:6)], pAUC005s))
    pAUC01s<-as.data.frame(cbind(res_cls001[,c(1:6)], pAUC01s))
    nMCCs001<-as.data.frame(cbind(res_cls001[,c(1:6)], nMCCs001))
    Gmeans001<-as.data.frame(cbind(res_cls001[,c(1:6)], Gmeans001))
    nMCCs005<-as.data.frame(cbind(res_cls001[,c(1:6)], nMCCs005))
    Gmeans005<-as.data.frame(cbind(res_cls001[,c(1:6)], Gmeans005))
    
    pAUC001s<-rank_metrics(pAUC001s)
    pAUC005s<-rank_metrics(pAUC005s)
    pAUC01s<-rank_metrics(pAUC01s)
    nMCCs001<-rank_metrics(nMCCs001)
    Gmeans001<-rank_metrics(Gmeans001)
    nMCCs005<-rank_metrics(nMCCs005)
    Gmeans005<-rank_metrics(Gmeans005)
    
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
    ranks_all<-as.data.frame(ranks_all)

    ranks_all$workflow<-paste0(ranks_all$ensemble_method,"||",ranks_all$cbn)
    # change order
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
      if(length(strs)>1 & length(grep('rank',strs))==0 & length(grep('ensemble',strs))==0){
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
    saveWorkbook(wb, file = paste0(res_root, 'combined_res/metrics_', platform, '_', acq, '_all_ensemble_',ensem_type,'.xlsx'), overwrite = TRUE)
    
    openxlsx::write.xlsx(as.data.frame(ranks_all), file=paste0(res_root, 'combined_res/ranks_all_', platform, '_', acq, '_ensemble_',ensem_type,'.xlsx'), sheetName = paste0('ranking_all'),append = T, rowNames = F)
    
  }
}



