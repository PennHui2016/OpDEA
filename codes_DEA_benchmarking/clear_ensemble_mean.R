############# clear ensemble inference results

data_fold<-'benchmark_res/ensemble/'

for (i in c(c(1:6))){
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
    for (j in c(1:length(datasets))){
      dt<-datasets[j]
      contrs<-contrss[j]
      cls001<-read.table(paste0(res_root, acq, '/', platform, '/', ensem_type,'/', dt, '_',
                                acq, '_', platform,'_cls_001_2.csv'), sep = ',')
      cls005<-read.table(paste0(res_root, acq, '/', platform, '/', ensem_type,'/', dt, '_',
                                acq, '_', platform,'_cls_005_2.csv'), sep = ',')
      paucs<-read.table(paste0(res_root, acq, '/', platform, '/', ensem_type,'/',dt, '_',
                               acq, '_', platform,'_paucs2.csv'), sep = ',')
      
      
      wf_strs_cls001<-paste0(cls001$V1,'|',cls001$V2,'|',cls001$V3,'|',cls001$V4,'|',cls001$V5)
      wf_strs_cls005<-paste0(cls005$V1,'|',cls005$V2,'|',cls005$V3,'|',cls005$V4,'|',cls005$V5)
      wf_strs_pauc<-paste0(paucs$V1,'|',paucs$V2,'|',paucs$V3,'|',paucs$V4,'|',paucs$V5)
      
      cls_metrics<-c('Acc', 'Prec', 'Rec', 'F1', 'F1w', 'Mcc', 'nMcc', 'geomean', 'tn', 'fp', 'fn', 'tp')
      cls_mt_all<-paste0(cls_metrics,'_all')
      uni_contr<-sort(unique(strsplit(gsub('condition', '',contrs), ';',fixed = T)[[1]]))
      for (contr in uni_contr) {
        cls_mt_all<-c(cls_mt_all, paste0(cls_metrics, '_',contr))
      }
      colname_cls<-c('ensemble_method','Platform','method','operation', 'cbn')
      
      pauc001_metrics<-c('pauc001_all',paste0('pauc001_', uni_contr))
      pauc005_metrics<-c('pauc005_all',paste0('pauc005_', uni_contr))
      pauc01_metrics<-c('pauc01_all',paste0('pauc01_', uni_contr))
      pauc001_metrics_lfq<-c('pauc001_all_lfq',paste0('pauc001_', uni_contr, '_lfq'))
      pauc005_metrics_lfq<-c('pauc005_all_lfq',paste0('pauc005_', uni_contr, '_lfq'))
      pauc01_metrics_lfq<-c('pauc01_all_lfq',paste0('pauc01_', uni_contr, '_lfq'))
      
      colname_cls_001<-c(colname_cls,'ensemble_info', cls_mt_all)
      colname_cls_005<-c(colname_cls, cls_mt_all)
      colname_pauc<-c(colname_cls, pauc001_metrics, pauc005_metrics, pauc01_metrics,
                      pauc001_metrics_lfq, pauc005_metrics_lfq, pauc01_metrics_lfq)
      
      colnames(cls001)<-colname_cls_001
      colnames(cls005)<-colname_cls_005
      colnames(paucs)<-colname_pauc
      
      ens_infos<-c()
      if(ensem_type=='mv'){
        for (ens_info in cls001$ensemble_info) {
          ens_infos<-c(ens_infos,strsplit(ens_info,'||', fixed = T)[[1]][1])
        }
      }else if(ensem_type=='topk'){
        for (ens_info in paucs$cbn) {
          ens_infos<-c(ens_infos,paste0('top',as.character(ens_info)))
        }
      }
      
      cls001$cbn<-ens_infos
      cls005$cbn<-ens_infos
      paucs$cbn<-ens_infos
      
      print(paste(as.character(i), '|', as.character(j), '|', ensem_type))
      write.table(cls001, paste0('benchmark_res/', 'clear_res/', dt,'_', platform, '_cls_001_ensemble_',ensem_type,'.csv'), sep = ',', col.names = T, row.names = F)
      write.table(cls005, paste0('benchmark_res/', 'clear_res/', dt,'_', platform,'_cls_005_ensemble_',ensem_type,'.csv'), sep = ',', col.names = T, row.names = F)
      write.table(paucs, paste0('benchmark_res/', 'clear_res/', dt,'_', platform,'_pAUCs_ensemble_',ensem_type,'.csv'), sep = ',', col.names = T, row.names = F)
    }
  }
}
#ensem_type = 'mv' #'topk'


