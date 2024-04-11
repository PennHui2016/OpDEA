res_root<-'benchmark_res/'
for (platform in c('spt')) { 

acq<-'DIA'

if (acq=='DIA' & platform=='DIANN'){
  dataset_info_file<-'data/dataset_info/DIA_DIANN.txt'
  intens_ty = c('top1', 'top3', 'LFQ', 'dlfq')
}else if(acq=='DIA' & platform=='spt'){
  dataset_info_file<-'data/dataset_info/DIA_spt.txt'
  intens_ty = c('top1', 'top3', 'LFQ', 'dlfq')
}

imputs<-c('MLE', 'missForest', 'blank', 'bpca', 'knn', 'QRILC', 'MinDet', 'MinProb', 'min', 'zero', 'nbavg', 'mice',
          "Impseq", "Impseqrob", "GMS", 'SeqKNN')
norms<-c('blank', "sum", "max", "center.mean", "center.median", "div.mean", "div.median", "quantiles",
         "quantiles.robust", "vsn", 'lossf', 'TIC', "Rlr", 'MBQN')

deas<-c('limma', 'ttest', 'DEP', 'ANOVA', 'SAM', 'ROTS', 'proDA')
dea_c<-c('edgeR', 'plgem', 'beta_binomial')
count_ty<-c('count')
avaliable_wfs<-vector()
for (dea in deas) {
  for (inten_ty in intens_ty) {
    for (imput in imputs) {
      for (norm in norms) {
        wf<-c(paste0(dea,'|',platform,'|',inten_ty,'|',imput, '|', norm), dea ,platform,inten_ty,imput, norm)
        avaliable_wfs<-rbind(avaliable_wfs, wf)
      }
    }
  }

}


for (dea in c('MSstats')) {
  for (inten_ty in c('all', 'top3')) {
    for (imput in c('TRUE', 'FALSE')) {
      for (norm in c('equalizeMedians', 'quantile', 'FALSE')) {
        wf<-c(paste0(dea,'|',platform,'|',inten_ty,'|',imput, '|', norm), dea ,platform,inten_ty,imput, norm)
        avaliable_wfs<-rbind(avaliable_wfs, wf)
      }
    }
  }
  
}

dataset_info<-read.table(dataset_info_file, header = T, sep = '\t')

datasets<-dataset_info$dataset
contrss<-dataset_info$test_contrs

#dt = datasets[1]

read_res_file<-function(filename){
  demo=file(filename,open="r")
  n=1
  ress<-vector()
  while ( TRUE ) {
    line = readLines(demo, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    res<-strsplit(line, ',', fixed = T)[[1]]
    if(length(grep('bug', res))==0){
      ress<-rbind(ress, res)
    }
    #cat(n,line,"\n")
    n=n+1
  }
  close(demo)
  return(as.data.frame(ress))
}


for (j in c(c(1:7))) {
  
  dt<-datasets[j]
  contrs<-contrss[j]
  cls001<-read_res_file(paste0(res_root, acq, '/', platform, dt, '_',
                            acq, '_', platform,'_cls_001_2.csv'))
  cls005<-read_res_file(paste0(res_root, acq, '/', platform, dt, '_',
                            acq, '_', platform,'_cls_005_2.csv'))
  paucs<-read_res_file(paste0(res_root, acq, '/', platform, dt, '_',
                            acq, '_', platform,'_paucs2.csv'))


  wf_strs_cls001<-paste0(cls001$V1,'|',cls001$V2,'|',cls001$V3,'|',cls001$V4,'|',cls001$V5)
  wf_strs_cls005<-paste0(cls005$V1,'|',cls005$V2,'|',cls005$V3,'|',cls005$V4,'|',cls005$V5)
  wf_strs_pauc<-paste0(paucs$V1,'|',paucs$V2,'|',paucs$V3,'|',paucs$V4,'|',paucs$V5)

cls_metrics<-c('Acc', 'Prec', 'Rec', 'F1', 'F1w', 'Mcc', 'nMcc', 'geomean', 'tn', 'fp', 'fn', 'tp')
cls_mt_all<-paste0(cls_metrics,'_all')
uni_contr<-sort(unique(strsplit(gsub('condition', '',contrs), ';',fixed = T)[[1]]))
for (contr in uni_contr) {
  cls_mt_all<-c(cls_mt_all, paste0(cls_metrics, '_',contr))
}
colname_cls<-c('DEA','Platform','Matrix','Imput', 'normalization')

pauc001_metrics<-c('pauc001_all',paste0('pauc001_', uni_contr))
pauc005_metrics<-c('pauc005_all',paste0('pauc005_', uni_contr))
pauc01_metrics<-c('pauc01_all',paste0('pauc01_', uni_contr))
pauc001_metrics_lfq<-c('pauc001_all_lfq',paste0('pauc001_', uni_contr, '_lfq'))
pauc005_metrics_lfq<-c('pauc005_all_lfq',paste0('pauc005_', uni_contr, '_lfq'))
pauc01_metrics_lfq<-c('pauc01_all_lfq',paste0('pauc01_', uni_contr, '_lfq'))

colname_cls_001<-c(colname_cls,'runTime', cls_mt_all)
colname_cls_005<-c(colname_cls, cls_mt_all)
colname_pauc<-c(colname_cls, pauc001_metrics, pauc005_metrics, pauc01_metrics,
                pauc001_metrics_lfq, pauc005_metrics_lfq, pauc01_metrics_lfq)

colnames(cls001)<-colname_cls_001
colnames(cls005)<-colname_cls_005
colnames(paucs)<-colname_pauc



cls_001_all<-vector()
cls_005_all<-vector()
pauc_all<-vector()
for (i in 1:length(avaliable_wfs[,1])) {
  print(paste0(as.character(i),'|',as.character(j)))
  if(i==6272){
    a=0
  }
  idx<-which(wf_strs_cls001==avaliable_wfs[i, 1])
  if(length(idx)==1){
    cls_001_all<-rbind(cls_001_all, cls001[idx,])
  }else if(length(idx)>1){
    cls_001_all<-rbind(cls_001_all, cls001[idx[length(idx)],])
  }else{
    cls_001_all<-rbind(cls_001_all, c(avaliable_wfs[i,c(2:6)],rep(0, length(cls001[1,])-5)))
  }
  colnames(cls_001_all)<-colname_cls_001

  idx<-which(wf_strs_cls005==avaliable_wfs[i, 1])
  if(length(idx)==1){
    cls_005_all<-rbind(cls_005_all, cls005[idx,])
  }else if(length(idx)>1){
    cls_005_all<-rbind(cls_005_all, cls005[idx[length(idx)],])
  }else{
    cls_005_all<-rbind(cls_005_all, c(avaliable_wfs[i,c(2:6)],rep(0, length(cls005[1,])-5)))
  }
  colnames(cls_005_all)<-colname_cls_005

  idx<-which(wf_strs_pauc==avaliable_wfs[i, 1])
  if(length(idx)==1){
    pauc_all<-rbind(pauc_all, paucs[idx,])
  }else if(length(idx)>1){
    pauc_all<-rbind(pauc_all, paucs[idx[length(idx)],])
  }else{
    pauc_all<-rbind(pauc_all, c(avaliable_wfs[i,c(2:6)],rep(0, length(paucs[1,])-5)))
  }
  colnames(pauc_all)<-colname_pauc
}
write.table(cls_001_all, paste0(res_root, 'clear_res/', dt,'_', platform, '_cls_001.csv'), sep = ',', col.names = T, row.names = F)
write.table(cls_005_all, paste0(res_root, 'clear_res/', dt,'_', platform,'_cls_005.csv'), sep = ',', col.names = T, row.names = F)
write.table(pauc_all, paste0(res_root, 'clear_res/', dt,'_', platform,'_pAUCs.csv'), sep = ',', col.names = T, row.names = F)
}
}


