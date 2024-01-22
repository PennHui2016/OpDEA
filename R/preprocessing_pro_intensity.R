
preprocessing_raw<-function(mat, design, imput='MinProb',normal='lossf', log2=T){
  library(MSnbase)
  library(GMSimpute)
  library(SeqKnn)
  library(rrcovNA)
  library(NormalyzerDE)
  library(limma)

  set.seed(123)
  pro_file<-mat
  intens<-pro_file[, 3:length(colnames(pro_file))]
  intens[intens==0]<-NA
  row.names(intens) <-pro_file[,1]

  colnames(intens)<-gsub('LFQ.intensity.','',colnames(intens))
  colnames(intens)<-gsub('.Spectral.Count','',colnames(intens))

  colnames(intens)<-gsub('.MaxLFQ.Intensity','',colnames(intens))
  colnames(intens)<-gsub('Intensity.','',colnames(intens))
  colnames(intens)<-gsub('.Intensity','',colnames(intens))
  colnames(intens)<-gsub('Top3.','',colnames(intens))

  idx_retain<-which(apply(intens,1,function(x) sum(is.na(x)))<length(intens[1,])*0.8)
  intens=intens[idx_retain,]
  pro_file<-pro_file[idx_retain,]

  intens_n<-as.data.frame(apply(intens, 2, function(x) as.numeric(x)))
  row.names(intens_n)<-row.names(intens)
  colnames(intens_n)<-colnames(intens)
  intens<-intens_n

  if (log2==T){
    intens_tran<-log2(as.matrix(intens))
  }else{
    intens_tran = as.matrix(intens)
  }

  fd <- data.frame(row.names(intens))
  row.names(fd)<-row.names(intens)
  pd <- data.frame(colnames(intens))
  row.names(pd) = colnames(intens)
  intens_imp<-MSnSet(intens_tran, fd, pd)

  #### normalization
  if(normal!='' & normal!='None'){
    if(normal=='Rlr'){
      intens_norm = performGlobalRLRNormalization(intens_imp@assayData[["exprs"]], noLogTransform = T)
      intens_norm<-MSnSet(exprs = intens_norm)
    }else if(normal=='lossf'){
      intens_norm = normalizeCyclicLoess(intens_imp@assayData[["exprs"]], method = "fast")
      intens_norm<-MSnSet(exprs = intens_norm)
    }else if(normal=='MBQN'){
      library(MBQN)
      mtx <- as.matrix(intens_imp@assayData[["exprs"]])
      mtx.trqn <- mbqn(mtx, FUN = mean)
      row.names(mtx.trqn) <- row.names(intens_imp@assayData[["exprs"]])
      colnames(mtx.trqn) <- colnames(intens_imp@assayData[["exprs"]])
      intens_norm <- mtx.trqn
      intens_norm<-MSnSet(exprs = intens_norm)
    }else if(normal=='TIC'){
      exp_mat = intens_imp@assayData[["exprs"]]
      exp_mat[which(is.na(exp_mat))]=0
      tics = colSums(exp_mat)
      intens_norm=vector()
      for (i in 1:length(tics)) {
        intens_norm = cbind(intens_norm, exp_mat[,i]/tics[i])
      }
      intens_norm = intens_norm * median(tics)
      intens_norm[intens_norm==0]=NA
      colnames(intens_norm)=colnames(exp_mat)
      intens_norm<-MSnSet(exprs = intens_norm)
    }else{
      intens_norm<-normalise(intens_imp, normal)
    }
  }else{
    intens_norm<-intens_imp
  }

  exprs_norm = intens_norm@assayData[["exprs"]]
  exprs_norm[is.infinite(exprs_norm)]=NA
  intens_norm<-MSnSet(exprs = exprs_norm)

  #### imputation
  if(imput!='' & imput!='None'){
    if(imput %in% c('bpca', 'knn', 'QRILC', 'MLE', 'MinDet', 'MinProb', 'min', 'zero', 'mixed', 'nbavg', 'with', 'None')){
      if(imput=='bpca'){
        intens_tran = exprs(intens_norm)
        idx_bpca<-which(apply(intens_tran,1,function(x) sum(is.na(x)))<length(intens_tran[1,]))
        intens_imp<-intens_tran[idx_bpca,]
        intens_imp<-MSnSet(exprs = intens_imp)
      }else if(imput=='knn'){
        intens_tran = exprs(intens_norm)
        idx_knn<-which(apply(intens_tran,1,function(x) sum(is.na(x)))<length(intens_tran[1,])*0.8)
        intens_imp<-intens_tran[idx_knn,]
        intens_imp<-MSnSet(exprs = intens_imp)
      }else{
        intens_imp=intens_norm
      }
      intens_imp<-MSnbase::impute(intens_imp, imput)

    }else if(imput=='mice'){
      library(mice)
      intens_tran = exprs(intens_norm)
      intens_imp<-mice(intens_tran)
      completeData <- complete(intens_imp,2)
      completeData<-as.matrix(completeData)
      intens_imp<-MSnSet(exprs = as.matrix(completeData))
    }else if(imput=='missForest'){
      library(missForest)
      intens_tran = exprs(intens_norm)
      intens_imp<-missForest(intens_tran)
      intens_imp<-MSnSet(exprs = intens_imp$ximp)
    }else if(imput=='Impseq' | imput=='Impseqrob'){
      #library(mi)
      intens_tran = exprs(intens_norm)
      #intens_imp<-missing_data.frame(intens_tran)
      mat = matrix(NA, nrow=length(intens_tran[,1]),ncol=length(intens_tran[1,]))
      for (i in 1:length(intens_tran[1,])) {
        mat[,i]=as.numeric(intens_tran[,i])
      }
      if(imput=='Impseq'){
        mat<-impSeq(mat)
      }else if(imput=='Impseqrob'){
        mat<-impSeqRob(mat)$x
      }

      intens_imp<-intens_tran
      for (i in 1:length(intens_tran[1,])) {
        intens_imp[,i]=mat[,i]
      }
      intens_imp<-MSnSet(exprs = intens_imp)
    }else if(imput=='GMS'){
      intens_tran = exprs(intens_norm)
      intens_imp<-GMS.Lasso(as.data.frame(intens_tran),nfolds=3, log.scale=F, TS.Lasso=TRUE)
      intens_imp<-MSnSet(exprs = as.matrix(intens_imp))
    }else if(imput=='GRR'){
      intens_tran = exprs(intens_norm)
      library(DreamAI)
      intens_imp<-impute.RegImpute(data=as.matrix(intens_tran), fillmethod = "row_mean", maxiter_RegImpute = 10,conv_nrmse = 1e-03)
      intens_imp<-MSnSet(exprs = intens_imp)
    }else if(imput=='SeqKNN'){
      intens_tran = exprs(intens_norm)
      intens_imp<-SeqKNN(intens_tran, k=10)
      intens_imp<-MSnSet(exprs = intens_imp)
    }
  }else{
    intens_imp = intens_norm
  }

  process_data<-pro_file

  if(imput=='bpca'){
    process_data<-pro_file[idx_bpca,]
  }else if(imput=='knn'){
    process_data<-pro_file[idx_knn,]
  }else{
    process_data<-pro_file
  }
  return(list(filtered=process_data, normed=intens_imp@assayData[["exprs"]]))
}

#prepro_res = preprocessing_raw(file_name, designs, platform, inten_type, imput = imput, normal = normal, log2=logT)
