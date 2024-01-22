#########################
### three types of input: 
###     1.fragpipe output --- combine_protein.tsv
###     2. maxquant output --- proteinGroup.txt
###     3. pure intensity matrix with protein id in the first column and sample names in the columan name
### preprocessing include:
###     1. remove decoy and contaminant
###     2. replace 0 with 'NA' if no NA included as missing values
###     #3. filter out proteins with only one non-NA value in a condition
###     4. if missing value imputation is required then conduct imputation with selected imputation method (KNN default)
###        include: 'bpca', 'knn', 'QRILC', 'MLE', 'MinDet', 'MinProb', 'min', 'zero', 'mixed', 'nbavg', 'with', 'none', 'mice', 'missForest', 'mi'
###     5. log2 transform and normailization (quantile or selected method)
###        include: quantiles, quantiles.robust, vsn, max, sum    
library(MSnbase)
library(GMSimpute)
library(SeqKnn)
library(rrcovNA)
library(NormalyzerDE)
library(limma)
set.seed(123)
preprocessing_raw_count<-function(file_name, designs, in_type, imput='knn',normal='quantiles', log2=F){
  pro_file<-read.table(file_name, sep = '\t', header = T)
  counts<-pro_file[, 3:length(colnames(pro_file))]
  counts[counts==0]<-NA
  row.names(counts) <-pro_file[,1]
  
  colnames(counts)<-gsub('.Spectral.Count','',colnames(counts))
  colnames(counts)<-gsub('MS.MS.count.','',colnames(counts))
  
  idx_retain<-which(apply(counts,1,function(x) sum(is.na(x)))<length(counts[1,])*0.8)
  counts=counts[idx_retain,]
  pro_file<-pro_file[idx_retain,]
  
  fd <- data.frame(row.names(counts))
  row.names(fd)<-row.names(counts)
  pd <- data.frame(colnames(counts))
  row.names(pd) = colnames(counts)
  intens_imp<-MSnSet(as.matrix(counts), fd, pd)
  
  #### normalization
  if(normal!=''){
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
      #intens_norm = normalizeCyclicLoess(intens_imp@assayData[["exprs"]], method = "fast")
      #intens_norm<-MSnSet(exprs = intens_norm)
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
  
  #intens_norm_m<-intens_norm@assayData[["exprs"]]
  exprs_norm = intens_norm@assayData[["exprs"]]
  exprs_norm[is.infinite(exprs_norm)]=NA
  intens_norm<-MSnSet(exprs = exprs_norm)
  
  #### imputation
  if(imput!=''){
    if(imput %in% c('bpca', 'knn', 'QRILC', 'MLE', 'MinDet', 'MinProb', 'min', 'zero', 'mixed', 'nbavg', 'with', 'none')){
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
      #library(mi)
      intens_tran = exprs(intens_norm)
      #idx_gms<-which(apply(intens_tran,1,function(x) sum(is.na(x)))<length(intens_tran[1,])*0.8)
      #intens_tran<-intens_tran[idx_gms,]
      #intens_imp<-missing_data.frame(intens_tran)
      intens_imp<-GMS.Lasso(as.data.frame(intens_tran),nfolds=3, log.scale=F, TS.Lasso=TRUE)
      intens_imp<-MSnSet(exprs = as.matrix(intens_imp))
    }else if(imput=='GRR'){
      intens_tran = exprs(intens_norm)
      #idx_GRR<-which(apply(intens_tran,1,function(x) sum(is.na(x)))<length(intens_tran[1,]))
      #intens_tran<-intens_tran[idx_GRR,]
      library(DreamAI)
      intens_imp<-impute.RegImpute(data=as.matrix(intens_tran), fillmethod = "row_mean", maxiter_RegImpute = 10,conv_nrmse = 1e-03)
      intens_imp<-MSnSet(exprs = intens_imp)
    }else if(imput=='SeqKNN'){
      #library(mi)
      intens_tran = exprs(intens_norm)
      #intens_imp<-missing_data.frame(intens_tran)
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

#prepro_res = preprocessing_raw_count(file_name, designs, platform, inten_type, imput = imput, normal = normal)