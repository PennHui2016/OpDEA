#########################
### three types of input: 
###     1.fragpipe output --- combine_protein.tsv
###     2. maxquant output --- proteinGroup.txt
###     3. pure intensity matrix with protein id in the first column and sample names in the columan name
### preprocessing include:
###     1. remove decoy and contaminant
###     2. replace 0 with 'NA' if no NA included as missing values
###     #3. filter out proteins with only one non-NA value in a condition
###     
###     4. log2 transform and normailization (quantile or selected method)
###        include: quantiles, quantiles.robust, vsn, max, sum  
###     5. if missing value imputation is required then conduct imputation with selected imputation method (KNN default)
###        include: 'bpca', 'knn', 'QRILC', 'MLE', 'MinDet', 'MinProb', 'min', 'zero', 'mixed', 'nbavg', 'with', 'none', 'mice', 'missForest', 'mi'
library(MSnbase)
#library(diann)
library(GMSimpute)
library(SeqKnn)
library(rrcovNA)
library(NormalyzerDE)
library(limma)

set.seed(123)
preprocessing_raw<-function(file_name, count_file, design, in_type, inten_type, imput='knn',normal='quantiles', log2=T){
  # in_type = 'Maxquant' | 'FragPipe' | 'DIANN' | 'spectronaut
  # inten_type = 'top0' | 'top3' | 'maxlfq' | 'directlfq'
  # file_name is the full address of the expression matrix
  # design_file is the uploaded file containing the experiment design containing at least: sample name/condition
  # design_file = '/home/hui/Documents/VM_share/DIA_data/yeast/design_yeast.txt'
  #raw_data = read.table(file_name, quote = "", sep = '\t', header = TRUE)
  # design = read.table(design_file, quote = "", sep = '\t', header = TRUE)
  
  pro_file<-read.table(file_name, sep = '\t', header = T)
  ct_file<-read.table(count_file, sep = '\t', header = T)
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
  pro_count<-ct_file[,1]
  idx_count<-match(pro_file[,1],pro_count)
  ct_file<-ct_file[idx_count,]
  ct_file$Protein<-pro_file$Protein
  ct_file$Organism<-pro_file$Organism
  ct_file[is.na(ct_file)]<-0

  colnames(ct_file)<-gsub('.Spectral.Count','',colnames(ct_file))
  colnames(ct_file)<-gsub('MS.MS.count.','',colnames(ct_file))
  
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

  #### imputation
  if(imput!=''){
    if(imput %in% c('bpca', 'knn', 'QRILC', 'MLE', 'MinDet', 'MinProb', 'min', 'zero', 'mixed', 'nbavg', 'with', 'none')){
      if(imput=='bpca'){
        intens_tran = exprs(intens_norm)
        #idx_bpca<-which(apply(intens_tran,1,function(x) sum(is.na(x)))<length(intens_tran[1,]))
        intens_imp<-intens_tran#[idx_bpca,]
        intens_imp<-MSnSet(exprs = intens_imp)
      }else if(imput=='knn'){
        intens_tran = exprs(intens_norm)
        #idx_knn<-which(apply(intens_tran,1,function(x) sum(is.na(x)))<length(intens_tran[1,])*0.8)
        intens_imp<-intens_tran#[idx_knn,]
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
  
    # if(imput=='bpca'){
    #   process_data<-pro_file[idx_bpca,]
    # }else if(imput=='knn'){
    #   process_data<-pro_file[idx_knn,]
    # }else{
    #   process_data<-pro_file
    # }
  return(list(filtered=process_data, counts=ct_file, normed=intens_imp@assayData[["exprs"]]))
}

#prepro_res = preprocessing_raw(file_name, count_file, designs, platform, inten_type, imput = imput, normal = normal)