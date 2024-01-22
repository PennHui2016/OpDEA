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
preprocessing_raw<-function(file_name, design, in_type, inten_type, imput='knn',normal='quantiles'){
  # in_type = 'Maxquant' | 'FragPipe' | 'DIANN' | 'spectronaut
  # inten_type = 'top0' | 'top3' | 'maxlfq' | 'directlfq'
  # file_name is the full address of the expression matrix
  # design_file is the uploaded file containing the experiment design containing at least: sample name/condition
  # design_file = '/home/hui/Documents/VM_share/DIA_data/yeast/design_yeast.txt'
  #raw_data = read.table(file_name, quote = "", sep = '\t', header = TRUE)
  # design = read.table(design_file, quote = "", sep = '\t', header = TRUE)
  if(in_type=='maxquant'){
    raw_data = read.table(file_name, quote = "", sep = '\t', header = TRUE)
    # remove decoy matches and matches to contaminant
    process_data = raw_data[!raw_data$Reverse=="+",]
    process_data = process_data[!process_data$Potential.contaminant=="+",]
    
    if(inten_type=='LFQ'){
      intens = process_data[,grep('LFQ.intensity.', colnames(process_data))]
      intens[intens==0] <- NA
    }else{
      intens = process_data[,grep('Intensity.', colnames(process_data))]
      intens[intens==0] <- NA
    }
    row.names(intens)=process_data$Majority.protein.IDs
    colnames(intens)=design$sample.name
  }else if(in_type=='fragpipe'){
    raw_data = read.table(file_name, quote = "", sep = '\t', header = TRUE)
    # fragpipe already remove decoy and contaminant
    idx_non_contam<-which(!grepl('contam_', raw_data$Protein))
    process_data = raw_data[idx_non_contam,]
    idx_non_decoy<- which(!grepl('rev_', process_data$Protein))
    process_data<-process_data[idx_non_decoy,]
    if(inten_type=='MaxLFQ'){
      intens = process_data[,grep('.MaxLFQ.Intensity', colnames(process_data))]
      intens[intens==0] <- NA
    }else{
      idxs<-c()
      for (i in 1:length(design$sample.name)) {
        idx<-which(colnames(process_data)==paste0(design$sample.name[i],'.Intensity'))
        idxs<-c(idxs, idx)
      }
      
      intens = process_data[,idxs]
      intens[intens==0] <- NA
    }
    row.names(intens)=process_data$Protein
    colnames(intens)=design$sample.name
  }else if(in_type=='DIANN'){
    
    if(inten_type=='raw'){
      #proteingroup<-diann_matrix(df, id.header="Protein.Group", quantity.header ='PG.Quantity', pg.q = 0.01)
      proteingroup<-read.table(gsub('.tsv', '.csv', gsub('report', 'raw', file_name)), sep = ',', header = TRUE)
      #row.names(proteingroup)<-proteingroup[,1]
      #proteingroup<-proteingroup[,2:length(proteingroup[1,])]
      intens<-proteingroup
      intens[intens==0] <- NA
    }else if(inten_type=='norm_raw'){
      #proteingroup<-diann_matrix(df, id.header="Protein.Group", quantity.header ='PG.Normalised', pg.q = 0.01)
      proteingroup<-read.table(gsub('.tsv', '.csv', gsub('report', 'norm_raw', file_name)), sep = ',', header = TRUE)
      #row.names(proteingroup)<-proteingroup[,1]
      #proteingroup<-proteingroup[,2:length(proteingroup[1,])]
      intens<-proteingroup
      intens[intens==0] <- NA
    }else if(inten_type=='MaxLFQ_raw'){
      #proteingroup<-diann_matrix(df, id.header="Protein.Group", quantity.header ='PG.MaxLFQ', pg.q = 0.01)
      proteingroup<-read.table(gsub('.tsv', '.csv', gsub('report', 'MaxLFQ_raw', file_name)), sep = ',', header = TRUE)
      #row.names(proteingroup)<-proteingroup[,1]
      #proteingroup<-proteingroup[,2:length(proteingroup[1,])]
      intens<-proteingroup
      intens[intens==0] <- NA
    }else if(inten_type=='MaxLFQ'){
      #proteingroup<-diann_maxlfq(df[df$Q.Value <= 0.01 & df$PG.Q.Value <= 0.01,], group.header="Protein.Group", id.header = "Precursor.Id", quantity.header = "Precursor.Normalised")
      proteingroup<-read.table(gsub('.tsv', '.csv', gsub('report', 'MaxLFQ', file_name)), sep = ',', header = TRUE)
      #row.names(proteingroup)<-proteingroup[,1]
      #proteingroup<-proteingroup[,2:length(proteingroup[1,])]
      intens<-proteingroup
      intens[intens==0] <- NA
    }else if(inten_type=='norm'){
      raw_data = read.table(file_name, quote = "", sep = '\t', header = TRUE)
      idx_non_contam<-which(!grepl('contam_', raw_data$Protein.Group))
      process_data = raw_data[idx_non_contam,]
      idx_non_decoy<- which(!grepl('rev_', process_data$Protein.Group))
      process_data<-process_data[idx_non_decoy,]
      df<-process_data
      proteingroup<-process_data
      idxs<-c()
      for (i in 1:length(design$sample.name)) {
        idx<-which(colnames(proteingroup)==design$sample.name[i])
        idxs<-c(idxs, idx)
      }
      
      intens = proteingroup[,idxs]
      intens[intens==0] <- NA
        row.names(intens)=proteingroup$Protein.Group
    }
    #row.names(intens)=proteingroup$Protein.Group
    cni<-c()
    for (q in 1:length(colnames(intens))) {
      idx_q = which(design$sample.name==colnames(intens)[q])
      cni<-c(cni, design$replicate[q])
    }
    colnames(intens)=cni
    #colnames(intens)=design$sample.name
  }
  
  intens_tran<-log2(as.matrix(intens))
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
    }else if(normal=='Loess'){
      intens_norm = normalizeCyclicLoess(intens_imp@assayData[["exprs"]], method = "fast")
      intens_norm<-MSnSet(exprs = intens_norm)
    }else if(normal=='TIC'){
      exp_mat = intens_imp@assayData[["exprs"]]
      tics = colSums(exp_mat)
      intens_norm=vector()
      for (i in 1:length(tics)) {
        intes_norm = cbind(intes_norm, exp_mat[,i]/tics[i])
      }
      intes_norm = intes_norm * median(tics)
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
    }else if(imput=='Impseq'){
      library(mi)
      intens_tran = exprs(intens_norm)
      #intens_imp<-missing_data.frame(intens_tran)
      intens_imp<-impSeq(intens_imp)
      intens_imp<-MSnSet(exprs = intens_imp)
    }else if(imput=='Impseqrob'){
      library(mi)
      intens_tran = exprs(intens_norm)
      #intens_imp<-missing_data.frame(intens_tran)
      intens_imp<-impSeqRob(intens_imp)
      intens_imp<-MSnSet(exprs = intens_imp)
    }else if(imput=='GMS'){
      library(mi)
      intens_tran = exprs(intens_norm)
      #intens_imp<-missing_data.frame(intens_tran)
      intens_imp<-GMS.Lasso(intens_imp,log.scale=TRUE,TS.Lasso=TRUE)
      intens_imp<-MSnSet(exprs = intens_imp)
    }else if(imput=='SeqKNN'){
      library(mi)
      intens_tran = exprs(intens_norm)
      #intens_imp<-missing_data.frame(intens_tran)
      intens_imp<-SeqKNN(intens_imp, k=10)
      intens_imp<-MSnSet(exprs = intens_imp)
    }
  }else{
    intens_imp = intens_norm
  }
  
  if(in_type=='DIANN'){
    if(length(colnames(proteingroup))==length(design$condition)){
      proteingroup<-cbind(row.names(proteingroup), proteingroup)
      colnames(proteingroup)[1]<-'Protein.ID'
      proteingroup<-as.data.frame(proteingroup)
    }
    if(imput=='bpca'){
      proteingroup<-proteingroup[idx_bpca,]
    }else if(imput=='knn'){
      proteingroup<-proteingroup[idx_knn,]
    }
    return(list(filtered=proteingroup, normed=intens_imp@assayData[["exprs"]]))
  }else{
    if(imput=='bpca'){
      process_data<-process_data[idx_bpca,]
    }else if(imput=='knn'){
      process_data<-process_data[idx_knn,]
    }
  return(list(filtered=process_data, normed=intens_imp@assayData[["exprs"]]))
  }
}