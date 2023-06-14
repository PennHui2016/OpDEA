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
set.seed(123)
preprocessing_raw_count<-function(file_name, design_file, in_type, imput='knn',normal='quantiles'){
  # in_type = 'maxquant' | 'fragpipe' | 'maxtrix'
  #in_type = 'maxquant'
  # file_name is the full address of the uploaded file
  # file_name = '/home/hui/Documents/VM_share/DIA_data/yeast/maxquant/txt/proteinGroups.txt'
  # design_file is the uploaded file containing the experiment design containing at least: sample name/condition
  # design_file = '/home/hui/Documents/VM_share/DIA_data/yeast/design_yeast.txt'
  raw_data = read.table(file_name, quote = "", sep = '\t', header = TRUE)
  design = read.table(design_file, quote = "", sep = '\t', header = TRUE)
  if(in_type=='maxquant'){
    # remove decoy matches and matches to contaminant
    process_data = raw_data[!raw_data$Reverse=="+",]
    process_data = process_data[!process_data$Potential.contaminant=="+",]
  
    counts = process_data[,grep('MS.MS.count.', colnames(process_data))]
    counts[counts==0] <- NA

    row.names(counts)=process_data$Majority.protein.IDs
    colnames(counts)=design$sample.name
  }else if(in_type=='fragpipe'){
    # fragpipe already remove decoy and contaminant
    idx_non_contam<-which(!grepl('contam_', raw_data$Protein))
    process_data = raw_data[idx_non_contam,]
    idx_non_decoy<- which(!grepl('rev_', process_data$Protein))
    process_data<-process_data[idx_non_decoy,]
    
    idxs<-c()
    for (i in 1:length(design$sample.name)) {
      idx<-which(colnames(process_data)==paste0(design$sample.name[i],'.Total.Spectral.Count'))
      idxs<-c(idxs, idx)
    }
      
    counts = process_data[,idxs]
    counts[counts==0] <- NA
    
    row.names(counts)=process_data$Protein
    colnames(counts)=design$sample.name
  }
  
  counts_tran<-log2(as.matrix(counts))
  #count<-apply(counts, 2, as.numeric)
  #count<-sapply(count, as.numeric)
   #count <- lapply(data.matrix(counts), as.numeric)
  #count<-2**(counts_tran)
  counts_imp<-MSnSet(exprs = as.matrix(counts_tran))
  
  #### normalization
  if(normal!=''){
    counts_norm<-normalise(counts_imp, normal)
  }else{
    counts_norm<-counts_imp
  }
  #counts_norm_m<-counts_norm@assayData[["exprs"]]
  
  #### imputation
  if(imput!=''){
    if(imput %in% c('bpca', 'knn', 'QRILC', 'MLE', 'MinDet', 'MinProb', 'min', 'zero', 'mixed', 'nbavg', 'with', 'none')){
      if(imput=='bpca'){
        intens_tran = exprs(counts_norm)
        intens_imp<-intens_tran[which(apply(intens_tran,1,function(x) sum(is.na(x)))<length(intens_tran[1,])),]
        counts_imp<-MSnSet(exprs = intens_imp)
      }else if(imput=='knn'){
        intens_tran = exprs(counts_norm)
        intens_imp<-intens_tran[which(apply(intens_tran,1,function(x) sum(is.na(x)))<length(intens_tran[1,])*0.8),]
        counts_imp<-MSnSet(exprs = intens_imp)
      }else{
        counts_imp<-counts_norm
        }
      counts_imp<-MSnbase::impute(counts_imp, imput)
    }else if(imput=='mice'){
      library(mice)
      intens_tran = exprs(counts_norm)
      intens_imp<-mice(intens_tran)
      completeData <- complete(intens_imp,2)
      completeData<-as.matrix(completeData)
      counts_imp<-MSnSet(exprs = as.matrix(completeData))
    }else if(imput=='missForest'){
      library(missForest)
      intens_tran = exprs(counts_norm)
      intens_imp<-missForest(intens_tran)
      counts_imp<-MSnSet(exprs = intens_imp$ximp)
    }else if(imput=='mi'){
      library(mi)
      intens_tran = exprs(counts_norm)
      intens_imp<-missing_data.frame(intens_tran)
      intens_imp<-mi(intens_imp)
      counts_imp<-MSnSet(exprs = intens_imp)
    }
  }else{
    counts_imp = counts_norm
  }
  
  #counts_imp<-2^counts_imp
  
  
  return(list(filtered=process_data, normed=2^counts_imp@assayData[["exprs"]]))
}