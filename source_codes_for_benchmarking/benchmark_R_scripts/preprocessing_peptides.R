#########################
### three types of input: 
###     1.fragpipe output --- combine_peptide.tsv
###     2. maxquant output --- peptides.txt
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
#file_name = '/home/hui/Documents/VM_share/DIA_data/yeast/maxquant/txt/peptides.txt'
#design_file = '/home/hui/Documents/VM_share/DIA_data/yeast/design_yeast.txt'
#in_type = 'maxquant'
#inten_type = ''
#imput = ''
#normal = ''
set.seed(123)
preprocessing_raw_pep<-function(file_name, design_file, in_type, inten_type, imput='knn',normal='quantiles'){
  #normal=''
  # in_type = 'maxquant' | 'fragpipe' | 'maxtrix'
  #in_type = 'maxquant'
  # inten_type = 'basic' | 'LFQ' for maxquant and 'basic' | 'MaxLFQ' for fragpipe (basic is default)
  # inten_type = '' 
  # file_name is the full address of the uploaded file
  
  # design_file is the uploaded file containing the experiment design containing at least: sample name/condition
  # design_file = '/home/hui/Documents/VM_share/DIA_data/yeast/design_yeast.txt'
  raw_data = read.table(file_name, quote = "", sep = '\t', header = TRUE)
  design = read.table(design_file, quote = "", sep = '\t', header = TRUE)
  col_idx = 0
  if(in_type=='maxquant'){
    # remove decoy matches and matches to contaminant
    process_data = raw_data[!raw_data$Reverse=="+",]
    process_data = process_data[!process_data$Potential.contaminant=="+",]
    
    if(inten_type=='LFQ'){
      col_idx = grep('LFQ.intensity.', colnames(process_data))
      intens = process_data[,col_idx]
      intens[intens==0] <- NA
    }else{
      col_idx = grep('Intensity.', colnames(process_data))
      intens = process_data[,col_idx]
      intens[intens==0] <- NA
    }
    row.names(intens)=process_data$Sequence
    colnames(intens)=design$sample.name
  }else if(in_type=='fragpipe'){
    # fragpipe already remove decoy and contaminant
    idx_non_contam<-which(!grepl('contam_', raw_data$Protein))
    process_data = raw_data[idx_non_contam,]
    idx_non_decoy<- which(!grepl('rev_', process_data$Protein))
    process_data<-process_data[idx_non_decoy,]
    if(inten_type=='MaxLFQ'){
      col_idx = grep('.MaxLFQ.Intensity', colnames(process_data))
      intens = process_data[,col_idx]
      intens[intens==0] <- NA
    }else{
      idxs<-c()
      for (i in 1:length(design$sample.name)) {
        idx<-which(colnames(process_data)==paste0(design$sample.name[i],'.Intensity'))
        idxs<-c(idxs, idx)
      }
      
      intens = process_data[,idxs]
      intens[intens==0] <- NA
       col_idx = idxs
    }
    row.names(intens)=process_data$Peptide.Sequence
    colnames(intens)=design$sample.name
  }
  
  #intens_tran<-as.matrix(intens)
  intens_tran<-log2(as.matrix(intens))
  intens_imp<-MSnSet(exprs = intens_tran)
  
  #### normalization
  if(normal!=''){
    intens_norm<-normalise(intens_imp, normal)
  }else{
    intens_norm<-intens_imp
  }
  
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
    }else if(imput=='mi'){
      library(mi)
      intens_tran = exprs(intens_norm)
      intens_imp<-missing_data.frame(intens_tran)
      intens_imp<-mi(intens_imp)
      intens_imp<-MSnSet(exprs = intens_imp)
    }
  }else{
    intens_imp = intens_norm
  }
  
  processed_replace = process_data
  value = 2^(intens_imp@assayData[["exprs"]])
  value[is.na(value)] <- 0
  if(imput=='bpca'){
    processed_replace = processed_replace[idx_bpca,]
    processed_replace[,col_idx]=value
  }else if(imput=='knn'){
    processed_replace = processed_replace[idx_knn,]
    processed_replace[,col_idx]=value
  }else{
    processed_replace[,col_idx]=value
  }
  return(list(filtered=process_data, normed=intens_imp@assayData[["exprs"]], processed = processed_replace, col=col_idx))
}

#prepro_res = preprocessing_raw_pep(file_name, design_file, in_type, inten_type, imput = imput, normal = normal)
