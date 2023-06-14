library(diann)
root_fold<-'/home/hui/Documents/VM_share/DIA_data/ecoli/DIANN/'
modes<-c('Variable', 'Wide', 'Narrow', 'Overlap')
inten_types<-c('raw', 'norm_raw', 'MaxLFQ', 'MaxLFQ_raw')

for (mode in modes) {
  fold<-paste0(root_fold, mode, '/report.tsv')
  raw_data = read.table(fold, quote = "", sep = '\t', header = TRUE)
  idx_non_contam<-which(!grepl('contam_', raw_data$Protein.Group))
  process_data = raw_data[idx_non_contam,]
  idx_non_decoy<- which(!grepl('rev_', process_data$Protein.Group))
  process_data<-process_data[idx_non_decoy,]
  df<-process_data
  for (inten_type in inten_types) {
    if(inten_type=='raw'){
      proteingroup<-diann_matrix(df, id.header="Protein.Group", quantity.header ='PG.Quantity', pg.q = 0.01)
      write.table(proteingroup, paste0(root_fold, mode, '/raw.csv'), sep = ',', col.names = TRUE, row.names = TRUE)
      #intens<-proteingroup
      #intens[intens==0] <- NA
    }else if(inten_type=='norm_raw'){
      proteingroup<-diann_matrix(df, id.header="Protein.Group", quantity.header ='PG.Normalised', pg.q = 0.01)
      write.table(proteingroup, paste0(root_fold, mode, '/norm_raw.csv'), sep = ',', col.names = TRUE, row.names = TRUE)
    }else if(inten_type=='MaxLFQ_raw'){
      proteingroup<-diann_matrix(df, id.header="Protein.Group", quantity.header ='PG.MaxLFQ', pg.q = 0.01)
      write.table(proteingroup, paste0(root_fold, mode, '/MaxLFQ_raw.csv'), sep = ',', col.names = TRUE, row.names = TRUE)
    }else if(inten_type=='MaxLFQ'){
      proteingroup<-diann_maxlfq(df[df$Q.Value <= 0.01 & df$PG.Q.Value <= 0.01,], group.header="Protein.Group", id.header = "Precursor.Id", quantity.header = "Precursor.Normalised")
      write.table(proteingroup, paste0(root_fold, mode, '/MaxLFQ.csv'), sep = ',', col.names = TRUE, row.names = TRUE)
    }
  }
}