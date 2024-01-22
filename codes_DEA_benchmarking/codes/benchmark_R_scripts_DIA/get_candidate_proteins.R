library(MSnbase)
library(diann)
set.seed(123)

args <- commandArgs(trailingOnly = TRUE)
in_type = args[1]
file_name<-args[2]
save_folder<-args[3]

#file_name<-'/home/hui/Documents/VM_share/DIA_data/yeast/frager/'
#file_name<-'/home/hui/Documents/VM_share/DIA_data/yeast/maxquant/txt/'
#file_name<-'/home/hui/Documents/VM_share/DIA_data/ecoli/DIANN/Variable/'
#file_name<-'/home/hui/Documents/VM_share/DIA_data/ecoli/DIANN/Variable/report.tsv'
#save_folder = '/home/hui/PycharmProjects/DL_proteomics/benchmark_res/metrics_new/test.csv'
#in_type = 'maxquant'
#in_type = 'fragpipe'
#in_type = 'DIANN'

if(in_type=='maxquant'){
  file_name <-paste0(file_name, 'proteinGroups.txt')
  raw_data = read.table(file_name, quote = "", sep = '\t', header = TRUE)
  # remove decoy matches and matches to contaminant
  process_data = raw_data[!raw_data$Reverse=="+",]
  process_data = process_data[!process_data$Potential.contaminant=="+",]
  candidate_proteins<-process_data$Protein.IDs
}else if(in_type=='fragpipe'){
  file_name <-paste0(file_name, 'combined_protein.tsv')
  raw_data = read.table(file_name, quote = "", sep = '\t', header = TRUE)
  # fragpipe already remove decoy and contaminant
  idx_non_contam<-which(!grepl('contam_', raw_data$Protein))
  process_data = raw_data[idx_non_contam,]
  idx_non_decoy<- which(!grepl('rev_', process_data$Protein))
  process_data<-process_data[idx_non_decoy,]
  candidate_proteins<-process_data$Protein
}else if(in_type=='DIANN'){
  file_name1 <-paste0(file_name, 'report.pg_matrix.tsv')
  raw_data1 = read.table(file_name1, quote = "", sep = '\t', header = TRUE)
  
  file_name2 <-paste0(file_name, 'report.tsv')
  raw_data2 = read.table(file_name2, quote = "", sep = '\t', header = TRUE)
  
  idx_non_contam1<-which(!grepl('contam_', raw_data1$Protein.Group))
  process_data1 = raw_data1[idx_non_contam1,]
  idx_non_decoy1<- which(!grepl('rev_', process_data1$Protein.Group))
  process_data1<-process_data1[idx_non_decoy1,]
  
  idx_non_contam2<-which(!grepl('contam_', raw_data2$Protein.Group))
  process_data2 = raw_data2[idx_non_contam2,]
  idx_non_decoy2<- which(!grepl('rev_', process_data2$Protein.Group))
  process_data2<-process_data2[idx_non_decoy2,]
  
  df<-process_data2
 
    proteingroup1<-diann_matrix(df, id.header="Protein.Group", quantity.header ='PG.Quantity', pg.q = 0.01)
    candidate_proteins1 = row.names(proteingroup1)
  
    proteingroup2<-diann_maxlfq(df[df$Q.Value <= 0.01 & df$PG.Q.Value <= 0.01,], group.header="Protein.Group", id.header = "Precursor.Id", quantity.header = "Precursor.Normalised")
    candidate_proteins2 = row.names(proteingroup2)

    proteingroup3<-process_data1
    candidate_proteins3 = proteingroup3$Protein.Ids
    candidate_proteins = unique(c(candidate_proteins1, candidate_proteins2, candidate_proteins3))
  #process_data<-proteingroup
}

write.table(as.data.frame(candidate_proteins), save_folder, sep=',', col.names = FALSE, row.names = FALSE)