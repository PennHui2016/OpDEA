library(diann)
fold_diann<-'/home/hui/Documents/VM_share/DIA_data/he_diann/'
save_fold<-'/data/res_DE_files/he_dia_data/'

pg_matrix<-read.table(paste0(fold_diann, 'report.pg_matrix.tsv'), quote = "", sep = '\t', header = TRUE)
report<-read.table(paste0(fold_diann, 'report.tsv'), quote = "", sep = '\t', header = TRUE)

c1A<-c('_05', '_09', '_13', '_21', '_25', '_29', '_33', '_41')
c1B<-c('_06', '_10', '_14', '_18', '_26', '_30', '_34', '_38')

getcomb<-function(n, m){
  as<-combn(c(1:8), n)
  idxs<-vector()
  idxs<-rbind(idxs, as[,1])
  for (i in 2:length(as[1,])) {
    ta<-as[,i]
    flag=0
    for (k in 1:length(idxs[,1])) {
      com<-intersect(ta, idxs[k,])
      if(length(com)<m){
          flag=flag+1
        }
      }
      if(flag==length(idxs[,1])){
        idxs<-rbind(idxs,ta)
      }
    
    }
  return(idxs)
}

cbns<-vector()
for (n in 3:7) {
  idxs<-getcomb(n,n-1)
  for (i in 1:length(idxs[,1])) {
    cbn = c(rep(0, 8))
    cbn[1:length(idxs[i,])] = idxs[i,]
    cbns<-rbind(cbns, cbn)
  }
}
cbns<-rbind(cbns, c(1:8))

set.seed(123)

rnd1<-sample(c(1:length(cbns[, 1])), length(cbns[,1]))
rnd2<-sample(c(1:length(cbns[, 1])), length(cbns[,1]))

rnd1<-cbns[rnd1,]
rnd2<-cbns[rnd2,]

write.table(cbind(rnd1,rnd2), paste0(save_fold, '/select_reps.csv'), sep = ',', row.names = FALSE, col.names = FALSE,quote = FALSE)

for (i in 1:length(rnd1[,1])) {
#for (i in 1:1) {
  idx1<-rnd1[i,]
  idx2<-rnd2[i,]
  idx1<-idx1[which(idx1!=0)]
  idx2<-idx2[which(idx2!=0)]
  repAs<-setdiff(c1A, c1A[idx1])
  repBs<-setdiff(c1B, c1B[idx2])
  name_pg<-colnames(pg_matrix)
  name_report<-colnames(report)
  col_report<-report$Run
  
  idx_crep_proA<-c()
  for (rep in repAs) {
    idx_crep_proA<-c(idx_crep_proA, grep(rep, name_pg))
  }
  
  idx_crep_proB<-c()
  for (rep in repBs) {
    idx_crep_proB<-c(idx_crep_proB, grep(rep, name_pg))
  }
  
  idx_remian_pro<-setdiff(c(1:length(name_pg)), c(idx_crep_proA, idx_crep_proB))
  out_pg<-pg_matrix[, idx_remian_pro]
  
  idx_report_A<-c()
  for (rep in c1A[idx1]) {
    idx_report_A<-c(idx_report_A, grep(rep, col_report))
  }
  idx_report_B<-c()
  for (rep in c1B[idx2]) {
    idx_report_B<-c(idx_report_B, grep(rep, col_report))
  }
  out_report<-report[c(idx_report_A, idx_report_B),]
  
  inten_types<-c('raw', 'norm_raw', 'MaxLFQ', 'MaxLFQ_raw')
  
  idx_non_contam<-which(!grepl('contam_', out_report$Protein.Group))
  process_data = out_report[idx_non_contam,]
  idx_non_decoy<- which(!grepl('rev_', process_data$Protein.Group))
  process_data<-process_data[idx_non_decoy,]
  df<-process_data
  for (inten_type in inten_types) {
    if(inten_type=='raw'){
      proteingroup<-diann_matrix(df, id.header="Protein.Names", quantity.header ='PG.Quantity', pg.q = 0.01)
      write.table(proteingroup, paste0(save_fold, 'raw', as.character(i), '.csv'), sep = ',', col.names = TRUE, row.names = TRUE)
      #intens<-proteingroup
      #intens[intens==0] <- NA
    }else if(inten_type=='norm_raw'){
      proteingroup<-diann_matrix(df, id.header="Protein.Names", quantity.header ='PG.Normalised', pg.q = 0.01)
      write.table(proteingroup, paste0(save_fold, 'norm_raw', as.character(i), '.csv'), sep = ',', col.names = TRUE, row.names = TRUE)
    }else if(inten_type=='MaxLFQ_raw'){
      proteingroup<-diann_matrix(df, id.header="Protein.Names", quantity.header ='PG.MaxLFQ', pg.q = 0.01)
      write.table(proteingroup, paste0(save_fold, 'MaxLFQ_raw', as.character(i), '.csv'), sep = ',', col.names = TRUE, row.names = TRUE)
    }else if(inten_type=='MaxLFQ'){
      proteingroup<-diann_maxlfq(df[df$Q.Value <= 0.01 & df$PG.Q.Value <= 0.01,], group.header="Protein.Names", id.header = "Precursor.Id", quantity.header = "Precursor.Normalised")
      write.table(proteingroup, paste0(save_fold, 'MaxLFQ', as.character(i), '.csv'), sep = ',', col.names = TRUE, row.names = TRUE)
    }
  }
  
  write.table(out_pg, paste0(save_fold, 'report.pg_matrix', as.character(i), '.tsv'), sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)
  write.table(out_report, paste0(save_fold, 'report', as.character(i), '.tsv'), sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)
  #write.table(out_mss, paste0(fold_diann, 'MSstats', as.character(i), '.csv'), sep = ',', row.names = FALSE, col.names = TRUE, quote = FALSE)
}

design_all<-read.table('/data/res_DE_files/he_dia_data/design_he_dia.txt', sep = '\t', header = TRUE)

for (i in 1:length(rnd1[,1])) {
  idx1<-rnd1[i,]
  idx2<-rnd2[i,]
  idx1<-idx1[which(idx1!=0)]
  idx2<-idx2[which(idx2!=0)]
  repAs<-c1A[idx1]
  repBs<-c1B[idx2]
  sample_name<-c(repAs, repBs)
  design<-vector()
  for (j in sample_name) {
    idx<-grep(j, design_all$sample.name)
    design<-rbind(design, design_all[idx,])
  }
  colnames(design)<-c('sample name', 'condition', 'replicate')
  write.table(design, paste0(save_fold, 'design', as.character(i), '.txt'), sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
}