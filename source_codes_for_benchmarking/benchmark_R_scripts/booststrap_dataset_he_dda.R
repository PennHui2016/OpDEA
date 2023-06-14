fold_frag<-'/home/hui/Documents/VM_share/DIA_data/human_ecoli/frag/'
fold_mq<-'/home/hui/Documents/VM_share/DIA_data/human_ecoli/combined/txt/'

frag_pro<-read.table(paste0(fold_frag, 'combined_protein.tsv'), quote = "", sep = '\t', header = TRUE)
frag_pep<-read.table(paste0(fold_frag, 'combined_peptide.tsv'), quote = "", sep = '\t', header = TRUE)
frag_mss<-read.table(paste0(fold_frag, 'MSstats.csv'), quote = "", sep = ',', header = TRUE)


mq_pro<-read.table(paste0(fold_mq, 'proteinGroups.txt'), quote = "", sep = '\t', header = TRUE)
mq_pep<-read.table(paste0(fold_mq, 'peptides.txt'), quote = "", sep = '\t', header = TRUE)
mq_mss<-read.table(paste0(fold_mq, 'evidence.txt'), quote = "", sep = '\t', header = TRUE)

c1A<-c('A_1', 'A_2', 'A_3', 'A_4', 'A_5', 'A_6', 'A_7', 'A_8')
c1B<-c('B_1', 'B_2', 'B_3', 'B_4', 'B_5', 'B_6', 'B_7', 'B_8')

c2A<-c('A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'A7', 'A8')
c2B<-c('B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8')

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

write.table(cbind(rnd1,rnd2), '/home/hui/Documents/VM_share/DIA_data/human_ecoli/select_reps.csv', sep = ',', row.names = FALSE, col.names = FALSE,quote = FALSE)

for (i in 1:length(rnd1[,1])) {
#for (i in 1:1) {
  idx1<-rnd1[i,]
  idx2<-rnd2[i,]
  idx1<-idx1[which(idx1!=0)]
  idx2<-idx2[which(idx2!=0)]
  repAs<-setdiff(c1A, c1A[idx1])
  repBs<-setdiff(c1B, c1B[idx2])
  name_pro<-colnames(frag_pro)
  name_pep<-colnames(frag_pep)
  col_mss<-paste0(frag_mss$Condition, '_', frag_mss$BioReplicate)
  
  idx_crep_proA<-c()
  for (rep in repAs) {
    idx_crep_proA<-c(idx_crep_proA, grep(rep, name_pro))
  }
  
  idx_crep_proB<-c()
  for (rep in repBs) {
    idx_crep_proB<-c(idx_crep_proB, grep(rep, name_pro))
  }
  
  idx_remian_pro<-setdiff(c(1:length(name_pro)), c(idx_crep_proA, idx_crep_proB))
  out_pro<-frag_pro[, idx_remian_pro]
  
  idx_crep_pepA<-c()
  for (rep in repAs) {
    idx_crep_pepA<-c(idx_crep_pepA, grep(rep, name_pep))
  }
  
  idx_crep_pepB<-c()
  for (rep in repBs) {
    idx_crep_pepB<-c(idx_crep_pepB, grep(rep, name_pep))
  }
  idx_remian_pep<-setdiff(c(1:length(name_pep)), c(idx_crep_pepA, idx_crep_pepB))
  out_pep<-frag_pep[, idx_remian_pep]
  
  idx_mss_A<-c()
  for (rep in c1A[idx1]) {
    idx_mss_A<-c(idx_mss_A, grep(rep, col_mss))
  }
  idx_mss_B<-c()
  for (rep in c1B[idx2]) {
    idx_mss_B<-c(idx_mss_B, grep(rep, col_mss))
  }
  out_mss<-frag_mss[c(idx_mss_A, idx_mss_B),]
  
  write.table(out_pro, paste0(fold_frag, 'combined_protein', as.character(i), '.tsv'), sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)
  write.table(out_pep, paste0(fold_frag, 'combined_peptide', as.character(i), '.tsv'), sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)
  write.table(out_mss, paste0(fold_frag, 'MSstats', as.character(i), '.csv'), sep = ',', row.names = FALSE, col.names = TRUE, quote = FALSE)
}

for (i in 1:length(rnd1[,1])) {
  idx1<-rnd1[i,]
  idx2<-rnd2[i,]
  idx1<-idx1[which(idx1!=0)]
  idx2<-idx2[which(idx2!=0)]
  repAs<-setdiff(c2A, c2A[idx1])
  repBs<-setdiff(c2B, c2B[idx2])
  name_pro<-colnames(mq_pro)
  name_pep<-colnames(mq_pep)
  col_mss<-mq_mss$Experiment #paste0(mq_mss$Condition, '_', frag_mss$BioReplicate)
  
  idx_crep_proA<-c()
  for (rep in repAs) {
    idx_crep_proA<-c(idx_crep_proA, grep(rep, name_pro))
  }
  
  idx_crep_proB<-c()
  for (rep in repBs) {
    idx_crep_proB<-c(idx_crep_proB, grep(rep, name_pro))
  }
  
  idx_remian_pro<-setdiff(c(1:length(name_pro)), c(idx_crep_proA, idx_crep_proB))
  out_pro<-mq_pro[, idx_remian_pro]
  
  idx_crep_pepA<-c()
  for (rep in repAs) {
    idx_crep_pepA<-c(idx_crep_pepA, grep(rep, name_pep))
  }
  
  idx_crep_pepB<-c()
  for (rep in repBs) {
    idx_crep_pepB<-c(idx_crep_pepB, grep(rep, name_pep))
  }
  idx_remian_pep<-setdiff(c(1:length(name_pep)), c(idx_crep_pepA, idx_crep_pepB))
  out_pep<-mq_pep[, idx_remian_pep]
  
  idx_mss_A<-c()
  for (rep in c2A[idx1]) {
    idx_mss_A<-c(idx_mss_A, grep(rep, col_mss))
  }
  idx_mss_B<-c()
  for (rep in c2B[idx2]) {
    idx_mss_B<-c(idx_mss_B, grep(rep, col_mss))
  }
  out_mss<-mq_mss[c(idx_mss_A, idx_mss_B),]
  
  write.table(out_pro, paste0(fold_mq, 'proteinGroups', as.character(i), '.txt'), sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)
  write.table(out_pep, paste0(fold_mq, 'peptides', as.character(i), '.txt'), sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)
  write.table(out_mss, paste0(fold_mq, 'evidence', as.character(i), '.txt'), sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)
}

msstats_design<-read.table('/home/hui/Documents/VM_share/DIA_data/human_ecoli/msstats_design_all.txt', sep = '\t', header = TRUE)

for (i in 1:length(rnd1[,1])) {
  idx1<-rnd1[i,]
  idx2<-rnd2[i,]
  idx1<-idx1[which(idx1!=0)]
  idx2<-idx2[which(idx2!=0)]
  repAs<-c1A[idx1]
  repBs<-c1B[idx2]
  sample_name<-c(repAs, repBs)
  condition<-c()
  for (j in 1:length(sample_name)) {
    condition<-c(condition, substr(sample_name[j], 1,1))
  }
  #condition<-c(apply(sample_name, 1, function(x) substr(x, 1,1)))
  replicate<-sample_name
  out_design<-cbind(sample_name, condition,replicate)
  colnames(out_design)<-c('sample name', 'condition', 'replicate')
  write.table(out_design, paste0('/home/hui/Documents/VM_share/DIA_data/human_ecoli/', 'design_he_frag', as.character(i), '.txt'), sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
  
  repAs<-c2A[idx1]
  repBs<-c2B[idx2]
  sample_name<-c(repAs, repBs)
  condition<-c()
  for (j in 1:length(sample_name)) {
    condition<-c(condition, substr(sample_name[j], 1,1))
  }
  replicate<-sample_name
  out_design<-cbind(sample_name, condition,replicate)
  colnames(out_design)<-c('sample name', 'condition', 'replicate')
  write.table(out_design, paste0('/home/hui/Documents/VM_share/DIA_data/human_ecoli/', 'design_he_mq', as.character(i), '.txt'), sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
  
  idx_design<-c()
  for (j in 1:length(sample_name)) {
    idx_design<-c(idx_design, grep(sample_name[j], msstats_design$BioReplicate))
  }
  #idx_design<-apply(sample_name, 2, function(x) grep(x, msstats_design$BioReplicate))
  out_design<-msstats_design[idx_design,]
  write.table(out_design, paste0('/home/hui/Documents/VM_share/DIA_data/human_ecoli/', 'design_he', as.character(i), '_msstats.txt'), sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
  
}
