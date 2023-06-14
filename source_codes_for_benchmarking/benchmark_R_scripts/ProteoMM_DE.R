library(ProteoMM)
source('/home/hui/PycharmProjects/DL_proteomics/quant_funcs/preprocessing_peptides_ProteoMM.R')

args <- commandArgs(trailingOnly = TRUE)
#print(args)
protein_file<-args[1]
peptide_file<-args[2]
design_file<-args[3]
msstats_file<-args[4]
in_type<-args[5]
inten_type<-args[6]
imput<-args[7]
normal<-args[8]
save_fold<-args[9]

if(imput=='blank'){
  imput=''
}

if(normal=='blank'){
  normal=''
}

if(inten_type=='basic'){
  inten_type=''
}
set.seed(123)
#file_name = '/home/hui/Documents/VM_share/DIA_data/yeast/maxquant/txt/peptides.txt'
file_name = peptide_file #'/home/hui/Documents/VM_share/DIA_data/yeast/frager/combined_peptide.tsv'
#design_file = '/home/hui/Documents/VM_share/DIA_data/yeast/design_yeast.txt'
#save_fold = '/home/hui/PycharmProjects/DL_proteomics/benchmark_res/yeast_'
designs = read.table(design_file, quote = "", sep = '\t', header = TRUE)

#in_type = 'maxquant'
# in_type = 'fragpipe'
# inten_type = ''
# imput = ''
# ## imput = '', 'bpca', 'knn', 'QRILC', 'MLE', 'MinDet', 'MinProb', 'min', 'zero', 'mixed', 'nbavg', 'with', 'none', 'mice', 'missForest', 'mi'
# normal = ""
## normal = '', "sum", "max", "center.mean", "center.median", "div.mean", "div.median", "diff.meda", "quantiles", "quantiles.robust" or "vsn"

prepro_res = preprocessing_raw_pep_MM(file_name, design_file, in_type, inten_type, imput = imput, normal = normal)

all_conts<-function(design_data){
  conds<-design_data$condition
  uni_con<-unique(conds)
  conts<-c()
  groups<-vector()
  #levels<-c()
  levels<-c(paste0('condition', uni_con))
  for (i in 1:(length(uni_con)-1)) {
    for(j in (i+1):length(uni_con)){
      conts<-c(conts, paste0('condition', uni_con[j], '-', 'condition', uni_con[i]))
      groups<-rbind(groups, c(uni_con[i], uni_con[j]))
    }
  }
  return(list(conts=conts, groups=groups, levels=levels))
}

consts<-all_conts(designs)


res_all<-vector()
for (k in 1:length(consts$conts)) {
  g1 = consts$groups[k,][1]
  g2 = consts$groups[k,][2]
  idx_g1 = grep(g1, designs$condition)
  idx_g2 = grep(g2, designs$condition)
  test_col = prepro_res$col[c(idx_g1, idx_g2)]
  intsCols = test_col
  if(in_type=='maxquant'){
    metaCols = c(grep('Sequence', colnames(prepro_res$processed)), which(colnames(prepro_res$processed)=='Proteins'))
  }else if(in_type=='fragpipe'){
    metaCols = c(grep('Peptide.Sequence', colnames(prepro_res$processed)), which(colnames(prepro_res$processed)=='Protein'))
  }
  
  m_logInts = make_intencities(prepro_res$processed, intsCols)
  m_prot.info = make_meta(prepro_res$processed, metaCols)
  m_logInts = convert_log2(m_logInts)
  
  m_logInts$na_g1 = apply(m_logInts[,c(grep(g1, colnames(m_logInts)))],1,function(x) sum(is.na(x)))
  m_logInts$na_g2 = apply(m_logInts[,c(grep(g2, colnames(m_logInts)))],1,function(x) sum(is.na(x)))
  # Filter protein table. DEqMS require minimum two values for each group.
  filter_idx = which(m_logInts$na_g1<length(idx_g1) & m_logInts$na_g2<length(idx_g2))
  m_logInts = m_logInts[filter_idx,][,1:(length(m_logInts[1,])-2)]
  m_prot.info = m_prot.info[filter_idx,]
  
  grps = as.factor(designs$condition[c(idx_g1, idx_g2)])
  
  set.seed(135) # results rarely vary due to the random seed for EigenMS
  
  if(normal!=''){
    # already normed
    mm_m_ints_norm = list(normalized=cbind(m_prot.info, m_logInts))
  }else{
    mm_m_ints_eig1 = eig_norm1(m=m_logInts,treatment=grps,prot.info=m_prot.info)
    mm_m_ints_eig1$h.c # check the number of bias trends detected
    mm_m_ints_norm = eig_norm2(rv=mm_m_ints_eig1)
  }
  
  if(imput!=''){
    #already imputed
    data_imp = mm_m_ints_norm$normalized[,3:length(mm_m_ints_norm$normalized[1,])]
    data_info = mm_m_ints_norm$normalized[,1:2]
  }else{
    mm_prot.info = mm_m_ints_norm$normalized[,1:2]
    mm_norm_m =  mm_m_ints_norm$normalized[,3:length(mm_m_ints_norm$normalized[1,])]
    
    set.seed(131) # important to reproduce the results later
    imp_mm = MBimpute(mm_norm_m, grps, prot.info=mm_prot.info,
                      pr_ppos=2, my.pi=0.05,
                      compute_pi=FALSE)
    data_imp = imp_mm$y_imputed
    data_info = mm_prot.info
  }
  
  res.ProteoMM = peptideLevel_DE(data_imp,
                           grps, data_info,
                           pr_ppos=2)
  
  res.ProteoMM$contrast<-consts$conts[k]
  res_all<-rbind(res_all, res.ProteoMM)
}

logFC<-c()
for (i in 1:length(res_all$FC)) {
  if (res_all$FC[i]<0){
    logFC<-c(logFC, -abs(log2(abs(res_all$FC[i]))))
  }else{
    logFC<-c(logFC, abs(log2(abs(res_all$FC[i]))))
  }
}

res_all$logFC<-logFC
colnames(res_all)[c(1,3,4)]<-c('protein','pvalue','adj.pvalue')
res_all<-res_all[order(res_all$pvalue),]
write.table(res_all, paste0(save_fold, 'ProteoMM_', in_type, '_', inten_type, '_', imput, '_', normal, '.csv'), sep = ',', row.names = TRUE, col.names = TRUE)

