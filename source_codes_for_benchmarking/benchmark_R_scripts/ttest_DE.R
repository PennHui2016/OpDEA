library(limma)
source('/home/hui/PycharmProjects/DL_proteomics/quant_funcs/preprocessing.R')

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
# design_file = '/home/hui/Documents/VM_share/DIA_data/yeast/design_yeast.txt'
# file_name = '/home/hui/Documents/VM_share/DIA_data/yeast/maxquant/txt/proteinGroups.txt'
file_name = protein_file #'/home/hui/Documents/VM_share/DIA_data/yeast/frager/combined_protein.tsv'
#save_fold = '/home/hui/PycharmProjects/DL_proteomics/benchmark_res/yeast_'
designs = read.table(design_file, quote = "", sep = '\t', header = TRUE)

# in_type = 'maxquant'
# # in_type = 'fragpipe'
# inten_type = ''
# imput = ''
# normal = ''
prepro_res = preprocessing_raw(file_name, design_file, in_type, inten_type, imput = imput, normal = normal)

#condition = as.factor(designs$condition)
#design = model.matrix(~0+condition) # fitting without intercept

all_conts<-function(design_data){
  conds<-design_data$condition
  uni_con<-unique(conds)
  conts<-c()
  groups<-vector()
  for (i in 1:(length(uni_con)-1)) {
    for(j in (i+1):length(uni_con)){
      conts<-c(conts, paste0('condition', uni_con[j], '-', 'condition', uni_con[i]))
      groups<-rbind(groups, c(uni_con[i], uni_con[j]))
    }
  }
  return(list(conts=conts, groups=groups))
}

consts<-all_conts(designs)

compute_logFC<-function(group1, group2){
  logFC<-c()
  for (i in 1:length(group1[,1])) {
    logfc<-rowMeans(group2[i,][which(group2[i,]!=0 & !is.na(group2[i,]))])-rowMeans(group1[i,][which(group1[i,]!=0 & !is.na(group1[i,]))])
    logFC<-c(logFC, logfc)
  }
  return(logFC)
}

res_all<-vector()
for (i in 1:length(consts$conts)) {
  groups<-consts$groups
  g1<-groups[i,][1]
  g2<-groups[i,][2]
  
  intens<-prepro_res$normed[,c(grep(g1, designs$condition), grep(g2, designs$condition))]
  condition = as.factor(designs[c(grep(g1, designs$condition), grep(g2, designs$condition)),]$condition)
  design = model.matrix(~0+condition) # fitting without intercept
  intens<-as.data.frame(intens)
  
  intens$na_g1 = apply(intens[,c(grep(g1, colnames(intens)))],1,function(x) sum(is.na(x)))
  intens$na_g2 = apply(intens[,c(grep(g2, colnames(intens)))],1,function(x) sum(is.na(x)))
  # Filter protein table. DEqMS require minimum two values for each group.
  filter_idx = which(intens$na_g1<2 & intens$na_g2<2)
  intens.filter = intens[filter_idx,][,1:(length(intens)-2)]
  
  pval<-c()
  for(j in 1:length(intens.filter[,1])){
    #print(i)
    x=intens.filter[j,]
    if(sd(x)==0){
      p=1
      pval<-c(pval, p)
    }else{
      p = t.test(as.numeric(x[c(grep(g1, condition))]), as.numeric(x[c(grep(g2, condition))]))$p.value
      pval<-c(pval, p)}
  }
  # pval = apply(intens.filter, 1, function(x)
  #   t.test(as.numeric(x[c(grep(g1, condition))]), as.numeric(x[c(grep(g2, condition))]))$p.value)
  
  logFC = compute_logFC(intens.filter[,c(grep(g1, condition))], intens.filter[,c(grep(g2, condition))])
  #logFC =rowMeans(intens.filter[,c(grep(g2, condition))])-rowMeans(intens.filter[,c(grep(g1, condition))])
  ttest.results = data.frame(protein=rownames(intens.filter),
                             logFC=logFC,P.Value = pval, 
                             adj.pval = p.adjust(pval,method = "BH"))
  #anova.results$PSMcount = psm.count.table[anova.results$gene,"count"]
  #ttest.results = anova.results[with(anova.results,order(P.Value)),]
  ttest.results$contrast<-consts$conts[i]
  res_all<-rbind(res_all, ttest.results)
}
colnames(res_all)[3:4]<-c('pvalue', 'adj.pvalue')
res_all<-res_all[order(res_all$pvalue),]
write.table(res_all, paste0(save_fold, 'ttest_', in_type, '_', inten_type, '_', imput, '_', normal, '.csv'), sep = ',', row.names = TRUE, col.names = TRUE)