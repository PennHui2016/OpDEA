library(DEqMS)
source('/home/hui/PycharmProjects/DL_proteomics/quant_funcs/preprocessing.R')
library(matrixStats)

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

#design_file = '/home/hui/Documents/VM_share/DIA_data/yeast/design_yeast.txt'
#file_name = '/home/hui/Documents/VM_share/DIA_data/yeast/maxquant/txt/proteinGroups.txt'
file_name = protein_file #'/home/hui/Documents/VM_share/DIA_data/yeast/frager/combined_protein.tsv'
#save_fold = '/home/hui/PycharmProjects/DL_proteomics/benchmark_res/yeast_'
designs = read.table(design_file, quote = "", sep = '\t', header = TRUE)
set.seed(123)
#in_type = 'maxquant'
# in_type = 'fragpipe'
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
  
  # use minimum peptide count among  samples
  # count unique+razor peptides used for quantification
  
  if(in_type=='maxquant'){
    col_idx_pep_c<-c(grep(paste0('Razor...unique.peptides.', groups[i,][1]), colnames(prepro_res$filtered)), grep(paste0('Razor...unique.peptides.', groups[i,][2]), colnames(prepro_res$filtered)))
    pep.count.table = data.frame(count = rowMins(as.matrix(prepro_res$filtered[filter_idx,][,col_idx_pep_c])),
                                 row.names = prepro_res$filtered$Majority.protein.IDs[filter_idx])
    # Minimum peptide count of some proteins can be 0
    # add pseudocount 1 to all proteins
    pep.count.table$count = pep.count.table$count+1
  }else if(in_type=='fragpipe'){
    spc_idx<-c()
    for (j in 1:length(designs$sample.name)) {
      idx<-which(colnames(prepro_res$filtered)==paste0(designs$sample.name[j],'.Total.Spectral.Count'))
      spc_idx<-c(spc_idx, idx)
    }
    counts<-prepro_res$filtered[filter_idx,][,spc_idx]
    col_idx_pep_c<-c(grep(groups[i,][1], colnames(counts)), grep(groups[i,][2], colnames(counts)))
    pep.count.table = data.frame(count = rowMins(as.matrix(counts[,col_idx_pep_c])),
                                 row.names = prepro_res$filtered$Protein[filter_idx])
    # Minimum peptide count of some proteins can be 0
    # add pseudocount 1 to all proteins
    pep.count.table$count = pep.count.table$count+1
  }
  fit1 = lmFit(intens.filter,design = design)
  cont <- makeContrasts(contrasts =  consts$conts[i], levels = design)
  fit2 = contrasts.fit(fit1,contrasts = cont)
  fit3 <- eBayes(fit2)
  fit3$count = pep.count.table[rownames(fit3$coefficients),"count"]
  fit3$sigma[which(fit3$sigma==0)]=1e-14
  fit4 = spectraCounteBayes(fit3)
  DEqMS.results = outputResult(fit4, coef_col = 1)
  #limma.results = topTable(fit3, adjust="BH", sort.by = 'logFC', n=Inf)
  DEqMS.results<-cbind(DEqMS.results, rep(consts$conts[i], length(DEqMS.results[,1])))
  res_all<-rbind(res_all, DEqMS.results)
}
colnames(res_all)[length(res_all[1,])]<-'contrast'
colnames(res_all)[c(7, 10,11)]<-c('protein', 'pvalue','adj.pvalue')
res_all<-res_all[order(res_all$pvalue),]
write.table(res_all, paste0(save_fold, 'DEqMS_', in_type, '_', inten_type, '_', imput, '_', normal, '.csv'), sep = ',', row.names = TRUE, col.names = TRUE)

