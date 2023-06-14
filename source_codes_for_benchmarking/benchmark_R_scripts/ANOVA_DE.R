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

# design_file = '/home/hui/Documents/VM_share/DIA_data/yeast/design_yeast.txt'
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
  
  fit1 = lmFit(intens.filter,design = design)
  ord.t = fit1$coefficients[, 1]/fit1$sigma/fit1$stdev.unscaled[, 1]
  ord.p = 2*pt(abs(ord.t), fit1$df.residual, lower.tail = FALSE)
  ord.q = p.adjust(ord.p,method = "BH")
  log2FC = fit1$coefficients[,2]-fit1$coefficients[,1]
  anova.results = data.frame(protein=row.names(intens.filter),
                             logFC=log2FC,
                             t=ord.t, 
                             P.Value=ord.p, 
                             adj.P.Val = ord.q)
  
  #anova.results$PSMcount = psm.count.table[anova.results$gene,"count"]
  #anova.results = anova.results[with(anova.results,order(P.Value)),]
  anova.results$contrast<-consts$conts[i]
  res_all<-rbind(res_all, anova.results)
}
colnames(res_all)[4:5]<-c('pvalue', 'adj.pvalue')
res_all<-res_all[order(res_all$pvalue),]
write.table(res_all, paste0(save_fold, 'ANOVA_', in_type, '_', inten_type, '_', imput, '_', normal, '.csv'), sep = ',', row.names = TRUE, col.names = TRUE)