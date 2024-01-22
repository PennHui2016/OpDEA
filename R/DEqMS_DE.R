

res_DEqMS<-function(mat,count, normal, imput, logT, designs){
  library(DEqMS)
  library(matrixStats)
  #source('./R/preprocessing_pro_intensity_DEqMS.R')

  if(imput=='No_MVI'){
    imput=''
  }

  if(normal=='No_normalization'){
    normal=''
  }

  if(logT=='T'){
    logT=T
  }else if(logT=='F'){
    logT=F
  }
  set.seed(123)

  prepro_res = preprocessing_raw_deqms(mat, count, designs, imput = imput, normal = normal, log2=logT)
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

    sample_names1<-designs$sample_name[grep(g1,designs$condition)]
    sample_names2<-designs$sample_name[grep(g2,designs$condition)]
    idx_g1=match(sample_names1,colnames(prepro_res$normed))
    idx_g2=match(sample_names2,colnames(prepro_res$normed))

    intens<-prepro_res$normed[,c(idx_g1, idx_g2)]
    condition = as.factor(designs[c(grep(g1, designs$condition), grep(g2, designs$condition)),]$condition)
    design = model.matrix(~0+condition) # fitting without intercept
    intens<-as.data.frame(intens)

    intens$na_g1 = apply(intens[,grep(g1,condition)],1,function(x) sum(is.na(x)))
    intens$na_g2 = apply(intens[,grep(g2,condition)],1,function(x) sum(is.na(x)))
    # Filter protein table. DEqMS require minimum two values for each group.
    #filter_idx = which(intens$na_g1<2 & intens$na_g2<2)
    filter_idx = which((length(idx_g1)-intens$na_g1)>=2 & (length(idx_g2)-intens$na_g2)>=2)
    intens.filter = intens[filter_idx,][,1:(length(intens)-2)]

    # use minimum peptide count among  samples
    # count unique+razor peptides used for quantification
    counts<-prepro_res$counts

    col_idx_pep_c<-c(match(sample_names1,colnames(counts)), match(sample_names2,colnames(counts)))
    count.table = data.frame(count = rowMins(as.matrix(counts[filter_idx,][,col_idx_pep_c])),
                             row.names = counts[,1][filter_idx])
    count.table$count = count.table$count+1

    fit1 = lmFit(intens.filter,design = design)
    cont <- makeContrasts(contrasts =  consts$conts[i], levels = design)
    fit2 = contrasts.fit(fit1,contrasts = cont)
    fit3 <- eBayes(fit2)
    fit3$count = count.table[rownames(fit3$coefficients),"count"]
    fit3$sigma[which(fit3$sigma==0)]=1e-14
    fit4 = spectraCounteBayes(fit3)
    DEqMS.results = outputResult(fit4, coef_col = 1)
    DEqMS.results<-cbind(DEqMS.results, rep(consts$conts[i], length(DEqMS.results[,1])))
    res_all<-rbind(res_all, DEqMS.results)
  }

  colnames(res_all)[length(res_all[1,])]<-'contrast'
  colnames(res_all)[c(7, 10,11)]<-c('protein', 'pvalue','adj.pvalue')
  res_all<-res_all[order(res_all$pvalue),]
  return(list(dea=res_all, processed=prepro_res$normed))
}

