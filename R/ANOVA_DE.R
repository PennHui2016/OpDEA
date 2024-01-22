

res_ANOVA<-function(mat, normal, imput, logT, designs){
  library(limma)
  #source('./R/preprocessing_pro_intensity.R')
  if(imput=='None'){
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
  prepro_res = preprocessing_raw(mat, designs, imput = imput, normal = normal, log2=logT)
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

    ## codes from:
    #https://www.bioconductor.org/packages/release/bioc/vignettes/DEqMS/inst/doc/DEqMS-package-vignette.html#compare-p-values-from-deqms-to-ordinary-t-test-anova-and-limma
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

    anova.results$contrast<-consts$conts[i]

    res_all<-rbind(res_all, anova.results)
  }

  colnames(res_all)[4:5]<-c('pvalue', 'adj.pvalue')
  res_all<-res_all[order(res_all$pvalue),]
  return(list(dea=res_all, processed=prepro_res$normed))
}

