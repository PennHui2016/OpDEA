

res_edgeR<-function(mat, normal, imput, logT, designs){
  library(edgeR)
  library(dplyr)
  #source('./R/preprocessing_pro_counts.R')

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

  prepro_res = preprocessing_raw_count(mat, designs, imput = imput, normal = normal)
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

    counts<-prepro_res$normed[,c(idx_g1, idx_g2)]
    condition = as.factor(designs[c(grep(g1, designs$condition), grep(g2, designs$condition)),]$condition)

    counts<-as.data.frame(counts)

    counts$na_g1 = apply(counts[,grep(g1,condition)],1,function(x) sum(is.na(x)))
    counts$na_g2 = apply(counts[,grep(g2,condition)],1,function(x) sum(is.na(x)))
    # Filter protein table. DEqMS require minimum two values for each group.
    #filter_idx = which(counts$na_g1<2 & counts$na_g2<2)
    filter_idx = which((length(idx_g1)-counts$na_g1)>=2 & (length(idx_g2)-counts$na_g2)>=2)
    counts.filter = counts[filter_idx,][,1:(length(counts)-2)]
    ## edgeR can't process NA, use 0 instead
    counts.filter<-as.matrix(counts.filter)
    counts.filter[is.na(counts.filter)] = 0
    counts.filter[which(counts.filter<0)] = 0
    counts.filter<-counts.filter+1

    idx<-which(apply(counts.filter, 1, function(x) sd(x))!=0)
    counts.filter<-counts.filter[idx,]

    y <- DGEList(counts = counts.filter, group = condition)
    y <- calcNormFactors(y)
    design <- model.matrix(~0+condition)
    y <- estimateDisp(y, design)
    fit <- glmQLFit(y, design)

    cont <- makeContrasts(contrasts =  consts$conts[i], levels = design)

    qlf <- glmQLFTest(fit, contrast=cont)
    res.edgeR <- topTags(qlf, n = Inf)$table
    res.edgeR$contrast <- consts$conts[i]
    res.edgeR <- res.edgeR %>% cbind(protein = rownames(.),.)

    res_all<-rbind(res_all, res.edgeR)
  }

  colnames(res_all)[length(res_all[1,])]<-'contrast'
  colnames(res_all)[c(5,6)]<-c('pvalue', 'adj.pvalue')
  res_all<-res_all[order(res_all$pvalue),]
  return(list(dea=res_all, processed=prepro_res$normed))
}

