

res_plgem<-function(mat, normal, imput, logT, designs){
  library(plgem)
  library(MSnbase)
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
    #condition = as.factor(designs[c(grep(g1, designs$condition), grep(g2, designs$condition)),]$condition)
    condition = data.frame(as.factor(designs[c(grep(g1, designs$condition), grep(g2, designs$condition)),]$condition))
    counts<-as.data.frame(counts)

    counts$na_g1 = apply(counts[,grep(g1,condition)],1,function(x) sum(is.na(x)))
    counts$na_g2 = apply(counts[,grep(g2,condition)],1,function(x) sum(is.na(x)))
    # Filter protein table. DEqMS require minimum two values for each group.
    #filter_idx = which(counts$na_g1<2 & counts$na_g2<2)
    filter_idx = which((length(idx_g1)-counts$na_g1)>=2 & (length(idx_g2)-counts$na_g2)>=2)
    counts.filter = counts[filter_idx,][,1:(length(counts)-2)]

    counts.filter<-as.matrix(counts.filter)
    counts.filter[is.na(counts.filter)] = 0
    counts.filter[which(counts.filter<0)] = 0
    counts.filter<-counts.filter+1

    colnames(condition)<-'conditionName'
    row.names(condition)<-colnames(counts.filter)
    data_test = ExpressionSet(as.matrix(counts.filter),
                              phenoData = AnnotatedDataFrame(data=condition)
    )
    set.seed(123)
    LPSdegList <- run.plgem(esdata=data_test)

    res.plgem<-as.data.frame(row.names(counts.filter))
    res.plgem$STN<-LPSdegList$PLGEM.STN
    res.plgem$logFC<-LPSdegList$PLGEM.STN

    res.plgem$pvalue<-LPSdegList$p.value
    res.plgem$adj.pvalue<-p.adjust(LPSdegList$p.value, method = 'BH')
    res.plgem$contrast<-consts$conts[i]


    res_all<-rbind(res_all, res.plgem)
  }

  colnames(res_all)[1:5]<-c('protein', 'STN', 'logFC', 'pvalue', 'adj.pvalue')
  res_all<-res_all[order(res_all$pvalue),]
  return(list(dea=res_all, processed=prepro_res$normed))
}

