
ens_multi_quant<-function(deas, ens_mth, ens_oper){
  library(dplyr)
  res1 = deas[[1]]

  #res1$protein<-apply(as.array(res1$protein), 1, function(x) dnum(x))
  base.vars <- c("protein", "contrast")
  b1<-data.frame(cbind(res1$protein, res1$contrast))
  colnames(b1)<-c('protein', 'contrast')
  res.Base<-b1

  ress<-list()

  for (i in 2:length(deas)) {
    res<-deas[[i]]
    #res$protein<-apply(as.array(res$protein), 1, function(x) dnum(x))
    ress[[i]]<-res
    #b<-data.frame(cbind(apply(as.array(res$protein), 1, function(x) dnum(x)), res$contrast))
    b<-data.frame(cbind(res$protein, res$contrast))
    colnames(b)<-c('protein', 'contrast')
    res.Base<-full_join(res.Base,b,base.vars)

  }

 if(ens_mth == 'hurdle'){

   res11 <- full_join(res.Base, res1)
   z11<- -qnorm(as.numeric(res11$pvalue)/2)*sign(as.numeric(res11$logFC))

   res.Hurdle <- res.Base

   logFCs<-as.numeric(res11$logFC)
   dfs<-(!is.na(z11))

   all_z<-z11^2
   for (i in 2:length(deas)) {
     resi <- full_join(res.Base, ress[[i]],base.vars)
     zi<- -qnorm(as.numeric(resi$pvalue)/2)*sign(as.numeric(resi$logFC))
     all_z<-cbind(all_z,zi^2)
     logFCs<-cbind(logFCs, as.numeric(resi$logFC))
     dfs = dfs + (!is.na(zi))
   }

   logFCs[is.na(logFCs)]<-0

   res.Hurdle$chisq <- all_z %>% rowSums(na.rm = TRUE)
   res.Hurdle$pensemble <- 1 - pchisq(res.Hurdle$chisq, df = dfs)

   res.Hurdle$logFC_avg<-apply(logFCs,1,function(x) sign(x[which(abs(x[which(!is.na(x) & !is.infinite(x))])==max(abs(x[which(!is.na(x) & !is.infinite(x))])))[1]]) * max(abs(x[which(!is.na(x) & !is.infinite(x))]))[1])

   res.Hurdle <- res.Hurdle %>% group_by(contrast) %>% mutate(qensemble = pensemble %>%  p.adjust(method = "BH")) %>% ungroup()
   res.Hurdle <- res.Hurdle %>% arrange(pensemble)

   colnames(res.Hurdle)[c(4,5,6)]<-c('pvalue', 'logFC', 'adj.pvalue')
   res_all = res.Hurdle


 }else{
   res11 <- full_join(res.Base, res1)
   res.Hurdle <- res.Base
   logFCs<-as.numeric(res11$logFC)
   pvalues<-as.numeric(res11$pvalue)

   for (i in 2:length(deas)) {
     resi <- full_join(res.Base, ress[[i]])
     logFCs<-cbind(logFCs, as.numeric(resi$logFC))
     pvalues<-cbind(pvalues,as.numeric(resi$pvalue))
   }

   logFCs[is.na(logFCs)]<-0
   pvalues[is.na(pvalues)]<-1

   if(ens_mth=='fisher'){

   res.Hurdle$pensemble <- apply(pvalues,1, function(x) aggregation::fisher(x))

 }else if(ens_th=='set'){

   if(operation=='min'){
     res.Hurdle$pensemble <- apply(pvalues,1, function(x) min(x[which(!is.na(x))]))
   }else if(operation=='median'){
     res.Hurdle$pensemble <- apply(pvalues,1, function(x) median(x[which(!is.na(x))]))
   }else if(operation=='max'){
     res.Hurdle$pensemble <- apply(pvalues,1, function(x) max(x[which(!is.na(x))]))
   }

 }
   res.Hurdle$logFC_avg<-apply(logFCs,1,function(x) sign(x[which(abs(x[which(!is.na(x) & !is.infinite(x))])==max(abs(x[which(!is.na(x) & !is.infinite(x))])))[1]]) * max(abs(x[which(!is.na(x) & !is.infinite(x))]))[1])

   res.Hurdle <- res.Hurdle %>% group_by(contrast) %>% mutate(qensemble = pensemble %>%  p.adjust(method = "BH")) %>% ungroup()
   res.Hurdle <- res.Hurdle %>% arrange(pensemble)

   colnames(res.Hurdle)[c(3,4,5)]<-c('pvalue', 'logFC', 'adj.pvalue')
   res_all = res.Hurdle
   res_all$logFC[is.na(res_all$logFC)]=0
 }
  return(res_all)
}


