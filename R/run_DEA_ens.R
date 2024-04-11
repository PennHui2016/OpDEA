run_DEA_fg_ens<-function(wfs, ens_mth, ens_opr, raw, evid, design, logFC_t, qval_t, python_path){
  ##source('./R/run_DEA_single.R')
  ##source('./R/ensemble.R')
  root_tmp = strsplit(raw,'0.tsv', fixed = T)
  temp_res_dir = paste0(root_tmp[[1]], 'res/')
  if(!dir.exists(temp_res_dir)){
    dir.create(temp_res_dir)
  }
  temp_fold = temp_res_dir
  designs = read.table(design, sep = '\t', header = TRUE)
  deas = list()
  ori_mats = list()
  processed_mats = list()
  ps=list()
  DEPs=list()
  for (i in 1:length(wfs[,1])) {
    exp = wfs[i,][1]
    norm = wfs[i,][2]
    imp = wfs[i,][3]
    dea = wfs[i,][4]

    if(imp=='missForest' | imp=='MLE'){
      imp='MinProb'
    }

    logT = 'T'
    mats = get_expression_matrix_fg_dda(raw, evid, design, exp, python_path)

    ori_mats[[i]]=mats$exp_mat

    if(dea=='limma'){
      #source('./R/limma_DE.R')
      res = res_limma(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='ROTS'){
      #source('./R/ROTS_DE.R')
      res = res_ROTS(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='DEP'){
      #source('./R/DEP_DE.R')
      res = res_DEP(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='proDA'){
      #source('./R/proDA_DE.R')
      res = res_proDA(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='DEqMS'){
      #source('./R/DEqMS_DE.R')
      counts = get_expression_matrix_fg_dda(raw, evid, design, 'count', python_path)
      res = res_DEqMS(mats$exp_mat, counts$exp_mat, norm, imp, logT, designs)
    }else if(dea=='plgem'){
      #source('./R/plgem_DE.R')
      res = res_plgem(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='edgeR'){
      #source('./R/edgeR.R')
      res = res_edgeR(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='beta_binomial'){
      #source('./R/beta_binomial.R')
      res = res_beta_binomial(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='ANOVA'){
      #source('./R/ANOVA.R')
      res = res_ANOVA(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='SAM'){
      #source('./R/SAM.R')
      res = res_SAM(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='ttest'){
      #source('./R/ttest.R')
      res = res_ttest(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='MSstats'){
      #source('./R/MSstats.R')
      res = res_MSstats(platform, raw, evid, acq, normal, imput, logT, designs)
    }
    DEA= res$dea
    deas[[i]]<-DEA
    #processed_mats[i] = res$processed
    p<-gen_figs(NULL, res$dea, logFC_t, qval_t, temp_fold, flag=paste0('E_', exp,'+N_', norm ,'+I_', imp, '+D_', dea))
    ps[[i]]<-p
    dep<-res$dea$protein[which(abs(as.numeric(res$dea$logFC))>=as.numeric(logFC_t) & as.numeric(res$dea$`adj.pvalue`)<=as.numeric(qval_t))]
    DEPs[[i]]<-dep
    write.table(res$dea, paste0(temp_fold, exp,'_DEA_res.csv'), sep = ',', col.names = T)
    write.table(mats$exp_mat, paste0(temp_fold, exp,'_original_matrix.csv'), sep = ',', col.names = T)
    write.table(res$processed, paste0(temp_fold, exp,'_preprocessed_matrix.csv'), sep = ',', col.names = T)
  }

  ens_res<-ens_multi_quant(deas, ens_mth, ens_oper)
  write.table(ens_res, paste0(temp_fold, exp,'_DEA_res.csv'), sep = ',', col.names = T)
  p_ens<-gen_figs(NULL, ens_res, logFC_t, qval_t, temp_fold, flag = '_ens_multi-quant')

  DEPs[[length(wfs[,1])+1]]=ens_res$protein[which(as.numeric(abs(ens_res$logFC))>=as.numeric(logFC_t) & as.numeric(ens_res$`adj.pvalue`)<=as.numeric(qval_t))]
  names(DEPs)<-c(wfs[,1],'ens_multi-quant')
  library(ggVennDiagram)
  library(patchwork)

  pg<-ggVennDiagram(DEPs) + scale_fill_gradient(low="gray",high = "red")
  p_all<-p_ens | pg
  zip(zipfile = paste0(temp_fold, "DEA_res.zip"), files = temp_fold, flags = '-r9Xj')
  #return(list(p=p_all, zipfile = paste0(temp_fold, "DEA_res.zip")))
  return(list(p=p_all, zipfile = paste0(temp_fold)))
}


run_DEA_mq_ens<-function(wfs, ens_mth, ens_opr, raw, evid, design, logFC_t, qval_t, python_path){
  #source('./R/run_DEA_single.R')
  #source('./R/ensemble.R')
  root_tmp = strsplit(raw,'0.txt', fixed = T)
  temp_res_dir = paste0(root_tmp[[1]], 'res/')
  if(!dir.exists(temp_res_dir)){
    dir.create(temp_res_dir)
  }
  temp_fold = temp_res_dir
  designs = read.table(design, sep = '\t', header = TRUE)
  deas = list()
  ori_mats = list()
  processed_mats = list()
  ps=list()
  DEPs=list()
  for (i in 1:length(wfs[,1])) {
    exp = wfs[i,][1]
    norm = wfs[i,][2]
    imp = wfs[i,][3]
    dea = wfs[i,][4]

    if(imp=='missForest' | imp=='MLE'){
      imp='MinProb'
    }

    logT = 'T'
    mats = get_expression_matrix_mq_dda(raw, evid, design, exp, python_path)

    ori_mats[[i]]=mats$exp_mat

    if(dea=='limma'){
      #source('./R/limma_DE.R')
      res = res_limma(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='ROTS'){
      #source('./R/ROTS_DE.R')
      res = res_ROTS(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='DEP'){
      #source('./R/DEP_DE.R')
      res = res_DEP(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='proDA'){
      #source('./R/proDA_DE.R')
      res = res_proDA(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='DEqMS'){
      #source('./R/DEqMS_DE.R')
      counts = get_expression_matrix_mq_dda(raw, evid, design, 'count', python_path)
      res = res_DEqMS(mats$exp_mat, counts$exp_mat, norm, imp, logT, designs)
    }else if(dea=='plgem'){
      #source('./R/plgem_DE.R')
      res = res_plgem(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='edgeR'){
      #source('./R/edgeR.R')
      res = res_edgeR(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='beta_binomial'){
      #source('./R/beta_binomial.R')
      res = res_beta_binomial(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='ANOVA'){
      #source('./R/ANOVA.R')
      res = res_ANOVA(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='SAM'){
      #source('./R/SAM.R')
      res = res_SAM(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='ttest'){
      #source('./R/ttest.R')
      res = res_ttest(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='MSstats'){
      #source('./R/MSstats.R')
      res = res_MSstats(platform, raw, evid, acq, normal, imput, logT, designs)
    }
    DEA= res$dea
    deas[[i]]<-DEA
    #processed_mats[i] = res$processed
    p<-gen_figs(NULL, res$dea, logFC_t, qval_t, temp_fold, flag=paste0('E_', exp,'+N_', norm ,'+I_', imp, '+D_', dea))
    ps[[i]]<-p
    dep<-res$dea$protein[which(abs(as.numeric(res$dea$logFC))>=as.numeric(logFC_t) & as.numeric(res$dea$`adj.pvalue`)<=as.numeric(qval_t))]
    DEPs[[i]]<-dep
    write.table(res$dea, paste0(temp_fold, exp,'_DEA_res.csv'), sep = ',', col.names = T)
    write.table(mats$exp_mat, paste0(temp_fold, exp,'_original_matrix.csv'), sep = ',', col.names = T)
    write.table(res$processed, paste0(temp_fold, exp,'_preprocessed_matrix.csv'), sep = ',', col.names = T)
  }

  ens_res<-ens_multi_quant(deas, ens_mth, ens_oper)
  write.table(ens_res, paste0(temp_fold, exp,'_DEA_res.csv'), sep = ',', col.names = T)
  p_ens<-gen_figs(NULL, ens_res, logFC_t, qval_t, temp_fold, flag = '_ens_multi-quant')

  DEPs[[length(wfs[,1])+1]]=ens_res$protein[which(as.numeric(abs(ens_res$logFC))>=as.numeric(logFC_t) & as.numeric(ens_res$`adj.pvalue`)<=as.numeric(qval_t))]
  names(DEPs)<-c(wfs[,1],'ens_multi-quant')
  library(ggVennDiagram)
  library(patchwork)

  pg<-ggVennDiagram(DEPs)+ scale_fill_gradient(low="gray",high = "red")
  p_all<-p_ens | pg
  zip(zipfile = paste0(temp_fold, "DEA_res.zip"), files = temp_fold, flags = '-r9Xj')
  #return(list(p=p_all, zipfile = paste0(temp_fold, "DEA_res.zip")))
  return(list(p=p_all, zipfile = paste0(temp_fold)))
}

run_DEA_diann_ens<-function(wfs, ens_mth, ens_opr, raw, evid, design, logFC_t, qval_t, python_path_diann_ens){
  #source('./R/run_DEA_single.R')
  #source('./R/ensemble.R')
  root_tmp = strsplit(evid,'0.txt', fixed = T)
  temp_res_dir = paste0(root_tmp[[1]], 'res/')
  if(!dir.exists(temp_res_dir)){
    dir.create(temp_res_dir)
  }
  temp_fold = temp_res_dir
  designs = read.table(design, sep = '\t', header = TRUE)
  deas = list()
  ori_mats = list()
  processed_mats = list()
  ps=list()
  DEPs=list()
  for (i in 1:length(wfs[,1])) {
    exp = wfs[i,][1]
    norm = wfs[i,][2]
    imp = wfs[i,][3]
    dea = wfs[i,][4]

    if(imp=='missForest' | imp=='MLE'){
      imp='MinProb'
    }

    if (exp=='dlfq'){
      logT = 'T'
    }else{
      logT = 'F'
    }
    #browser()
    mats = get_expression_matrix_diann_dia(raw, evid, design, exp, python_path_diann_ens)

    ori_mats[[i]]=mats$exp_mat

    if(dea=='limma'){
      #source('./R/limma_DE.R')
      res = res_limma(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='ROTS'){
      #source('./R/ROTS_DE.R')
      res = res_ROTS(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='DEP'){
      #source('./R/DEP_DE.R')
      res = res_DEP(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='proDA'){
      #source('./R/proDA_DE.R')
      res = res_proDA(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='DEqMS'){
      #source('./R/DEqMS_DE.R')
      counts = get_expression_matrix_fg_dda(raw, evid, design, 'count', python_path_diann_ens)
      res = res_DEqMS(mats$exp_mat, counts$exp_mat, norm, imp, logT, designs)
    }else if(dea=='plgem'){
      #source('./R/plgem_DE.R')
      res = res_plgem(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='edgeR'){
      #source('./R/edgeR.R')
      res = res_edgeR(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='beta_binomial'){
      #source('./R/beta_binomial.R')
      res = res_beta_binomial(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='ANOVA'){
      #source('./R/ANOVA.R')
      res = res_ANOVA(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='SAM'){
      #source('./R/SAM.R')
      res = res_SAM(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='ttest'){
      #source('./R/ttest.R')
      res = res_ttest(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='MSstats'){
      #source('./R/MSstats.R')
      res = res_MSstats(platform, raw, evid, acq, normal, imput, logT, designs)
    }
    DEA= res$dea
    deas[[i]]<-DEA
    #processed_mats[i] = res$processed
    p<-gen_figs(NULL, res$dea, logFC_t, qval_t, temp_fold, flag=paste0('E_', exp,'+N_', norm ,'+I_', imp, '+D_', dea))
    ps[[i]]<-p
    dep<-res$dea$protein[which(abs(as.numeric(res$dea$logFC))>=as.numeric(logFC_t) & as.numeric(res$dea$`adj.pvalue`)<=as.numeric(qval_t))]
    DEPs[[i]]<-dep
    write.table(res$dea, paste0(temp_fold, exp,'_DEA_res.csv'), sep = ',', col.names = T)
    write.table(mats$exp_mat, paste0(temp_fold, exp,'_original_matrix.csv'), sep = ',', col.names = T)
    write.table(res$processed, paste0(temp_fold, exp,'_preprocessed_matrix.csv'), sep = ',', col.names = T)
  }

  ens_res<-ens_multi_quant(deas, ens_mth, ens_oper)
  write.table(ens_res, paste0(temp_fold, exp,'_DEA_res.csv'), sep = ',', col.names = T)
  p_ens<-gen_figs(NULL, ens_res, logFC_t, qval_t, temp_fold, flag = '_ens_multi-quant')

  DEPs[[length(wfs[,1])+1]]=ens_res$protein[which(as.numeric(abs(ens_res$logFC))>=as.numeric(logFC_t) & as.numeric(ens_res$`adj.pvalue`)<=as.numeric(qval_t))]
  names(DEPs)<-c(wfs[,1],'ens_multi-quant')
  library(ggVennDiagram)
  library(patchwork)

  pg<-ggVennDiagram(DEPs)+ scale_fill_gradient(low="gray",high = "red")
  p_all<-p_ens | pg
  zip(zipfile = paste0(temp_fold, "DEA_res.zip"), files = temp_fold, flags = '-r9Xj')
  #return(list(p=p_all, zipfile = paste0(temp_fold, "DEA_res.zip")))
  return(list(p=p_all, zipfile = paste0(temp_fold)))
}


run_DEA_spt_ens<-function(wfs, ens_mth, ens_opr, raw, evid, design, logFC_t, qval_t, python_path_spt_ens){
  #source('./R/run_DEA_single.R')
  #source('./R/ensemble.R')
  root_tmp = strsplit(evid,'0.txt', fixed = T)
  temp_res_dir = paste0(root_tmp[[1]], 'res/')
  if(!dir.exists(temp_res_dir)){
    dir.create(temp_res_dir)
  }
  temp_fold = temp_res_dir
  designs = read.table(design, sep = '\t', header = TRUE)
  deas = list()
  ori_mats = list()
  processed_mats = list()
  ps=list()
  DEPs=list()
  for (i in 1:length(wfs[,1])) {
    exp = wfs[i,][1]
    norm = wfs[i,][2]
    imp = wfs[i,][3]
    dea = wfs[i,][4]

    if(imp=='missForest' | imp=='MLE'){
      imp='MinProb'
    }

    if (exp=='dlfq'){
      logT = 'T'
    }else{
      logT = 'F'
    }

    mats = get_expression_matrix_spt_dia(raw, evid, design, exp, python_path_spt_ens)

    ori_mats[[i]]=mats$exp_mat

    if(dea=='limma'){
      #source('./R/limma_DE.R')
      res = res_limma(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='ROTS'){
      #source('./R/ROTS_DE.R')
      res = res_ROTS(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='DEP'){
      #source('./R/DEP_DE.R')
      res = res_DEP(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='proDA'){
      #source('./R/proDA_DE.R')
      res = res_proDA(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='DEqMS'){
      #source('./R/DEqMS_DE.R')
      counts = get_expression_matrix_fg_dda(raw, evid, design, 'count', python_path_spt_ens)
      res = res_DEqMS(mats$exp_mat, counts$exp_mat, norm, imp, logT, designs)
    }else if(dea=='plgem'){
      #source('./R/plgem_DE.R')
      res = res_plgem(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='edgeR'){
      #source('./R/edgeR.R')
      res = res_edgeR(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='beta_binomial'){
      #source('./R/beta_binomial.R')
      res = res_beta_binomial(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='ANOVA'){
      #source('./R/ANOVA.R')
      res = res_ANOVA(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='SAM'){
      #source('./R/SAM.R')
      res = res_SAM(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='ttest'){
      #source('./R/ttest.R')
      res = res_ttest(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='MSstats'){
      #source('./R/MSstats.R')
      res = res_MSstats(platform, raw, evid, acq, normal, imput, logT, designs)
    }
    DEA= res$dea
    deas[[i]]<-DEA
    #processed_mats[i] = res$processed
    p<-gen_figs(NULL, res$dea, logFC_t, qval_t, temp_fold, flag=paste0('E_', exp,'+N_', norm ,'+I_', imp, '+D_', dea))
    ps[[i]]<-p
    dep<-res$dea$protein[which(abs(as.numeric(res$dea$logFC))>=as.numeric(logFC_t) & as.numeric(res$dea$`adj.pvalue`)<=as.numeric(qval_t))]
    DEPs[[i]]<-dep
    write.table(res$dea, paste0(temp_fold, exp,'_DEA_res.csv'), sep = ',', col.names = T)
    write.table(mats$exp_mat, paste0(temp_fold, exp,'_original_matrix.csv'), sep = ',', col.names = T)
    write.table(res$processed, paste0(temp_fold, exp,'_preprocessed_matrix.csv'), sep = ',', col.names = T)
  }

  ens_res<-ens_multi_quant(deas, ens_mth, ens_oper)
  write.table(ens_res, paste0(temp_fold, exp,'_DEA_res.csv'), sep = ',', col.names = T)
  p_ens<-gen_figs(NULL, ens_res, logFC_t, qval_t, temp_fold, flag = '_ens_multi-quant')

  DEPs[[length(wfs[,1])+1]]=ens_res$protein[which(as.numeric(abs(ens_res$logFC))>=as.numeric(logFC_t) & as.numeric(ens_res$`adj.pvalue`)<=as.numeric(qval_t))]
  names(DEPs)<-c(wfs[,1],'ens_multi-quant')
  library(ggVennDiagram)
  library(patchwork)

  pg<-ggVennDiagram(DEPs)+ scale_fill_gradient(low="gray",high = "red")
  p_all<-p_ens | pg
  zip(zipfile = paste0(temp_fold, "DEA_res.zip"), files = temp_fold, flags = '-r9Xj')
  #return(list(p=p_all, zipfile = paste0(temp_fold, "DEA_res.zip")))
  return(list(p=p_all, zipfile = paste0(temp_fold)))
}
