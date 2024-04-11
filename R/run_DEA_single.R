run_dlfq<-function(raw, evid, platform, python_path){
    #root_fold<-''#'inst/app/www/'
    dlfq_path = system.file("app", "www", "run_dlfq.py", package = "OpDEA")
    python_path = system.file("app", "www", "directlfq/python.exe", package = "OpDEA")
    #system(paste0('./directlfq/python', ' ', root_fold,'www/run_dlfq.py ', platform, ' ', raw, ' ', evid))
    #system(paste0('./directlfq/python', ' ', dlfq_path, ' ', platform, ' ', raw, ' ', evid))
    #system(paste0(python_path, ' ','./R/run_dlfq.py ', platform, ' ', raw, ' ', evid))
    #system(paste0(python_path, ' ',dlfq_path, ' ', platform, ' ', raw, ' ', evid))
    print(paste0(python_path, ' ', dlfq_path, ' ', platform, ' ', raw, ' ', evid))
    system(paste0(python_path, ' ', dlfq_path, ' ', platform, ' ', raw, ' ', evid))
    return(paste0(evid, '.protein_intensities.tsv'))
}

get_organism_mq<-function(protein_ids){
  Organs<-c()
  for (i in 1:length(protein_ids)) {
    if(protein_ids[i]==''){
      maj_organ='NULL'
    }else{
      proteins<-strsplit(protein_ids[i], ';',fixed = T)[[1]]
      orgs=c()
      for (pro in proteins) {

        strs<-strsplit(pro,'_', fixed = T)[[1]]
        Organ<-strs[length(strs)]
        if(strs[1]=='rev' | strs[1]=='REV'){
          Organ='decoy'
        }
        if(strs[1]=='con' | strs[1]=='CON'){
          Organ='contam'
        }

        orgs=c(orgs, Organ)
      }
      uni_og = unique(orgs)
      if(length(uni_og)>1){
        maj_organ = 'mixed'
      }else{
        maj_organ = uni_og
      }
    }
    Organs<-c(Organs,maj_organ)
  }
  return(Organs)
}

get_organism_frag<-function(protein_ids, indist_pro){
  Organs<-c()
  for (i in 1:length(protein_ids)) {
    if(length(indist_pro[i])>0){
      proteins<-c(protein_ids[i], strsplit(indist_pro[i], ';',fixed = T)[[1]])
    }else{
      proteins<-c(protein_ids[i])
    }
    orgs=c()
    for (pro in proteins) {
      strs<-strsplit(pro,'_', fixed = T)[[1]]
      Organ<-strs[length(strs)]
      if(strs[1]=='rev' | strs[1]=='REV'){
        Organ='decoy'
      }
      if(strs[1]=='con' | strs[1]=='CON'){
        Organ='contam'
      }
      orgs=c(orgs, Organ)
    }
    uni_og = unique(orgs)
    if(length(uni_og)>1){
      maj_organ = 'mixed'
    }else{
      maj_organ = uni_og
    }
    Organs<-c(Organs,maj_organ)
  }
  return(Organs)
}

get_organism_dia<-function(protein_group, protein_name){
  proteins<-c()
  for (i in 1:length(protein_group)) {
    pro = strsplit(protein_group[i], ';', fixed=T)[[1]]
    pro_name = strsplit(protein_name[i], ';', fixed=T)[[1]]
    pro_sp<-''
    for (j in 1:length(pro)) {
      if(j==1){
        if(length(grep('_UPS', pro))==1){
          pro_sp=paste0('sp|',pro[j])
        }else{
          pro_sp=paste0('sp|',pro[j],'|',pro_name[j])}
      }else{
        pro_spj=paste0('sp|',pro[j],'|',pro_name[j])
        pro_sp<-paste0(pro_sp, ';', pro_spj)
      }

    }
    proteins<-c(proteins, pro_sp)
  }

  organisms<-get_organism_mq(proteins)
  return(list(sp=proteins, orga=organisms))
}

get_expression_matrix_fg_dda<-function(raw, evid, design, exp, python_path){
  platform = 'FragPipe'
  out_exp = ''
  root_tmp = strsplit(raw,'0.tsv', fixed = T)
  temp_res_dir = paste0(root_tmp[[1]], 'res/')
  dir.create(temp_res_dir)
  if(exp=='dlfq'){
    dlfq_path<-run_dlfq(raw, evid, platform, python_path)
    dlfq_table<-read.table(dlfq_path, header = T, sep = '\t', quote = "")
    idx<-setdiff(c(1:length(colnames(dlfq_table))),
                 c(grep('protein',colnames(dlfq_table)),grep('Protein',colnames(dlfq_table))))

    organims_all<-get_organism_mq(dlfq_table$protein)
    out_dlfq<-cbind(dlfq_table$protein, organims_all,
                    dlfq_table[, idx])
    colnames(out_dlfq)[1:2]<-c('Protein','Organism')
    out_exp = out_dlfq
  }else{
  protein_table<-read.table(raw, header = T, sep = '\t', quote = "")
  inten_idx<-grep('Intensity',colnames(protein_table))
  lfq_idx<-grep('MaxLFQ.Intensity',colnames(protein_table))
  frag_inten_idx<-setdiff(inten_idx, lfq_idx)
  count_idx_a<-grep('.Spectral.Count',colnames(protein_table))
  uni_c_idx<-grep('.Unique.Spectral.Count',colnames(protein_table))
  total_c_idx<-grep('.Total.Spectral.Count',colnames(protein_table))
  cbt_c_idx<-grep('Combined.Spectral.Count',colnames(protein_table))
  count_idx<-setdiff(count_idx_a, c(uni_c_idx, total_c_idx, cbt_c_idx))

  organims_all<-get_organism_frag(protein_table$Protein,
                                  protein_table$Indistinguishable.Proteins)

  out_count_frag<-cbind(protein_table$Protein, organims_all,
                        protein_table[,count_idx])
  out_iten_frag<-cbind(protein_table$Protein, organims_all,
                       protein_table[,frag_inten_idx])
  out_inten_maxlfq<-cbind(protein_table$Protein, organims_all,
                          protein_table[,lfq_idx])
  if(exp=='count'){
    out_exp = out_count_frag
  }else if(exp=='top0'){
    out_exp = out_iten_frag
  }else if(exp=='top3'){
    out_exp = out_iten_frag
  }else if(exp=='LFQ'){
    out_exp = out_inten_maxlfq
  }
  }

  colnames(out_exp)[1:2]<-c('Protein','Organism')
  return(list(exp_mat=out_exp, res_temp=temp_res_dir))
}

get_expression_matrix_fg_tmt<-function(raw, evid, design, exp){
  platform = 'FragPipe'
  out_exp = ''
  root_tmp = strsplit(raw,'0.tsv', fixed = T)
  temp_res_dir = paste0(root_tmp[[1]], 'res/')
  dir.create(temp_res_dir)
  if(exp=='Philosopher'){

    pro_table<-read.table(raw, header = T, sep = '\t', quote = "")
    designs<-read.table(design, header = T, sep = '\t', quote = "")

    sample = designs$sample
    condition = designs$condition
    replicate = designs$replicate
    colns<-c()
    idx_col<-c()
    for (samp_i in 1:length(sample)) {
      idx_samp = grep(sample[samp_i], colnames(pro_table))
      idx_col<-c(idx_col, idx_samp)
      colns<-c(colns, sample[samp_i])
    }
    organims_all<-get_organism_frag(pro_table$Protein,
                                    pro_table$Indistinguishable.Proteins)

    out_phi<-cbind(pro_table$Protein, organims_all,
                   pro_table[,idx_col])
    colnames(out_phi)<-c('Protein', 'Organism', colns)
    out_exp = out_phi
  }else{
    pro_table<-read.table(raw, header = T, sep = '\t', quote = "")
    designs<-read.table(design, header = T, sep = '\t')

    sample = designs$sample
    condition = designs$condition
    replicate = designs$replicate
    colns<-c()
    idx_col<-c()
    for (samp_i in 1:length(sample)) {
      idx_samp = grep(sample[samp_i], colnames(pro_table))
      idx_col<-c(idx_col, idx_samp)
      colns<-c(colns, sample[samp_i])
    }
    out_abd<-cbind(pro_table$Index, pro_table$Gene, pro_table[idx_col])
    colnames(out_abd)<-c('Protein', 'Gene', colns)
    out_exp = out_abd

  }
  return(list(exp_mat=out_exp, res_temp=temp_res_dir))
}

get_expression_matrix_mq_dda<-function(raw, evid, design, exp, python_path){
  platform = 'Maxquant'
  out_exp = ''
  root_tmp = strsplit(raw,'0.txt', fixed = T)
  temp_res_dir = paste0(root_tmp[[1]], 'res/')
  dir.create(temp_res_dir)
  if(exp=='dlfq'){
    dlfq_path<-run_dlfq(raw, evid, platform, python_path)
    dlfq_table<-read.table(dlfq_path, header = T, sep = '\t', quote = "")
    idx<-setdiff(c(1:length(colnames(dlfq_table))),
                 c(grep('protein',colnames(dlfq_table)),grep('Protein',colnames(dlfq_table))))

    organims_all<-get_organism_mq(dlfq_table$protein)
    out_dlfq<-cbind(dlfq_table$protein, organims_all,
                    dlfq_table[, idx])
    colnames(out_dlfq)[1:2]<-c('Protein','Organism')
    out_exp = out_dlfq
  }else{
    protein_table = read.table(raw, header = T, sep = '\t', quote = "")
    protein_table<-protein_table[which(protein_table$Reverse!="+" & protein_table$Potential.contaminant!="+"),]

    #spectral counts
    idx_count<-grep('MS.MS.count.',colnames(protein_table))
    idx_mq_inten<-grep('Intensity.',colnames(protein_table))
    idx_lfq<-grep('LFQ.intensity.',colnames(protein_table))
    idx_top3<-grep('Top3.',colnames(protein_table))

    #organs<-get_organism_mq(protein_table$Protein.IDs)
    organims_all<-get_organism_mq(protein_table$Protein.IDs)

    out_mq_count<-cbind(protein_table$Protein.IDs,
                        organims_all,
                        protein_table[,idx_count])
    colnames(out_mq_count)[c(1:2)]<-c('Protein', 'Organism')

    out_mq_inten<-cbind(protein_table$Protein.IDs,
                        organims_all,
                        protein_table[,idx_mq_inten])
    colnames(out_mq_inten)[c(1:2)]<-c('Protein', 'Organism')

    if(length(idx_top3)>0){
    out_mq_top3<-cbind(protein_table$Protein.IDs,
                       organims_all,
                       protein_table[,idx_top3])
    colnames(out_mq_top3)[c(1:2)]<-c('Protein', 'Organism')
    }

    out_mq_lfq<-cbind(protein_table$Protein.IDs,
                      organims_all,
                      protein_table[,idx_lfq])
    colnames(out_mq_lfq)[c(1:2)]<-c('Protein', 'Organism')

    if(exp=='count'){
      out_exp = out_mq_count
    }else if(exp=='top0'){
      out_exp = out_mq_inten
    }else if(exp=='top3'){
      out_exp = out_mq_top3
    }else if(exp=='LFQ'){
      out_exp = out_mq_lfq
    }
  }
  colnames(out_exp)[1:2]<-c('Protein','Organism')
  return(list(exp_mat=out_exp, res_temp=temp_res_dir))
}

get_expression_matrix_mq_tmt<-function(raw, evid, design, exp){
  platform = 'Maxquant'
  out_exp = ''
  root_tmp = strsplit(raw,'0.tsv', fixed = T)
  temp_res_dir = paste0(root_tmp[[1]], 'res/')
  dir.create(temp_res_dir)
  designs<-read.table(design, header = T, sep = '\t', quote = "")

  sample = designs$sample
  condition = designs$condition
  replicate = designs$replicate

  # tmtIntegrator_abundance log2
  pro_table<-read.table(raw, header = T, sep = '\t', quote = "")

  idx_int<-grep('Reporter.intensity.corrected.', colnames(pro_table))

  organims_all<-get_organism_mq(pro_table$Protein.IDs)

  out_pro<-cbind(pro_table$Protein.IDs,
                 organims_all,
                 pro_table[,idx_int])
  colnames(out_pro)[c(1:2)]<-c('Protein', 'Organism')
  out_exp = out_pro
  return(list(exp_mat=out_exp, res_temp=temp_res_dir))
}

get_expression_matrix_diann_dia<-function(raw, evid, design, exp, python_path){
  platform = 'DIANN'
  out_exp = ''
  root_tmp = strsplit(evid,'0.tsv', fixed = T)
  temp_res_dir = paste0(root_tmp[[1]], 'res/')
  dir.create(temp_res_dir)
  designs<-read.table(design, header = T, sep = '\t')
  if(exp=='dlfq'){
    #browser()
    dlfq_path<-run_dlfq(raw, evid, platform, python_path)
    dlfq_table<-read.table(dlfq_path, header = T, sep = '\t', quote = "")
    colnames(dlfq_table)<-gsub('X', '', colnames(dlfq_table))

    sp_orga_dlfq = get_organism_dia(dlfq_table$Protein.Group, dlfq_table$Protein.Names)
    out_table = data.frame(Protein=sp_orga_dlfq$sp, Organism=sp_orga_dlfq$orga)
    sample_idx<-c()
    sample_name<-c()
    for (hn in 1:length(colnames(dlfq_table))) {
      idx<-grep(colnames(dlfq_table)[hn], designs$file)
      if(length(idx)==1){
        sample_idx<-c(sample_idx, hn)
        sample_name<-c(sample_name, designs$sample[idx])
      }
    }
    colnames(dlfq_table)[sample_idx]<-sample_name
    out_table<-cbind(out_table, dlfq_table[,sample_idx])

    out_exp = out_table
  }else{
    library(iq)
    #source('./R/iq-fast_MOD.R')
    if(exp!='LFQ'){
      if(exp=='top1'){
        N=1
      }else if(exp=='top3'){
        N=3
      }
      tab=process_long_format(evid,
                              output_filename = paste0(temp_res_dir, 'DIANN_', exp,'.tsv'),
                              annotation_col = c("Protein.Names", "Genes"),
                              normalization = "median",
                              filter_double_less = c("Global.Q.Value" = "0.01", "Global.PG.Q.Value" = "0.01"),
                              method='topN', N=N)
    }else{
      tab=process_long_format(evid,
                              output_filename = paste0(temp_res_dir, 'DIANN_', 'maxlfq','.tsv'),
                              annotation_col = c("Protein.Names", "Genes"),
                              normalization = "median",
                              filter_double_less = c("Global.Q.Value" = "0.01", "Global.PG.Q.Value" = "0.01"),
                              method='maxlfq', N=NULL)
    }
  tab<-as.data.frame(tab)
  sp_orga = get_organism_dia(tab$Protein.Group, tab$Protein.Names)
  out_table<-data.frame(Protein=sp_orga$sp, Organism=sp_orga$orga)

  out_data <- cbind(out_table, subset(tab, select = -c(Protein.Group, Protein.Names, Genes)))

  com<-intersect(designs$file, colnames(out_data))
  idx_col<-match(com, colnames(out_data))
  idx_des<-match(com, designs$file)
  colnames(out_data)[idx_col]<-designs$sample[idx_des]
  out_exp = out_data
  }
  colnames(out_exp)[1:2]<-c('Protein','Organism')
  return(list(exp_mat=out_exp, res_temp=temp_res_dir))
}

get_expression_matrix_spt_dia<-function(raw, evid, design, exp, python_path){
  platform = 'spt'
  out_exp = ''
  root_tmp = strsplit(evid,'0.tsv', fixed = T)
  temp_res_dir = paste0(root_tmp[[1]], 'res/')
  dir.create(temp_res_dir)
  designs<-read.table(design, header = T, sep = '\t')
  if(exp=='dlfq'){
    dlfq_path<-run_dlfq(raw, evid, platform, python_path)
    dlfq_table<-read.table(dlfq_path, header = T, sep = '\t', quote = "")
    colnames(dlfq_table)<-gsub('X', '', colnames(dlfq_table))

    out_table = data.frame(Protein=dlfq_table$protein, Organism=dlfq_table$`PG.Genes`)
    sample_idx<-c()
    sample_name<-c()
    for (hn in 1:length(colnames(dlfq_table))) {
      idx<-grep(colnames(dlfq_table)[hn], designs$file)
      if(length(idx)==1){
        sample_idx<-c(sample_idx, hn)
        sample_name<-c(sample_name, designs$sample[idx])
      }
    }
    colnames(dlfq_table)[sample_idx]<-sample_name
    out_table<-cbind(out_table, dlfq_table[,sample_idx])

    out_exp = out_table
  }else{
    library(iq)
    #source('./R/iq-fast_MOD.R')
    if(exp!='LFQ'){
      if(exp=='top1'){
        N=1
      }else if(exp=='top3'){
        N=3
      }
      tab=process_long_format(evid,
                              output_filename = paste0(temp_res_dir, 'spt_', 'topN','.tsv'),
                              sample_id  = "R.FileName",
                              primary_id = "PG.ProteinGroups",
                              secondary_id = c("EG.Library", "FG.Id", "FG.Charge", "F.FrgIon", "F.Charge", "F.FrgLossType"),
                              intensity_col = "F.PeakArea",
                              annotation_col = c("PG.Genes", "PG.ProteinNames", "PG.FastaFiles"),
                              filter_string_equal = c("F.ExcludedFromQuantification" = "False"),
                              filter_double_less = c("PG.Qvalue" = "0.01", "EG.Qvalue" = "0.01"),
                              log2_intensity_cutoff = 0,
                              normalization = "median",
                              method='topN', N=N)
    }else{
      tab=process_long_format(evid,
                              output_filename = paste0(temp_res_dir, 'spt_', 'maxlfq','.tsv'),
                              sample_id  = "R.FileName",
                              primary_id = "PG.ProteinGroups",
                              secondary_id = c("EG.Library", "FG.Id", "FG.Charge", "F.FrgIon", "F.Charge", "F.FrgLossType"),
                              intensity_col = "F.PeakArea",
                              annotation_col = c("PG.Genes", "PG.ProteinNames", "PG.FastaFiles"),
                              filter_string_equal = c("F.ExcludedFromQuantification" = "False"),
                              filter_double_less = c("PG.Qvalue" = "0.01", "EG.Qvalue" = "0.01"),
                              log2_intensity_cutoff = 0,
                              normalization = "median",
                              method='maxlfq', N=NULL)
    }
    tab<-as.data.frame(tab)
    sp_orga = get_organism_dia(tab$PG.ProteinGroups, tab$PG.ProteinNames)
    out_table<-data.frame(Protein=sp_orga$sp, Organism=sp_orga$orga)

    out_data <- cbind(out_table, subset(tab, select = -c(PG.ProteinGroups, PG.ProteinNames, PG.Genes, PG.FastaFiles)))

    com<-intersect(designs$file, colnames(out_data))
    idx_col<-match(com, colnames(out_data))
    idx_des<-match(com, designs$file)
    colnames(out_data)[idx_col]<-designs$sample[idx_des]
    out_exp = out_data
  }
  colnames(out_exp)[1:2]<-c('Protein','Organism')
  return(list(exp_mat=out_exp, res_temp=temp_res_dir))
}

gen_figs<-function(exp_matrix, dea_res, logFC_t, qval, temp_fold, flag=''){
  library(ggplot2)
  library(ggrepel)
  library(reshape2)
  library(ggpubr)
  logFC_t<-as.numeric(logFC_t)
  qval<-as.numeric(qval)
  melt_dt<-as.data.frame(dea_res)
  #colnames(melt_anova_ps_dat)<-c('workflow','pvalue','rank','setting','method')
  melt_dt$logFC<-as.numeric(melt_dt$logFC)
  melt_dt$log10q<--log10(as.numeric(melt_dt$`adj.pvalue`))
  melt_dt$stable<-as.factor(ifelse((as.numeric(melt_dt$`adj.pvalue`) <= qval & abs(as.numeric(melt_dt$logFC)) >= logFC_t),
                                             'DEP',
                                             'non-DEP'))
  melt_dt$label<-melt_dt$protein
  melt_dt$label[which(melt_dt$stable=='non-DEP')]=''


  p<-ggplot(melt_dt, aes(x=logFC, y=log10q,color=stable)) +
    geom_point(alpha=0.5, size=2) +
    theme_bw(base_size = 12) +
    #scale_shape_manual(values=c(1, 3,4, 8, 13, 24))+
    xlab("log2(fold change)") +
    ylab("-Log10(adj.pvalue)") +
    #theme(plot.title = element_text(size=15,hjust = 0.5)) +
    scale_colour_manual(values = c("purple",'gray')) +
    geom_hline(yintercept = -log10(qval), lty = 4) +
    geom_vline(xintercept = c(-logFC_t, logFC_t), lty = 4)+
    labs(title = flag)+
    geom_label_repel(data = melt_dt, aes(label = label),
                     size = 3,box.padding = unit(0.5, "lines"),
                     point.padding = unit(0.8, "lines"),
                     segment.color = "black",
                     show.legend = FALSE, max.overlaps = 10)+
    theme(legend.title=element_text(size=12),legend.text=element_text(size=12), legend.position = "top")+
    guides(color=guide_legend(nrow=3), shape=guide_legend(nrow=3))

  p
  ggsave(paste0(temp_fold, "volcano", flag,".pdf"), width = 10, height = 10, units = "cm")
  return(p)
}

#platform, acq, raw_file, evid_path, design_file, mat, norm, imp, dea
run_DEA<-function(platform, acq, raw, evid, design, exp, norm, imp, dea, logFC, qval, python_path){
  write.table(data.frame(a=c('3','run_DEA')), "log.txt", sep="\t", col.names=FALSE, row.names=FALSE, append=TRUE)
  if(platform=='FragPipe' & acq=='DDA'){
    logT = 'T'
    mats = get_expression_matrix_fg_dda(raw, evid, design, exp, python_path)
    designs = read.table(design, sep = '\t', header = TRUE)
    if(dea=='limma'){
      #source('./R/limma_DE.R')
      write.table(data.frame(a=c('4','res_limma')), "log.txt", sep="\t", col.names=FALSE, row.names=FALSE, append=TRUE)
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
      #source('./R/ANOVA_DE.R')
      res = res_ANOVA(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='SAM'){
      #source('./R/SAM_DE.R')
      res = res_SAM(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='ttest'){
      #source('./R/ttest_DE.R')
      res = res_ttest(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='MSstats'){
      #source('./R/MSstats_DE.R')
      res = res_MSstats(platform, raw, evid, acq, normal, imput, logT, designs)
    }
    write.table(res$dea, paste0(mats$res_temp,'DEA_res.csv'), sep = ',', col.names = T)
    write.table(mats$exp_mat, paste0(mats$res_temp,'original_matrix.csv'), sep = ',', col.names = T)
    write.table(res$processed, paste0(mats$res_temp,'preprocessed_matrix.csv'), sep = ',', col.names = T)
    p<-gen_figs(mats$exp_mat, res$dea, logFC, qval, mats$res_temp, flag=paste0('E_', exp,'+N_', norm ,'+I_', imp, '+D_', dea))
    zip(zipfile = paste0(mats$res_temp, "DEA_res.zip"), files = mats$res_temp, flags = '-r9Xj')
    #return(list(p=p, zipfile = paste0(mats$res_temp, "DEA_res.zip")))
    return(list(p=p, zipfile = paste0(mats$res_temp)))

  }else if(platform=='Maxquant' & acq=='DDA'){
    logT = 'T'
    mats = get_expression_matrix_mq_dda(raw, evid, design, exp, python_path)
    designs = read.table(design, sep = '\t', header = TRUE)
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
      #source('./R/ANOVA_DE.R')
      res = res_ANOVA(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='SAM'){
      #source('./R/SAM_DE.R')
      res = res_SAM(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='ttest'){
      #source('./R/ttest_DE.R')
      res = res_ttest(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='MSstats'){
      #source('./R/MSstats_DE.R')
      res = res_MSstats(platform, raw, evid, acq, normal, imput, logT, designs)
    }
    write.table(res$dea, paste0(mats$res_temp,'DEA_res.csv'), sep = ',', col.names = T)
    write.table(mats$exp_mat, paste0(mats$res_temp,'original_matrix.csv'), sep = ',', col.names = T)
    write.table(res$processed, paste0(mats$res_temp,'preprocessed_matrix.csv'), sep = ',', col.names = T)
    p<-gen_figs(mats$exp_mat, res$dea, logFC, qval, mats$res_temp, flag=paste0('E_', exp,'+N_', norm ,'+I_', imp, '+D_', dea))
    zip(zipfile = paste0(mats$res_temp, "DEA_res.zip"), files = mats$res_temp, flags = '-r9Xj')
    #return(list(p=p, zipfile = paste0(mats$res_temp, "DEA_res.zip")))
    return(list(p=p, zipfile = paste0(mats$res_temp)))
  }else if(platform=='DIANN' & acq=='DIA'){
    if(exp=='dlfq'){
      logT = 'T'
    }else{
      logT = 'F'
    }

    mats = get_expression_matrix_diann_dia(raw, evid, design, exp, python_path)
    designs = read.table(design, sep = '\t', header = TRUE)
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
    }else if(dea=='ANOVA'){
      #source('./R/ANOVA_DE.R')
      res = res_ANOVA(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='SAM'){
      #source('./R/SAM_DE.R')
      res = res_SAM(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='ttest'){
      #source('./R/ttest_DE.R')
      res = res_ttest(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='MSstats'){
      #source('./R/MSstats_DE.R')
      res = res_MSstats(platform, raw, evid, acq, normal, imput, logT, designs)
    }
    write.table(res$dea, paste0(mats$res_temp,'DEA_res.csv'), sep = ',', col.names = T)
    write.table(mats$exp_mat, paste0(mats$res_temp,'original_matrix.csv'), sep = ',', col.names = T)
    write.table(res$processed, paste0(mats$res_temp,'preprocessed_matrix.csv'), sep = ',', col.names = T)
    p<-gen_figs(mats$exp_mat, res$dea, logFC, qval, mats$res_temp, flag=paste0('E_', exp,'+N_', norm ,'+I_', imp, '+D_', dea))
    zip(zipfile = paste0(mats$res_temp, "DEA_res.zip"), files = mats$res_temp, flags = '-r9Xj')
    #return(list(p=p, zipfile = paste0(mats$res_temp, "DEA_res.zip")))
    return(list(p=p, zipfile = paste0(mats$res_temp)))

  }else if(platform=='spt' & acq=='DIA'){
    if(exp=='dlfq'){
      logT = 'T'
    }else{
      logT = 'F'
    }
    mats = get_expression_matrix_spt_dia(raw, evid, design, exp, python_path)
    designs = read.table(design, sep = '\t', header = TRUE)
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
    }else if(dea=='ANOVA'){
      #source('./R/ANOVA_DE.R')
      res = res_ANOVA(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='SAM'){
      #source('./R/SAM_DE.R')
      res = res_SAM(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='ttest'){
      #source('./R/ttest_DE.R')
      res = res_ttest(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='MSstats'){
      #source('./R/MSstats_DE.R')
      res = res_MSstats(platform, raw, evid, acq, normal, imput, logT, designs)
    }
    write.table(res$dea, paste0(mats$res_temp,'DEA_res.csv'), sep = ',', col.names = T)
    write.table(mats$exp_mat, paste0(mats$res_temp,'original_matrix.csv'), sep = ',', col.names = T)
    write.table(res$processed, paste0(mats$res_temp,'preprocessed_matrix.csv'), sep = ',', col.names = T)
    p<-gen_figs(mats$exp_mat, res$dea, logFC, qval, mats$res_temp, flag=paste0('E_', exp,'+N_', norm ,'+I_', imp, '+D_', dea))
    zip(zipfile = paste0(mats$res_temp, "DEA_res.zip"), files = mats$res_temp, flags = '-r9Xj')
    #return(list(p=p, zipfile = paste0(mats$res_temp, "DEA_res.zip")))
    return(list(p=p, zipfile = paste0(mats$res_temp)))

  }else if(platform=='FragPipe' & acq=='TMT'){
    logT = 'T'
    mats = get_expression_matrix_fg_tmt(raw, evid, design, exp)
    designs = read.table(design, sep = '\t', header = TRUE)
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
    }else if(dea=='ANOVA'){
      #source('./R/ANOVA_DE.R')
      res = res_ANOVA(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='SAM'){
      #source('./R/SAM_DE.R')
      res = res_SAM(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='ttest'){
      #source('./R/ttest_DE.R')
      res = res_ttest(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='MSstats'){
      #source('./R/MSstats_DE.R')
      res = res_MSstats(platform, raw, evid, acq, normal, imput, logT, designs)
    }
    write.table(res$dea, paste0(mats$res_temp,'DEA_res.csv'), sep = ',', col.names = T)
    write.table(mats$exp_mat, paste0(mats$res_temp,'original_matrix.csv'), sep = ',', col.names = T)
    write.table(res$processed, paste0(mats$res_temp,'preprocessed_matrix.csv'), sep = ',', col.names = T)
    p<-gen_figs(mats$exp_mat, res$dea, logFC, qval, mats$res_temp, flag=paste0('E_', exp,'+N_', norm ,'+I_', imp, '+D_', dea))
    zip(zipfile = paste0(mats$res_temp, "DEA_res.zip"), files = mats$res_temp, flags = '-r9Xj')
    #return(list(p=p, zipfile = paste0(mats$res_temp, "DEA_res.zip")))
    return(list(p=p, zipfile = paste0(mats$res_temp)))

  }else if(platform=='Maxquant' & acq=='TMT'){
    logT = 'T'
    mats = get_expression_matrix_mq_tmt(raw, evid, design, exp)
    designs = read.table(design, sep = '\t', header = TRUE)
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
    }else if(dea=='ANOVA'){
      #source('./R/ANOVA_DE.R')
      res = res_ANOVA(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='SAM'){
      #source('./R/SAM_DE.R')
      res = res_SAM(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='ttest'){
      #source('./R/ttest_DE.R')
      res = res_ttest(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='MSstats'){
      #source('./R/MSstats_DE.R')
      res = res_MSstats(platform, raw, evid, acq, normal, imput, logT, designs)
    }
    write.table(res$dea, paste0(mats$res_temp,'DEA_res.csv'), sep = ',', col.names = T)
    write.table(mats$exp_mat, paste0(mats$res_temp,'original_matrix.csv'), sep = ',', col.names = T)
    write.table(res$processed, paste0(mats$res_temp,'preprocessed_matrix.csv'), sep = ',', col.names = T)
    p<-gen_figs(mats$exp_mat, res$dea, logFC, qval, mats$res_temp, flag=paste0('E_', exp,'+N_', norm ,'+I_', imp, '+D_', dea))
    zip(zipfile = paste0(mats$res_temp, "DEA_res.zip"), files = mats$res_temp, flags = '-r9Xj')
    #return(list(p=p, zipfile = paste0(mats$res_temp, "DEA_res.zip")))
    return(list(p=p, zipfile = paste0(mats$res_temp)))
  }

}

run_DEA_TMT_fg<-function(platform, acq, abd, ratio, phi, design, exp, norm, imp, dea, logFC, qval){
  if(platform=='FragPipe' & acq=='TMT'){
    logT = 'T'
    if(exp=='abd'){
      raw=abd
      evid=phi
      logT = 'F'
    }else if(exp=='ratio'){
      raw=ratio
      evid=phi
      logT = 'F'
    }
    mats = get_expression_matrix_fg_tmt(raw, evid, design, exp)

    designs = read.table(design, sep = '\t', header = TRUE)
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
    }else if(dea=='ANOVA'){
      #source('./R/ANOVA_DE.R')
      res = res_ANOVA(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='SAM'){
      #source('./R/SAM_DE.R')
      res = res_SAM(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='ttest'){
      #source('./R/ttest_DE.R')
      res = res_ttest(mats$exp_mat, norm, imp, logT, designs)
    }else if(dea=='MSstats'){
      #source('./R/MSstats_DE.R')
      res = res_MSstats(platform, raw, evid, acq, normal, imput, logT, designs)
    }
    write.table(res$dea, paste0(mats$res_temp,'DEA_res.csv'), sep = ',', col.names = T)
    write.table(mats$exp_mat, paste0(mats$res_temp,'original_matrix.csv'), sep = ',', col.names = T)
    write.table(res$processed, paste0(mats$res_temp,'preprocessed_matrix.csv'), sep = ',', col.names = T)
    p<-gen_figs(mats$exp_mat, res$dea, logFC, qval, mats$res_temp, flag=paste0('E_', exp,'+N_', norm ,'+I_', imp, '+D_', dea))
    zip(zipfile = paste0(mats$res_temp, "DEA_res.zip"), files = mats$res_temp, flags = '-r9Xj')
    #return(list(p=p, zipfile = paste0(mats$res_temp, "DEA_res.zip")))
    return(list(p=p, zipfile = paste0(mats$res_temp)))
  }
}

run_DEA_fg_DDA<-function(input, mat, norm, imp, dea, raw_file, evid_path, design_file, logFC, qval,python_path){
  platform = 'FragPipe'
  acq = 'DDA'

  # mat = data_fg_sug_prefer()$expression_matrix[1]
  # norm = data_fg_sug_prefer()$normalization[1]
  # imp = data_fg_sug_prefer()$imputation[1]
  # dea = data_fg_sug_prefer()$DEA_tool[1]
  # logFC=input$logFC_fg1
  # qval = input$adjp_fg1
  DEA_res<-run_DEA(platform, acq, raw_file, evid_path, design_file, mat, norm, imp, dea, logFC, qval, python_path)
  return(DEA_res)
}

run_DEA_mq_DDA<-function(input, mat, norm, imp, dea, raw_file, evid_path, design_file, logFC, qval, python_path_mq){
  platform = 'Maxquant'
  acq = 'DDA'
  # raw_file<-input$upload_mq_raw$datapath
  # evid_path<-input$upload_mq_evid$datapath
  # design_file<-input$upload_mq_dg$datapath
  # mat = data_mq_sug_prefer()$expression_matrix[1]
  # norm = data_mq_sug_prefer()$normalization[1]
  # imp = data_mq_sug_prefer()$imputation[1]
  # dea = data_mq_sug_prefer()$DEA_tool[1]
  DEA_res<-run_DEA(platform, acq, raw_file, evid_path, design_file, mat, norm, imp, dea, logFC, qval, python_path_mq)
  return(DEA_res)
}

run_DEA_diann_DIA<-function(input, mat, norm, imp, dea, raw_file, evid_path, design_file, logFC, qval, python_path_diann){
  platform = 'DIANN'
  acq = 'DIA'
  # raw_file<-''
  # evid_path<-input$upload_diann_evid$datapath
  # design_file<-input$upload_diann_dg$datapath
  # mat = data_diann_sug_prefer()$expression_matrix[1]
  # norm = data_diann_sug_prefer()$normalization[1]
  # imp = data_diann_sug_prefer()$imputation[1]
  # dea = data_diann_sug_prefer()$DEA_tool[1]
  DEA_res<-run_DEA(platform, acq, raw_file, evid_path, design_file, mat, norm, imp, dea, logFC, qval, python_path_diann)
  return(DEA_res)
}

run_DEA_spt_DIA<-function(input, mat, norm, imp, dea, raw_file, evid_path, design_file, logFC, qval, python_path_spt){
  platform = 'spt'
  acq = 'DIA'
  # raw_file<-''
  # evid_path<-input$upload_spt_evid$datapath
  # design_file<-input$upload_spt_dg$datapath
  # mat = data_spt_sug_prefer()$expression_matrix[1]
  # norm = data_spt_sug_prefer()$normalization[1]
  # imp = data_spt_sug_prefer()$imputation[1]
  # dea = data_spt_sug_prefer()$DEA_tool[1]
  DEA_res<-run_DEA(platform, acq, raw_file, evid_path, design_file, mat, norm, imp, dea, logFC, qval, python_path_spt)
  return(DEA_res)
}

run_DEA_fg_TMT<-function(input, mat, norm, imp, dea, abd_file, ratio_path, phi_path, design_file, logFC, qval){
  platform = 'FragPipe'
  acq = 'TMT'
  # abd_file<-input$upload_fg_tmt_abd$datapath
  # ratio_path<-input$upload_fg_tmt_rat$datapath
  # phi_path<-input$upload_fg_tmt_phi$datapath
  # design_file<-input$upload_fg_tmt_dg$datapath
  # mat = data_spt_sug_prefer()$expression_matrix[1]
  # norm = data_spt_sug_prefer()$normalization[1]
  # imp = data_spt_sug_prefer()$imputation[1]
  # dea = data_spt_sug_prefer()$DEA_tool[1]
  DEA_res<-run_DEA_TMT_fg(platform, acq, abd_file, ratio_path, phi_path, design_file, mat, norm, imp, dea, logFC, qval)
  return(DEA_res)
}

run_DEA_mq_TMT<-function(input, mat, norm, imp, dea, raw_file, evid_path, design_file, logFC, qval){
  platform = 'Maxquant'
  acq = 'TMT'
  # raw_file<-input$upload_tmt_mq_raw$datapath
  # evid_path<-input$upload_tmt_mq_evid$datapath
  # design_file<-input$upload_tmt_mq_dg$datapath
  # mat = data_spt_sug_prefer()$expression_matrix[1]
  # norm = data_spt_sug_prefer()$normalization[1]
  # imp = data_spt_sug_prefer()$imputation[1]
  # dea = data_spt_sug_prefer()$DEA_tool[1]
  DEA_res<-run_DEA(platform, acq, raw_file, evid_path, design_file, mat, norm, imp, dea, logFC, qval)
  return(DEA_res)
}
