###########
##  process outputs from different quntification platforms
######

## 1.software specific
## fragpipe:: intensity (topN 0); maxquant: intensity (topN 0); DIANN: topN1; spectranaut: 


## 2. top3:: fragpipe rerun; maxquant: output; DIANN: iq-mod; spectranaut: iq-mod


## 3. maxlfq:: fragpipe: output; maxquant: output; DIANN: iq-mod; spectranaut: iq-mod


## 4. directlfq::  fragpipe: directlfq; maxquant: directlfq; DIANN: directlfq; spectranaut: directlfq


## 5. spectral counts:: fragpipe: output; maxquant: output;

save_fold_DDA_frag<-'D:/data/benchmark/data/DDA/FragPipe/'
save_fold_DDA_mq<-'D:/data/benchmark/data/DDA/Maxquant/'
save_fold_DIA_DIANN<-'D:/data/benchmark/data/DIA/DIANN/'
save_fold_DIA_spt<-'D:/data/benchmark/data/DIA/Spectronaunt/'
save_fold_tmt_frag<-'D:/data/benchmark/data/TMT/FragPipe/'
save_fold_tmt_mq<-'D:/data/benchmark/data/TMT/Maxquant/'
dataset_info_path<-'D:/data/benchmark/data/dataset_info/'

get_uniport_id_UPS<-function(){
  library(seqinr)
  ups_pro<-read.fasta('E:/MS_data/PXD002099/ups1-ups2-sequences.fasta',as.string = TRUE, seqtype = "AA")
  all_ups_pro<-names(ups_pro)
  uniport_id<-c()
  for (ups in all_ups_pro) {
    uniport_id<-c(uniport_id,strsplit(ups,'|',fixed = T)[[1]][1])
  }
  uniport_id<-gsub('ups','',uniport_id)
  return(uniport_id)
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


get_idx_True_Organism<-function(organisms, true_organism){
  idx=c()
  for (org in true_organism) {
    idx_org<-grep(org, organisms)
    idx<-c(idx, idx_org)
  }
  return(idx)
}

get_allType_matrix_frag_dda<-function(dataset_info){
  for(i in 1:length(dataset_info$dataset)){
    datasetname<-dataset_info$dataset[i]
    true_organism <- strsplit(dataset_info$organism[i], ';', fixed = T)[[1]]
    de_organism<- strsplit(dataset_info$trueDE[i], ';', fixed = T)[[1]]
  protein_file<-paste0(dataset_info$output_folder[i],'TOP0/', dataset_info$main_output_protein_file[i])
  protein_table<-read.table(protein_file, header = T, sep = '\t', quote = "")
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
  idx_true_organisms<-get_idx_True_Organism(organims_all, true_organism)
  out_count_frag<-cbind(protein_table$Protein[idx_true_organisms], organims_all[idx_true_organisms], 
                        protein_table[,count_idx][idx_true_organisms,])
  out_iten_frag<-cbind(protein_table$Protein[idx_true_organisms], organims_all[idx_true_organisms], 
                       protein_table[,frag_inten_idx][idx_true_organisms,])
  out_inten_maxlfq<-cbind(protein_table$Protein[idx_true_organisms], organims_all[idx_true_organisms],
                          protein_table[,lfq_idx][idx_true_organisms,])
  
  colnames(out_iten_frag)[c(1,2)]<-c('Protein', 'Organism')
  colnames(out_inten_maxlfq)[c(1,2)]<-c('Protein', 'Organism')
  colnames(out_count_frag)[c(1,2)]<-c('Protein', 'Organism')
  write.table(out_iten_frag, paste0(save_fold_DDA_frag, datasetname,'_FragPipe_pro_intensity.tsv'),col.names = T,row.names = F, sep = '\t')
  write.table(out_inten_maxlfq, paste0(save_fold_DDA_frag, datasetname,'_FragPipe_pro_maxlfq.tsv'),col.names = T,row.names = F, sep = '\t')
  write.table(out_count_frag, paste0(save_fold_DDA_frag, datasetname,'_FragPipe_pro_count.tsv'),col.names = T,row.names = F, sep = '\t')
  
  #top3
  protein_file<-paste0(dataset_info$output_folder[i],'TOP3/', dataset_info$main_output_protein_file[i])
  protein_table<-read.table(protein_file, header = T, sep = '\t', quote = "")
  inten_idx<-grep('Intensity',colnames(protein_table))
  lfq_idx<-grep('MaxLFQ.Intensity',colnames(protein_table))
  frag_inten_idx<-setdiff(inten_idx, lfq_idx)
  organims_all<-get_organism_frag(protein_table$Protein,
                                  protein_table$Indistinguishable.Proteins)
  idx_true_organisms<-get_idx_True_Organism(organims_all, true_organism)
  out_iten_frag<-cbind(protein_table$Protein[idx_true_organisms], organims_all[idx_true_organisms],
                       protein_table[,frag_inten_idx][idx_true_organisms,])
  colnames(out_iten_frag)[c(1,2)]<-c('Protein', 'Organism')
  write.table(out_iten_frag, paste0(save_fold_DDA_frag, datasetname,'_FragPipe_top3_pro_intensity.tsv'),col.names = T,row.names = F, sep = '\t')
  
  #directlfq
  dlfq_file<-paste0(dataset_info$output_folder[i],'TOP3/', dataset_info$directlfq_file[i])
  dlfq_table<-read.table(dlfq_file, header = T, sep = '\t', quote = "")
  
  idx<-setdiff(c(1:length(colnames(dlfq_table))), 
               c(grep('protein',colnames(dlfq_table)),grep('Protein',colnames(dlfq_table))))
  
  organims_all<-get_organism_mq(dlfq_table$protein)
  idx_true_organisms<-get_idx_True_Organism(organims_all, true_organism)
  
  out_dlfq<-cbind(dlfq_table$protein[idx_true_organisms], organims_all[idx_true_organisms],
                  dlfq_table[, idx][idx_true_organisms,])
  colnames(out_dlfq)[c(1,2)]<-c('Protein', 'Organism')
  write.table(out_dlfq, paste0(save_fold_DDA_frag, datasetname,'_FragPipe_dlfq_pro_intensity.tsv'),col.names = T,row.names = F, sep = '\t')
  
  design_file<-paste0(dataset_info$output_folder[i], 'TOP3/',dataset_info$design_file[i])
  design_table<-read.table(design_file, sep = '\t', header = T, quote = "")
  write.table(design_table, paste0(save_fold_DDA_frag, datasetname,'_FragPipe_design.tsv'),col.names = T,row.names = F, sep = '\t')
  
  mstats_design<-data.frame(cbind(design_table$condition, design_table$replicate))
  colnames(mstats_design)<-c('Condition','BioReplicate')
  mstats_design$IsotopeLabelType<-c(rep('L',length(design_table$condition)))
  
  file_names<-c()
  for (fl in design_table$file) {
    flsp=strsplit(fl,'\\', fixed=T)[[1]]
    flc=flsp[length(flsp)]
    flc=strsplit(flc,'.',fixed = T)[[1]][1]
    file_names<-c(file_names,flc)
  }
  mstats_design$Raw.file=file_names
  mstats_design$Run=file_names
  write.table(mstats_design, paste0(save_fold_DDA_frag, datasetname,'_FragPipe_design_msstats.tsv'),col.names = T,row.names = F, sep = '\t')

  
  peptide_file<-paste0(dataset_info$output_folder[i],'TOP0/', dataset_info$main_output_peptide_file[i])
  peptide_table<-read.table(peptide_file, header = T, sep = '\t', quote = "")
  inten_idx<-grep('.Intensity',colnames(peptide_table))
  lfq_idx<-grep('.MaxLFQ.Intensity',colnames(peptide_table))
  frag_inten_idx<-setdiff(inten_idx, lfq_idx)
  
  organims_all<-get_organism_mq(peptide_table$Protein)
  idx_true_organisms<-get_idx_True_Organism(organims_all, true_organism)
  out_iten_frag<-cbind(peptide_table$Peptide.Sequence[idx_true_organisms], 
                       peptide_table$Protein[idx_true_organisms],
                       organims_all[idx_true_organisms],
                       peptide_table[,frag_inten_idx][idx_true_organisms,])
  out_inten_maxlfq<-cbind(peptide_table$Peptide.Sequence[idx_true_organisms], 
                          peptide_table$Protein[idx_true_organisms],
                          organims_all[idx_true_organisms],
                          peptide_table[,lfq_idx][idx_true_organisms,])
  
  colnames(out_iten_frag)[c(1,2,3)]<-c('Peptide Sequence', 'Protein', 'Organism')
  colnames(out_inten_maxlfq)[c(1,2,3)]<-c('Peptide Sequence', 'Protein', 'Organism')

  write.table(out_iten_frag, paste0(save_fold_DDA_frag, datasetname,'_FragPipe_pep_intensity.tsv'),col.names = T,row.names = F, sep = '\t')
  write.table(out_inten_maxlfq, paste0(save_fold_DDA_frag, datasetname,'_FragPipe_pep_maxlfq.tsv'),col.names = T,row.names = F, sep = '\t')
  
  all_proteins<-rbind(cbind(protein_table$Protein, get_organism_frag(protein_table$Protein,
                                                                     protein_table$Indistinguishable.Proteins)),
                      cbind(dlfq_table$protein, get_organism_mq(dlfq_table$protein)),
                      cbind(peptide_table$Protein, get_organism_mq(peptide_table$Protein)))
  idx_pro<-c()
  for (org in true_organism) {
    idx_org = grep(org, all_proteins[,2])
    idx_pro<-c(idx_pro, idx_org)
  }
  all_proteins<-all_proteins[idx_pro,]
  labs<-c(rep(0, length(all_proteins[,1])))
  idx_pos<-c()
  for (de_org in de_organism) {
    idx_org = grep(org, all_proteins[,2])
    if(de_org=='UPS'){#contam added by FragPipe are overlapped with human
      idx_human = grep("HUMAN", all_proteins[,2])
      uniport_UPS<-get_uniport_id_UPS()
      for (ups in uniport_UPS) {
        idx_ups<-grep(ups,all_proteins[,1])
        idx_org<-c(idx_org, idx_ups)
      }
    }
    idx_pos<-c(idx_pos, idx_org)
  }
  labs[idx_pos] = 1
  all_proteins = cbind(all_proteins, labs)
  
  colnames(all_proteins)<-c('Protein', 'Organism', 'DEP')
  
  uni_pro = unique(all_proteins[,1])
  idx_uni<-match(uni_pro, all_proteins[,1])
  
  all_proteins<-all_proteins[idx_uni,]

  #all_proteins<-as.data.frame(all_proteins)
  write.table(all_proteins, paste0(save_fold_DDA_frag, datasetname,'_FragPipe_all_proteins.tsv'),col.names = T,row.names = F, sep = '\t')
  }
}

#dataset_info<-read.table(paste0(dataset_info_path, 'DDA_frag.txt'),
#                         header = T, sep = '\t')

#get_allType_matrix_frag_dda(dataset_info)

get_allType_matrix_mq_dda<-function(dataset_info){
  for(i in 1:length(dataset_info$dataset)){
    datasetname<-dataset_info$dataset[i]
    true_organism <- strsplit(dataset_info$organism[i], ';', fixed = T)[[1]]
    de_organism<- strsplit(dataset_info$trueDE[i], ';', fixed = T)[[1]]
    protein_file<-paste0(dataset_info$output_folder[i], dataset_info$main_output_protein_file[i])
    protein_table<-read.table(protein_file, header = T, sep = '\t', quote = "")
    protein_table<-protein_table[which(protein_table$Reverse!="+" & protein_table$Potential.contaminant!="+"),]
    
    #spectral counts
    idx_count<-grep('MS.MS.count.',colnames(protein_table))
    idx_mq_inten<-grep('Intensity.',colnames(protein_table))
    idx_lfq<-grep('LFQ.intensity.',colnames(protein_table))
    idx_top3<-grep('Top3.',colnames(protein_table))
    
    #organs<-get_organism_mq(protein_table$Protein.IDs)
    organims_all<-get_organism_mq(protein_table$Protein.IDs)
    idx_true_organisms<-get_idx_True_Organism(organims_all, true_organism)
    #protein_table<-protein_table[which(organs!='decoy' & organs!='contam'),]
    
    out_mq_count<-cbind(protein_table$Protein.IDs[idx_true_organisms], 
                        organims_all[idx_true_organisms],
                        protein_table[,idx_count][idx_true_organisms,])
    colnames(out_mq_count)[c(1:2)]<-c('Protein', 'Organism')
    write.table(out_mq_count, paste0(save_fold_DDA_mq, datasetname, '_Maxquant_pro_count.tsv'), sep = '\t', col.names = T, row.names = F)
    
    out_mq_inten<-cbind(protein_table$Protein.IDs[idx_true_organisms], 
                        organims_all[idx_true_organisms], 
                        protein_table[,idx_mq_inten][idx_true_organisms,])
    colnames(out_mq_inten)[c(1:2)]<-c('Protein', 'Organism')
    write.table(out_mq_inten, paste0(save_fold_DDA_mq, datasetname, '_Maxquant_pro_intensity.tsv'), sep = '\t', col.names = T, row.names = F)
    
    out_mq_top3<-cbind(protein_table$Protein.IDs[idx_true_organisms], 
                       organims_all[idx_true_organisms],
                       protein_table[,idx_top3][idx_true_organisms,])
    colnames(out_mq_top3)[c(1:2)]<-c('Protein', 'Organism')
    write.table(out_mq_top3, paste0(save_fold_DDA_mq, datasetname, '_Maxquant_top3_pro_intensity.tsv'), sep = '\t', col.names = T, row.names = F)
    
    out_mq_lfq<-cbind(protein_table$Protein.IDs[idx_true_organisms], 
                      organims_all[idx_true_organisms],
                      protein_table[,idx_lfq][idx_true_organisms,])
    colnames(out_mq_lfq)[c(1:2)]<-c('Protein', 'Organism')
    write.table(out_mq_lfq, paste0(save_fold_DDA_mq, datasetname, '_Maxquant_pro_maxlfq.tsv'), sep = '\t', col.names = T, row.names = F)
    
    dlfq_file<-paste0(dataset_info$output_folder[i], dataset_info$directlfq_file[i])
    dlfq_table<-read.table(dlfq_file, header = T, sep = '\t', quote = "")
    #organism<-get_organism_mq(dlfq_table$protein)
    
    #dlfq_table<-dlfq_table[which(organism!='decoy' & organism!='contam'),]
    organims_all<-get_organism_mq(dlfq_table$protein)
    idx_true_organisms<-get_idx_True_Organism(organims_all, true_organism)
    
    idx<-setdiff(c(1:length(colnames(dlfq_table))), 
                 c(grep('protein',colnames(dlfq_table)),grep('Leading.razor.protein',colnames(dlfq_table))))
    out_dlfq<-cbind(dlfq_table$protein[idx_true_organisms], 
                    organims_all[idx_true_organisms], 
                    dlfq_table[, idx][idx_true_organisms,])
    colnames(out_dlfq)[c(1,2)]<-c('Protein', 'Organism')
    write.table(out_dlfq, paste0(save_fold_DDA_mq, datasetname,'_Maxquant_dlfq_pro_intensity.tsv'),col.names = T,row.names = F, sep = '\t')
    
    design_file<-paste0(dataset_info$design_file[i])
    design_table<-read.table(design_file, sep = '\t', header = T, quote = "")
    write.table(design_table, paste0(save_fold_DDA_mq, datasetname,'_Maxquant_design.tsv'),col.names = T,row.names = F, sep = '\t')
    
    mstats_design<-data.frame(cbind(design_table$condition, design_table$replicate))
    colnames(mstats_design)<-c('Condition','BioReplicate')
    mstats_design$IsotopeLabelType<-c(rep('L',length(design_table$condition)))
    
    file_names<-c()
    for (fl in design_table$file) {
      flsp=strsplit(fl,'\\', fixed=T)[[1]]
      flc=flsp[length(flsp)]
      flc=strsplit(flc,'.',fixed = T)[[1]][1]
      file_names<-c(file_names,flc)
    }
    mstats_design$Raw.file=file_names
    mstats_design$Run=file_names
    write.table(mstats_design, paste0(save_fold_DDA_mq, datasetname,'_Maxquant_design_msstats.tsv'),col.names = T,row.names = F, sep = '\t')
    
    
    peptide_file<-paste0(dataset_info$output_folder[i], dataset_info$main_output_peptide_file[i])
    peptide_table<-read.table(peptide_file, header = T, sep = '\t', quote = "")
    idx_mq_inten<-grep('Intensity.',colnames(peptide_table))
    idx_lfq<-grep('LFQ.intensity.',colnames(peptide_table))
    
    organims_all<-get_organism_mq(peptide_table$Proteins)
    idx_true_organisms<-get_idx_True_Organism(organims_all, true_organism)
    
    out_mq_inten<-cbind(peptide_table$Sequence[idx_true_organisms],
                        peptide_table$Proteins[idx_true_organisms],
                        organims_all[idx_true_organisms],
                        peptide_table[,idx_mq_inten][idx_true_organisms,])
    colnames(out_mq_inten)[c(1:3)]<-c('Sequence', 'Protein', 'Organism')
    write.table(out_mq_inten, paste0(save_fold_DDA_mq, datasetname, '_Maxquant_pep_intensity.tsv'), sep = '\t', col.names = T, row.names = F)
    
    out_mq_lfq<-cbind(peptide_table$Sequence[idx_true_organisms],
                      peptide_table$Proteins[idx_true_organisms],
                      organims_all[idx_true_organisms],
                      peptide_table[,idx_lfq][idx_true_organisms,])
    colnames(out_mq_lfq)[c(1:3)]<-c('Sequence', 'Protein', 'Organism')
    write.table(out_mq_lfq, paste0(save_fold_DDA_mq, datasetname, '_Maxquant_pep_maxlfq.tsv'), sep = '\t', col.names = T, row.names = F)
    
    all_proteins<-rbind(cbind(protein_table$Protein.IDs, get_organism_mq(protein_table$Protein.IDs)),
                        cbind(dlfq_table$protein, get_organism_mq(dlfq_table$protein)),
                        cbind(peptide_table$Proteins, get_organism_mq(peptide_table$Proteins)))
    
    idx_pro<-c()
    for (org in true_organism) {
      idx_org = grep(org, all_proteins[,2])
      idx_pro<-c(idx_pro, idx_org)
    }
    all_proteins<-all_proteins[idx_pro,]
    labs<-c(rep(0, length(all_proteins[,1])))
    idx_pos<-c()
    for (de_org in de_organism) {
      idx_org = grep(org, all_proteins[,2])
      if(de_org=='UPS'){#contam added by FragPipe are overlapped with human
        idx_human = grep("HUMAN", all_proteins[,2])
        uniport_UPS<-get_uniport_id_UPS()
        for (ups in uniport_UPS) {
          idx_ups<-grep(ups,all_proteins[,1])
          idx_org<-c(idx_org, idx_ups)
        }
      }
      idx_pos<-c(idx_pos, idx_org)
    }
    labs[idx_pos] = 1
    all_proteins = cbind(all_proteins, labs)
    
    colnames(all_proteins)<-c('Protein', 'Organism', 'DEP')
    
    uni_pro = unique(all_proteins[,1])
    idx_uni<-match(uni_pro, all_proteins[,1])
    
    all_proteins<-all_proteins[idx_uni,]
    write.table(all_proteins, paste0(save_fold_DDA_mq, datasetname,'_Maxquant_all_proteins.tsv'),col.names = T,row.names = F, sep = '\t')
    
  }
}

#dataset_info<-read.table(paste0(dataset_info_path, 'DDA_mq.txt'),
#                         header = T, sep = '\t')

#get_allType_matrix_mq_dda(dataset_info)


###################### DIA

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

prepare_diann_out<-function(report_file, datasetname, true_organism, method, norm, N, out_name, designs){
  library(iq)
  tab=process_long_format(report_file,
                          output_filename = paste0(save_fold_DIA_DIANN, datasetname, '_DIANN_', method,'.tsv'),
                          annotation_col = c("Protein.Names", "Genes"),
                          normalization = norm,#"median",
                          filter_double_less = c("Global.Q.Value" = "0.01", "Global.PG.Q.Value" = "0.01"),
                          method=method, N=N)
  
  tab<-as.data.frame(tab)
  sp_orga = get_organism_dia(tab$Protein.Group, tab$Protein.Names)
  out_table<-data.frame(Protein=sp_orga$sp, Organism=sp_orga$orga)
  #tab$Protein = sp_orga$sp
  #tab$Organism = sp_orga$orga
  idx_true_organisms<-get_idx_True_Organism(sp_orga$orga, true_organism)
  out_data <- cbind(out_table, subset(tab, select = -c(Protein.Group, Protein.Names, Genes)))
  
  com<-intersect(designs$file, colnames(out_data))
  idx_col<-match(com, colnames(out_data))
  idx_des<-match(com, designs$file)
  colnames(out_data)[idx_col]<-designs$sample[idx_des]
  
  write.table(out_data[idx_true_organisms,], out_name, sep = '\t', col.names = T, row.names = F)
  return(sp_orga)
}

get_allType_matrix_diann_dia<-function(dataset_info){
  for(i in 1:length(dataset_info$dataset)){
    datasetname<-dataset_info$dataset[i]
    true_organism <- strsplit(dataset_info$organism[i], ';', fixed = T)[[1]]
    de_organism<- strsplit(dataset_info$trueDE[i], ';', fixed = T)[[1]]
    report_file<-paste0(dataset_info$output_folder[i], dataset_info$main_output_evidence[i])
    
    design_file<-dataset_info$design_file[i]
    design_msstats_file<-dataset_info$msstat_design_path[i]
    designs<-read.table(design_file, header = T, sep = '\t')
    msstats_designs<-read.table(design_msstats_file, header = T, sep = '\t')
    #designs$file<-gsub(':', '\\.', gsub("\\\\", '\\.', designs$file))
    source('D:/data/benchmark/iq-master/R/iq-fast_MOD.R')
    #top1(maxInt)
    #no normalization ('none')
    out_name = paste0(save_fold_DIA_DIANN, datasetname, '_DIANN_top1.tsv')
    sp_org_top1<-prepare_diann_out(report_file, datasetname, true_organism, 'topN', 'none', 1, out_name, designs)
    # top1_table<-read.table(out_name, header = T, sep = '\t')
    # com<-intersect(designs$file, colnames(top1_table))
    # idx_col<-match(com, colnames(top1_table))
    # idx_des<-match(com, designs$file)
    # colnames(top1_table)[idx_col]<-designs$sample[idx_des]
    # write.table(top1_table, out_name, sep = '\t', col.names = T, row.names = F)
    #top3
    #no normalization ('none')
    out_name = paste0(save_fold_DIA_DIANN, datasetname, '_DIANN_top3.tsv')
    sp_org_top3<-prepare_diann_out(report_file, datasetname, true_organism, 'topN', 'none', 3, out_name, designs)
    
    #MaxLFQ
    #no normalization ("none")
    out_name = paste0(save_fold_DIA_DIANN, datasetname, '_DIANN_maxlfq_nonorm.tsv')
    sp_org_maxlfq_nonorm<-prepare_diann_out(report_file, datasetname, true_organism, "maxlfq", 'none', NULL, out_name, designs)
    
    #MaxLFQ
    #normalization ("median")
    out_name = paste0(save_fold_DIA_DIANN, datasetname, '_DIANN_maxlfq.tsv')
    sp_org_maxlfq<-prepare_diann_out(report_file, datasetname, true_organism, "maxlfq", 'median', NULL, out_name, designs)
    
    #median_polish
    #no normalization ("none")
    out_name = paste0(save_fold_DIA_DIANN, datasetname, '_DIANN_median_polish.tsv')
    sp_org_median_polish<-prepare_diann_out(report_file, datasetname, true_organism, "median_polish", 'none', NULL, out_name, designs)
    
    #sum
    #out_name = paste0(save_fold_DIA_DIANN, datasetname, '_DIANN_sum.tsv')
    #sp_org_sum<-prepare_diann_out(report_file, datasetname, true_organism, "sum", 'none', NULL, out_name)
    
    #directLFQ
    dlfq_file<-paste0(dataset_info$output_folder[i], dataset_info$directlfq_file[i])
    dlfq_table<-read.table(dlfq_file, header = T, sep = '\t', quote = "")
    colnames(dlfq_table)<-gsub('X', '', colnames(dlfq_table))
    
    sp_orga_dlfq = get_organism_dia(dlfq_table$Protein.Group, dlfq_table$Protein.Names)
    out_table = data.frame(Protein=sp_orga_dlfq$sp, Organism=sp_orga_dlfq$orga)
    # dlfq_table$Protein = sp_orga_dlfq$sp
    # dlfq_table$Organism = sp_orga_dlfq$orga
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
    
    idx_true_organisms<-get_idx_True_Organism(sp_orga_dlfq$orga, true_organism)
    out_name = paste0(save_fold_DIA_DIANN, datasetname, '_DIANN_dlfq.tsv')
    write.table(out_table[idx_true_organisms,], out_name, sep = '\t', col.names = T, row.names = F)
    
    all_proteins<-rbind(cbind(sp_org_top1$sp, sp_org_top1$orga),
                        cbind(sp_org_top3$sp, sp_org_top3$orga),
                        cbind(sp_org_maxlfq$sp, sp_org_maxlfq$orga),
                        cbind(sp_org_maxlfq_nonorm$sp, sp_org_maxlfq_nonorm$orga),
                        cbind(sp_org_median_polish$sp, sp_org_median_polish$orga),
                        # cbind(sp_org_sum$sp, sp_org_sum$orga),
                        cbind(sp_orga_dlfq$sp, sp_orga_dlfq$orga))

    idx_pro<-c()
    for (org in true_organism) {
      idx_org = grep(org, all_proteins[,2])
      idx_pro<-c(idx_pro, idx_org)
    }
    all_proteins<-all_proteins[idx_pro,]
    labs<-c(rep(0, length(all_proteins[,1])))
    idx_pos<-c()
    for (de_org in de_organism) {
      idx_org = grep(org, all_proteins[,2])
      if(de_org=='UPS'){#contam added by FragPipe are overlapped with human
        idx_human = grep("HUMAN", all_proteins[,2])
        uniport_UPS<-get_uniport_id_UPS()
        for (ups in uniport_UPS) {
          idx_ups<-grep(ups,all_proteins[,1])
          idx_org<-c(idx_org, idx_ups)
        }
      }
      idx_pos<-c(idx_pos, idx_org)
    }
    labs[idx_pos] = 1
    all_proteins = cbind(all_proteins, labs)

    colnames(all_proteins)<-c('Protein', 'Organism', 'DEP')

    uni_pro = unique(all_proteins[,1])
    idx_uni<-match(uni_pro, all_proteins[,1])

    all_proteins<-all_proteins[idx_uni,]
    write.table(all_proteins, paste0(save_fold_DIA_DIANN, datasetname,'_DIANN_all_proteins.tsv'),col.names = T,row.names = F, sep = '\t')
    write.table(designs, paste0(save_fold_DIA_DIANN, datasetname,'_DIANN_design.tsv'),col.names = T,row.names = F, sep = '\t')
    write.table(msstats_designs, paste0(save_fold_DIA_DIANN, datasetname,'_DIANN_design_msstats.tsv'),col.names = T,row.names = F, sep = '\t')
    
  }
}

dataset_info<-read.table(paste0(dataset_info_path, 'DIA_DIANN.txt'),
                         header = T, sep = '\t')

get_allType_matrix_diann_dia(dataset_info[2:7,])

get_sp_pro_from_fasta<-function(proteinlist, fasta_file){
  library(seqinr)
  db_pro<-read.fasta(fasta_file,as.string = TRUE, seqtype = "AA")
  all_db_pro<-as.vector(unlist(names(db_pro)))
  sp_pros<-c()
  for (i in 1:length(proteinlist)) {
    proi=strsplit(proteinlist[i],';',fixed = T)[[1]]
    sp_pro<-''
    for (j in 1:length(proi)) {
      idx_pro_db<-grep(proi[j],all_db_pro)
      if(length(idx_pro_db)>0){
      if(j==1){
        sp_pro=all_db_pro[idx_pro_db]
      }else if(j>1){
        sp_pro=paste0(sp_pro,';',all_db_pro[idx_pro_db])
      }
      }
    }
    sp_pros<-c(sp_pros,sp_pro)
  }
  organisms<-get_organism_mq(sp_pros)
  return(list(sp=sp_pros, orga=organisms))
}

prepare_spt_out<-function(report_file, datasetname, true_organism, method, norm, N, out_name, designs){
  library(iq)
  tab=process_long_format(report_file,
                          output_filename = paste0(save_fold_DIA_spt, datasetname, '_spt_', method,'.tsv'),
                          sample_id  = "R.FileName",
                          primary_id = "PG.ProteinGroups",
                          secondary_id = c("EG.Library", "FG.Id", "FG.Charge", "F.FrgIon", "F.Charge", "F.FrgLossType"),
                          intensity_col = "F.PeakArea",
                          annotation_col = c("PG.Genes", "PG.ProteinNames", "PG.FastaFiles"),
                          filter_string_equal = c("F.ExcludedFromQuantification" = "False"),
                          filter_double_less = c("PG.Qvalue" = "0.01", "EG.Qvalue" = "0.01"),
                          log2_intensity_cutoff = 0,
                          normalization = norm,#"median",
                          method=method, N=N)
 
   tab<-as.data.frame(tab)
  sp_orga = get_organism_dia(tab$PG.ProteinGroups, tab$PG.ProteinNames)
  #tab$Protein = sp_orga$sp
  #tab$Organism = sp_orga$orga
  
  out_table<-data.frame(Protein=sp_orga$sp, Organism=sp_orga$orga)
  #tab$Protein = sp_orga$sp
  #tab$Organism = sp_orga$orga
  idx_true_organisms<-get_idx_True_Organism(sp_orga$orga, true_organism)
  out_data <- cbind(out_table, subset(tab, select = -c(PG.ProteinGroups, PG.ProteinNames, PG.Genes, PG.FastaFiles)))
  
  com<-intersect(designs$file, colnames(out_data))
  idx_col<-match(com, colnames(out_data))
  idx_des<-match(com, designs$file)
  colnames(out_data)[idx_col]<-designs$sample[idx_des]
  
  write.table(out_data[idx_true_organisms,], out_name, sep = '\t', col.names = T, row.names = F)
  
  #idx_true_organisms<-get_idx_True_Organism(sp_orga$orga, true_organism)
  
  #write.table(tab[idx_true_organisms,], out_name, sep = '\t', col.names = T, row.names = F)
  return(sp_orga)
}

get_allType_matrix_spt_dia<-function(dataset_info){
  for(i in 1:length(dataset_info$dataset)){
    datasetname<-dataset_info$dataset[i]
    true_organism <- strsplit(dataset_info$organism[i], ';', fixed = T)[[1]]
    de_organism<- strsplit(dataset_info$trueDE[i], ';', fixed = T)[[1]]
    report_file<-paste0(dataset_info$output_folder[i], dataset_info$main_output_evidence[i])
    
    design_file<-dataset_info$design_file[i]
    design_msstats_file<-dataset_info$msstat_design_path[i]
    designs<-read.table(design_file, header = T, sep = '\t')
    msstats_designs<-read.table(design_msstats_file, header = T, sep = '\t')
    source('D:/data/benchmark/iq-master/R/iq-fast_MOD.R')
    #top1(maxInt)
    #no normalization ('none')
    out_name = paste0(save_fold_DIA_spt, datasetname, '_spt_top1.tsv')
    #sp_org_top1<-prepare_spt_out(report_file, datasetname, true_organism, 'topN', 'none', 1, out_name, designs)
    
    #top3
    #no normalization ('none')
    out_name = paste0(save_fold_DIA_spt, datasetname, '_spt_top3.tsv')
    #sp_org_top3<-prepare_spt_out(report_file, datasetname, true_organism, 'topN', 'none', 3, out_name, designs)
    
    #MaxLFQ
    #no normalization ("none")
    out_name = paste0(save_fold_DIA_spt, datasetname, '_spt_maxlfq_nonorm.tsv')
    #sp_org_maxlfq_nonorm<-prepare_spt_out(report_file, datasetname, true_organism, "maxlfq", 'none', NULL, out_name, designs)
    
    #MaxLFQ
    #normalization ("median")
    out_name = paste0(save_fold_DIA_spt, datasetname, '_spt_maxlfq.tsv')
    #sp_org_maxlfq<-prepare_spt_out(report_file, datasetname, true_organism, "maxlfq", 'median', NULL, out_name, designs)
    
    #median_polish
    #no normalization ("none")
    out_name = paste0(save_fold_DIA_spt, datasetname, '_spt_median_polish.tsv')
    #sp_org_median_polish<-prepare_spt_out(report_file, datasetname, true_organism, "median_polish", 'none', NULL, out_name, designs)
    
    #sum
    #out_name = paste0(save_fold_DIA_spt, datasetname, '_spt_sum.tsv')
    #sp_org_sum<-prepare_diann_out(report_file, datasetname, true_organism, "sum", 'none', NULL, out_name)
    
    #directLFQ
    dlfq_file<-paste0(dataset_info$output_folder[i], dataset_info$directlfq_file[i])
    dlfq_table<-read.table(dlfq_file, header = T, sep = '\t', quote = "")
    colnames(dlfq_table)<-gsub('X', '', colnames(dlfq_table))
    sp_orga_dlfq = get_sp_pro_from_fasta(dlfq_table$PG.ProteinGroups, dataset_info$db_file[i])
    out_table = data.frame(Protein=sp_orga_dlfq$sp, Organism=sp_orga_dlfq$orga)
    # dlfq_table$Protein = sp_orga_dlfq$sp
    # dlfq_table$Organism = sp_orga_dlfq$orga
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
    
    idx_true_organisms<-get_idx_True_Organism(sp_orga_dlfq$orga, true_organism)
    out_name = paste0(save_fold_DIA_spt, datasetname, '_spt_dlfq.tsv')
    write.table(out_table[idx_true_organisms,], out_name, sep = '\t', col.names = T, row.names = F)
    
    #  all_proteins<-rbind(cbind(sp_org_top1$sp, sp_org_top1$orga),
    #                     cbind(sp_org_top3$sp, sp_org_top3$orga),
    #                     cbind(sp_org_maxlfq$sp, sp_org_maxlfq$orga),
    #                     cbind(sp_org_maxlfq_nonorm$sp, sp_org_maxlfq_nonorm$orga),
    #                     cbind(sp_org_median_polish$sp, sp_org_median_polish$orga),
    #                     #cbind(sp_org_sum$sp, sp_org_sum$orga),
    #                     cbind(sp_orga_dlfq$sp, sp_orga_dlfq$orga))
    # 
    # idx_pro<-c()
    # for (org in true_organism) {
    #   idx_org = grep(org, all_proteins[,2])
    #   idx_pro<-c(idx_pro, idx_org)
    # }
    # all_proteins<-all_proteins[idx_pro,]
    # labs<-c(rep(0, length(all_proteins[,1])))
    # idx_pos<-c()
    # for (de_org in de_organism) {
    #   idx_org = grep(org, all_proteins[,2])
    #   if(de_org=='UPS'){#contam added by FragPipe are overlapped with human
    #     idx_human = grep("HUMAN", all_proteins[,2])
    #     uniport_UPS<-get_uniport_id_UPS()
    #     for (ups in uniport_UPS) {
    #       idx_ups<-grep(ups,all_proteins[,1])
    #       idx_org<-c(idx_org, idx_ups)
    #     }
    #   }
    #   idx_pos<-c(idx_pos, idx_org)
    # }
    # labs[idx_pos] = 1
    # all_proteins = cbind(all_proteins, labs)
    # 
    # colnames(all_proteins)<-c('Protein', 'Organism', 'DEP')
    # 
    # uni_pro = unique(all_proteins[,1])
    # idx_uni<-match(uni_pro, all_proteins[,1])
    # 
    # all_proteins<-all_proteins[idx_uni,]
    # write.table(all_proteins, paste0(save_fold_DIA_spt, datasetname,'_spt_all_proteins.tsv'),col.names = T,row.names = F, sep = '\t')
    # write.table(designs, paste0(save_fold_DIA_spt, datasetname,'_spt_design.tsv'),col.names = T,row.names = F, sep = '\t')
    # write.table(msstats_designs, paste0(save_fold_DIA_spt, datasetname,'_spt_design_msstats.tsv'),col.names = T,row.names = F, sep = '\t')
  }
}
dataset_info<-read.table(paste0(dataset_info_path, 'DIA_spt.txt'),
                         header = T, sep = '\t')

get_allType_matrix_spt_dia(dataset_info)




###################### TMT
# for fragpipe, philosiper intensity,protein abundance and ratio were used, for maxquant, only reporter ion corrected intensity were used
get_allType_matrix_frag_tmt<-function(dataset_info){
  for(i in 3:length(dataset_info$dataset)){
    print(i)
    
    datasetname<-dataset_info$dataset[i]
    true_organism <- strsplit(dataset_info$organism[i], ';', fixed = T)[[1]]
    de_organism<- strsplit(dataset_info$trueDE[i], ';', fixed = T)[[1]]
    pro_abun_file<-paste0(dataset_info$output_folder[i], dataset_info$main_output_abundance_file[i])
    
    design_file<-paste0(dataset_info$design_file[i])
    designs<-read.table(design_file, header = T, sep = '\t', quote = "")
    
    sample = designs$sample
    condition = designs$condition
    replicate = designs$replicate
    
    # tmtIntegrator_abundance log2
    pro_table<-read.table(pro_abun_file, header = T, sep = '\t', quote = "")
    
    colns<-c()
    idx_col<-c()
    for (samp_i in 1:length(sample)) {
      idx_samp = grep(sample[samp_i], colnames(pro_table))
      idx_col<-c(idx_col, idx_samp)
      colns<-c(colns, sample[samp_i])
    }
    sp_orga_abd = get_sp_pro_from_fasta(pro_table$Index, dataset_info$db_file[i])
    #sp_orga_dlfq = get_organism_dia(dlfq_table$PG.ProteinGroups, dlfq_table$PG.ProteinNames)
    pro_table$Protein = sp_orga_abd$sp
    pro_table$Organism = sp_orga_abd$orga
    idx_true_organisms<-get_idx_True_Organism(sp_orga_abd$orga, true_organism)
    pro_table = pro_table[idx_true_organisms,]
    
    if(i==2){
      a=0
    }
    
    out_abd<-cbind(pro_table$Protein, pro_table$Organism, pro_table[idx_col])
    colnames(out_abd)<-c('Protein', 'Organism', colns)
    out_name = paste0(save_fold_tmt_frag, datasetname, '_FragPipe_tmt_abd.tsv')
    write.table(out_abd, out_name, sep = '\t', col.names = T, row.names = F)
    
    #tmtIntegrator_ratio log2-log2_ref
    pro_ratio_file<-paste0(dataset_info$output_folder[i], dataset_info$main_output_ratio_file[i])
    pro_table<-read.table(pro_ratio_file, header = T, sep = '\t', quote = "")
    
    colns<-c()
    idx_col<-c()
    for (samp_i in 1:length(sample)) {
      idx_samp = grep(sample[samp_i], colnames(pro_table))
      idx_col<-c(idx_col, idx_samp)
      colns<-c(colns, sample[samp_i])
    }
    sp_orga_rat = get_sp_pro_from_fasta(pro_table$Index, dataset_info$db_file[i])
    #sp_orga_dlfq = get_organism_dia(dlfq_table$PG.ProteinGroups, dlfq_table$PG.ProteinNames)
    pro_table$Protein = sp_orga_rat$sp
    pro_table$Organism = sp_orga_rat$orga
    idx_true_organisms<-get_idx_True_Organism(sp_orga_rat$orga, true_organism)
    
    pro_table = pro_table[idx_true_organisms,]
  
    out_rat<-cbind(pro_table$Protein, pro_table$Organism, pro_table[idx_col])
    colnames(out_rat)<-c('Protein', 'Organism', colns)
    out_name = paste0(save_fold_tmt_frag, datasetname, '_FragPipe_tmt_ratio.tsv')
    write.table(out_rat, out_name, sep = '\t', col.names = T, row.names = F)
    
    # philosopher intensity
    pro_phi_file<-paste0(dataset_info$output_Philosopher_folder[i], 'protein.tsv')
    pro_table<-read.table(pro_phi_file, header = T, sep = '\t', quote = "")
    colns<-c()
    idx_col<-c()
    for (samp_i in 1:length(sample)) {
      idx_samp = grep(sample[samp_i], colnames(pro_table))
      idx_col<-c(idx_col, idx_samp)
      colns<-c(colns, sample[samp_i])
    }
    
    organims_all<-get_organism_frag(pro_table$Protein,
                                    pro_table$Indistinguishable.Proteins)

    idx_true_organisms<-get_idx_True_Organism(organims_all, true_organism)
    
    out_phi<-cbind(pro_table$Protein[idx_true_organisms], organims_all[idx_true_organisms], 
                   pro_table[,idx_col][idx_true_organisms,])
    colnames(out_phi)<-c('Protein', 'Organism', colns)
    out_name = paste0(save_fold_tmt_frag, datasetname, '_FragPipe_tmt_phi.tsv')
    write.table(out_phi, out_name, sep = '\t', col.names = T, row.names = F)
    
    all_proteins<-rbind(cbind(sp_orga_abd$sp, sp_orga_abd$orga),
                        cbind(sp_orga_rat$sp, sp_orga_rat$orga),
                        cbind(out_phi$Protein, out_phi$Organism))
    
    idx_pro<-c()
    for (org in true_organism) {
      idx_org = grep(org, all_proteins[,2])
      idx_pro<-c(idx_pro, idx_org)
    }
    all_proteins<-all_proteins[idx_pro,]
    labs<-c(rep(0, length(all_proteins[,1])))
    idx_pos<-c()
    for (de_org in de_organism) {
      idx_org = grep(org, all_proteins[,2])
      if(de_org=='UPS'){#contam added by FragPipe are overlapped with human
        idx_human = grep("HUMAN", all_proteins[,2])
        uniport_UPS<-get_uniport_id_UPS()
        for (ups in uniport_UPS) {
          idx_ups<-grep(ups,all_proteins[,1])
          idx_org<-c(idx_org, idx_ups)
        }
      }
      idx_pos<-c(idx_pos, idx_org)
    }
    labs[idx_pos] = 1
    all_proteins = cbind(all_proteins, labs)
    
    colnames(all_proteins)<-c('Protein', 'Organism', 'DEP')
    
    uni_pro = unique(all_proteins[,1])
    idx_uni<-match(uni_pro, all_proteins[,1]) 
    
    all_proteins<-all_proteins[idx_uni,]
    write.table(all_proteins, paste0(save_fold_tmt_frag, datasetname,'_FragPipe_tmt_all_proteins.tsv'),col.names = T,row.names = F, sep = '\t')
    
    designs$sample_name<-designs$sample
    write.table(designs, paste0(save_fold_tmt_frag, datasetname,'_FragPipe_tmt_design.tsv'),col.names = T,row.names = F, sep = '\t')
    
    MSstat_design_file<-paste0(dataset_info$msstat_design_path[i])
    MSstat_design<-read.table(MSstat_design_file, header = T, sep = ',', quote = "")
    MSstat_condi<-MSstat_design$Condition
    for (samp_i in 1:length(designs$sample)) {
      idx_cd<-grep(designs$channel[samp_i], MSstat_design$Channel)
      MSstat_condi[idx_cd]=designs$condition[samp_i]
    }
    MSstat_design$Condition<-MSstat_condi
    write.table(designs, paste0(save_fold_tmt_frag, datasetname,'_FragPipe_tmt_MSstats_design.tsv'),col.names = T,row.names = F, sep = '\t')
  }
}

# dataset_info<-read.table(paste0(dataset_info_path, 'TMT_frag.txt'),
#                                                     header = T, sep = '\t')
#get_allType_matrix_frag_tmt(dataset_info)


get_allType_matrix_mq_tmt<-function(dataset_info){
  for(i in 1:length(dataset_info$dataset)){
    datasetname<-dataset_info$dataset[i]
    true_organism <- strsplit(dataset_info$organism[i], ';', fixed = T)[[1]]
    de_organism<- strsplit(dataset_info$trueDE[i], ';', fixed = T)[[1]]
    pro_file<-paste0(dataset_info$output_folder[i], dataset_info$main_output_protein_file[i])
    
    design_file<-paste0(dataset_info$design_file[i])
    designs<-read.table(design_file, header = T, sep = '\t', quote = "")
    
    sample = designs$sample
    condition = designs$condition
    replicate = designs$replicate
    
    # tmtIntegrator_abundance log2
    pro_table<-read.table(pro_file, header = T, sep = '\t', quote = "")
    
    idx_int<-grep('Reporter.intensity.corrected.', colnames(pro_table))
    
    organims_all<-get_organism_mq(pro_table$Protein.IDs)
    idx_true_organisms<-get_idx_True_Organism(organims_all, true_organism)
    #protein_table<-protein_table[which(organs!='decoy' & organs!='contam'),]
    
    out_pro<-cbind(pro_table$Protein.IDs[idx_true_organisms], 
                        organims_all[idx_true_organisms],
                   pro_table[,idx_int][idx_true_organisms,])
    colnames(out_pro)[c(1:2)]<-c('Protein', 'Organism')
    out_name = paste0(save_fold_tmt_mq, datasetname, '_Maxquant_tmt_intensity.tsv')
    write.table(out_pro, out_name, sep = '\t', col.names = T, row.names = F)
    
    
    all_proteins<-rbind(#cbind(sp_orga_abd$sp, sp_orga_abd$orga),
      cbind(out_pro$Protein, out_pro$Organism))
    
    idx_pro<-c()
    for (org in true_organism) {
      idx_org = grep(org, all_proteins[,2])
      idx_pro<-c(idx_pro, idx_org)
    }
    all_proteins<-all_proteins[idx_pro,]
    labs<-c(rep(0, length(all_proteins[,1])))
    idx_pos<-c()
    for (de_org in de_organism) {
      idx_org = grep(org, all_proteins[,2])
      if(de_org=='UPS'){#contam added by FragPipe are overlapped with human
        idx_human = grep("HUMAN", all_proteins[,2])
        uniport_UPS<-get_uniport_id_UPS()
        for (ups in uniport_UPS) {
          idx_ups<-grep(ups,all_proteins[,1])
          idx_org<-c(idx_org, idx_ups)
        }
      }
      idx_pos<-c(idx_pos, idx_org)
    }
    labs[idx_pos] = 1
    all_proteins = cbind(all_proteins, labs)
    
    colnames(all_proteins)<-c('Protein', 'Organism', 'DEP')
    
    uni_pro = unique(all_proteins[,1])
    idx_uni<-match(uni_pro, all_proteins[,1]) 
    
    all_proteins<-all_proteins[idx_uni,]
    write.table(all_proteins, paste0(save_fold_tmt_mq, datasetname,'_Maxquant_tmt_all_proteins.tsv'),col.names = T,row.names = F, sep = '\t')
    
    designs$sample_name<-designs$sample
    write.table(designs, paste0(save_fold_tmt_mq, datasetname,'_Maxquant_tmt_design.tsv'),col.names = T,row.names = F, sep = '\t')
    
    MSstat_design_file<-paste0(dataset_info$msstat_design_path[i])
    MSstat_design<-read.table(MSstat_design_file, header = T, sep = ',', quote = "")
    MSstat_condi<-MSstat_design$Condition
    for (samp_i in 1:length(designs$sample)) {
      idx_cd<-grep(designs$channel[samp_i], MSstat_design$Channel)
      MSstat_condi[idx_cd]=designs$condition[samp_i]
    }
    MSstat_design$Condition<-MSstat_condi
    write.table(designs, paste0(save_fold_tmt_mq, datasetname,'_Maxquant_tmt_MSstats_design.tsv'),col.names = T,row.names = F, sep = '\t')
  }
}

# dataset_info<-read.table(paste0(dataset_info_path, 'TMT_mq.txt'),
#                          header = T, sep = '\t')
#get_allType_matrix_mq_tmt(dataset_info)
