###########
##  process outputs from different quntification platforms
######

## 1.software specific
## fragpipe:: intensity (topN 0); maxquant: intensity (topN 0); DIANN: topN1; spectranaut: 


## 2. top3:: fragpipe rerun; maxquant: output; DIANN: iq-mod; spectranaut: iq-mod


## 3. maxlfq:: fragpipe: output; maxquant: output; DIANN: iq-mod; spectranaut: iq-mod


## 4. directlfq::  fragpipe: directlfq; maxquant: directlfq; DIANN: directlfq; spectranaut: directlfq


## 5. spectral counts:: fragpipe: output; maxquant: output;

save_fold_DDA_frag<-'E:/benchmarking/data/DDA/FragPipe/'
save_fold_DDA_mq<-'E:/benchmarking/data/DDA/maxquant/'
save_fold_DIA_DIANN<-'D:/data/benchmark/data/DIA/DIANN/'
save_fold_DIA_spt<-'D:/data/benchmark/data/DIA/spectraunt/'
save_fold_tmt_frag<-'D:/data/benchmark/data/TMT/FragPipe/'
save_fold_tmt_mq<-'D:/data/benchmark/data/TMT/Maxquant/'

get_organism<-function(protein_ids){
  Organs<-c()
  for (i in 1:length(protein_ids)) {
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
    num_uni_og<-c()
    for (j in uni_og) {
      num_uni_og<-c(num_uni_og, length(which(orgs==j)))
    }
    idx=which(num_uni_og==max(num_uni_og))
    if(length(idx)>1){
      maj_organ = uni_og[idx[1]]
    }else{
      maj_organ = uni_og[idx]
    }
    Organs<-c(Organs,maj_organ)
  }
  return(Organs)
}

get_allType_matrix_frag_dda<-function(dataset_info){
  for(i in 1:length(dataset_info$dataset)){
    datasetname<-dataset_info$dataset[i]
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
  
  out_count_frag<-cbind(protein_table$Protein, get_organism(protein_table$Protein), protein_table[,count_idx])
  out_iten_frag<-cbind(protein_table$Protein, get_organism(protein_table$Protein), protein_table[,frag_inten_idx])
  out_inten_maxlfq<-cbind(protein_table$Protein, get_organism(protein_table$Protein), protein_table[,lfq_idx])
  
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
  out_iten_frag<-cbind(protein_table$Protein, get_organism(protein_table$Protein), protein_table[,frag_inten_idx])
  colnames(out_iten_frag)[c(1,2)]<-c('Protein', 'Organism')
  write.table(out_iten_frag, paste0(save_fold_DDA_frag, datasetname,'_FragPipe_top3_pro_intensity.tsv'),col.names = T,row.names = F, sep = '\t')
  
  #directlfq
  dlfq_file<-paste0(dataset_info$output_folder[i],'TOP3/', dataset_info$directlfq_file[i])
  dlfq_table<-read.table(dlfq_file, header = T, sep = '\t', quote = "")
  idx<-setdiff(c(1:length(colnames(dlfq_table))), 
               c(grep('protein',colnames(dlfq_table)),grep('Protein',colnames(dlfq_table))))
  out_dlfq<-cbind(dlfq_table$protein, get_organism(dlfq_table$protein), dlfq_table[, idx])
  colnames(out_dlfq)[c(1,2)]<-c('Protein', 'Organism')
  write.table(out_dlfq, paste0(save_fold_DDA_frag, datasetname,'_FragPipe_dlfq_pro_intensity.tsv'),col.names = T,row.names = F, sep = '\t')
  
  design_file<-paste0(dataset_info$output_folder[i], 'TOP3/',dataset_info$design_file[i])
  design_table<-read.table(design_file, sep = '\t', header = T, quote = "")
  write.table(design_table, paste0(save_fold_DDA_frag, datasetname,'_FragPipe_design.tsv'),col.names = T,row.names = F, sep = '\t')
  
  peptide_file<-paste0(dataset_info$output_folder[i],'TOP0/', dataset_info$main_output_peptide_file[i])
  peptide_table<-read.table(peptide_file, header = T, sep = '\t', quote = "")
  inten_idx<-grep('.Intensity',colnames(peptide_table))
  lfq_idx<-grep('.MaxLFQ.Intensity',colnames(peptide_table))
  frag_inten_idx<-setdiff(inten_idx, lfq_idx)
  out_iten_frag<-cbind(peptide_table$Peptide.Sequence, peptide_table$Protein, get_organism(peptide_table$Protein), peptide_table[,frag_inten_idx])
  out_inten_maxlfq<-cbind(peptide_table$Peptide.Sequence, peptide_table$Protein, get_organism(peptide_table$Protein), peptide_table[,lfq_idx])
  
  colnames(out_iten_frag)[c(1,2,3)]<-c('Peptide Sequence', 'Protein', 'Organism')
  colnames(out_inten_maxlfq)[c(1,2,3)]<-c('Peptide Sequence', 'Protein', 'Organism')

  write.table(out_iten_frag, paste0(save_fold_DDA_frag, datasetname,'_FragPipe_pep_intensity.tsv'),col.names = T,row.names = F, sep = '\t')
  write.table(out_inten_maxlfq, paste0(save_fold_DDA_frag, datasetname,'_FragPipe_pep_maxlfq.tsv'),col.names = T,row.names = F, sep = '\t')
  
  }
}

dataset_info<-read.table('E:/benchmarking/data/dataset_info/DDA_frag.txt',
                         header = T, sep = '\t')

#get_allType_matrix_frag_dda(dataset_info)

get_allType_matrix_mq_dda<-function(dataset_info){
  for(i in 1:length(dataset_info$dataset)){
    datasetname<-dataset_info$dataset[i]
    datasetname<-dataset_info$dataset[i]
    protein_file<-paste0(dataset_info$output_folder[i], dataset_info$main_output_protein_file[i])
    protein_table<-read.table(protein_file, header = T, sep = '\t', quote = "")
    protein_table<-protein_table[which(protein_table$Reverse!="+" & protein_table$Potential.contaminant!="+"),]
    
    #spectral counts
    idx_count<-grep('MS.MS.count.',colnames(protein_table))
    idx_mq_inten<-grep('Intensity.',colnames(protein_table))
    idx_lfq<-grep('LFQ.intensity.',colnames(protein_table))
    idx_top3<-grep('Top3.',colnames(protein_table))
    
    organs<-get_organism(protein_table$Majority.protein.IDs)
    protein_table<-protein_table[which(organs!='decoy' & organs!='contam'),]
    
    out_mq_count<-cbind(protein_table$Majority.protein.IDs, organs, protein_table[,idx_count])
    colnames(out_mq_count)[c(1:2)]<-c('Protein', 'Organism')
    write.table(out_mq_count, paste0(save_fold_DDA_mq, datasetname, '_Maxquant_counts.tsv'), sep = '\t', col.names = T, row.names = F)
    
    out_mq_inten<-cbind(protein_table$Majority.protein.IDs, organs, protein_table[,idx_mq_inten])
    colnames(out_mq_inten)[c(1:2)]<-c('Protein', 'Organism')
    write.table(out_mq_inten, paste0(save_fold_DDA_mq, datasetname, '_Maxquant_intensity.tsv'), sep = '\t', col.names = T, row.names = F)
    
    out_mq_top3<-cbind(protein_table$Majority.protein.IDs, organs, protein_table[,idx_top3])
    colnames(out_mq_top3)[c(1:2)]<-c('Protein', 'Organism')
    write.table(out_mq_top3, paste0(save_fold_DDA_mq, datasetname, '_Maxquant_top3_intensity.tsv'), sep = '\t', col.names = T, row.names = F)
    
    out_mq_lfq<-cbind(protein_table$Majority.protein.IDs, organs, protein_table[,idx_lfq])
    colnames(out_mq_lfq)[c(1:2)]<-c('Protein', 'Organism')
    write.table(out_mq_lfq, paste0(save_fold_DDA_mq, datasetname, '_Maxquant_maxlfq_intensity.tsv'), sep = '\t', col.names = T, row.names = F)
    
    dlfq_file<-paste0(dataset_info$output_folder[i], dataset_info$directlfq_file[i])
    dlfq_table<-read.table(dlfq_file, header = T, sep = '\t', quote = "")
    organism<-get_organism(dlfq_table$Leading.razor.protein)
    dlfq_table<-dlfq_table[which(organism!='decoy' & organism!='contam'),]
    organism<-organism[which(organism!='decoy' & organism!='contam')]
    idx<-setdiff(c(1:length(colnames(dlfq_table))), 
                 c(grep('protein',colnames(dlfq_table)),grep('Leading.razor.protein',colnames(dlfq_table))))
    out_dlfq<-cbind(dlfq_table$protein, organism, dlfq_table[, idx])
    colnames(out_dlfq)[c(1,2)]<-c('Protein', 'Organism')
    write.table(out_dlfq, paste0(save_fold_DDA_mq, datasetname,'_Maxquant_dlfq_intensity.tsv'),col.names = T,row.names = F, sep = '\t')
    
    design_file<-paste0(dataset_info$design_file[i])
    design_table<-read.table(design_file, sep = '\t', header = T, quote = "")
    write.table(design_table, paste0(save_fold_DDA_mq, datasetname,'_Maxquant_design.tsv'),col.names = T,row.names = F, sep = '\t')
    
    peptide_file<-paste0(dataset_info$output_folder[i], dataset_info$main_output_peptide_file[i])
    peptide_table<-read.table(peptide_file, header = T, sep = '\t', quote = "")
    idx_mq_inten<-grep('Intensity.',colnames(peptide_table))
    idx_lfq<-grep('LFQ.intensity.',colnames(peptide_table))
    
    organs<-get_organism(peptide_table$Leading.razor.protein)
    out_mq_inten<-cbind(peptide_table$Sequence, peptide_table$Leading.razor.protein, organs, peptide_table[,idx_mq_inten])
    colnames(out_mq_inten)[c(1:3)]<-c('Sequence', 'Protein', 'Organism')
    write.table(out_mq_inten, paste0(save_fold_DDA_mq, datasetname, '_Maxquant_pep_intensity.tsv'), sep = '\t', col.names = T, row.names = F)
    
    out_mq_lfq<-cbind(peptide_table$Sequence, peptide_table$Leading.razor.protein, organs, peptide_table[,idx_lfq])
    colnames(out_mq_lfq)[c(1:3)]<-c('Sequence', 'Protein', 'Organism')
    write.table(out_mq_lfq, paste0(save_fold_DDA_mq, datasetname, '_Maxquant_maxlfq_pep_intensity.tsv'), sep = '\t', col.names = T, row.names = F)
    
  }
}

dataset_info<-read.table('E:/benchmarking/data/dataset_info/DDA_mq.txt',
                         header = T, sep = '\t')

get_allType_matrix_mq_dda(dataset_info)



protein_file<-paset0()
protein_table

dataset_type = 'DDA_Frag'
if(dataset_type == 'DDA_Frag'){
  dataset_info<-read.table('D:/data/benchmark/data/dataset_info/DDA_Frag.txt',
                           header = T, sep = '\t')
  
}
DDA_dataset=c('HYE6600735_LFQ', 'HYEqe735_LFQ','HYEtims735_LFQ',
              'HYtims134_LFQ','HEtims425_LFQ','YUltq099_LFQ',
              'YUltq819_LFQ','HEqe408_LFQ','HYqfl683_LFQ',
              'HYEtims777_LFQ')

source_folder = ''

dataset = ''
save_file=''
source_folder = ''
main_output= ''
design_file = ''