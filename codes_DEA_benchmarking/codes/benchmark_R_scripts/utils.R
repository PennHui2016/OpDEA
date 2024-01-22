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
      proteins<-c(protein_ids[i], strsplit(indist_pro[i], ',',fixed = T)[[1]])
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
  
  get_labels<-function(DEA_res, true_organism, de_organism, true_lgfc){
    proteins<-as.character(DEA_res$protein)
    contrast<-DEA_res$contrast
    uni_cont<-unique(contrast)
    organs<-get_organism_mq(proteins)
    true_lgfcs<-read.table(true_lgfc,sep = '\t', header = T)
    all_proteins<-cbind(proteins, organs)
    # idx_pro<-c()
    # for (org in true_organism) {
    #   idx_org = grep(org, all_proteins[,2])
    #   idx_pro<-c(idx_pro, idx_org)
    # }
    # all_proteins<-all_proteins[idx_pro,]
    labs<-c(rep(0, length(all_proteins[,1])))
    idx_pos<-c()
    for (de_org in de_organism) {
      idx_org = grep(de_org, all_proteins[,2])
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
    
    TlogFC<-c(rep(0,length(all_proteins[,1])))
    for (cont in uni_cont) {
      cont_c<-gsub('condition','',cont)
      conts<-strsplit(cont_c,'-',fixed = T)[[1]]
      ci1<-which(true_lgfcs$condition==conts[1])
      ci2<-which(true_lgfcs$condition==conts[2])
      for (deo in de_organism) {
        idx<-which(colnames(true_lgfcs)==deo)
        tlfc=log2(true_lgfcs[,idx][ci1]/true_lgfcs[,idx][ci2])
        TlogFC[which(contrast==cont & organs==deo)] = tlfc
        if(deo=='UPS'){
          TlogFC[which(contrast==cont & labs==1)] = tlfc
        }
      }
    }
    
    all_proteins = cbind(all_proteins, TlogFC)
    colnames(all_proteins)<-c('Protein', 'Organism', 'DEP', 'TlogFC')
    
    # uni_pro = unique(all_proteins[,1])
    # idx_uni<-match(uni_pro, all_proteins[,1])
    # 
    all_proteins<-as.data.frame(all_proteins)
    
    return(all_proteins)
  }
  #labs<-get_labels(res_all, true_organisms, DE_organisms, true_lgfc)
  #labs<-get_labels(res_all$protein, true_organism, DE_organism)
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