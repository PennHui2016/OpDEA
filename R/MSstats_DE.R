

res_ttest<-function(platform, raw, evid, acq, normal, imput, logT, designs){
  library(MSstats)
  library(MSnbase)
  library(readr)
  library(tidyverse)

  set.seed(123)

  if(platform=='FragPipe'){
    evidence_file_name = ''
    protein_file_name = evid
  }else if(platform=='Maxquant'){
    evidence_file_name = evid
    protein_file_name = raw
  }else if(platform=='DIANN' | platform=='Spectronaut'){
    evidence_file_name = evid
    infile <- read.table(evidence_file_name, quote = "", sep = '\t', header = TRUE)
  }

  #designs = read.table(msstats_design_file, sep = '\t', header = TRUE)

  if(normal=='FALSE'){
    normal=FALSE
  }

  if(imput=='FALSE'){
    imput=FALSE
  }else if(imput=='TRUE'){
    imput=TRUE
  }

  all_conts<-function(design_data){
    conds<-design_data$Condition
    uni_con<-unique(conds)
    conts<-c()
    groups<-vector()
    #levels<-c()
    levels<-c(paste0('condition', uni_con))
    code_cont<-matrix(data=0, nrow = round(length(uni_con)*(length(uni_con)-1)/2), ncol = length(uni_con))
    flag=1
    for (i in 1:(length(uni_con)-1)) {
      for(j in (i+1):length(uni_con)){
        conts<-c(conts, paste0('condition', uni_con[j], '-', 'condition', uni_con[i]))
        groups<-rbind(groups, c(uni_con[i], uni_con[j]))
        code_cont[flag,i]=-1
        code_cont[flag,j]=1
        flag=flag+1
      }
    }
    colnames(code_cont)=uni_con
    row.names(code_cont) = conts
    return(list(conts=conts, groups=groups, levels=levels, code_cont=code_cont))
  }

  consts<-all_conts(designs)

  if(platform=='Maxquant' & acq == 'DDA'){

    proteinGroups <- read.table(protein_file_name, quote = "", sep = '\t', header = TRUE)
    infile <- read.table(evidence_file_name, quote = "", sep = '\t', header = TRUE)
    annot<-designs

    raw <- MaxQtoMSstatsFormat(evidence = infile,
                               annotation = annot,
                               proteinGroups = proteinGroups)

    raw$Run<-gsub('\\+', '', raw$Run)
    designs$Run<-gsub('\\+', '', designs$Run)
    for (run in designs$Run) {
      raw$BioReplicate[grep(run, raw$Run)]=designs$BioReplicate[grep(run, designs$Run)]
      raw$Condition[grep(run, raw$Run)]=designs$Condition[grep(run, designs$Run)]
    }

    QuantData <- dataProcess(raw, censoredInt="NA", normalization=normal,
                             MBimpute = imput,
                             featureSubset = inten_type)
  }else if (platform=='FragPipe' & acq == 'DDA'){
    raw <- read_csv(protein_file_name, na = c("", "NA", "0"))
    raw$BioReplicate[is.na(raw$BioReplicate)]=0
    raw$ProteinName <- factor(raw$ProteinName)
    raw$PeptideSequence <- factor(raw$PeptideSequence)

    raw$Run<-gsub('\\+', '', raw$Run)
    designs$Run<-gsub('\\+', '', designs$Run)
    for (run in designs$Run) {
      raw$BioReplicate[grep(run, raw$Run)]=designs$BioReplicate[grep(run, designs$Run)]
      raw$Condition[grep(run, raw$Run)]=designs$Condition[grep(run, designs$Run)]
    }

    QuantData <- dataProcess(raw, censoredInt="NA", normalization=normal,
                             MBimpute = imput,
                             featureSubset = inten_type)
  }else if(platform=='DIANN' | platform=='spt' ){

    annot<-designs

    if(platform=='DIANN'){
      raw <- DIANNtoMSstatsFormat(infile,
                                  annotation = annot
      )
    }else if(platform=='spt'){
      raw <- SpectronauttoMSstatsFormat(infile,
                                        annotation = annot#,
      )
    }

    QuantData <- dataProcess(raw, censoredInt="NA", normalization=normal,
                             MBimpute = imput,
                             featureSubset = inten_type)
  }else if(platform=='Maxquant'){
    # Read in MaxQuant files
    proteinGroups <- read.table(protein_file_name, quote = "", sep = '\t', header = TRUE)
    infile <- read.table(evidence_file_name, quote = "", sep = '\t', header = TRUE)
    annot<-designs

    raw <- MaxQtoMSstatsTMTFormat(evidence = infile,
                                  annotation = annot,
                                  proteinGroups = proteinGroups)

    QuantData<- proteinSummarization(
      raw,
      method = inten_type, #"msstats"(default), "MedianPolish", "Median", "LogSum"
      global_norm = normal,
      reference_norm = reference_norm,

      MBimpute = imput,
    )
  }else if (platform=='FragPipe'){
    frag_msstats<-read_csv(protein_file_name, na = c("", "NA", "0"))
    raw <- PhilosophertoMSstatsTMTFormat(
      frag_msstats,
      annotation=designs)

    QuantData<- proteinSummarization(
      raw,
      method = inten_type, #"msstats"(default), "MedianPolish", "Median", "LogSum"
      global_norm = normal,
      reference_norm = reference_norm,

      MBimpute = imput,

    )
  }
  levels(QuantData$ProteinLevelData$GROUP)

  if (acq=='TMT'){
    testResultOneComparison <- groupComparisonTMT(contrast.matrix = consts$code_cont, data = QuantData)
  }else{
    testResultOneComparison <- groupComparison(contrast.matrix = consts$code_cont, data = QuantData)
  }
  res.MSstats<-testResultOneComparison$ComparisonResult%>% arrange(pvalue)
  colnames(res.MSstats)[1:3]<-c('protein','contrast', 'logFC')

  if(platform=='DIANN'){
    pro_strs<-paste0(infile$Protein.Group, '|', infile$Protein.Ids, '|',
                     infile$Protein.Names)
    uni_pro_strs<-unique(pro_strs)
    idx_uni<-match(uni_pro_strs, pro_strs)

    map_table<-infile[idx_uni, ]
    map_table<-cbind(map_table$Protein.Group, map_table$Protein.Ids,
                     map_table$Protein.Names)

    map_str<-paste0(map_table[,1], '|', map_table[,2], '|',
                    map_table[,3])
    pro_replace<-c()
    for (i in 1:length(res.MSstats$protein)) {
      pro<-res.MSstats$protein[i]
      idx1<-which(map_table[,3]==pro)
      if(length(idx1)>0){
        #pro_split<-strsplit(pro, ';', fixed = T)[[1]]
        if(length(grep('UPS', map_str[idx1[1]]))>0){
          sep_ups<-strsplit(map_table[,1][idx1[1]],';', fixed=t)[[1]]
          cat_prostr=''
          for (j in 1:length(sep_ups)) {
            if(cat_prostr==''){
              cat_prostr<-paste(cat_prostr, paste0('sp|',sep_ups[j]))
            }else{
              cat_prostr<-paste(cat_prostr,';', paste0('sp|',sep_ups[j]))
            }

          }
          pro_replace<-c(pro_replace, cat_prostr)
        }else{
          uniports<-strsplit(map_table[,1][idx1[1]], ';', fixed = T)[[1]]
          prons<-strsplit(map_table[,3][idx1[1]], ';', fixed = T)[[1]]
          cat_prostr=''
          for (j in 1:length(uniports)) {
            if(cat_prostr==''){
              cat_prostr<-paste0(cat_prostr, paste0('sp|',uniports[j],'|',prons[j]))
            }else{
              cat_prostr<-paste0(cat_prostr,';', paste0('sp|',uniports[j],'|',prons[j]))
            }

          }
          pro_replace<-c(pro_replace, cat_prostr)
        }

      }
    }
    res.MSstats$protein<-pro_replace
  }else if(platform=='spt'){
    all_ups<-get_UPS_info()
    pro_strs<-paste0(infile$PG.ProteinGroups, '|', infile$PG.ProteinNames)
    uni_pro_strs<-unique(pro_strs)
    idx_uni<-match(uni_pro_strs, pro_strs)

    map_table<-infile[idx_uni, ]
    map_table<-cbind(map_table$PG.ProteinGroups, map_table$PG.ProteinNames)

    map_str<-paste0(map_table[,1], '|', map_table[,2])
    pro_replace<-c()
    for (i in 1:length(res.MSstats$protein)) {
      pro<-as.character(res.MSstats$protein[i])
      idx1<-which(map_table[,1]==pro)

      if(length(idx1)>0){
        pro_split<-strsplit(pro, ';', fixed = T)[[1]]
        if(length(grep('UPS', map_str[idx1[1]]))>0){
          sep_ups<-strsplit(map_table[,1][idx1[1]],';', fixed=t)[[1]]
          cat_prostr=''
          for (j in 1:length(sep_ups)) {
            idx_ups<-grep(sep_ups[j], all_ups)
            if(cat_prostr==''){
              cat_prostr<-paste(cat_prostr, paste0('sp|',all_ups[idx_ups]))
            }else{
              cat_prostr<-paste(cat_prostr,';', paste0('sp|',all_ups[idx_ups]))
            }

          }
          pro_replace<-c(pro_replace, cat_prostr)
        }else{
          uniports<-strsplit(map_table[,1][idx1[1]], ';', fixed = T)[[1]]
          prons<-strsplit(map_table[,2][idx1[1]], ';', fixed = T)[[1]]
          cat_prostr=''
          for (j in 1:length(uniports)) {
            if(cat_prostr==''){
              cat_prostr<-paste0(cat_prostr, paste0('sp|',uniports[j],'|',prons[j]))
            }else{
              cat_prostr<-paste0(cat_prostr,';', paste0('sp|',uniports[j],'|',prons[j]))
            }

          }
          pro_replace<-c(pro_replace, cat_prostr)
        }

      }
    }
    res.MSstats$protein<-pro_replace
  }

  res_all<-res.MSstats
  res_all<-res_all[order(res_all$pvalue),]
  return(list(dea=res_all, processed=prepro_res$normed))
}

