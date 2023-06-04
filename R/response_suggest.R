getFilsug_fg<-function(input){
  
  data<-wf_frag
  
  if (!is.null(input$rb11)) {
    data <- data[data$expression_matrix %in% input$rb11,]
  }
  if (!is.null(input$rb12)) {
    data <- data[data$normalization %in% input$rb12,]
  }
  if (!is.null(input$rb13)) {
    data <- data[data$imputation %in% input$rb13,]
  }
  if (!is.null(input$rb14)) {
    data <- data[data$DEA_tool %in% input$rb14,]
  }
  #print(data[1,][1])
  #wf = data$workflow
  #rwf = metric_frag$pauc001[which(metric_frag$pauc001$workflow==wf),]
  return(data)
}

getFilsug_mq<-function(input){
  
  data<-wf_mq
  
  if (!is.null(input$rb21)) {
    data <- data[data$expression_matrix %in% input$rb21,]
  }
  if (!is.null(input$rb22)) {
    data <- data[data$normalization %in% input$rb22,]
  }
  if (!is.null(input$rb23)) {
    data <- data[data$imputation %in% input$rb23,]
  }
  if (!is.null(input$rb24)) {
    data <- data[data$DEA_tool %in% input$rb24,]
  }

  return(data)
}

getFilsug_dia<-function(input){
  
  data<-wf_mq
  
  if (!is.null(input$rb31)) {
    data <- data[data$expression_matrix %in% input$rb31,]
  }
  if (!is.null(input$rb32)) {
    data <- data[data$normalization %in% input$rb32,]
  }
  if (!is.null(input$rb33)) {
    data <- data[data$imputation %in% input$rb33,]
  }
  if (!is.null(input$rb34)) {
    data <- data[data$DEA_tool %in% input$rb34,]
  }

  return(data)
}


get_sug_ens_fg<-function(input){
  c1=0
  c2=0
  c3=0
  wf_bs_pro = 'DEqMS|fragpipe|MaxLFQ|MinProb|' 
  wf_bs_pep = 'ProteoMM|fragpipe|MaxLFQ|MinDet|max'
  wf_bs_sc = 'plgem|fragpipe||QRILC|div.mean'
 
  idx = ''
  
  if(length(intersect(input$cgFg1, 'pro'))==1){
    c1=1
  }
  if(length(intersect(input$cgFg1, 'pep'))==1){
    c2=1
  }
  if(length(intersect(input$cgFg1, 'sc'))==1){
    c3=1
  }
  if(c1==1 & c2==1 & c3==1){
    idx = which(wf_frag$workflow==wf_bs_pro | wf_frag$workflow==wf_bs_pep | wf_frag$workflow==wf_bs_sc)
  }else if(c1==1 & c2==1 & c3!=1){
    idx = which(wf_frag$workflow==wf_bs_pro | wf_frag$workflow==wf_bs_pep)
  }
  if(length(idx)>1){
  ens_wf = wf_frag[idx,]
  out_texts = ''
  for (i in 1:length(idx)) {
    out_text = paste0(ens_wf$workflow[i],':',
    '        expression matrix:',ens_wf$expression_matrix[i],
    '        normalization:',ens_wf$normalization[i],
    '        imputation:',ens_wf$imputation[i],
    '        DEA tool:', ens_wf$DEA_tool[i])
    if(i==1){
      out_texts = paste0(out_texts, out_text)
    }else{
      out_texts = paste0(out_texts, "<br>", out_text)
    }
  }
  out_texts = paste0(out_texts, "<br>", " The model used for p-value integration should be: ", "<br>", "    Hurdle")
  }else{
    out_texts = 'The ensemle inference may not improve DEA performance, please use single optimal workflow instead!!!'
  }
  return(out_texts)
}

get_sug_ens_mq<-function(input){
  c1=0
  c2=0
  c3=0
  wf_bs_pro = 'DEqMS|maxquant||MLE|center.median' 
  wf_bs_pep = 'ProteoMM|maxquant|LFQ|MinDet|max'
  wf_bs_sc = 'plgem|maxquant||QRILC|div.mean'
  
  idx = ''
  
  if(length(intersect(input$cgMq1, 'pro'))==1){
    c1=1
  }
  if(length(intersect(input$cgMq1, 'pep'))==1){
    c2=1
  }
  if(length(intersect(input$cgMq1, 'sc'))==1){
    c3=1
  }
  if(c1==1 & c2==1 & c3==1){
    idx = which(wf_mq$workflow==wf_bs_pro | wf_mq$workflow==wf_bs_pep | wf_mq$workflow==wf_bs_sc)
  }else if(c1==1 & c2==1 & c3!=1){
    idx = which(wf_mq$workflow==wf_bs_pro | wf_mq$workflow==wf_bs_pep)
  }
  if(length(idx)>1){
    ens_wf = wf_mq[idx,]
    out_texts = ''
    for (i in 1:length(idx)) {
      out_text = paste0(ens_wf$workflow[i],':',
                        '        expression matrix:',ens_wf$expression_matrix[i],
                        '        normalization:',ens_wf$normalization[i],
                        '        imputation:',ens_wf$imputation[i],
                        '        DEA tool:', ens_wf$DEA_tool[i])
      if(i==1){
        out_texts = paste0(out_texts, out_text)
      }else{
        out_texts = paste0(out_texts, "<br>", out_text)
      }
    }
    out_texts = paste0(out_texts, "<br>", " The model used for p-value integration should be: ", "<br>", "    Hurdle")
  }else{
    out_texts = 'The ensemle inference may not improve DEA performance, please use single optimal workflow instead!!!'
  }
  return(out_texts)
}

get_sug_ens_dia<-function(input){
  
  idx = c(1:20)
  if(length(idx)>1){
    ens_wf = wf_dia[idx,]
    out_texts = ''
    for (i in 1:length(idx)) {
      out_text = paste0(ens_wf$workflow[i],':',
                        '        expression matrix:',ens_wf$expression_matrix[i],
                        '        normalization:',ens_wf$normalization[i],
                        '        imputation:',ens_wf$imputation[i],
                        '        DEA tool:', ens_wf$DEA_tool[i])
      if(i==1){
        out_texts = paste0(out_texts, out_text)
      }else{
        out_texts = paste0(out_texts, "<br>", out_text)
      }
    }
    out_texts = paste0(out_texts, "<br>", " The model used for p-value integration should be: ", "<br>", "    set minimum")
  }else{
    out_texts = 'The ensemle inference may not improve DEA performance, please use single optimal workflow instead!!!'
  }
  return(out_texts)
}