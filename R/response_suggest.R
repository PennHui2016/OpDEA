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

getFilsug_diann<-function(input){

  data<-wf_dia

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

getFilsug_spt<-function(input){

  data<-wf_spt

  if (!is.null(input$rb41)) {
    data <- data[data$expression_matrix %in% input$rb41,]
  }
  if (!is.null(input$rb42)) {
    data <- data[data$normalization %in% input$rb42,]
  }
  if (!is.null(input$rb43)) {
    data <- data[data$imputation %in% input$rb43,]
  }
  if (!is.null(input$rb44)) {
    data <- data[data$DEA_tool %in% input$rb44,]
  }

  return(data)
}

getFilsug_tmt_fg<-function(input){

  data<-wf_fg_tmt

  if (!is.null(input$rb51)) {
    data <- data[data$expression_matrix %in% input$rb51,]
  }
  if (!is.null(input$rb52)) {
    data <- data[data$normalization %in% input$rb52,]
  }
  if (!is.null(input$rb53)) {
    data <- data[data$imputation %in% input$rb53,]
  }
  if (!is.null(input$rb54)) {
    data <- data[data$DEA_tool %in% input$rb54,]
  }

  return(data)
}

getFilsug_tmt_mq<-function(input){

  data<-wf_mq_tmt

  if (!is.null(input$rb61)) {
    data <- data[data$expression_matrix %in% input$rb61,]
  }
  if (!is.null(input$rb62)) {
    data <- data[data$normalization %in% input$rb62,]
  }
  if (!is.null(input$rb63)) {
    data <- data[data$imputation %in% input$rb63,]
  }
  if (!is.null(input$rb64)) {
    data <- data[data$DEA_tool %in% input$rb64,]
  }

  return(data)
}

get_sug_ens_fg<-function(input){

  ens_cbn<-wf_fg_dda_ens$cbn[1]
  exps=strsplit(ens_cbn, '|', fixed = T)[[1]]
  ens_mth = wf_fg_dda_ens$method[1]

  out_texts = ''
  wfs=vector()
  for (i in 1:length(exps)) {
    idx=which(wf_frag$expression_matrix==exps[i])[1]
    wf=c(wf_frag$expression_matrix[idx], wf_frag$normalization[idx], wf_frag$imputation[idx], wf_frag$DEA_tool[idx])
    wfs<-rbind(wfs, wf)
    out_text = paste0(wf_frag$workflow[idx],':',
    '        expression matrix:',wf_frag$expression_matrix[idx],
    '        normalization:',wf_frag$normalization[idx],
    '        imputation:',wf_frag$imputation[idx],
    '        DEA tool:', wf_frag$DEA_tool[idx])
    if(i==1){
      out_texts = paste0(out_texts, out_text)
    }else{
      out_texts = paste0(out_texts, "<br>", out_text)
    }
  }
  out_texts = paste0(out_texts, "<br>", " The model used for p-value integration should be: ", "<br>", "    Hurdle")

  return(list(text=out_texts, wfs=wfs, ens_mth=ens_mth))
}

get_sug_ens_mq<-function(input){
  ens_cbn<-wf_mq_dda_ens$cbn[1]
  exps=strsplit(ens_cbn, '|', fixed = T)[[1]]
  ens_mth = wf_mq_dda_ens$method[1]

  out_texts = ''
  wfs=vector()
  for (i in 1:length(exps)) {
    idx=which(wf_mq$expression_matrix==exps[i])[1]
    wf=c(wf_mq$expression_matrix[idx], wf_mq$normalization[idx], wf_mq$imputation[idx], wf_mq$DEA_tool[idx])
    wfs<-rbind(wfs, wf)
    out_text = paste0(wf_mq$workflow[idx],':',
                      '        expression matrix:',wf_mq$expression_matrix[idx],
                      '        normalization:',wf_mq$normalization[idx],
                      '        imputation:',wf_mq$imputation[idx],
                      '        DEA tool:', wf_mq$DEA_tool[idx])
    if(i==1){
      out_texts = paste0(out_texts, out_text)
    }else{
      out_texts = paste0(out_texts, "<br>", out_text)
    }
  }
  out_texts = paste0(out_texts, "<br>", " The model used for p-value integration should be: ", "<br>", "    Hurdle")

  return(list(text=out_texts, wfs=wfs, ens_mth=ens_mth))
}

get_sug_ens_diann<-function(input){

  ens_cbn<-wf_diann_dia_ens$cbn[1]
  exps=strsplit(ens_cbn, '|', fixed = T)[[1]]
  ens_mth = wf_diann_dia_ens$method[1]

  out_texts = ''
  wfs=vector()
  for (i in 1:length(exps)) {
    idx=which(wf_dia$expression_matrix==exps[i])[1]
    wf=c(wf_dia$expression_matrix[idx], wf_dia$normalization[idx], wf_dia$imputation[idx], wf_dia$DEA_tool[idx])
    wfs<-rbind(wfs, wf)
    out_text = paste0(wf_dia$workflow[idx],':',
                      '        expression matrix:',wf_dia$expression_matrix[idx],
                      '        normalization:',wf_dia$normalization[idx],
                      '        imputation:',wf_dia$imputation[idx],
                      '        DEA tool:', wf_dia$DEA_tool[idx])
    if(i==1){
      out_texts = paste0(out_texts, out_text)
    }else{
      out_texts = paste0(out_texts, "<br>", out_text)
    }
  }
  out_texts = paste0(out_texts, "<br>", " The model used for p-value integration should be: ", "<br>", "    Hurdle")

  return(list(text=out_texts, wfs=wfs, ens_mth=ens_mth))
}


get_sug_ens_spt<-function(input){

  ens_cbn<-wf_spt_dia_ens$cbn[1]
  exps=strsplit(ens_cbn, '|', fixed = T)[[1]]
  ens_mth = wf_diann_dia_ens$method[1]

  out_texts = ''
  wfs=vector()
  for (i in 1:length(exps)) {
    idx=which(wf_spt$expression_matrix==exps[i])[1]
    wf=c(wf_spt$expression_matrix[idx], wf_spt$normalization[idx], wf_spt$imputation[idx], wf_spt$DEA_tool[idx])
    wfs<-rbind(wfs, wf)
    out_text = paste0(wf_spt$workflow[idx],':',
                      '        expression matrix:',wf_spt$expression_matrix[idx],
                      '        normalization:',wf_spt$normalization[idx],
                      '        imputation:',wf_spt$imputation[idx],
                      '        DEA tool:', wf_spt$DEA_tool[idx])
    if(i==1){
      out_texts = paste0(out_texts, out_text)
    }else{
      out_texts = paste0(out_texts, "<br>", out_text)
    }
  }
  out_texts = paste0(out_texts, "<br>", " The model used for p-value integration should be: ", "<br>", "    Hurdle")

  return(list(text=out_texts, wfs=wfs, ens_mth=ens_mth))
}
