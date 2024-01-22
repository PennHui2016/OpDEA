#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
library(shinydashboard)
library(shiny)
library(threejs)
library(DT)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(ggsci)
library("readxl")
library(ggalluvial)
#library(shinyFiles)

app_server <- function(input, output, session) {
  output$menu <- renderMenu({
    sidebarMenu(id = "smenu",
                menuItem("Introduction", tabName = "front"),
                menuItem("Benchmarking",
                         helpText("ranking workflows"),
                         menuSubItem("DDA_LFQ-FragPipe", "models1"),
                         menuSubItem("DDA_LFQ-Maxquant", "models2"),
                         menuSubItem("DIA_LFQ-DIA-NN", "models3"),
                         menuSubItem("DIA_LFQ-Spectronaut", "models4"),
                         menuSubItem("TMT-FragPipe", "models5"),
                         menuSubItem("TMT-Maxquant", "models6")
                ),
                menuItem("Suggestion & DEA",
                         helpText("recommend workflows & DEA"),
                         menuSubItem("DDA", "sug1"),
                         menuSubItem("DIA", "sug2"),
                         menuSubItem("TMT", "sug3")

                ),
                menuItem("Data", tabName = "accData"),
                menuItem("Help", helpText("tourial"),
                         menuSubItem("Toolkit introduction", "Help0"),
                         menuSubItem("Benchmarking", "Help1"),
                         menuSubItem("Suggestion", "Help2"),
                         menuSubItem("Contact", "Help3"),
                         tabName = "Help")
    )
  })

  data1 <- reactive({
    getTabFil_fg(input)
  })
  data2 <- reactive({
    getTabFil_mq(input)
  })
  data3 <- reactive({
    getTabFil_dia(input)
  })

  data4 <- reactive({
    getTabFil_spt(input)
  })
  data5 <- reactive({
    getTabFil_fg_tmt(input)
  })
  data6 <- reactive({
    getTabFil_mq_tmt(input)
  })

  output$tabFrag <- DT::renderDataTable({
    data<- wf_frag
    if (!is.null(input$filter1)) {
      data <- data[data$expression_matrix %in% input$filter1,]
    }
    if (!is.null(input$filter2)) {
      data <- data[data$normalization %in% input$filter2,]
    }
    if (!is.null(input$filter3)) {
      data <- data[data$imputation %in% input$filter3,]
    }
    if (!is.null(input$filter4)) {
      data <- data[data$DEA_tool %in% input$filter4,]
    }
    DT::datatable(data)
  },options = list(paging = FALSE))

  output$boxplot1 <- renderPlot({
    row<-input$tabFrag_rows_selected
    if (length(row) == 1) {
      pt<-plot_metrics(data1()[row,][1][[1]], 'FG_DDA')
      pt
    }
  })

  output$tabMq <- DT::renderDataTable({
    data<- wf_mq
    if (!is.null(input$filter5)) {
      data <- data[data$expression_matrix %in% input$filter5,]
    }
    if (!is.null(input$filter6)) {
      data <- data[data$normalization %in% input$filter6,]
    }
    if (!is.null(input$filter7)) {
      data <- data[data$imputation %in% input$filter7,]
    }
    if (!is.null(input$filter8)) {
      data <- data[data$DEA_tool %in% input$filter8,]
    }
    DT::datatable(data)
  },options = list(paging = FALSE))

  output$boxplot2 <- renderPlot({
    row<-input$tabMq_rows_selected
    if (length(row) == 1) {
      pt<-plot_metrics(data2()[row,][1][[1]], 'MQ_DDA')
      pt
    }
  })

  output$tabDia <- DT::renderDataTable({
    data<- wf_spt
    if (!is.null(input$filter9)) {
      data <- data[data$expression_matrix %in% input$filter9,]
    }
    if (!is.null(input$filter10)) {
      data <- data[data$normalization %in% input$filter10,]
    }
    if (!is.null(input$filter11)) {
      data <- data[data$imputation %in% input$filter11,]
    }
    if (!is.null(input$filter12)) {
      data <- data[data$DEA_tool %in% input$filter12,]
    }
    DT::datatable(data)
  },options = list(paging = FALSE))

  output$boxplot3 <- renderPlot({
    row<-input$tabDia_rows_selected
    if (length(row) == 1) {
      pt<-plot_metrics(data3()[row,][1][[1]], 'DIANN_DIA')
      pt
    }
  })

  output$tabDia1 <- DT::renderDataTable({
    data<- wf_spt
    if (!is.null(input$filter13)) {
      data <- data[data$expression_matrix %in% input$filter13,]
    }
    if (!is.null(input$filter14)) {
      data <- data[data$normalization %in% input$filter14,]
    }
    if (!is.null(input$filter15)) {
      data <- data[data$imputation %in% input$filter15,]
    }
    if (!is.null(input$filter16)) {
      data <- data[data$DEA_tool %in% input$filter16,]
    }
    DT::datatable(data)
  },options = list(paging = FALSE))

  output$boxplot4 <- renderPlot({
    row<-input$tabDia1_rows_selected
    if (length(row) == 1) {
      pt<-plot_metrics(data4()[row,][1][[1]], 'spt_DIA')
      pt
    }
  })


  output$tabTMT1 <- DT::renderDataTable({
    data<- wf_fg_tmt
    if (!is.null(input$filter17)) {
      data <- data[data$expression_matrix %in% input$filter17,]
    }
    if (!is.null(input$filter18)) {
      data <- data[data$normalization %in% input$filter18,]
    }
    if (!is.null(input$filter19)) {
      data <- data[data$imputation %in% input$filter19,]
    }
    if (!is.null(input$filter20)) {
      data <- data[data$DEA_tool %in% input$filter20,]
    }
    DT::datatable(data)
  },options = list(paging = FALSE))

  output$boxplot5 <- renderPlot({
    row<-input$tabTMT1_rows_selected
    if (length(row) == 1) {
      pt<-plot_metrics(data5()[row,][1][[1]], 'FG_TMT')
      pt
    }
  })

  output$tabTMT2 <- DT::renderDataTable({
    data<- wf_mq_tmt
    if (!is.null(input$filter21)) {
      data <- data[data$expression_matrix %in% input$filter21,]
    }
    if (!is.null(input$filter22)) {
      data <- data[data$normalization %in% input$filter22,]
    }
    if (!is.null(input$filter23)) {
      data <- data[data$imputation %in% input$filter23,]
    }
    if (!is.null(input$filter24)) {
      data <- data[data$DEA_tool %in% input$filter24,]
    }
    DT::datatable(data)
  },options = list(paging = FALSE))

  output$boxplot6 <- renderPlot({
    row<-input$tabTMT2_rows_selected
    if (length(row) == 1) {
      pt<-plot_metrics(data6()[row,][1][[1]], 'MQ_TMT')
      pt
    }
  })

  output$fgcv <- renderPlot({
    lopocv_fg_fig
  })
  output$mqcv <- renderPlot({
    lopocv_mq_fig
  })
  output$diacv <- renderPlot({
    lopocv_dia_fig
  })
  output$cvCls <- renderPlot({
    fig_catboost
  })
  output$feaImp <- renderPlot({
    imp_catboost
  })

  output$dD_fg1 <- downloadHandler(
    filename = 'metric_fragpipe_pAUC0.01.csv',
    content = function(file){write.csv(metric_frag$pauc001, file, row.names = FALSE)}
  )
  output$dD_fg2 <- downloadHandler(
    filename = 'metric_fragpipe_pAUC0.05.csv',
    content = function(file){write.csv(metric_frag$pauc005, file, row.names = FALSE)}
  )
  output$dD_fg3 <- downloadHandler(
    filename = 'metric_fragpipe_pAUC0.1.csv',
    content = function(file){write.csv(metric_frag$pauc01, file, row.names = FALSE)}
  )
  output$dD_fg4 <- downloadHandler(
    filename = 'metric_fragpipe_nMCC.csv',
    content = function(file){write.csv(metric_frag$nMCC, file, row.names = FALSE)}
  )
  output$dD_fg5 <- downloadHandler(
    filename = 'metric_fragpipe_G-mean.csv',
    content = function(file){write.csv(metric_frag$Gmean, file, row.names = FALSE)}
  )
  output$dD_fg6 <- downloadHandler(
    filename = 'fragpipe_ranks.csv',
    content = function(file){write.csv(wf_frag, file, row.names = FALSE)}
  )

  output$dD_mq1 <- downloadHandler(
    filename = 'metric_maxquant_pAUC0.01.csv',
    content = function(file){write.csv(metric_mq$pauc001, file, row.names = FALSE)}
  )
  output$dD_mq2 <- downloadHandler(
    filename = 'metric_maxquant_pAUC0.05.csv',
    content = function(file){write.csv(metric_mq$pauc005, file, row.names = FALSE)}
  )
  output$dD_mq3 <- downloadHandler(
    filename = 'metric_maxquant_pAUC0.1.csv',
    content = function(file){write.csv(metric_mq$pauc01, file, row.names = FALSE)}
  )
  output$dD_mq4 <- downloadHandler(
    filename = 'metric_maxquant_nMCC.csv',
    content = function(file){write.csv(metric_mq$nMCC, file, row.names = FALSE)}
  )
  output$dD_mq5 <- downloadHandler(
    filename = 'metric_maxquant_G-mean.csv',
    content = function(file){write.csv(metric_mq$Gmean, file, row.names = FALSE)}
  )
  output$dD_mq6 <- downloadHandler(
    filename = 'maxquant_ranks.csv',
    content = function(file){write.csv(wf_mq, file, row.names = FALSE)}
  )

  output$dD_dia1 <- downloadHandler(
    filename = 'metric_DIANN_pAUC0.01.csv',
    content = function(file){write.csv(metric_dia$pauc001, file, row.names = FALSE)}
  )
  output$dD_dia2 <- downloadHandler(
    filename = 'metric_DIANN_pAUC0.05.csv',
    content = function(file){write.csv(metric_dia$pauc005, file, row.names = FALSE)}
  )
  output$dD_dia3 <- downloadHandler(
    filename = 'metric_DIANN_pAUC0.1.csv',
    content = function(file){write.csv(metric_dia$pauc01, file, row.names = FALSE)}
  )
  output$dD_dia4 <- downloadHandler(
    filename = 'metric_DIANN_nMCC.csv',
    content = function(file){write.csv(metric_dia$nMCC, file, row.names = FALSE)}
  )
  output$dD_dia5 <- downloadHandler(
    filename = 'metric_DIANN_G-mean.csv',
    content = function(file){write.csv(metric_dia$Gmean, file, row.names = FALSE)}
  )
  output$dD_dia6 <- downloadHandler(
    filename = 'DIANN_ranks.csv',
    content = function(file){write.csv(wf_dia, file, row.names = FALSE)}
  )

  data_fg_sug_prefer <- reactive({
    getFilsug_fg(input)
  })

  output$fg_sug_prefer <- renderUI({

    HTML(paste0("<pre>",data_fg_sug_prefer()$workflow[1],':',
                '        expression matrix:',data_fg_sug_prefer()$expression_matrix[1],
                '        normalization:',data_fg_sug_prefer()$normalization[1],
                '        imputation:',data_fg_sug_prefer()$imputation[1],
                '        DEA tool:', data_fg_sug_prefer()$DEA_tool[1],"<pre>"))#,


  })

  output$volc_fg1<-renderUI({
    if (input$runDEA_fg1==0) return()
    #source('./R/run_DEA_single.R')
    if(!is.null(input$rb11)){
      mat = input$rb11
    }else{
      mat =data_fg_sug_prefer()$expression_matrix[1]
    }
    if(!is.null(input$rb12)){
      norm = input$rb12
    }else{
      norm =data_fg_sug_prefer()$normalization[1]
    }
    if(!is.null(input$rb13)){
      imp = input$rb13
    }else{
      imp =data_fg_sug_prefer()$imputation[1]
    }
    if(!is.null(input$rb14)){
      dea = input$rb14
    }else{
      dea =data_fg_sug_prefer()$DEA_tool[1]
    }

    raw_file<-gsub("\\\\", "/",input$upload_fg_raw)
    evid_path<-gsub("\\\\", "/",input$upload_fg_ion)
    design_file<-gsub("\\\\", "/", input$upload_fg_dg)

    logFC=input$logFC_fg1
    qval = input$adjp_fg1

    python_path = input$python_path_fg1

    t1 = proc.time()[3]
    #tic()

    DEA_res = run_DEA_fg_DDA(input, mat, norm, imp, dea, raw_file, evid_path, design_file, logFC, qval, python_path)
    #toc()
    t2 = proc.time()[3]

    output$volc_fg1 <- renderPlot({
      DEA_res$p
    })

    output$DEA_res_fg1 <- renderText(paste0('Finished with: ', as.character(round(t2-t1, 2)), 's \n;', 'Result files are stored at: ', DEA_res$zipfile))
    #   downloadHandler(
    #   filename <- function() {
    #     paste("output", "zip", sep=".")
    #   },
    #
    #   content <- function(file) {
    #     file.copy(DEA_res$zipfile, file)
    #
    #   },
    #   contentType = "application/zip"
    # )

  })

  output$volc_fg1_1<-renderUI({
    if (input$runDEA_fg1_1==0) return()
    #source('./R/run_DEA_single.R')
    if(!is.null(input$rb11)){
      mat_fg1 = input$rb11
    }else{
      mat_fg1 =data_fg_sug_prefer()$expression_matrix[1]
    }
    if(!is.null(input$rb12)){
      norm_fg1 = input$rb12
    }else{
      norm_fg1 =data_fg_sug_prefer()$normalization[1]
    }
    if(!is.null(input$rb13)){
      imp_fg1 = input$rb13
    }else{
      imp_fg1 =data_fg_sug_prefer()$imputation[1]
    }
    if(!is.null(input$rb14)){
      dea_fg1 = input$rb14
    }else{
      dea_fg1 =data_fg_sug_prefer()$DEA_tool[1]
    }

    raw_file_fg1<-input$upload_fg_raw1
    evid_path_fg1<-input$upload_fg_ion1
    design_file_fg1<-input$upload_fg_dg1
    logFC_fg1=input$logFC_fg1_1
    qval_fg1 = input$adjp_fg1_1
    python_path1 = ''#input$python_path_fg1_1
    write.table(data.frame(a=c('1','start')), "log.txt", sep="\t", col.names=FALSE, row.names=FALSE, append=TRUE)
    t1 = proc.time()[3]
    DEA_res_fg1 = run_DEA_fg_DDA(input, mat_fg1, norm_fg1, imp_fg1, dea_fg1, raw_file_fg1, evid_path_fg1, design_file_fg1, logFC_fg1, qval_fg1, python_path1)
    t2 = proc.time()[3]
    write.table(data.frame(a=c('2','end')), "log.txt", sep="\t", col.names=FALSE, row.names=FALSE, append=TRUE)
    output$volc_fg1_1 <- renderPlot({
      DEA_res_fg1$p
    })

    output$DEA_res_fg1_1 <- renderText(paste0('Finished with: ', as.character(round(t2-t1, 2)), 's \n;', 'Result files are stored at: ', DEA_res_fg1$zipfile))

    #   downloadHandler(
    #   filename <- function() {
    #     paste("output", "zip", sep=".")
    #   },
    #
    #   content <- function(file) {
    #     file.copy(DEA_res_fg1$zipfile, file)
    #
    #   },
    #   contentType = "application/zip"
    # )

  })

  output$fg_sug_ens <- renderUI({
    ens_sug = get_sug_ens_fg(input)

    if (input$suggest2==0) return()
    HTML(paste0("<pre>",ens_sug$text,"<pre>"))


  })

  output$volc_fg_ens<-renderUI({

    ens_sug = get_sug_ens_fg(input)
    raw_file_fg_ens<-input$upload_fg_raw_ens
    evid_path_fg_ens<-input$upload_fg_ion_ens
    design_file_fg_ens<-input$upload_fg_dg_ens
    logFC_fg_ens=as.numeric(input$logFC_fg_ens)
    qval_fg_ens = as.numeric(input$adjp_fg_ens)

    python_path_ens = ''#input$python_path_fg_ens

    t1 = proc.time()[3]

    des_res_ens<-run_DEA_fg_ens(ens_sug$wfs, ens_sug$ens_mth, ens_sug$opr, raw_file_fg_ens, evid_path_fg_ens, design_file_fg_ens, logFC_fg_ens, qval_fg_ens, python_path_ens)

    t2 = proc.time()[3]

    output$volc_fg_ens <- renderPlot({
      des_res_ens$p
    })

    output$DEA_res_fg_ens <- renderText(paste0('Finished with: ', as.character(round(t2-t1, 2)), 's \n;', 'Result files are stored at: ', des_res_ens$zipfile))

    #   downloadHandler(
    #   filename <- function() {
    #     paste("output", "zip", sep=".")
    #   },
    #
    #   content <- function(file) {
    #     file.copy(des_res_ens$zipfile, file)
    #
    #   },
    #   contentType = "application/zip"
    # )
  })


  data_mq_sug_prefer <- reactive({
    getFilsug_mq(input)
  })

  output$mq_sug_prefer <- renderUI({
    if (input$rb2=='N' | input$suggest3==0) return()
    HTML(paste0("<pre>",data_mq_sug_prefer()$workflow[1],':',
                '        expression matrix:',data_mq_sug_prefer()$expression_matrix[1],
                '        normalization:',data_mq_sug_prefer()$normalization[1],
                '        imputation:',data_mq_sug_prefer()$imputation[1],
                '        DEA tool:', data_mq_sug_prefer()$DEA_tool[1],"<pre>"))
  })

  output$volc_mq<-renderUI({
    if (input$runDEA_mq==0) return()
    #source('./R/run_DEA_single.R')
    if(!is.null(input$rb21)){
      mat = input$rb21
    }else{
      mat =data_mq_sug_prefer()$expression_matrix[1]
    }
    if(!is.null(input$rb22)){
      norm = input$rb22
    }else{
      norm =data_mq_sug_prefer()$normalization[1]
    }
    if(!is.null(input$rb23)){
      imp = input$rb23
    }else{
      imp =data_mq_sug_prefer()$imputation[1]
    }
    if(!is.null(input$rb24)){
      dea = input$rb24
    }else{
      dea =data_mq_sug_prefer()$DEA_tool[1]
    }

    raw_file<-input$upload_mq_pro
    evid_path<-input$upload_mq_evi
    design_file<-input$upload_mq_dg

    logFC=input$logFC_mq
    qval = input$adjp_mq

    python_path_mq=''#input$python_path_mq

    t1 = proc.time()[3]

    DEA_res_mq0 = run_DEA_mq_DDA(input, mat, norm, imp, dea, raw_file, evid_path, design_file, logFC, qval, python_path_mq)

    t2 = proc.time()[3]

    output$volc_mq <- renderPlot({
      DEA_res_mq0$p
    })

    output$DEA_res_mq <- renderText(paste0('Finished with: ', as.character(round(t2-t1, 2)), 's \n;', 'Result files are stored at: ', DEA_res_mq0$zipfile))
    #   downloadHandler(
    #   filename <- function() {
    #     paste("output", "zip", sep=".")
    #   },
    #
    #   content <- function(file) {
    #     file.copy(DEA_res_mq0$zipfile, file)
    #
    #   },
    #   contentType = "application/zip"
    # )

  })

  output$volc_mq1<-renderUI({
    if (input$runDEA_mq1==0) return()
    #source('./R/run_DEA_single.R')
    if(!is.null(input$rb21)){
      mat_mq1 = input$rb21
    }else{
      mat_mq1 =data_mq_sug_prefer()$expression_matrix[1]
    }
    if(!is.null(input$rb22)){
      norm_mq1 = input$rb22
    }else{
      norm_mq1 =data_mq_sug_prefer()$normalization[1]
    }
    if(!is.null(input$rb23)){
      imp_mq1 = input$rb23
    }else{
      imp_mq1 =data_mq_sug_prefer()$imputation[1]
    }
    if(!is.null(input$rb24)){
      dea_mq1 = input$rb24
    }else{
      dea_mq1 =data_mq_sug_prefer()$DEA_tool[1]
    }

    raw_file_mq1<-input$upload_mq_pro1
    evid_path_mq1<-input$upload_mq_evi1
    design_file_mq1<-input$upload_mq_dg1
    logFC_mq1=input$logFC_mq1
    qval_mq1 = input$adjp_mq1
    python_path_mq1=''#input$python_path_mq1

    t1 = proc.time()[3]
    DEA_res_mq01 = run_DEA_mq_DDA(input, mat_mq1, norm_mq1, imp_mq1, dea_mq1, raw_file_mq1, evid_path_mq1, design_file_mq1, logFC_mq1, qval_mq1, python_path_mq1)
    t2 = proc.time()[3]

    output$volc_mq1 <- renderPlot({
      DEA_res_mq01$p
    })

    output$DEA_res_mq1 <- renderText(paste0('Finished with: ', as.character(round(t2-t1, 2)), 's \n;', 'Result files are stored at: ', DEA_res_mq01$zipfile))

    #   downloadHandler(
    #   filename <- function() {
    #     paste("output", "zip", sep=".")
    #   },
    #
    #   content <- function(file) {
    #     file.copy(DEA_res_mq01$zipfile, file)
    #
    #   },
    #   contentType = "application/zip"
    # )

  })

  output$mq_sug_ens <- renderUI({
    ens_sug_mq = get_sug_ens_mq(input)
    if (input$suggest4==0) return()
    HTML(paste0("<pre>",ens_sug_mq$text,"<pre>"))
  })

  output$volc_mq_ens<-renderUI({
    ens_sug_mq = get_sug_ens_mq(input)
    raw_file_mq_ens<-input$upload_mq_pro_ens
    evid_path_mq_ens<-input$upload_mq_evi_ens
    design_file_mq_ens<-input$upload_mq_dg_ens
    logFC_mq_ens=as.numeric(input$logFC_mq_ens)
    qval_mq_ens = as.numeric(input$adjp_mq_ens)

    python_path_mq_ens=''#input$python_path_mq_ens

    t1 = proc.time()[3]
    des_res_ens_mq<-run_DEA_mq_ens(ens_sug_mq$wfs, ens_sug_mq$ens_mth, ens_sug_mq$opr, raw_file_mq_ens, evid_path_mq_ens, design_file_mq_ens, logFC_mq_ens, qval_mq_ens, python_path_mq_ens)
    t2 = proc.time()[3]

    output$volc_mq_ens <- renderPlot({
      des_res_ens_mq$p
    })

    output$DEA_res_mq_ens <- renderText(paste0('Finished with: ', as.character(round(t2-t1, 2)), 's \n;', 'Result files are stored at: ', des_res_ens_mq$zipfile))
    #   downloadHandler(
    #   filename <- function() {
    #     paste("output", "zip", sep=".")
    #   },
    #
    #   content <- function(file) {
    #     file.copy(des_res_ens_mq$zipfile, file)
    #
    #   },
    #   contentType = "application/zip"
    # )
  })

  data_diann_sug_prefer <- reactive({
    getFilsug_diann(input)
  })
  output$diann_sug_prefer <- renderUI({

    if (input$rb3=='N' | input$suggest5==0) return()
    HTML(paste0("<pre>",data_diann_sug_prefer()$workflow[1],':',
                '        expression matrix:',data_diann_sug_prefer()$expression_matrix[1],
                '        normalization:',data_diann_sug_prefer()$normalization[1],
                '        imputation:',data_diann_sug_prefer()$imputation[1],
                '        DEA tool:', data_diann_sug_prefer()$DEA_tool[1],"<pre>"))
  })

  output$volc_diann<-renderUI({
    if (input$runDEA_diann==0) return()
    #source('./R/run_DEA_single.R')
    if(!is.null(input$rb31)){
      mat = input$rb31
    }else{
      mat =data_diann_sug_prefer()$expression_matrix[1]
    }
    if(!is.null(input$rb32)){
      norm = input$rb32
    }else{
      norm =data_diann_sug_prefer()$normalization[1]
    }
    if(!is.null(input$rb33)){
      imp = input$rb33
    }else{
      imp =data_diann_sug_prefer()$imputation[1]
    }
    if(!is.null(input$rb34)){
      dea = input$rb34
    }else{
      dea =data_diann_sug_prefer()$DEA_tool[1]
    }

    raw_file<-'NULL'
    evid_path<-input$upload_diann_report
    design_file<-input$upload_diann_dg

    logFC = input$logFC_diann
    qval = input$adjp_diann

    python_path_diann=''#input$python_path_diann

    t1 = proc.time()[3]
    DEA_res_diann0 = run_DEA_diann_DIA(input, mat, norm, imp, dea, raw_file, evid_path, design_file, logFC, qval, python_path_diann)
    t2 = proc.time()[3]

    output$volc_diann <- renderPlot({
      DEA_res_diann0$p
    })

    output$DEA_res_diann <- renderText(paste0('Finished with: ', as.character(round(t2-t1, 2)), 's \n;', 'Result files are stored at: ', DEA_res_diann0$zipfile))
    #   downloadHandler(
    #   filename <- function() {
    #     paste("output", "zip", sep=".")
    #   },
    #
    #   content <- function(file) {
    #     file.copy(DEA_res_diann0$zipfile, file)
    #
    #   },
    #   contentType = "application/zip"
    # )

  })

  output$volc_diann1<-renderUI({

    if (input$runDEA_diann1==0) return()
    #source('./R/run_DEA_single.R')
    if(!is.null(input$rb31)){
      mat_diann1 = input$rb31
    }else{
      mat_diann1 =data_diann_sug_prefer()$expression_matrix[1]
    }
    if(!is.null(input$rb32)){
      norm_diann1 = input$rb32
    }else{
      norm_diann1 =data_diann_sug_prefer()$normalization[1]
    }
    if(!is.null(input$rb33)){
      imp_diann1 = input$rb33
    }else{
      imp_diann1 =data_diann_sug_prefer()$imputation[1]
    }
    if(!is.null(input$rb34)){
      dea_diann1 = input$rb34
    }else{
      dea_diann1 =data_diann_sug_prefer()$DEA_tool[1]
    }

    raw_file_diann1<-'NULL'
    evid_path_diann1<-input$upload_diann_report1
    design_file_diann1<-input$upload_diann_dg1
    logFC_diann1=input$logFC_diann1
    qval_diann1 = input$adjp_diann1
    python_path_diann1=''#input$python_path_diann1

    t1 = proc.time()[3]
    DEA_res_diann01 = run_DEA_diann_DIA(input, mat_diann1, norm_diann1, imp_diann1, dea_diann1, raw_file_diann1, evid_path_diann1, design_file_diann1, logFC_diann1, qval_diann1, python_path_diann1)
    t2 = proc.time()[3]

    output$volc_diann1 <- renderPlot({
      DEA_res_diann01$p
    })

    output$DEA_res_diann1 <- renderText(paste0('Finished with: ', as.character(round(t2-t1, 2)), 's \n;', 'Result files are stored at: ', DEA_res_diann01$zipfile))
    #   downloadHandler(
    #   filename <- function() {
    #     paste("output", "zip", sep=".")
    #   },
    #
    #   content <- function(file) {
    #     file.copy(DEA_res_diann01$zipfile, file)
    #
    #   },
    #   contentType = "application/zip"
    # )

  })

  output$diann_sug_ens <- renderUI({
    ens_sug_diann = get_sug_ens_diann(input)
    if (input$suggest6==0) return()
    HTML(paste0("<pre>",ens_sug_diann$text,"<pre>"))
  })

  output$volc_diann_ens<-renderUI({

    ens_sug_diann = get_sug_ens_diann(input)
    raw_file_diann_ens<-'NULL'
    evid_path_diann_ens<-input$upload_diann_ens
    design_file_diann_ens<-input$upload_diann_dg_ens
    logFC_diann_ens=as.numeric(input$logFC_diann_ens)
    qval_diann_ens = as.numeric(input$adjp_diann_ens)
    python_path_diann_ens=''#input$python_path_diann_ens

    t1 = proc.time()[3]
    des_res_ens_diann<-run_DEA_diann_ens(ens_sug_diann$wfs, ens_sug_diann$ens_mth, ens_sug_diann$opr, raw_file_diann_ens, evid_path_diann_ens, design_file_diann_ens, logFC_diann_ens, qval_diann_ens, python_path_diann_ens)
    t2 = proc.time()[3]

    output$volc_diann_ens <- renderPlot({
      des_res_ens_diann$p
    })

    output$DEA_res_diann_ens <- renderText(paste0('Finished with: ', as.character(round(t2-t1, 2)), 's \n;', ' Result files are stored at:', des_res_ens_diann$zipfile))
    #   downloadHandler(
    #   filename <- function() {
    #     paste("output", "zip", sep=".")
    #   },
    #
    #   content <- function(file) {
    #     file.copy(des_res_ens_diann$zipfile, file)
    #
    #   },
    #   contentType = "application/zip"
    # )
  })

  data_spt_sug_prefer <- reactive({
    getFilsug_spt(input)
  })
  output$spt_sug_prefer <- renderUI({

    if (input$rb4=='N' | input$suggest7==0) return()
    HTML(paste0("<pre>",data_spt_sug_prefer()$workflow[1],':',
                '        expression matrix:',data_spt_sug_prefer()$expression_matrix[1],
                '        normalization:',data_spt_sug_prefer()$normalization[1],
                '        imputation:',data_spt_sug_prefer()$imputation[1],
                '        DEA tool:', data_spt_sug_prefer()$DEA_tool[1],"<pre>"))
  })

  output$volc_spt<-renderUI({
    if (input$runDEA_spt==0) return()
    #source('./R/run_DEA_single.R')
    if(!is.null(input$rb41)){
      mat = input$rb41
    }else{
      mat =data_spt_sug_prefer()$expression_matrix[1]
    }
    if(!is.null(input$rb42)){
      norm = input$rb42
    }else{
      norm =data_spt_sug_prefer()$normalization[1]
    }
    if(!is.null(input$rb43)){
      imp = input$rb43
    }else{
      imp =data_spt_sug_prefer()$imputation[1]
    }
    if(!is.null(input$rb44)){
      dea = input$rb44
    }else{
      dea =data_spt_sug_prefer()$DEA_tool[1]
    }

    raw_file<-'NULL'
    evid_path<-input$upload_spt_report
    design_file<-input$upload_spt_dg

    logFC = input$logFC_spt
    qval = input$adjp_spt

    python_path_spt=''#input$python_path_spt

    t1 = proc.time()[3]
    DEA_res_spt0 = run_DEA_spt_DIA(input, mat, norm, imp, dea, raw_file, evid_path, design_file, logFC, qval, python_path_spt)
    t2 = proc.time()[3]

    output$volc_spt <- renderPlot({
      DEA_res_spt0$p
    })

    output$DEA_res_spt <- renderText(paste0('Finished with: ', as.character(round(t2-t1, 2)), 's \n;', ' Result files are stored at: ', DEA_res_spt0$zipfile))
    #   downloadHandler(
    #   filename <- function() {
    #     paste("output", "zip", sep=".")
    #   },
    #
    #   content <- function(file) {
    #     file.copy(DEA_res_spt0$zipfile, file)
    #
    #   },
    #   contentType = "application/zip"
    # )

  })

  output$volc_spt1<-renderUI({

    if (input$runDEA_spt1==0) return()
    #source('./R/run_DEA_single.R')
    if(!is.null(input$rb41)){
      mat_spt1 = input$rb41
    }else{
      mat_spt1 =data_spt_sug_prefer()$expression_matrix[1]
    }
    if(!is.null(input$rb42)){
      norm_spt1 = input$rb42
    }else{
      norm_spt1 =data_spt_sug_prefer()$normalization[1]
    }
    if(!is.null(input$rb43)){
      imp_spt1 = input$rb43
    }else{
      imp_spt1 =data_spt_sug_prefer()$imputation[1]
    }
    if(!is.null(input$rb44)){
      dea_spt1 = input$rb44
    }else{
      dea_spt1 =data_spt_sug_prefer()$DEA_tool[1]
    }

    raw_file_spt1<-'NULL'
    evid_path_spt1<-input$upload_spt_report1
    design_file_spt1<-input$upload_spt_dg1
    logFC_spt1=input$logFC_spt1
    qval_spt1 = input$adjp_spt1
    python_path_spt1=''#input$python_path_spt1

    t1 = proc.time()[3]
    DEA_res_spt01 = run_DEA_spt_DIA(input, mat_spt1, norm_spt1, imp_spt1, dea_spt1, raw_file_spt1, evid_path_spt1, design_file_spt1, logFC_spt1, qval_spt1, python_path_spt1)
    t2 = proc.time()[3]

    output$volc_spt1 <- renderPlot({
      DEA_res_spt01$p
    })

    output$DEA_res_spt1 <- renderText(paste0('Finished with: ', as.character(round(t2-t1, 2)), 's \n;', ' Result files are stored at: ', DEA_res_spt01$zipfile))
    #   downloadHandler(
    #   filename <- function() {
    #     paste("output", "zip", sep=".")
    #   },
    #
    #   content <- function(file) {
    #     file.copy(DEA_res_spt01$zipfile, file)
    #
    #   },
    #   contentType = "application/zip"
    # )

  })

  output$spt_sug_ens <- renderUI({
    ens_sug_spt = get_sug_ens_spt(input)
    if (input$suggest6==0) return()
    HTML(paste0("<pre>",ens_sug_spt$text,"<pre>"))
  })

  output$volc_spt_ens<-renderUI({

    ens_sug_spt = get_sug_ens_spt(input)
    raw_file_spt_ens<-'NULL'
    evid_path_spt_ens<-input$upload_spt_ens
    design_file_spt_ens<-input$upload_spt_dg_ens
    logFC_spt_ens=as.numeric(input$logFC_spt_ens)
    qval_spt_ens = as.numeric(input$adjp_spt_ens)
    python_path_spt_ens=''#input$python_path_spt_ens

    t1 = proc.time()[3]
    des_res_ens_spt<-run_DEA_spt_ens(ens_sug_spt$wfs, ens_sug_spt$ens_mth, ens_sug_spt$opr, raw_file_spt_ens, evid_path_spt_ens, design_file_spt_ens, logFC_spt_ens, qval_spt_ens, python_path_spt_ens)
    t2 = proc.time()[3]

    output$volc_spt_ens <- renderPlot({
      des_res_ens_spt$p
    })

    output$DEA_res_spt_ens <- renderText(paste0('Finished with: ', as.character(round(t2-t1, 2)), 's \n;', ' Result files are stored at: ', des_res_ens_spt$zipfile))
    #   downloadHandler(
    #   filename <- function() {
    #     paste("output", "zip", sep=".")
    #   },
    #
    #   content <- function(file) {
    #     file.copy(des_res_ens_spt$zipfile, file)
    #
    #   },
    #   contentType = "application/zip"
    # )
  })

  ######### FG_TMT DEA
  data_tmt_fg_sug_prefer <- reactive({
    getFilsug_tmt_fg(input)
  })
  output$fg_tmt_sug_prefer <- renderUI({

    if (input$rb5=='N' | input$suggest9==0) return()
    HTML(paste0("<pre>",data_tmt_fg_sug_prefer()$workflow[1],':',
                '        expression matrix:',data_tmt_fg_sug_prefer()$expression_matrix[1],
                '        normalization:',data_tmt_fg_sug_prefer()$normalization[1],
                '        imputation:',data_tmt_fg_sug_prefer()$imputation[1],
                '        DEA tool:', data_tmt_fg_sug_prefer()$DEA_tool[1],"<pre>"))
  })

  output$volc_tmt_fg<-renderUI({
    if (input$runDEA_tmt_fg==0) return()
    #source('./R/run_DEA_single.R')
    if(!is.null(input$rb51)){
      mat = input$rb51
    }else{
      mat =data_tmt_fg_sug_prefer()$expression_matrix[1]
    }
    if(!is.null(input$rb52)){
      norm = input$rb52
    }else{
      norm =data_tmt_fg_sug_prefer()$normalization[1]
    }
    if(!is.null(input$rb53)){
      imp = input$rb53
    }else{
      imp =data_tmt_fg_sug_prefer()$imputation[1]
    }
    if(!is.null(input$rb54)){
      dea = input$rb54
    }else{
      dea =data_tmt_fg_sug_prefer()$DEA_tool[1]
    }

    #raw_file<-'NULL'
    abd_file<-input$upload_abd
    ratio_path<-input$upload_ratio
    phi_path<-input$upload_phi
    design_file<-input$upload_tmt_fg_dg

    logFC = input$logFC_tmt_fg
    qval = input$adjp_tmt_fg

    #python_path_tmt_fg=input$python_path_tmt_fg
    t1 = proc.time()[3]
    DEA_res_tmt_fg0 = run_DEA_fg_TMT(input, mat, norm, imp, dea, abd_file, ratio_path, phi_path, design_file, logFC, qval)
    t2 = proc.time()[3]

    output$volc_tmt_fg <- renderPlot({
      DEA_res_tmt_fg0$p
    })

    output$DEA_res_tmt_fg <- renderText(paste0('Finished with: ', as.character(round(t2-t1, 2)), 's \n;', ' Result files are stored at: ', DEA_res_tmt_fg0$zipfile))
    #   downloadHandler(
    #   filename <- function() {
    #     paste("output", "zip", sep=".")
    #   },
    #
    #   content <- function(file) {
    #     file.copy(DEA_res_tmt_fg0$zipfile, file)
    #
    #   },
    #   contentType = "application/zip"
    # )

  })

  output$volc_tmt_fg1<-renderUI({

    if (input$runDEA_tmt_fg1==0) return()
    #source('./R/run_DEA_single.R')
    if(!is.null(input$rb51)){
      mat_tmt_fg1 = input$rb51
    }else{
      mat_tmt_fg1 =data_tmt_fg_sug_prefer()$expression_matrix[1]
    }
    if(!is.null(input$rb52)){
      norm_tmt_fg1 = input$rb52
    }else{
      norm_tmt_fg1 =data_tmt_fg_sug_prefer()$normalization[1]
    }
    if(!is.null(input$rb53)){
      imp_tmt_fg1 = input$rb53
    }else{
      imp_tmt_fg1 =data_tmt_fg_sug_prefer()$imputation[1]
    }
    if(!is.null(input$rb54)){
      dea_tmt_fg1 = input$rb54
    }else{
      dea_tmt_fg1 =data_tmt_fg_sug_prefer()$DEA_tool[1]
    }

    abd_file1<-input$upload_abd1
    ratio_path1<-input$upload_ratio1
    phi_path1<-input$upload_phi1
    design_file1<-input$upload_tmt_fg_dg1

    logFC1 = input$logFC_tmt_fg1
    qval1 = input$adjp_tmt_fg1

    t1 = proc.time()[3]
    DEA_res_tmt_fg01 = run_DEA_fg_TMT(input, mat_tmt_fg1, norm_tmt_fg1, imp_tmt_fg1, dea_tmt_fg1, abd_file1, ratio_path1, phi_path1, design_file1, logFC1, qval1)
    t2 = proc.time()[3]

    output$volc_tmt_fg1 <- renderPlot({
      DEA_res_tmt_fg01$p
    })

    output$DEA_res_tmt_fg1 <- renderText(paste0('Finished with: ', as.character(round(t2-t1, 2)), 's \n;', ' Result files are stored at: ', DEA_res_tmt_fg01$zipfile))
    #   downloadHandler(
    #   filename <- function() {
    #     paste("output", "zip", sep=".")
    #   },
    #
    #   content <- function(file) {
    #     file.copy(DEA_res_tmt_fg01$zipfile, file)
    #
    #   },
    #   contentType = "application/zip"
    # )

  })

  ######### MQ_TMT DEA
  data_tmt_mq_sug_prefer <- reactive({
    getFilsug_tmt_mq(input)
  })
  output$mq_tmt_sug_prefer <- renderUI({

    if (input$rb6=='N' | input$suggest11==0) return()
    HTML(paste0("<pre>",data_tmt_mq_sug_prefer()$workflow[1],':',
                '        expression matrix:',data_tmt_mq_sug_prefer()$expression_matrix[1],
                '        normalization:',data_tmt_mq_sug_prefer()$normalization[1],
                '        imputation:',data_tmt_mq_sug_prefer()$imputation[1],
                '        DEA tool:', data_tmt_mq_sug_prefer()$DEA_tool[1],"<pre>"))
  })

  output$volc_tmt_mq<-renderUI({
    if (input$runDEA_tmt_mq==0) return()
    #source('./R/run_DEA_single.R')

    mat = 'intensity'
    if(!is.null(input$rb62)){
      norm = input$rb62
    }else{
      norm =data_tmt_mq_sug_prefer()$normalization[1]
    }
    if(!is.null(input$rb63)){
      imp = input$rb63
    }else{
      imp =data_tmt_mq_sug_prefer()$imputation[1]
    }
    if(!is.null(input$rb64)){
      dea = input$rb64
    }else{
      dea =data_tmt_mq_sug_prefer()$DEA_tool[1]
    }

    raw_file_fg_tmt<-input$upload_tmt_mq_pro
    evid_path_fg_tmt<-input$upload_tmt_mq_evi

    design_file_fg_tmt<-input$upload_tmt_mq_dg

    logFC_fg_tmt = input$logFC_tmt_mq
    qval_fg_tmt= input$adjp_tmt_mq

    t1 = proc.time()[3]
    DEA_res_tmt_mq0 = run_DEA_mq_TMT(input, mat, norm, imp, dea, raw_file_fg_tmt, evid_path_fg_tmt, design_file_fg_tmt, logFC_fg_tmt, qval_fg_tmt)
    t2 = proc.time()[3]

    output$volc_tmt_mq <- renderPlot({
      DEA_res_tmt_mq0$p
    })

    output$DEA_res_tmt_mq <- renderText(paste0('Finished with: ', as.character(round(t2-t1, 2)), 's \n;', ' Result files are stored at: ', DEA_res_tmt_mq0$zipfile))
    #   downloadHandler(
    #   filename <- function() {
    #     paste("output", "zip", sep=".")
    #   },
    #
    #   content <- function(file) {
    #     file.copy(DEA_res_tmt_mq0$zipfile, file)
    #
    #   },
    #   contentType = "application/zip"
    # )

  })

  output$volc_tmt_mq1<-renderUI({

    if (input$runDEA_tmt_mq1==0) return()
    #source('./R/run_DEA_single.R')

    mat = 'intensity'
    if(!is.null(input$rb62)){
      norm_tmt_mq1 = input$rb62
    }else{
      norm_tmt_mq1 =data_tmt_mq_sug_prefer()$normalization[1]
    }
    if(!is.null(input$rb63)){
      imp_tmt_mq1 = input$rb63
    }else{
      imp_tmt_mq1 =data_tmt_mq_sug_prefer()$imputation[1]
    }
    if(!is.null(input$rb64)){
      dea_tmt_mq1 = input$rb64
    }else{
      dea_tmt_mq1 =data_tmt_mq_sug_prefer()$DEA_tool[1]
    }

    raw_file1<-input$upload_tmt_mq_pro1
    evid_path1<-input$upload_tmt_mq_evi1

    design_file1<-input$upload_tmt_mq_dg1

    logFC1 = input$logFC_tmt_mq1
    qval1 = input$adjp_tmt_mq1

    t1 = proc.time()[3]
    DEA_res_tmt_mq01 = run_DEA_mq_TMT(input, mat, norm_tmt_mq1, imp_tmt_mq1, dea_tmt_mq1, raw_file1, evid_path1, design_file1, logFC1, qval1)
    t2 = proc.time()[3]

    output$volc_tmt_mq1 <- renderPlot({
      DEA_res_tmt_mq01$p
    })

    output$DEA_res_tmt_mq1 <- renderText(paste0('Finished with: ', as.character(round(t2-t1, 2)), 's \n;', ' Result files are stored at: ', DEA_res_tmt_mq01$zipfile))
    #   downloadHandler(
    #   filename <- function() {
    #     paste("output", "zip", sep=".")
    #   },
    #
    #   content <- function(file) {
    #     file.copy(DEA_res_tmt_mq01$zipfile, file)
    #
    #   },
    #   contentType = "application/zip"
    # )

  })

  # output$DDA_exp <- downloadHandler(
  #   filename <- function() {
  #     paste("DDA_expression_matrices", "zip", sep=".")
  #   },
  #
  #   content <- function(file) {
  #     file.copy('./R/example/DDA.zip', file)
  #
  #   },
  #   contentType = "application/zip"
  # )
  #
  # output$DIA_exp <- downloadHandler(
  #   filename <- function() {
  #     paste("DIA_expression_matrices", "zip", sep=".")
  #   },
  #
  #   content <- function(file) {
  #     file.copy('./R/example/DIA.zip', file)
  #
  #   },
  #   contentType = "application/zip"
  # )
  #
  # output$TMT_exp <- downloadHandler(
  #   filename <- function() {
  #     paste("TMT_expression_matrices", "zip", sep=".")
  #   },
  #
  #   content <- function(file) {
  #     file.copy('./R/example/TMT.zip', file)
  #
  #   },
  #   contentType = "application/zip"
  # )
  #
  # output$DDA_met <- downloadHandler(
  #   filename <- function() {
  #     paste("DDA_performance_metrics", "zip", sep=".")
  #   },
  #
  #   content <- function(file) {
  #     file.copy('./R/example/metrics_DDA.zip', file)
  #
  #   },
  #   contentType = "application/zip"
  # )
  #
  # output$DIA_met <- downloadHandler(
  #   filename <- function() {
  #     paste("DIA_performance_metrics", "zip", sep=".")
  #   },
  #
  #   content <- function(file) {
  #     file.copy('./R/example/metrics_DIA.zip', file)
  #
  #   },
  #   contentType = "application/zip"
  # )
  #
  # output$TMT_met <- downloadHandler(
  #   filename <- function() {
  #     paste("TMT_performance_metrics", "zip", sep=".")
  #   },
  #
  #   content <- function(file) {
  #     file.copy('./R/example/metrics_TMT.zip', file)
  #
  #   },
  #   contentType = "application/zip"
  # )
  #
  # output$DDA_rk <- downloadHandler(
  #   filename <- function() {
  #     paste("DDA_workflow_ranks", "zip", sep=".")
  #   },
  #
  #   content <- function(file) {
  #     file.copy('./R/example/ranks_DDA.zip', file)
  #
  #   },
  #   contentType = "application/zip"
  # )
  #
  # output$DIA_rk <- downloadHandler(
  #   filename <- function() {
  #     paste("DIA_workflow_ranks", "zip", sep=".")
  #   },
  #
  #   content <- function(file) {
  #     file.copy('./R/example/ranks_DIA.zip', file)
  #
  #   },
  #   contentType = "application/zip"
  # )
  #
  # output$TMT_rk <- downloadHandler(
  #   filename <- function() {
  #     paste("TMT_workflow_ranks", "zip", sep=".")
  #   },
  #
  #   content <- function(file) {
  #     file.copy('./R/example/ranks_TMT.zip', file)
  #
  #   },
  #   contentType = "application/zip"
  # )
}
