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
app_server <- function(input, output, session) {
  output$menu <- renderMenu({
    sidebarMenu(id = "smenu",
                menuItem("Introduction", tabName = "front"),
                menuItem("Benchmarking",
                         helpText("ranking workflows"),
                         menuSubItem("FragPipe", "models1"),
                         menuSubItem("maxquant", "models2"),
                         menuSubItem("DIA-NN", "models3")
                ),
                menuItem("Suggestion",
                         helpText("recommend workflows"),
                         menuSubItem("DDA", "sug1"),
                         menuSubItem("DIA", "sug2")
                ),
                menuItem("Data", tabName = "accData"),
                menuItem("Help", helpText("tourial"),
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
      pt<-plot_metrics(data1()[row,][1][[1]], 'FragPipe')
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
      pt<-plot_metrics(data2()[row,][1][[1]], 'maxquant')
      pt
    }
  })

  output$tabDia <- DT::renderDataTable({
    data<- wf_dia
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
      pt<-plot_metrics(data3()[row,][1][[1]], 'DIANN')
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
    if (input$rb1=='N' | input$suggest1==0) return()
    HTML(paste0("<pre>",data_fg_sug_prefer()$workflow[1],':',
                '        expression matrix:',data_fg_sug_prefer()$expression_matrix[1],
                '        normalization:',data_fg_sug_prefer()$normalization[1],
                '        imputation:',data_fg_sug_prefer()$imputation[1],
                '        DEA tool:', data_fg_sug_prefer()$DEA_tool[1],"<pre>"))
  })

  data_fg_sug_ens <- reactive({
    get_sug_ens_fg(input)
  })

  output$fg_sug_ens <- renderUI({
    a=data_fg_sug_ens()
    if (input$suggest2==0) return()
    HTML(paste0("<pre>",data_fg_sug_ens()[1],"<pre>"))
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

  data_mq_sug_ens <- reactive({
    get_sug_ens_mq(input)
  })

  output$mq_sug_ens <- renderUI({
    #a=data_fg_sug_ens()
    if (input$suggest4==0) return()
    HTML(paste0("<pre>",data_mq_sug_ens()[1],"<pre>"))
  })

  data_dia_sug_prefer <- reactive({
    getFilsug_dia(input)
  })
  output$dia_sug_prefer <- renderUI({
    if (input$rb3=='N' | input$suggest5==0) return()
    HTML(paste0("<pre>",data_mq_sug_prefer()$workflow[1],':',
                '        expression matrix:',data_dia_sug_prefer()$expression_matrix[1],
                '        normalization:',data_dia_sug_prefer()$normalization[1],
                '        imputation:',data_dia_sug_prefer()$imputation[1],
                '        DEA tool:', data_dia_sug_prefer()$DEA_tool[1],"<pre>"))
  })

  data_dia_sug_ens <- reactive({
    get_sug_ens_dia(input)
  })

  output$dia_sug_ens <- renderUI({
    #a=data_fg_sug_ens()
    if (input$suggest6==0) return()
    HTML(paste0("<pre>",data_dia_sug_ens()[1],"<pre>"))
  })
}
