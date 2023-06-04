tabOptFg1<-fluidRow(
  shinydashboard::box(width = 6, height = 150,
      status = "primary", solidHeader = T,
      h3("Do you have any preferred selections?"),
      radioButtons("rb1", "preferred selections:",
                   list("None selected" = "", 'Yes'='Y',
                        'No' = 'N'), inline = T),
      actionButton("suggest1", "Suggest Workflow")
  ),
  shinydashboard::box(width = 6, height = 150, status = 'primary',solidHeader = T,
      h3("Do you have the following expression matrix types?"),
      checkboxGroupInput("cgFg1", "expression matrix types:",
                         c("protein" = "pro",
                           "peptide" = "pep",
                           "spectra counts" = "sc"), inline = TRUE),
      conditionalPanel(condition = "input.cgFg1.length>1",
                       actionButton("suggest2", "Clik to apply ensemble inference")
      )
  ),
  conditionalPanel(
    condition = "input.suggest1 % 2==1 & input.rb1.includes('N')",
    shinydashboard::box(
      width = 12, height = 200,status = "success", solidHeader = T,
      h5('We suggest to apply the FragPipe-specific top-ranked workflow:'),
      pre(paste0(wf_frag$workflow[1], ':',
               '        expression matrix:', wf_frag$expression_matrix[1],
               '        normalization:', wf_frag$normalization[1],
               '        MVI:', wf_frag$imputation[1],
               '        DEA tool:', wf_frag$DEA_tool[1])),
      h5('or your can apply one of the workflows including the following choices:'),
      pre(paste0('expression matrix: protein MaxLFQ intensity',
                 '          normalization: None',
                 '          MVI: [MinProb, MinDect]',
                 '          DEA tool: limma')),
      h5("You can view the details of your selected workflow in benchmarking-FragPipe page")
    )),
  conditionalPanel(
    condition = "input.suggest1 % 2==1 & input.rb1.includes('Y')",
    shinydashboard::box(
      width = 12, height = 200,status = "success", solidHeader = T,
      h5("Per your info, our suggestions are as follows:"),
      #textOutput("fg_sug_prefer"),
      htmlOutput("fg_sug_prefer"),
      br(),
      p("You can view the details of this workflow in Benchmarking-FragPipe page")
    )
  ),
conditionalPanel(
  condition = "input.suggest2>0",
  shinydashboard::box(width = 12, height = 200,status = "success", solidHeader = T,
      strong(h5("Per your info, we suggest to use the following workflows to conduct ensemble inference:")),
      htmlOutput("fg_sug_ens")#,
      #strong(pre("Model used for p-value integration should be:      Hurdle"))
  )
)
)

tabOptMq1<-fluidRow(
  shinydashboard::box(width = 6, height = 150,
      status = "primary", solidHeader = T,
      h3("Do you have any preferred selections?"),
      radioButtons("rb2", "preferred selections:",
                   list("None selected" = "", 'Yes'='Y',
                        'No' = 'N'), inline = T),
      actionButton("suggest3", "Suggest Workflow")
  ),
  shinydashboard::box(width = 6, height = 150, status = 'primary',solidHeader = T,
      h3("Do you have the following expression matrix types?"),
      checkboxGroupInput("cgMq1", "expression matrix types:",
                         c("protein" = "pro",
                           "peptide" = "pep",
                           "spectra counts" = "sc"), inline = TRUE),
      conditionalPanel(condition = "input.cgMq1.length>1",
                       actionButton("suggest4", "Clik to apply ensemble inference")
      )
  ),
  conditionalPanel(
    condition = "input.suggest3 % 2==1 & input.rb2.includes('N')",
    shinydashboard::box(
      width = 12, height = 200,status = "success", solidHeader = T,
      h5('We suggest to apply the maxquant-specific top-ranked workflow:'),
      pre(paste0(wf_mq$workflow[1], ':',
                 '        expression matrix:', wf_frag$expression_matrix[1],
                 '        normalization:', wf_frag$normalization[1],
                 '        MVI:', wf_frag$imputation[1],
                 '        DEA tool:', wf_frag$DEA_tool[1])),
      h5('or your can apply one of the workflows including the following choices:'),
      pre(paste0('expression matrix: protein peak intensity',
                 '          normalization: [None, center.mean, center.median]',
                 '          MVI: [MinProb, MinDect, missForest]',
                 '          DEA tool: limma, DEqMS')),
      h5("You can view the details of your selected workflow in benchmarking-FragPipe page")
    )),
  conditionalPanel(
    condition = "input.suggest3 % 2==1 & input.rb2.includes('Y')",
    shinydashboard::box(
      width = 12, height = 200,status = "success", solidHeader = T,
      h5("Per your info, our suggestions are as follows:"),
      #textOutput("fg_sug_prefer"),
      htmlOutput("mq_sug_prefer"),
      br(),
      p("You can view the details of this workflow in Benchmarking-maxquant page")
    )
  ),
  conditionalPanel(
    condition = "input.suggest4>0",
    shinydashboard::box(width = 12, height = 200,status = "success", solidHeader = T,
        strong(h5("Per your info, we suggest to use the following workflows to conduct ensemble inference:")),
        htmlOutput("mq_sug_ens")#,
        #strong(pre("Model used for p-value integration should be:      Hurdle"))
    )
  )
)

tabpanel1<-tabsetPanel(
        tabPanel("FragPipe",
                 tabOptFg1,

                 conditionalPanel(
                   condition = "input.rb1.includes('Y')",
                   shinydashboard::box(width=6, title = "please check in which step:", solidHeader = TRUE,
                       status = "warning",
                       checkboxGroupInput('sug1cg1',
                                          'Steps',
                                          c("expression matrix type" = "exp",
                                            "normalization" = "norm",
                                            "MVI" = "imp",
                                            "DEA" = "dea"), inline = TRUE))
                 ),
                 conditionalPanel(
                   condition = "input.sug1cg1.includes('exp') & input.rb1.includes('Y')",
                   shinydashboard::box(width=6, title = "expression matrix types", solidHeader = TRUE,
                       status = "warning",
                       checkboxGroupInput("rb11","Data acqusition type:",
                       c("protein intensity" = "protein intensity",
                            "protein MaxLFQ" = "protein MaxLFQ",
                            "peptide intensity" = "peptide intensity",
                            "peptide MaxLFQ" = "peptide MaxLFQ",
                            "spectra counts" = "spectra counts"),inline = T))
                 ),
                 conditionalPanel(
                   condition = "input.sug1cg1.includes('norm') & input.rb1.includes('Y')",
                   shinydashboard::box(width=6, title = "Normalization methods",  solidHeader = TRUE,
                       status = "warning",
                       checkboxGroupInput("rb12","Normalization:",
                                    c("No_normalization" = "No_normalization",
                                         "center.mean" = "center.mean",
                                         "center.median" = "center.median",
                                         "max" = "max",
                                         "sum" = "sum",
                                         "vsn" = "vsn",
                                         "quantiles" = "quantiles",
                                         "quantiles.robust" = "quantiles.robust",
                                         "div.mean" = "div.mean",
                                         "div.median" = "div.median"),inline = T))
                 ),
                 conditionalPanel(
                   condition = "input.sug1cg1.includes('imp') & input.rb1.includes('Y')",
                   shinydashboard::box(width=6, title = "MVI algorithms", solidHeader = TRUE,
                       status = "warning",
                       checkboxGroupInput("rb13","MVI:",
                                    c("MinProb" = "MinProb",
                                         "QRILC" = "QRILC",
                                         "MinDet" = "MinDet",
                                         "missForest" = "missForest",
                                         "nbavg" = "nbavg",
                                         "zero" = "zero",
                                         "bpca" = "bpca",
                                         "MLE" = "MLE",
                                         "knn" = "knn",
                                         "No_MVI" = "No_MVI",
                                         "min" = "min",
                                         "mice" = "mice"),inline = T))
                 ),
                 conditionalPanel(
                   condition = "input.sug1cg1.includes('dea') & input.rb1.includes('Y')",
                   shinydashboard::box(width=6, title = "DEA tools", solidHeader = TRUE,
                       status = "warning",
                       checkboxGroupInput("rb14","DEA:",
                                    list("ANOVA" = "ANOVA",
                                         "DEP" = "DEP",
                                         "DEqMS" = "DEqMS",
                                         'limma' = "limma",
                                         "proDA" = "proDA",
                                         "ROTS" = "ROTS",
                                         "SAM" = "SAM",
                                         "siggenes" = "siggenes",
                                         "ttest" = "ttest",
                                         "beta_binomial" = "beta_binomial",
                                         "edgeR" = "edgeR",
                                         "plgem" = "plgem",
                                         "msqrob2" = "msqrob2",
                                         "ProteoMM" = "ProteoMM"),inline = T))
                 )
                 ),
        tabPanel("maxquant",
                 tabOptMq1,
                 conditionalPanel(
                   condition = "input.rb2.includes('Y')",
                   shinydashboard::box(width=6, title = "please check in which step:", solidHeader = TRUE,
                       status = "warning",
                       checkboxGroupInput('sug1cg2',
                                          'Steps',
                                          c("expression matrix type" = "exp",
                                            "normalization" = "norm",
                                            "MVI" = "imp",
                                            "DEA" = "dea"), inline = TRUE))
                 ),
                 conditionalPanel(
                   condition = "input.sug1cg2.includes('exp') & input.rb2.includes('Y')",
                   shinydashboard::box(width=6, title = "expression matrix types", solidHeader = TRUE,
                       status = "warning",
                       checkboxGroupInput("rb21","Data acqusition type:",
                                    c("protein intensity" = "protein intensity",
                                      "protein MaxLFQ" = "protein MaxLFQ",
                                      "peptide intensity" = "peptide intensity",
                                      "peptide MaxLFQ" = "peptide MaxLFQ",
                                      "spectra counts" = "spectra counts"),inline = T))
                 ),
                 conditionalPanel(
                   condition = "input.sug1cg2.includes('norm') & input.rb2.includes('Y')",
                   shinydashboard::box(width=6, title = "Normalization methods", solidHeader = TRUE,
                       status = "warning",
                       checkboxGroupInput("rb22","Normalization:",
                                    c("None" = "None",
                                      "center.mean" = "center.mean",
                                      "center.median" = "center.median",
                                      "max" = "max",
                                      "sum" = "sum",
                                      "vsn" = "vsn",
                                      "quantiles" = "quantiles",
                                      "quantiles.robust" = "quantiles.robust",
                                      "div.mean" = "div.mean",
                                      "div.median" = "div.median"),inline = T))
                 ),
                 conditionalPanel(
                   condition = "input.sug1cg2.includes('imp') & input.rb2.includes('Y')",
                   shinydashboard::box(width=6, title = "MVI algorithms", solidHeader = TRUE,
                       status = "warning",
                       checkboxGroupInput("rb23","MVI:",
                                    c("MinProb" = "MinProb",
                                      "QRILC" = "QRILC",
                                      "MinDet" = "MinDet",
                                      "missForest" = "missForest",
                                      "nbavg" = "nbavg",
                                      "zero" = "zero",
                                      "bpca" = "bpca",
                                      "MLE" = "MLE",
                                      "knn" = "knn",
                                      "None" = "None",
                                      "min" = "min",
                                      "mice" = "mice"),inline = T))
                 ),
                 conditionalPanel(
                   condition = "input.sug1cg2.includes('dea') & input.rb2.includes('Y')",
                   shinydashboard::box(width=6, title = "DEA tools", solidHeader = TRUE,
                       status = "warning",
                       checkboxGroupInput("rb24","DEA:",
                                    c("ANOVA" = "ANOVA",
                                         "DEP" = "DEP",
                                         "DEqMS" = "DEqMS",
                                         'limma' = "limma",
                                         "proDA" = "proDA",
                                         "ROTS" = "ROTS",
                                         "SAM" = "SAM",
                                         "siggenes" = "siggenes",
                                         "ttest" = "ttest",
                                         "beta_binomial" = "beta_binomial",
                                         "edgeR" = "edgeR",
                                         "plgem" = "plgem",
                                         "msqrob2" = "msqrob2",
                                         "ProteoMM" = "ProteoMM"),inline = T))
                 )
        # ),
        # tabPanel("Non-specify", checkboxGroupInput('sugcg',
        #                                         'Do you have any preferred selections?',
        #                                         c("quantification platform" = "qua",
        #                                           "expression matrix type" = "exp",
        #                                           "normalization" = "norm",
        #                                           "MVI" = "imp",
        #                                           "DEA" = "dea"), inline = TRUE))
        )
)


tabOptDia1<-fluidRow(
  shinydashboard::box(width = 6, height = 200,
      status = "primary", solidHeader = T,
      h3("Do you have any preferred selections?"),
      radioButtons("rb3", "preferred selections:",
                   list("None selected" = "", 'Yes'='Y',
                        'No' = 'N'), inline = T),
      actionButton("suggest5", "Suggest Workflow")
  ),
  shinydashboard::box(width = 6, height = 200, status = 'primary',solidHeader = T,
      h3("Do you want to try ensemble inference (this may take more time)?"),
      radioButtons("cgDia1", "whether ensemble inference:",
                   list("None selected" = "",
                        "Yes" = "Y",
                        "No" = "N"), inline = TRUE),
      conditionalPanel(condition = "input.cgDia1.includes('Y')",
                       actionButton("suggest6", "Clik to apply ensemble inference")
      )
  ),
  conditionalPanel(
    condition = "input.suggest5 % 2==1 & input.rb3.includes('N')",
    shinydashboard::box(
      width = 12, height = 200,status = "success", solidHeader = T,
      h5('We suggest to apply the DIA-NN-specific top-ranked workflow:'),
      pre(paste0(wf_dia$workflow[1], ':',
                 '        expression matrix:', wf_dia$expression_matrix[1],
                 '        normalization:', wf_dia$normalization[1],
                 '        MVI:', wf_dia$imputation[1],
                 '        DEA tool:', wf_dia$DEA_tool[1])),
      h5('or your can apply one of the workflows including the following choices:'),
      pre(paste0('expression matrix: "MaxLFQ intensity"',
                 '          normalization: [None, center.mean, center.median]',
                 '          MVI: [MinProb, MinDect]',
                 '          DEA tool: [limma, ROTS]')),
      h5("You can view the details of your selected workflow in benchmarking-FragPipe page")
    )),
  conditionalPanel(
    condition = "input.suggest5 % 2==1 & input.rb3.includes('Y')",
    shinydashboard::box(
      width = 12, height = 200,status = "success", solidHeader = T,
      h5("Per your info, our suggestions are as follows:"),
      #textOutput("fg_sug_prefer"),
      htmlOutput("dia_sug_prefer"),
      br(),
      p("You can view the details of this workflow in Benchmarking-FragPipe page")
    )
  ),
  conditionalPanel(
    condition = "input.cgDia1.includes('Y') & input.suggest6>0",
    shinydashboard::box(width = 12, height = 500,status = "success", solidHeader = T,
        strong(h5("Per your info, we suggest to use the following workflows to conduct ensemble inference:")),
        htmlOutput("dia_sug_ens")#,
        #strong(pre("Model used for p-value integration should be:      Hurdle"))
    )
  )
)

tabpanel2<-tabsetPanel(
  tabPanel("DIA-NN",
           tabOptDia1,

           conditionalPanel(
             condition = "input.rb3.includes('Y')",
             shinydashboard::box(width=6, title = "please check in which step:", solidHeader = TRUE,
                 status = "warning",
                 checkboxGroupInput('sug1cg3',
                                    'Steps',
                                    c("expression matrix type" = "exp",
                                      "normalization" = "norm",
                                      "MVI" = "imp",
                                      "DEA" = "dea"), inline = TRUE))
           ),
           conditionalPanel(
             condition = "input.sug1cg3.includes('exp') & input.rb3.includes('Y')",
             shinydashboard::box(width=6, title = "expression matrix types", solidHeader = TRUE,
                 status = "warning",
                 checkboxGroupInput("rb31","Data acqusition type:",
                                    c("MaxLFQ" = "MaxLFQ",
                                        "norm" = "norm",
                                        "MaxLFQ_raw" = "MaxLFQ_raw",
                                        "norm_raw" = "norm_raw",
                                        "raw" = "raw"),inline = T))
           ),
           conditionalPanel(
             condition = "input.sug1cg3.includes('norm') & input.rb3.includes('Y')",
             shinydashboard::box(width=6, title = "Normalization methods",  solidHeader = TRUE,
                 status = "warning",
                 checkboxGroupInput("rb32","Normalization:",
                                    c("None" = "None",
                                      "center.mean" = "center.mean",
                                      "center.median" = "center.median",
                                      "max" = "max",
                                      "sum" = "sum",
                                      "vsn" = "vsn",
                                      "quantiles" = "quantiles",
                                      "quantiles.robust" = "quantiles.robust",
                                      "div.mean" = "div.mean",
                                      "div.median" = "div.median"),inline = T))
           ),
           conditionalPanel(
             condition = "input.sug1cg3.includes('imp') & input.rb3.includes('Y')",
             shinydashboard::box(width=6, title = "MVI algorithms", solidHeader = TRUE,
                 status = "warning",
                 checkboxGroupInput("rb33","MVI:",
                                    c("MinProb" = "MinProb",
                                      "QRILC" = "QRILC",
                                      "MinDet" = "MinDet",
                                      "missForest" = "missForest",
                                      "nbavg" = "nbavg",
                                      "zero" = "zero",
                                      "bpca" = "bpca",
                                      "MLE" = "MLE",
                                      "knn" = "knn",
                                      "No_MVI" = "No_MVI",
                                      "min" = "min",
                                      "mice" = "mice"),inline = T))
           ),
           conditionalPanel(
             condition = "input.sug1cg3.includes('dea') & input.rb3.includes('Y')",
             shinydashboard::box(width=6, title = "DEA tools", solidHeader = TRUE,
                 status = "warning",
                 checkboxGroupInput("rb34","DEA:",
                                    list("ANOVA" = "ANOVA",
                                         "DEP" = "DEP",
                                         'limma' = "limma",
                                         "proDA" = "proDA",
                                         "ROTS" = "ROTS",
                                         "SAM" = "SAM",
                                         "siggenes" = "siggenes",
                                         "ttest" = "ttest"),inline = T))
           )
  # ),
  #   tabPanel("Non-specify", checkboxGroupInput('sugcg',
  #                                            'Do you have any preferred selections?',
  #                                            c("quantification platform" = "qua",
  #                                              "expression matrix type" = "exp",
  #                                              "normalization" = "norm",
  #                                              "MVI" = "imp",
  #                                              "DEA" = "dea"), inline = TRUE)
  )
  )

