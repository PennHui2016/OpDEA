
tabOptFg1<-fluidRow(
  shinydashboard::box(width = 6, height = 200,
      status = "primary", solidHeader = T,
      h3("Do you have any preferred selections?"),
      radioButtons("rb1", "preferred selections:",
                   list("None selected" = "", 'Yes'='Y',
                        'No' = 'N'), inline = T),
      actionButton("suggest1", "Suggest Workflow"),
      splitLayout(p('click the botton again to hide !!!', style = "color:red")#,

                  )#,

      # column(3,p('click the botton again to hide !!!', style = "color:red")),
      # column(3,downloadLink("example_fg", "download_example_data_for_DEA.zip"))
  ),
  shinydashboard::box(width = 6, height = 200, status = 'primary',solidHeader = T,
      #h3("Do you have the following expression matrix types?"),
      # checkboxGroupInput("cgFg1", "expression matrix types:",
      #                    c("top0" = "top0",
      #                      "top3" = "top3",
      #                      "count" = "sc",
      #                        "MaxLFQ" = 'LFQ',
      #                      "directLFQ" = "dlfq"), inline = TRUE),
      #conditionalPanel(condition = "input.cgFg1.length>1",
      p('ensemble inference may help improve the true positive rates, Do you help try it?'),
      actionButton("suggest2", "Try ensemble inference !!!"),
      p('click the botton again to hide !!!', style = "color:red"),
      p("can obtain download_example_data_for_DEA.zip from http://www.ai4pro.tech:3838/")
      #)
  ),
  conditionalPanel(
    condition = "input.suggest1 % 2==1 & input.rb1.includes('N')",
    shinydashboard::box(
      width = 12, height = 1000,status = "success", solidHeader = T,
      h5('We suggest to apply the FragPipe-specific top-ranked workflow:'),
      pre(paste0(wf_frag$workflow[1], ':',
               '        expression matrix:', wf_frag$expression_matrix[1],
               '        normalization:', wf_frag$normalization[1],
               '        MVI:', wf_frag$imputation[1],
               '        DEA tool:', wf_frag$DEA_tool[1])),
      h5('or your can apply one of the workflows including the following choices:'),
      pre(paste0('expression matrix: [directLFQ] intensity',
                 '          normalization: None',
                 '          MVI: [SeqKNN]',
                 '          DEA tool: [limma, ROTS]')),
      h5("You can view the details of your selected workflow in benchmarking-DDA_LFQ-FragPipe page"),
      # shinyFilesButton("GetFile_fg1", "Choose a file" ,
      #                  title = "Please specify the combined_protein.tsv from FragPipe:", multiple = FALSE,
      #                  buttonType = "default", class = NULL),
      #
      # textOutput("file_fg1"),

      column(4,textInput("upload_fg_raw1", "paste path to combined_protein.tsv from FragPipe",width='100%',
                         value = 'E:/proteomics/maus2/test/OpDEA-master/R/example/FragPipe_DDA/combined_protein.tsv')),
      column(4,textInput("upload_fg_ion1", "paste path to combined_ion.tsv file from FragPipe",width='100%',
                         value = 'E:/proteomics/maus2/test/OpDEA-master/R/example/FragPipe_DDA/combined_ion.tsv')),
      column(4,textInput("upload_fg_dg1", "paste path to your design file with columns of sample name, replicate and condition",width='100%',
                         value = 'E:/proteomics/maus2/test/OpDEA-master/R/example/FragPipe_DDA/YUltq099_LFQ_FragPipe_design.tsv')),
      br(),
      p('thresholds:'),
      column(3,textInput("logFC_fg1_1", "logFC", value = '1')),
      column(3,textInput("adjp_fg1_1", "adj.p-value", value = '0.05')),
      br(),
      # column(12,textInput("python_path_fg1_1", "paste the python installation path with directLFQ python package (https://github.com/MannLabs/directlfq) installed",
      #                     value = 'D:/software/anaconda/envs/directlfq/python')),
      column(12,actionButton("runDEA_fg1_1", "DEA", class = "btn-success")),
      conditionalPanel(
        condition = "input.runDEA_fg1_1>0",
        shinydashboard::box(
          title = "Volcano plot of the differential expression analysis (DEA)",
          background = "light-blue",
          width = 12,
          height = 500,
          column(12,p("DEA with the suggested workflow, please wait for the results...")),
          column(12,span(textOutput("DEA_res_fg1_1"), style="color:white")),
          column(12,plotOutput("volc_fg1_1"))#,
          #downloadLink("DEA_res_fg1_1", "DEA results.zip")

        )
      )
    )),
  conditionalPanel(
    condition = "input.suggest1 % 2==1 & input.rb1.includes('Y')",
    shinydashboard::box(
      width = 12, height = 1000,status = "success", solidHeader = T,
      h5("Per your info, our suggestions are as follows:"),
      htmlOutput("fg_sug_prefer"),
      p("You can view the details of this workflow in benchmarking-DDA_LFQ-FragPipe page"),
      br(),
      column(4,textInput("upload_fg_raw", "paste path to combined_protein.tsv from FragPipe",width='100%',
                         value = 'E:/proteomics/maus2/test/OpDEA-master/R/example/FragPipe_DDA/combined_protein.tsv')),
      column(4,textInput("upload_fg_ion", "paste path to combined_ion.tsv file from FragPipe",width='100%',
                         value = 'E:/proteomics/maus2/test/OpDEA-master/R/example/FragPipe_DDA/combined_ion.tsv')),
      column(4,textInput("upload_fg_dg", "paste path to your design file with columns of sample name, replicate and condition",width='100%',
                         value = 'E:/proteomics/maus2/test/OpDEA-master/R/example/FragPipe_DDA/YUltq099_LFQ_FragPipe_design.tsv')),
      br(),
      p('thresholds:'),
      column(3,textInput("logFC_fg1", "logFC", value = '1')),
      column(3,textInput("adjp_fg1", "adj.p-value", value = '0.05')),
      br(),
      # column(12,textInput("python_path_fg1", "paste the python installation path with directLFQ python package (https://github.com/MannLabs/directlfq) installed",
      #                     value = 'D:/software/anaconda/envs/directlfq/python')),
      column(12,actionButton("runDEA_fg1", "DEA", class = "btn-success")),
      conditionalPanel(
        condition = "input.runDEA_fg1>0",
        shinydashboard::box(
          title = "Volcano plot of the differential expression analysis (DEA)",
          background = "light-blue",
          width = 12,
          height = 500,
          column(12,p("DEA with the suggested workflow, please wait for the results...")),
          column(12,span(textOutput("DEA_res_fg1"), style="color:white")),
          column(12,plotOutput("volc_fg1"))#,
          #downloadLink("DEA_res_fg1", "DEA results.zip")

        )
)

    )
  ),
conditionalPanel(
  condition = "input.suggest2>0 & input.suggest2 % 2 == 1",
  shinydashboard::box(width = 12, height = 1000,status = "success", solidHeader = T,
      strong(h5("We suggest to use the following workflows to conduct ensemble inference:")),
      htmlOutput("fg_sug_ens"),
      br(),
      column(4,textInput("upload_fg_raw_ens", "paste path to combined_protein.tsv from FragPipe",width='100%',
                         value = 'E:/proteomics/maus2/test/OpDEA-master/R/example/FragPipe_DDA/combined_protein.tsv')),
      column(4,textInput("upload_fg_ion_ens", "paste path to combined_ion.tsv file from FragPipe",width='100%',
                         value = 'E:/proteomics/maus2/test/OpDEA-master/R/example/FragPipe_DDA/combined_ion.tsv')),
      column(4,textInput("upload_fg_dg_ens", "paste path to your design file with columns of sample name, replicate and condition",width='100%',
                         value = 'E:/proteomics/maus2/test/OpDEA-master/R/example/FragPipe_DDA/YUltq099_LFQ_FragPipe_design.tsv')),
      br(),
      p('thresholds:'),
      column(3,textInput("logFC_fg_ens", "logFC", '1')),
      column(3,textInput("adjp_fg_ens", "adj.p-value", '0.05')),
      br(),
      # column(12,textInput("python_path_fg_ens", "paste python installation path with directLFQ python package (https://github.com/MannLabs/directlfq) installed",
      #                     value = 'D:/software/anaconda/envs/directlfq/python')),
      column(12,actionButton("runDEA_fg_ens", "DEA", class = "btn-success")),
      conditionalPanel(
        condition = "input.runDEA_fg_ens>0",
        shinydashboard::box(
          title = "Volcano plot of the differential expression analysis (DEA)",
          background = "light-blue",
          width = 12,
          height = 500,
          column(12,p("DEA with the suggested ensemble inference workflow, please wait for the results...")),
          column(12,textOutput("DEA_res_fg_ens")),
          column(12,plotOutput("volc_fg_ens"))#,
          #downloadLink("DEA_res_fg_ens", "DEA results.zip")

        )
      )
  )
)
)

tabOptMq1<-fluidRow(
  shinydashboard::box(width = 6, height = 200,
      status = "primary", solidHeader = T,
      h3("Do you have any preferred selections?"),
      radioButtons("rb2", "preferred selections:",
                   list("None selected" = "", 'Yes'='Y',
                        'No' = 'N'), inline = T),
      actionButton("suggest3", "Suggest Workflow"),
      splitLayout(p('click the botton again to hide !!!', style = "color:red")#,
      )#,

  ),
  shinydashboard::box(width = 6, height = 200, status = 'primary',solidHeader = T,

      p('ensemble inference may help improve the true positive rates, Do you help try it?'),
      actionButton("suggest4", "Try ensemble inference!!!"),
      p('click the botton again to hide !!!', style = "color:red"),
      p("Can obtain download_example_data_for_DEA.zip from http://www.ai4pro.tech:3838/")
      #)
  ),
  conditionalPanel(
    condition = "input.suggest3 % 2==1 & input.rb2.includes('N')",
    shinydashboard::box(
      width = 12, height = 1000,status = "success", solidHeader = T,
      h5('We suggest to apply the maxquant-specific top-ranked workflow:'),
      pre(paste0(wf_mq$workflow[1], ':',
                 '        expression matrix:', wf_mq$expression_matrix[1],
                 '        normalization:', wf_mq$normalization[1],
                 '        MVI:', wf_mq$imputation[1],
                 '        DEA tool:', wf_mq$DEA_tool[1])),
      h5('or your can apply one of the workflows including the following choices:'),
      pre(paste0('expression matrix: directLFQ intensity',
                 '          normalization: [None]',
                 '          MVI: [Impseq]',
                 '          DEA tool: limma')),
      h5("You can view the details of your selected workflow in benchmarking-DDA_LFQ-Maxquant page"),
      column(4,textInput("upload_mq_pro1", "paste path of proteinGroups.txt file from Maxquant",width='100%',
                         value = 'E:/proteomics/maus2/test/OpDEA-master/R/example/Maxquant_DDA/proteinGroups.txt')),
      column(4,textInput("upload_mq_evi1", "paste path of evidence.txt file from Maxquant",width='100%',
                         value = 'E:/proteomics/maus2/test/OpDEA-master/R/example/Maxquant_DDA/evidence.txt')),
      column(4,textInput("upload_mq_dg1", "paste path of design file with columns of sample name, replicate and condition",width='100%',
                         value = 'E:/proteomics/maus2/test/OpDEA-master/R/example/Maxquant_DDA/HEqe408_LFQ_Maxquant_design.tsv')),
      br(),
      p('thresholds:'),
      column(3,textInput("logFC_mq1", "logFC", value = '1')),
      column(3,textInput("adjp_mq1", "adj.p-value", value = '0.05')),
      br(),
      # column(12,textInput("python_path_mq1", "paste python installation path with directLFQ python package (https://github.com/MannLabs/directlfq) installed",
      #                     value = 'D:/software/anaconda/envs/directlfq/python')),
      column(12,actionButton("runDEA_mq1", "DEA", class = "btn-success")),
      conditionalPanel(
        condition = "input.runDEA_mq1>0",
        shinydashboard::box(
          title = "Volcano plot of the differential expression analysis (DEA)",
          background = "light-blue",
          width = 12,
          height = 500,
          column(12,p("DEA with the suggested workflow, please wait for the results...")),
          column(12,span(textOutput("DEA_res_mq1"), style="color:white")),
          column(12,plotOutput("volc_mq1"))#,
          #downloadLink("DEA_res_mq1", "DEA results.zip")
        )
      )
    )),
  conditionalPanel(
    condition = "input.suggest3 % 2==1 & input.rb2.includes('Y')",
    shinydashboard::box(
      width = 12, height = 1000,status = "success", solidHeader = T,
      h5("Per your info, our suggestions are as follows:"),
      #textOutput("fg_sug_prefer"),
      htmlOutput("mq_sug_prefer"),
      br(),
      p("You can view the details of this workflow in benchmarking-DDA_LFQ-Maxquant page"),

      column(4,textInput("upload_mq_pro", "paste path of proteinGroups.txt file from Maxquant",width='100%',
                         value = 'E:/proteomics/maus2/test/OpDEA-master/R/example/Maxquant_DDA/proteinGroups.txt')),
      column(4,textInput("upload_mq_evi", "paste path of evidence.txt file from Maxquant",width='100%',
                         value = 'E:/proteomics/maus2/test/OpDEA-master/R/example/Maxquant_DDA/evidence.txt')),
      column(4,textInput("upload_mq_dg", "paste path of design file with columns of sample name, replicate and condition",width='100%',
                         value = 'E:/proteomics/maus2/test/OpDEA-master/R/example/Maxquant_DDA/HEqe408_LFQ_Maxquant_design.tsv')),
      br(),
      p('thresholds:'),
      column(3,textInput("logFC_mq", "logFC", value = '1')),
      column(3,textInput("adjp_mq", "adj.p-value", value = '0.05')),
      br(),
      # column(12,textInput("python_path_mq", "paste python installation path with directLFQ python package (https://github.com/MannLabs/directlfq) installed",
      #                     value = 'D:/software/anaconda/envs/directlfq/python')),

      column(12,actionButton("runDEA_mq", "DEA", class = "btn-success")),
      conditionalPanel(
        condition = "input.runDEA_mq>0",
        shinydashboard::box(
          title = "Volcano plot of the differential expression analysis (DEA)",
          background = "light-blue",
          width = 12,
          height = 500,
          column(12,p("DEA with the suggested workflow, please wait for the results...")),
          column(12,textOutput("DEA_res_mq")),
          column(12,plotOutput("volc_mq"))#,
          #downloadLink("DEA_res_mq", "DEA results.zip")
        )
      )
    )
  ),
  conditionalPanel(
    condition = "input.suggest4>0 & input.suggest4 % 2 ==1",
    shinydashboard::box(width = 12, height = 1000,status = "success", solidHeader = T,
        strong(h5("We suggest to use the following workflows to conduct ensemble inference:")),
        htmlOutput("mq_sug_ens"),
        br(),
        column(4,textInput("upload_mq_pro_ens", "paste path of proteinGroups.txt file from Maxquant",width='100%',
                           value = 'E:/proteomics/maus2/test/OpDEA-master/R/example/Maxquant_DDA/proteinGroups.txt')),
        column(4,textInput("upload_mq_evi_ens", "paste path of evidence.txt file from Maxquant",width='100%',
                           value = 'E:/proteomics/maus2/test/OpDEA-master/R/example/Maxquant_DDA/evidence.txt')),
        column(4,textInput("upload_mq_dg_ens", "paste path of design file with columns of sample name, replicate and condition",width='100%',
                           value = 'E:/proteomics/maus2/test/OpDEA-master/R/example/Maxquant_DDA/HEqe408_LFQ_Maxquant_design.tsv')),
        br(),
        p('thresholds:'),
        column(3,textInput("logFC_mq_ens", "logFC", value = '1')),
        column(3,textInput("adjp_mq_ens", "adj.p-value", value = '0.05')),
        br(),
        # column(12,textInput("python_path_mq_ens", "paste python installation path with directLFQ python package (https://github.com/MannLabs/directlfq) installed",
        #                     value = 'D:/software/anaconda/envs/directlfq/python')),
        column(12,actionButton("runDEA_mq_ens", "DEA", class = "btn-success")),
        conditionalPanel(
          condition = "input.runDEA_mq_ens>0",
          shinydashboard::box(
            title = "Volcano plot of the differential expression analysis (DEA)",
            background = "light-blue",
            width = 12,
            height = 500,
            column(12,p("DEA with the suggested ensemble inference workflow, please wait for the results...")),
            column(12,textOutput("DEA_res_mq_ens")),
            column(12,plotOutput("volc_mq_ens"))#,
            #downloadLink("DEA_res_mq_ens", "DEA results.zip")
          )
        )
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
                       c("top0" = "top0",
                            "top3" = "top3",
                            "count" = "count",
                            "MaxLFQ" = "LFQ",
                            "directLFQ" = "dlfq"),inline = T))
                 ),
                 conditionalPanel(
                   condition = "input.sug1cg1.includes('norm') & input.rb1.includes('Y')",
                   shinydashboard::box(width=6, title = "Normalization methods",  solidHeader = TRUE,
                       status = "warning",
                       checkboxGroupInput("rb12","Normalization:",
                                    c("No_normalization" = "None",
                                      "center.mean" = "center.mean",
                                      "center.median" = "center.median",
                                      "max" = "max",
                                      "sum" = "sum",
                                      "vsn" = "vsn",
                                      "quantiles" = "quantiles",
                                      "quantiles.robust" = "quantiles.robust",
                                      "div.mean" = "div.mean",
                                      "div.median" = "div.median",
                                      'lossf' = 'lossf',
                                      'TIC' = 'TIC',
                                      "Rlr" = "Rlr",
                                      'MBQN' = 'MBQN'),inline = T))
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
                                      "No_MVI" = "None",
                                      "min" = "min",
                                      "mice" = "mice",
                                      "Impseq" = "Impseq",
                                      "Impseqrob" = "Impseqrob",
                                      "GMS" = "GMS",
                                      'SeqKNN' = 'SeqKNN'),inline = T))
                 ),
                 conditionalPanel(
                   condition = "input.sug1cg1.includes('dea') & input.rb1.includes('Y')",
                   shinydashboard::box(width=6, title = "DEA tools", solidHeader = TRUE,
                       status = "warning",
                       checkboxGroupInput("rb14","DEA:",
                                    c("ANOVA" = "ANOVA",
                                         "DEP" = "DEP",
                                         "DEqMS" = "DEqMS",
                                         'limma' = "limma",
                                         "proDA" = "proDA",
                                         "ROTS" = "ROTS",
                                         "SAM" = "SAM",
                                         "ttest" = "ttest",
                                         "beta_binomial" = "beta_binomial",
                                         "edgeR" = "edgeR",
                                         "plgem" = "plgem"#,
                                         #'MSstats' = 'MSstats'
                                      ),inline = T))
                 )
                 ),
        tabPanel("Maxquant",
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
                                          c("top0" = "top0",
                                            "top3" = "top3",
                                            "count" = "count",
                                            "LFQ" = "LFQ",
                                            "dlfq" = "dlfq"),inline = T))
                 ),
                 conditionalPanel(
                   condition = "input.sug1cg2.includes('norm') & input.rb2.includes('Y')",
                   shinydashboard::box(width=6, title = "Normalization methods", solidHeader = TRUE,
                       status = "warning",
                       checkboxGroupInput("rb22","Normalization:",
                                          c("No_normalization" = "None",
                                            "center.mean" = "center.mean",
                                            "center.median" = "center.median",
                                            "max" = "max",
                                            "sum" = "sum",
                                            "vsn" = "vsn",
                                            "quantiles" = "quantiles",
                                            "quantiles.robust" = "quantiles.robust",
                                            "div.mean" = "div.mean",
                                            "div.median" = "div.median",
                                            'lossf' = 'lossf',
                                            'TIC' = 'TIC',
                                            "Rlr" = "Rlr",
                                            'MBQN' = 'MBQN'),inline = T))
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
                                            "No_MVI" = "None",
                                            "min" = "min",
                                            "mice" = "mice",
                                            "Impseq" = "Impseq",
                                            "Impseqrob" = "Impseqrob",
                                            "GMS" = "GMS",
                                            'SeqKNN' = 'SeqKNN'),inline = T))
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
                                            "ttest" = "ttest",
                                            "beta_binomial" = "beta_binomial",
                                            "edgeR" = "edgeR",
                                            "plgem" = "plgem"#,
                                            #'MSstats' = 'MSstats'
                                            ),inline = T))
                 )

        )
)

tabOptDia1<-fluidRow(
  shinydashboard::box(width = 6, height = 200,
                      status = "primary", solidHeader = T,
                      h3("Do you have any preferred selections?"),
                      radioButtons("rb3", "preferred selections:",
                                   list("None selected" = "", 'Yes'='Y',
                                        'No' = 'N'), inline = T),
                      actionButton("suggest5", "Suggest Workflow"),
                      splitLayout(p('click the botton again to hide !!!', style = "color:red")#,
                      )#,

  ),
  shinydashboard::box(width = 6, height = 200, status = 'primary',solidHeader = T,
                      p('ensemble inference may help improve the true positive rates, Do you help try it?'),
                      actionButton("suggest6", "Try ensemble inference !!!"),
                      p('click the botton again to hide!!!', style = "color:red"),
                      p("Can obtain download_example_data_for_DEA.zip from http://www.ai4pro.tech:3838/")#,
  ),
  conditionalPanel(
    condition = "input.suggest5 % 2==1 & input.rb3.includes('N')",
    #suggest6=0,
    shinydashboard::box(
      width = 12, height = 1000,status = "success", solidHeader = T,
      h5('We suggest to apply the DIA-NN-specific top-ranked workflow:'),
      pre(paste0(wf_dia$workflow[1], ':',
                 '        expression matrix:', wf_dia$expression_matrix[1],
                 '        normalization:', wf_dia$normalization[1],
                 '        MVI:', wf_dia$imputation[1],
                 '        DEA tool:', wf_dia$DEA_tool[1])),
      h5('or your can apply one of the workflows including the following choices:'),
      pre(paste0('expression matrix: "directLFQ intensity"',
                 '          normalization: [None]',
                 '          MVI: [MinDect]',
                 '          DEA tool: [limma]')),
      h5("You can view the details of your selected workflow in benchmarking-DIA_LFQ-DIANN page"),

      column(4,textInput("upload_diann_report1", "paste path of report.tsv file from DIA-NN",width='100%',
                         value = 'E:/proteomics/maus2/test/OpDEA-master/R/example/DIANN_DIA/report.tsv')),
      column(4,textInput("upload_diann_dg1", "paste path of design file with columns of sample name, replicate and condition",width='100%',
                         value = 'E:/proteomics/maus2/test/OpDEA-master/R/example/DIANN_DIA/HEqe408_DIA_DIANN_design.tsv')),
      br(),
      p('thresholds:'),
      column(3,textInput("logFC_diann1", "logFC", value = '1')),
      column(3,textInput("adjp_diann1", "adj.p-value", value = '0.05')),
      br(),
      # column(12,textInput("python_path_diann1", "paste python installation path with directLFQ python package (https://github.com/MannLabs/directlfq) installed",
      #                     value = 'D:/software/anaconda/envs/directlfq/python')),

      column(12,actionButton("runDEA_diann1", "DEA", class = "btn-success")),
      conditionalPanel(
        condition = "input.runDEA_diann1>0",
        shinydashboard::box(
          title = "Volcano plot of the differential expression analysis (DEA)",
          background = "light-blue",
          width = 12,
          height = 500,
          column(12,p("DEA with the suggested workflow, please wait for the results...")),
          column(12,textOutput("DEA_res_diann1")),
          column(12,plotOutput("volc_diann1"))#,
          #downloadLink("DEA_res_diann1", "DEA results.zip")
        )
      )
    )),
  conditionalPanel(
    condition = "input.suggest5 % 2==1 & input.rb3.includes('Y')",
    shinydashboard::box(
      width = 12, height = 1000,status = "success", solidHeader = T,
      h5("Per your info, our suggestions are as follows:"),
      #textOutput("fg_sug_prefer"),
      htmlOutput("diann_sug_prefer"),
      br(),
      p("You can view the details of this workflow in benchmarking-DIA_LFQ-DIANN page"),

      column(4,textInput("upload_diann_report", "paste path of report.tsv file from DIA-NN",width='100%',
                         value = 'E:/proteomics/maus2/test/OpDEA-master/R/example/DIANN_DIA/report.tsv')),
      column(4,textInput("upload_diann_dg", "paste path of design file with columns of sample name, replicate and condition",width='100%',
                         value = 'E:/proteomics/maus2/test/OpDEA-master/R/example/DIANN_DIA/HEqe408_DIA_DIANN_design.tsv')),
      br(),
      p('thresholds:'),
      column(3,textInput("logFC_diann", "logFC", value = '1')),
      column(3,textInput("adjp_diann", "adj.p-value", value = '0.05')),
      br(),
      # column(12,textInput("python_path_diann", "paste python installation path with directLFQ python package (https://github.com/MannLabs/directlfq) installed",
      #                     value = 'D:/software/anaconda/envs/directlfq/python')),

      column(12,actionButton("runDEA_diann", "DEA", class = "btn-success")),
      conditionalPanel(
        condition = "input.runDEA_diann>0",
        shinydashboard::box(
          title = "Volcano plot of the differential expression analysis (DEA)",
          background = "light-blue",
          width = 12,
          height = 500,
          p("DEA with the suggested workflow, please wait for the results..."),
          column(12,textOutput("DEA_res_diann")),
          column(12,plotOutput("volc_diann"))#,
          #downloadLink("DEA_res_diann", "DEA results.zip")
        )
      )
    )
  ),
  conditionalPanel(
    condition = "input.suggest6>0 & input.suggest6 % 2==1",
    shinydashboard::box(width = 12, height = 500,status = "success", solidHeader = T,
                        strong(h5("We suggest to use the following workflows to conduct ensemble inference:")),
                        htmlOutput("diann_sug_ens"),
                        br(),

                        column(4,textInput("upload_diann_ens", "paste path of report.tsv file from DIA-NN",width='100%',
                                           value = 'E:/proteomics/maus2/test/OpDEA-master/R/example/DIANN_DIA/report.tsv')),
                        column(4,textInput("upload_diann_dg_ens", "paste path of design file with columns of sample name, replicate and condition",width='100%',
                                           value = 'E:/proteomics/maus2/test/OpDEA-master/R/example/DIANN_DIA/HEqe408_DIA_DIANN_design.tsv')),
                        br(),
                        p('thresholds:'),
                        column(3,textInput("logFC_diann_ens", "logFC", value = '1')),
                        column(3,textInput("adjp_diann_ens", "adj.p-value", value = '0.05')),
                        br(),
                        # column(12,textInput("python_path_diann_ens", "paste python installation path with directLFQ python package (https://github.com/MannLabs/directlfq) installed",
                        #                     value = 'D:/software/anaconda/envs/directlfq/python')),

                        column(12,actionButton("runDEA_diann_ens", "DEA", class = "btn-success")),
                        conditionalPanel(
                          condition = "input.runDEA_diann_ens>0",
                          shinydashboard::box(
                            title = "Volcano plot of the differential expression analysis (DEA)",
                            background = "light-blue",
                            width = 12,
                            height = 500,
                            column(12,p("DEA with the suggested ensemble inference workflow, please wait for the results...")),
                            column(12,textOutput("DEA_res_diann_ens")),
                            column(12,plotOutput("volc_diann_ens"))#,
                            #downloadLink("DEA_res_diann_ens", "DEA results.zip")
                          )
                        )
    )
  )
)

tabOptDia2<-fluidRow(
  shinydashboard::box(width = 6, height = 200,
                      status = "primary", solidHeader = T,
                      h3("Do you have any preferred selections?"),
                      radioButtons("rb4", "preferred selections:",
                                   list("None selected" = "", 'Yes'='Y',
                                        'No' = 'N'), inline = T),
                      actionButton("suggest7", "Suggest Workflow"),
                      splitLayout(p('click the botton again to hide!!!', style = "color:red")#,
                      )#,

  ),
  shinydashboard::box(width = 6, height = 200, status = 'primary',solidHeader = T,
                      p('ensemble inference may help improve the true positive rates, Do you help try it?'),
                      actionButton("suggest8", "Try ensemble inference !!!"),
                      p('click the botton again to hide!!!', style = "color:red"),
                      p("Can obtain download_example_data_for_DEA.zip from http://www.ai4pro.tech:3838/")

  ),
  conditionalPanel(
    condition = "input.suggest7 % 2==1 & input.rb4.includes('N')",
    shinydashboard::box(
      width = 12, height = 1000,status = "success", solidHeader = T,
      h5('We suggest to apply the Spectronaut-specific top-ranked workflow:'),
      pre(paste0(wf_dia$workflow[1], ':',
                 '        expression matrix:', wf_dia$expression_matrix[1],
                 '        normalization:', wf_dia$normalization[1],
                 '        MVI:', wf_dia$imputation[1],
                 '        DEA tool:', wf_dia$DEA_tool[1])),
      h5('or your can apply one of the workflows including the following choices:'),
      pre(paste0('expression matrix: "directLFQ intensity"',
                 '          normalization: [None]',
                 '          MVI: [Impseq]',
                 '          DEA tool: [ROTS, limma]')),
      h5("You can view the details of your selected workflow in benchmarking-DIA_LFQ-Spectronaut page"),

      column(4,textInput("upload_spt_report1", "paste path of report.tsv file from Spectronaut",width='100%',
                         value = 'E:/proteomics/maus2/test/OpDEA-master/R/example/Spectronaut_DIA/20231013_224410_EA100915_Report.tsv')),
      column(4,textInput("upload_spt_dg1", "paste path of design file with columns of sample name, replicate and condition",width='100%',
                         value = 'E:/proteomics/maus2/test/OpDEA-master/R/example/Spectronaut_DIA/HEqe408_DIA_spt_design.tsv')),
      br(),
      p('thresholds:'),
      column(3,textInput("logFC_spt1", "logFC", value = '1')),
      column(3,textInput("adjp_spt1", "adj.p-value", value = '0.05')),
      br(),
      # column(12,textInput("python_path_spt1", "paste python installation path with directLFQ python package (https://github.com/MannLabs/directlfq) installed",
      #                     value = 'D:/software/anaconda/envs/directlfq/python')),

      column(12,actionButton("runDEA_spt1", "DEA", class = "btn-success")),
      conditionalPanel(
        condition = "input.runDEA_spt1>0",
        shinydashboard::box(
          title = "Volcano plot of the differential expression analysis (DEA)",
          background = "light-blue",
          width = 12,
          height = 500,
          column(12,p("DEA with the suggested workflow, please wait for the results...")),
          column(12,textOutput("DEA_res_spt1")),
          column(12,plotOutput("volc_spt1"))#,
          #downloadLink("DEA_res_spt1", "DEA results.zip")
        )
      )
    )),
  conditionalPanel(
    condition = "input.suggest7 % 2==1 & input.rb4.includes('Y')",
    shinydashboard::box(
      width = 12, height = 1000,status = "success", solidHeader = T,
      h5("Per your info, our suggestions are as follows:"),
      #textOutput("fg_sug_prefer"),
      htmlOutput("spt_sug_prefer"),
      br(),
      p("You can view the details of this workflow in benchmarking-DIA_LFQ-Spectronaut page"),

      column(4,textInput("upload_spt_report", "paste path of report.tsv file from Spectronaut",width='100%',
                         value = 'E:/proteomics/maus2/test/OpDEA-master/R/example/Spectronaut_DIA/20231013_224410_EA100915_Report.tsv')),
      column(4,textInput("upload_spt_dg", "paste path of design file with columns of sample name, replicate and condition",width='100%',
                         value = 'E:/proteomics/maus2/test/OpDEA-master/R/example/Spectronaut_DIA/HEqe408_DIA_spt_design.tsv')),
      br(),
      p('thresholds:'),
      column(3,textInput("logFC_spt", "logFC", value = '1')),
      column(3,textInput("adjp_spt", "adj.p-value", value = '0.05')),
      br(),
      # column(12,textInput("python_path_spt", "paste python installation path with directLFQ python package (https://github.com/MannLabs/directlfq) installed",
      #                     value = 'D:/software/anaconda/envs/directlfq/python')),

      column(12,actionButton("runDEA_spt", "DEA", class = "btn-success")),
      conditionalPanel(
        condition = "input.runDEA_spt>0",
        shinydashboard::box(
          title = "Volcano plot of the differential expression analysis (DEA)",
          background = "light-blue",
          width = 12,
          height = 500,
          column(12,p("DEA with the suggested workflow, please wait for the results...")),
          column(12,textOutput("DEA_res_spt")),
          column(12,plotOutput("volc_spt"))#,
          #downloadLink("DEA_res_spt", "DEA results.zip")
        )
      )
    )
  ),
  conditionalPanel(
    condition = "input.suggest8 % 2 == 1 & input.suggest8>0",
    shinydashboard::box(width = 12, height = 1000,status = "success", solidHeader = T,
                        strong(h5("We suggest to use the following workflows to conduct ensemble inference:")),
                        htmlOutput("spt_sug_ens"),
                        br(),

                        column(4,textInput("upload_spt_ens", "paste path of report.tsv file from Spectronaut",width='100%',
                                           value = 'E:/proteomics/maus2/test/OpDEA-master/R/example/Spectronaut_DIA/20231013_224410_EA100915_Report.tsv')),
                        column(4,textInput("upload_spt_dg_ens", "paste path of design file with columns of sample name, replicate and condition",width='100%',
                                           value = 'E:/proteomics/maus2/test/OpDEA-master/R/example/Spectronaut_DIA/HEqe408_DIA_spt_design.tsv')),
                        br(),
                        p('thresholds:'),
                        column(3,textInput("logFC_spt_ens", "logFC", value = '1')),
                        column(3,textInput("adjp_spt_ens", "adj.p-value", value = '0.05')),
                        br(),
                        # column(12,textInput("python_path_spt", "paste python installation path with directLFQ python package (https://github.com/MannLabs/directlfq) installed",
                        #                     value = 'D:/software/anaconda/envs/directlfq/python')),

                        column(12,actionButton("runDEA_spt_ens", "DEA", class = "btn-success")),
                        conditionalPanel(
                          condition = "input.runDEA_spt_ens>0",
                          shinydashboard::box(
                            title = "Volcano plot of the differential expression analysis (DEA)",
                            background = "light-blue",
                            width = 12,
                            height = 500,
                            column(12,p("DEA with the suggested ensemble inference workflow, please wait for the results...")),
                            column(12,textOutput("DEA_res_spt_ens")),
                            column(12,plotOutput("volc_spt_ens"))#,
                            #downloadLink("DEA_res_spt_ens", "DEA results.zip")
                          )
                        )
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
                                    c("top1" = "top1",
                                      "top3" = "top3",
                                      "MaxLFQ" = 'LFQ',
                                      "directLFQ" = "dlfq"),inline = T))
           ),
           conditionalPanel(
             condition = "input.sug1cg3.includes('norm') & input.rb3.includes('Y')",
             shinydashboard::box(width=6, title = "Normalization methods",  solidHeader = TRUE,
                 status = "warning",
                 checkboxGroupInput("rb32","Normalization:",
                                    c("No_normalization" = "None",
                                      "center.mean" = "center.mean",
                                      "center.median" = "center.median",
                                      "max" = "max",
                                      "sum" = "sum",
                                      "vsn" = "vsn",
                                      "quantiles" = "quantiles",
                                      "quantiles.robust" = "quantiles.robust",
                                      "div.mean" = "div.mean",
                                      "div.median" = "div.median",
                                      'lossf' = 'lossf',
                                      'TIC' = 'TIC',
                                      "Rlr" = "Rlr",
                                      'MBQN' = 'MBQN'
                                    ),inline = T))
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
                                      "No_MVI" = "None",
                                      "min" = "min",
                                      "mice" = "mice",
                                      "Impseq" = "Impseq",
                                      "Impseqrob" = "Impseqrob",
                                      "GMS" = "GMS",
                                      'SeqKNN' = 'SeqKNN'),inline = T))
           ),
           conditionalPanel(
             condition = "input.sug1cg3.includes('dea') & input.rb3.includes('Y')",
             shinydashboard::box(width=6, title = "DEA tools", solidHeader = TRUE,
                 status = "warning",
                 checkboxGroupInput("rb34","DEA:",
                                    c("ANOVA" = "ANOVA",
                                      "DEP" = "DEP",
                                      'limma' = "limma",
                                      "proDA" = "proDA",
                                      "ROTS" = "ROTS",
                                      "SAM" = "SAM",
                                      "ttest" = "ttest"#,
                                      #'MSstats' = 'MSstats'
                                    ),inline = T))
           )
  ),
  tabPanel("Spectronaut",
           tabOptDia2,

           conditionalPanel(
             condition = "input.rb4.includes('Y')",
             shinydashboard::box(width=6, title = "please check in which step:", solidHeader = TRUE,
                                 status = "warning",
                                 checkboxGroupInput('sug1cg4',
                                                    'Steps',
                                                    c("expression matrix type" = "exp",
                                                      "normalization" = "norm",
                                                      "MVI" = "imp",
                                                      "DEA" = "dea"), inline = TRUE))
           ),
           conditionalPanel(
             condition = "input.sug1cg4.includes('exp') & input.rb4.includes('Y')",
             shinydashboard::box(width=6, title = "expression matrix types", solidHeader = TRUE,
                                 status = "warning",
                                 checkboxGroupInput("rb41","Data acqusition type:",
                                                    c("top1" = "top1",
                                      "top3" = "top3",
                                      "MaxLFQ" = 'LFQ',
                                      "directLFQ" = "dlfq"),inline = T))
           ),
           conditionalPanel(
             condition = "input.sug1cg4.includes('norm') & input.rb4.includes('Y')",
             shinydashboard::box(width=6, title = "Normalization methods",  solidHeader = TRUE,
                                 status = "warning",
                                 checkboxGroupInput("rb42","Normalization:",
                                                    c("No_normalization" = "None",
                                      "center.mean" = "center.mean",
                                      "center.median" = "center.median",
                                      "max" = "max",
                                      "sum" = "sum",
                                      "vsn" = "vsn",
                                      "quantiles" = "quantiles",
                                      "quantiles.robust" = "quantiles.robust",
                                      "div.mean" = "div.mean",
                                      "div.median" = "div.median",
                                      'lossf' = 'lossf',
                                      'TIC' = 'TIC',
                                      "Rlr" = "Rlr",
                                      'MBQN' = 'MBQN'
                                    ),inline = T))
           ),
           conditionalPanel(
             condition = "input.sug1cg4.includes('imp') & input.rb4.includes('Y')",
             shinydashboard::box(width=6, title = "MVI algorithms", solidHeader = TRUE,
                                 status = "warning",
                                 checkboxGroupInput("rb43","MVI:",
                                                    c("MinProb" = "MinProb",
                                      "QRILC" = "QRILC",
                                      "MinDet" = "MinDet",
                                      "missForest" = "missForest",
                                      "nbavg" = "nbavg",
                                      "zero" = "zero",
                                      "bpca" = "bpca",
                                      "MLE" = "MLE",
                                      "knn" = "knn",
                                      "No_MVI" = "None",
                                      "min" = "min",
                                      "mice" = "mice",
                                      "Impseq" = "Impseq",
                                      "Impseqrob" = "Impseqrob",
                                      "GMS" = "GMS",
                                      'SeqKNN' = 'SeqKNN'),inline = T))
           ),
           conditionalPanel(
             condition = "input.sug1cg4.includes('dea') & input.rb4.includes('Y')",
             shinydashboard::box(width=6, title = "DEA tools", solidHeader = TRUE,
                                 status = "warning",
                                 checkboxGroupInput("rb44","DEA:",
                                                    c("ANOVA" = "ANOVA",
                                      "DEP" = "DEP",
                                      'limma' = "limma",
                                      "proDA" = "proDA",
                                      "ROTS" = "ROTS",
                                      "SAM" = "SAM",
                                      "ttest" = "ttest"#,
                                      #'MSstats' = 'MSstats'
                                    ),inline = T))
           )

  )
  )

tabOptFg2<-fluidRow(
  shinydashboard::box(width = 12, height = 200,
                      status = "primary", solidHeader = T,
                      h3("Do you have any preferred selections?"),
                      radioButtons("rb5", "preferred selections:",
                                   list("None selected" = "", 'Yes'='Y',
                                        'No' = 'N'), inline = T),
                      actionButton("suggest9", "Suggest Workflow"),
                      splitLayout(p('click the botton again to hide !!!', style = "color:red"),
                                  p("Can obtain download_example_data_for_DEA.zip from http://www.ai4pro.tech:3838/")
                      ),
                      #,
  ),
  conditionalPanel(
    condition = "input.suggest9 % 2==1 & input.rb5.includes('N')",
    shinydashboard::box(
      width = 12, height = 1000,status = "success", solidHeader = T,
      h5('We suggest to apply the TMT-FragPipe-specific top-ranked workflow:'),
      pre(paste0(wf_fg_tmt$workflow[1], ':',
                 '        expression matrix:', wf_fg_tmt$expression_matrix[1],
                 '        normalization:', wf_fg_tmt$normalization[1],
                 '        MVI:', wf_fg_tmt$imputation[1],
                 '        DEA tool:', wf_fg_tmt$DEA_tool[1])),
      h5('or your can apply one of the workflows including the following choices:'),
      pre(paste0('expression matrix: TMT-Integrator abundance',
                 '          normalization: None',
                 '          MVI: [SeqKNN]',
                 '          DEA tool: limma')),
      h5("You can view the details of your selected workflow in benchmarking-TMT-FragPipe page"),
      br(),
      p("You can view the details of this workflow in benchmarking-DIA_LFQ-Spectronaut page"),
      column(3,textInput("upload_abd1", "paste path of tmt-report/abundance_protein_MD.tsv file from FragPipe",width='100%',
                         value = 'E:/proteomics/maus2/test/OpDEA-master/R/example/FragPipe_TMT/tmt-report/abundance_protein_MD.tsv')),
      column(3,textInput("upload_ratio1", "paste path of tmt-report/ratio_protein_MD file from FragPipe",width='100%',
                         value = 'E:/proteomics/maus2/test/OpDEA-master/R/example/FragPipe_TMT/tmt-report/ratio_protein_MD.tsv')),
      column(3,textInput("upload_phi1", "paste path of protein.tsv file from FragPipe Philosopher",width='100%',
                         value = 'E:/proteomics/maus2/test/OpDEA-master/R/example/FragPipe_TMT/protein.tsv')),
      column(3,textInput("upload_tmt_fg_dg1", "paste path of design file with columns of sample name, replicate and condition",width='100%',
                         value = 'E:/proteomics/maus2/test/OpDEA-master/R/example/FragPipe_TMT/HEqe277_TMT10_FragPipe_tmt_design.tsv')),
      br(),
      p('thresholds:'),
      column(3,textInput("logFC_tmt_fg1", "logFC", value = '1')),
      column(3,textInput("adjp_tmt_fg1", "adj.p-value", value = '0.05')),
      # br(),
      # column(12,textInput("python_path_tmt_fg1", "paste python installation path with directLFQ python package (https://github.com/MannLabs/directlfq) installed",
      #                     value = 'D:/software/anaconda/envs/directlfq/python')),

      column(12,actionButton("runDEA_tmt_fg1", "DEA", class = "btn-success")),
      conditionalPanel(
        condition = "input.runDEA_tmt_fg1>0",
        shinydashboard::box(
          title = "Volcano plot of the differential expression analysis (DEA)",
          background = "light-blue",
          width = 12,
          height = 500,
          column(12,p("DEA with the suggested workflow, please wait for the results...")),
          column(12,textOutput("DEA_res_tmt_fg1")),
          column(12,plotOutput("volc_tmt_fg1"))#,
          #downloadLink("DEA_res_tmt_fg1", "DEA results.zip")
        )
      )
    )),
  conditionalPanel(
    condition = "input.suggest9 % 2==1 & input.rb5.includes('Y')",
    shinydashboard::box(
      width = 12, height = 1000,status = "success", solidHeader = T,
      h5("Per your info, our suggestions are as follows:"),
      #textOutput("fg_sug_prefer"),
      htmlOutput("fg_tmt_sug_prefer"),
      # pre(paste0('expression matrix: protein MaxLFQ intensity',
      #            '          normalization: None',
      #            '          MVI: [MinProb, MinDect]',
      #            '          DEA tool: limma')),
      h5("You can view the details of your selected workflow in benchmarking-TMT-FragPipe page"),
      br(),
      p("You can view the details of this workflow in benchmarking-DIA_LFQ-Spectronaut page"),

      column(3,textInput("upload_abd", "paste path of tmt-report/abundance_protein_MD.tsv file from FragPipe",width='100%',
                         value = 'E:/proteomics/maus2/test/OpDEA-master/R/example/FragPipe_TMT/tmt-report/abundance_protein_MD.tsv')),
      column(3,textInput("upload_ratio", "paste path of tmt-report/ratio_protein_MD file from FragPipe",width='100%',
                         value = 'E:/proteomics/maus2/test/OpDEA-master/R/example/FragPipe_TMT/tmt-report/ratio_protein_MD.tsv')),
      column(3,textInput("upload_phi", "paste path of protein.tsv file from FragPipe Philosopher",width='100%',
                         value = 'E:/proteomics/maus2/test/OpDEA-master/R/example/FragPipe_TMT/protein.tsv')),
      column(3,textInput("upload_tmt_fg_dg", "paste path of design file with columns of sample name, replicate and condition",width='100%',
                         value = 'E:/proteomics/maus2/test/OpDEA-master/R/example/FragPipe_TMT/HEqe277_TMT10_FragPipe_tmt_design.tsv')),
      br(),
      p('thresholds:'),
      column(3,textInput("logFC_tmt_fg", "logFC", value = '1')),
      column(3,textInput("adjp_tmt_fg", "adj.p-value", value = '0.05')),
      # br(),
      # column(12,textInput("python_path_tmt_fg", "paste python installation path with directLFQ python package (https://github.com/MannLabs/directlfq) installed",
      #                     value = 'D:/software/anaconda/envs/directlfq/python')),

      column(12,actionButton("runDEA_tmt_fg", "DEA", class = "btn-success")),
      conditionalPanel(
        condition = "input.runDEA_tmt_fg>0",
        shinydashboard::box(
          title = "Volcano plot of the differential expression analysis (DEA)",
          background = "light-blue",
          width = 12,
          height = 500,
          column(12,p("DEA with the suggested workflow, please wait for the results...")),
          column(12,textOutput("DEA_res_tmt_fg")),
          column(12,plotOutput("volc_tmt_fg"))#,
          #downloadLink("DEA_res_tmt_fg", "DEA results.zip")
        )
      )
    )
  )
)

tabOptMq2<-fluidRow(
  shinydashboard::box(width = 12, height = 200,
                      status = "primary", solidHeader = T,
                      h3("Do you have any preferred selections?"),
                      radioButtons("rb6", "preferred selections:",
                                   list("None selected" = "", 'Yes'='Y',
                                        'No' = 'N'), inline = T),
                      actionButton("suggest11", "Suggest Workflow"),
                      splitLayout(p('click the botton again to hide !!!', style = "color:red"),
                                  p("Can obtain download_example_data_for_DEA.zip from http://www.ai4pro.tech:3838/")#,
                                  )#,

  ),
  conditionalPanel(
    condition = "input.suggest11 % 2==1 & input.rb6.includes('N')",
    shinydashboard::box(
      width = 12, height = 1000,status = "success", solidHeader = T,
      h5('We suggest to apply the TMT-Maxquant-specific top-ranked workflow:'),
      pre(paste0(wf_mq_tmt$workflow[1], ':',
                 '        expression matrix:', wf_mq_tmt$expression_matrix[1],
                 '        normalization:', wf_mq_tmt$normalization[1],
                 '        MVI:', wf_mq_tmt$imputation[1],
                 '        DEA tool:', wf_mq_tmt$DEA_tool[1])),
      h5('or your can apply one of the workflows including the following choices:'),
      pre(paste0('expression matrix: Reporter intensity',
                 '          normalization: None',
                 '          MVI: [bpca]',
                 '          DEA tool: [limma, ROTS]')),
      h5("You can view the details of your selected workflow in benchmarking-TMT-FragPipe page"),
      br(),
      p("You can view the details of this workflow in benchmarking-DIA_LFQ-Spectronaut page"),

      column(4,textInput("upload_tmt_mq_pro1", "paste path of proteinGroups.txt file from Maxquant",width='100%',
                         value = 'E:/proteomics/maus2/test/OpDEA-master/R/example/Maxquant_TMT/proteinGroups.txt')),
      column(4,textInput("upload_tmt_mq_evi1", "paste path of evidence.txt file from Maxquant",width='100%',
                         value = 'E:/proteomics/maus2/test/OpDEA-master/R/example/Maxquant_TMT/evidence.txt')),
      column(4,textInput("upload_tmt_mq_dg1", "paste path of design file with columns of sample name, replicate and condition",width='100%',
                         value = 'E:/proteomics/maus2/test/OpDEA-master/R/example/Maxquant_TMT/HEqe277_TMT10_Maxquant_tmt_design.tsv')),
      br(),
      p('thresholds:'),
      column(3,textInput("logFC_tmt_mq1", "logFC", value = '1')),
      column(3,textInput("adjp_tmt_mq1", "adj.p-value", value = '0.05')),
      # br(),
      # column(12,textInput("python_path_tmt_mq1", "paste python installation path with directLFQ python package (https://github.com/MannLabs/directlfq) installed",
      #                     value = 'D:/software/anaconda/envs/directlfq/python')),

      column(12,actionButton("runDEA_tmt_mq1", "DEA", class = "btn-success")),
      conditionalPanel(
        condition = "input.runDEA_tmt_mq1>0",
        shinydashboard::box(
          title = "Volcano plot of the differential expression analysis (DEA)",
          background = "light-blue",
          width = 12,
          height = 500,
          column(12,p("DEA with the suggested workflow, please wait for the results...")),
          column(12,textOutput("DEA_res_tmt_mq1")),
          column(12,plotOutput("volc_tmt_mq1"))#,
          #downloadLink("DEA_res_tmt_mq1", "DEA results.zip")
        )
      )
    )),
  conditionalPanel(
    condition = "input.suggest11 % 2==1 & input.rb6.includes('Y')",
    shinydashboard::box(
      width = 12, height = 1000,status = "success", solidHeader = T,
      h5("Per your info, our suggestions are as follows:"),
      #textOutput("fg_sug_prefer"),
      htmlOutput("mq_tmt_sug_prefer"),
      br(),
      p("You can view the details of this workflow in benchmarking-TMT-FragPipe page"),
      br(),
      p("You can view the details of this workflow in benchmarking-DIA_LFQ-Spectronaut page"),

      column(4,textInput("upload_tmt_mq_pro", "paste path of proteinGroups.txt file from Maxquant",width='100%',
                         value = 'E:/proteomics/maus2/test/OpDEA-master/R/example/Maxquant_TMT/proteinGroups.txt')),
      column(4,textInput("upload_tmt_mq_evi", "paste path of evidence.txt file from Maxquant",width='100%',
                         value = 'E:/proteomics/maus2/test/OpDEA-master/R/example/Maxquant_TMT/evidence.txt')),
      column(4,textInput("upload_tmt_mq_dg", "paste path of design file with columns of sample name, replicate and condition",width='100%',
                         value = 'E:/proteomics/maus2/test/OpDEA-master/R/example/Maxquant_TMT/HEqe277_TMT10_Maxquant_tmt_design.tsv')),
      br(),
      p('thresholds:'),
      column(3,textInput("logFC_tmt_mq", "logFC", value = '1')),
      column(3,textInput("adjp_tmt_mq", "adj.p-value", value = '0.05')),
      # br(),
      # column(12,textInput("python_path_tmt_mq", "paste python installation path with directLFQ python package (https://github.com/MannLabs/directlfq) installed",
      #                     value = 'D:/software/anaconda/envs/directlfq/python')),

      column(12,actionButton("runDEA_tmt_mq", "DEA", class = "btn-success")),
      conditionalPanel(
        condition = "input.runDEA_tmt_mq>0",
        shinydashboard::box(
          title = "Volcano plot of the differential expression analysis (DEA)",
          background = "light-blue",
          width = 12,
          height = 500,
          column(12,p("DEA with the suggested workflow, please wait for the results...")),
          column(12,textOutput("DEA_res_tmt_mq")),
          column(12,plotOutput("volc_tmt_mq"))#,
          #downloadLink("DEA_res_tmt_mq", "DEA results.zip")
        )
      )
    )
  )#,
)



tabpanel3<-tabsetPanel(
  tabPanel("FragPipe",
           tabOptFg2,

           conditionalPanel(
             condition = "input.rb5.includes('Y')",
             shinydashboard::box(width=6, title = "please check in which step:", solidHeader = TRUE,
                                 status = "warning",
                                 checkboxGroupInput('sug1cg5',
                                                    'Steps',
                                                    c("expression matrix type" = "exp",
                                                      "normalization" = "norm",
                                                      "MVI" = "imp",
                                                      "DEA" = "dea"), inline = TRUE))
           ),
           conditionalPanel(
             condition = "input.sug1cg5.includes('exp') & input.rb5.includes('Y')",
             shinydashboard::box(width=6, title = "expression matrix types", solidHeader = TRUE,
                                 status = "warning",
                                 checkboxGroupInput("rb51","Data acqusition type:",
                                                    c("abd" = "abd",
                                                      "ratio" = "ratio",
                                                      "phi" = "phi"
                                                      #"norm_raw" = "norm_raw",
                                                    ),inline = T))
           ),
           conditionalPanel(
             condition = "input.sug1cg5.includes('norm') & input.rb5.includes('Y')",
             shinydashboard::box(width=6, title = "Normalization methods",  solidHeader = TRUE,
                                 status = "warning",
                                 checkboxGroupInput("rb52","Normalization:",
                                                    c("No_normalization" = "None",
                                                      "center.mean" = "center.mean",
                                                      "center.median" = "center.median",
                                                      "max" = "max",
                                                      "sum" = "sum",
                                                      "vsn" = "vsn",
                                                      "quantiles" = "quantiles",
                                                      "quantiles.robust" = "quantiles.robust",
                                                      "div.mean" = "div.mean",
                                                      "div.median" = "div.median",
                                                      'lossf' = 'lossf',
                                                      'TIC' = 'TIC',
                                                      "Rlr" = "Rlr",
                                                      'MBQN' = 'MBQN'
                                                    ),inline = T))
           ),
           conditionalPanel(
             condition = "input.sug1cg5.includes('imp') & input.rb5.includes('Y')",
             shinydashboard::box(width=6, title = "MVI algorithms", solidHeader = TRUE,
                                 status = "warning",
                                 checkboxGroupInput("rb53","MVI:",
                                                    c("MinProb" = "MinProb",
                                                      "QRILC" = "QRILC",
                                                      "MinDet" = "MinDet",
                                                      "missForest" = "missForest",
                                                      "nbavg" = "nbavg",
                                                      "zero" = "zero",
                                                      "bpca" = "bpca",
                                                      "MLE" = "MLE",
                                                      "knn" = "knn",
                                                      "No_MVI" = "None",
                                                      "min" = "min",
                                                      "mice" = "mice",
                                                      "Impseq" = "Impseq",
                                                      "Impseqrob" = "Impseqrob",
                                                      "GMS" = "GMS",
                                                      'SeqKNN' = 'SeqKNN'),inline = T))
           ),
           conditionalPanel(
             condition = "input.sug1cg5.includes('dea') & input.rb5.includes('Y')",
             shinydashboard::box(width=6, title = "DEA tools", solidHeader = TRUE,
                                 status = "warning",
                                 checkboxGroupInput("rb54","DEA:",
                                                    c("ANOVA" = "ANOVA",
                                                      "DEP" = "DEP",
                                                      #"DEqMS" = "DEqMS",
                                                      'limma' = "limma",
                                                      "proDA" = "proDA",
                                                      "ROTS" = "ROTS",
                                                      "SAM" = "SAM",
                                                      "ttest" = "ttest"#,
                                                      #"beta_binomial" = "beta_binomial",
                                                      #"edgeR" = "edgeR",
                                                      #"plgem" = "plgem",
                                                      #'MSstats' = 'MSstats'
                                                    ),inline = T))
           )
  ),
  tabPanel("Maxquant",
           tabOptMq2,
           conditionalPanel(
             condition = "input.rb6.includes('Y')",
             shinydashboard::box(width=6, title = "please check in which step:", solidHeader = TRUE,
                                 status = "warning",
                                 checkboxGroupInput('sug1cg6',
                                                    'Steps',
                                                    c(#"expression matrix type" = "exp",
                                                      "normalization" = "norm",
                                                      "MVI" = "imp",
                                                      "DEA" = "dea"), inline = TRUE))
           ),
           # conditionalPanel(
           #   condition = "input.sug1cg6.includes('exp') & input.rb6.includes('Y')",
           #   shinydashboard::box(width=6, title = "expression matrix types", solidHeader = TRUE,
           #                       status = "warning",
           #                       checkboxGroupInput("rb61","Data acqusition type:",
           #                                          c("protein intensity" = "protein intensity",
           #                                            "protein MaxLFQ" = "protein MaxLFQ",
           #                                            "peptide intensity" = "peptide intensity",
           #                                            "peptide MaxLFQ" = "peptide MaxLFQ",
           #                                            "spectra counts" = "spectra counts"),inline = T))
           # ),
           conditionalPanel(
             condition = "input.sug1cg6.includes('norm') & input.rb6.includes('Y')",
             shinydashboard::box(width=6, title = "Normalization methods", solidHeader = TRUE,
                                 status = "warning",
                                 checkboxGroupInput("rb62","Normalization:",
                                                    c("No_normalization" = "None",
                                                      "center.mean" = "center.mean",
                                                      "center.median" = "center.median",
                                                      "max" = "max",
                                                      "sum" = "sum",
                                                      "vsn" = "vsn",
                                                      "quantiles" = "quantiles",
                                                      "quantiles.robust" = "quantiles.robust",
                                                      "div.mean" = "div.mean",
                                                      "div.median" = "div.median",
                                                      'lossf' = 'lossf',
                                                      'TIC' = 'TIC',
                                                      "Rlr" = "Rlr",
                                                      'MBQN' = 'MBQN'
                                                    ),inline = T))
           ),
           conditionalPanel(
             condition = "input.sug1cg6.includes('imp') & input.rb6.includes('Y')",
             shinydashboard::box(width=6, title = "MVI algorithms", solidHeader = TRUE,
                                 status = "warning",
                                 checkboxGroupInput("rb63","MVI:",
                                                    c("MinProb" = "MinProb",
                                                      "QRILC" = "QRILC",
                                                      "MinDet" = "MinDet",
                                                      "missForest" = "missForest",
                                                      "nbavg" = "nbavg",
                                                      "zero" = "zero",
                                                      "bpca" = "bpca",
                                                      "MLE" = "MLE",
                                                      "knn" = "knn",
                                                      "No_MVI" = "None",
                                                      "min" = "min",
                                                      "mice" = "mice",
                                                      "Impseq" = "Impseq",
                                                      "Impseqrob" = "Impseqrob",
                                                      "GMS" = "GMS",
                                                      'SeqKNN' = 'SeqKNN'),inline = T))
           ),
           conditionalPanel(
             condition = "input.sug1cg6.includes('dea') & input.rb6.includes('Y')",
             shinydashboard::box(width=6, title = "DEA tools", solidHeader = TRUE,
                                 status = "warning",
                                 checkboxGroupInput("rb64","DEA:",
                                                    c("ANOVA" = "ANOVA",
                                                      "DEP" = "DEP",
                                                      #"DEqMS" = "DEqMS",
                                                      'limma' = "limma",
                                                      "proDA" = "proDA",
                                                      "ROTS" = "ROTS",
                                                      "SAM" = "SAM",
                                                      "ttest" = "ttest"#,
                                                      #"beta_binomial" = "beta_binomial",
                                                      #"edgeR" = "edgeR",
                                                      #"plgem" = "plgem",
                                                      #'MSstats' = 'MSstats'
                                                    ),inline = T))
           )

  )
)
