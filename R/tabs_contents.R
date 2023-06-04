root_fold<-'inst/app/www/'

tab1 <- fluidRow(
  shinydashboard::box(width=6, height = 100, checkboxGroupInput("filter1", "expression type of DDA:",
                         c("protein intensity" = "protein intensity",
                           "protein MaxLFQ" = "protein MaxLFQ",
                           "peptide intensity" = "peptide intensity",
                           "peptide MaxLFQ" = "peptide MaxLFQ",
                           "spectra counts" = "spectra counts"), inline = TRUE)),
  shinydashboard::box(width=6, height = 100,checkboxGroupInput("filter2", "normalization:",
                         c("No_normalization" = "No_normalization",
                           "center.mean" = "center.mean",
                           "center.median" = "center.median",
                           "max" = "max",
                           "sum" = "sum",
                           "vsn" = "vsn",
                           "quantiles" = "quantiles",
                           "quantiles.robust" = "quantiles.robust",
                           "div.mean" = "div.mean",
                           "div.median" = "div.median"
                         ), inline = TRUE)),
  shinydashboard::box(width=6, height = 100,checkboxGroupInput("filter3", "MVI:",
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
                           "mice" = "mice"), inline = TRUE)),
  shinydashboard::box(width=6, height = 100,checkboxGroupInput("filter4", "DEA tool:",
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
                           "ProteoMM" = "ProteoMM"
                         ), inline = TRUE)),
  shinydashboard::box(width = 6, height = 600,
      "workflow benchmarking",
      br(),
      column(width = 12,DT::dataTableOutput("tabFrag"), style = "height:500px; overflow-x: scroll;")
  ),
  shinydashboard::box(width = 6, height = 600,
      "performance distributions",
      br(),
      plotOutput("boxplot1", height = 500, width = 500)
  )
)

tab2 <- fluidRow(
  shinydashboard::box(width=6, height = 100,checkboxGroupInput("filter5", "expression type of DDA:",
                         c("protein intensity" = "protein intensity",
                           "protein MaxLFQ" = "protein MaxLFQ",
                           "peptide intensity" = "peptide intensity",
                           "peptide MaxLFQ" = "peptide MaxLFQ",
                           "spectra counts" = "spectra counts"), inline = TRUE)),
  shinydashboard::box(width=6, height = 100,checkboxGroupInput("filter6", "normalization:",
                         c("None" = "None",
                           "center.mean" = "center.mean",
                           "center.median" = "center.median",
                           "max" = "max",
                           "sum" = "sum",
                           "vsn" = "vsn",
                           "quantiles" = "quantiles",
                           "quantiles.robust" = "quantiles.robust",
                           "div.mean" = "div.mean",
                           "div.median" = "div.median"
                         ), inline = TRUE)),
  shinydashboard::box(width=6, height = 100,checkboxGroupInput("filter7", "MVI:",
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
                           "mice" = "mice"), inline = TRUE)),
  shinydashboard::box(width=6, height = 100,checkboxGroupInput("filter8", "DEA tool:",
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
                           "ProteoMM" = "ProteoMM"
                         ), inline = TRUE)),
  shinydashboard::box(width = 6, height = 600,
      "workflow benchmarking",
      br(),
      column(width = 12,DT::dataTableOutput("tabMq"), style = "height:500px; overflow-x: scroll;")
  ),
  shinydashboard::box(width = 6, height = 600,
      "performance distributions",
      br(),
      plotOutput("boxplot2", height = 500, width = 500)
  )
)

tab3 <- fluidRow(
  shinydashboard::box(width=6, height = 100,checkboxGroupInput("filter9", "expression type of DIA:",
                         c("MaxLFQ" = "MaxLFQ",
                           "norm" = "norm",
                           "MaxLFQ_raw" = "MaxLFQ_raw",
                           "norm_raw" = "norm_raw",
                           "raw" = "raw"), inline = TRUE)),
  shinydashboard::box(width=6, height = 100,checkboxGroupInput("filter10", "normalization:",
                         c("None" = "None",
                           "center.mean" = "center.mean",
                           "center.median" = "center.median",
                           "max" = "max",
                           "sum" = "sum",
                           "vsn" = "vsn",
                           "quantiles" = "quantiles",
                           "quantiles.robust" = "quantiles.robust",
                           "div.mean" = "div.mean",
                           "div.median" = "div.median"
                         ), inline = TRUE)),
  shinydashboard::box(width=6, height = 100,checkboxGroupInput("filter11", "MVI:",
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
                           "mice" = "mice"), inline = TRUE)),
  shinydashboard::box(width=6, height = 100,checkboxGroupInput("filter12", "DEA tool:",
                         c("ANOVA" = "ANOVA",
                           "DEP" = "DEP",
                           'limma' = "limma",
                           "proDA" = "proDA",
                           "ROTS" = "ROTS",
                           "SAM" = "SAM",
                           "siggenes" = "siggenes",
                           "ttest" = "ttest"
                         ), inline = TRUE)),
  shinydashboard::box(width = 6, height = 600,
      "workflow benchmarking",
      br(),
      column(width = 12,DT::dataTableOutput("tabDia"), style = "height:500px; overflow-x: scroll;")
  ),
  shinydashboard::box(width = 6, height = 600,
      "performance distributions",
      br(),
      plotOutput("boxplot3", height = 500, width = 500)
  )
)

tabIntro<-fluidRow(
  shinydashboard::box(
    title = "Abstract",
    background = "green",
    width = 12,
    p("The process of identifying phenotype-specific or differentially expressed proteins from proteomic data typically involves five steps: quantifying raw data, constructing an expression matrix, normalizing the matrix, imputing missing data, and conducting differential expression analysis. However, since there are multiple options available for each step, selecting an inappropriate or suboptimal combination of options may not yield optimal results. To learn optimal combinations to create high performing workflows, we conducted 10,808 experiments to assess the performance of exhaustive combinations of options for each step across 12 gold standard spike-in datasets. We ranked these workflows based on five performance indicators and found that the top-performing workflows were unique to the quantitative platform being used (FragPipe, MaxQuant, or DIA-NN). On average, leave-one-project-out cross-validation (LOPOCV)  Spearman correlations for benchmarking FragPipe and MaxQuant workflows were higher than 0.7, demonstrating their superior ability to recommend optimal workflows for new data. We also found that workflows could be accurately classified into performance levels with average F1 scores and Matthew's correlation coefficients no less than 0.79 in 10-fold cross-validations for all three platforms. The selection of DEA tools and normalization methods was more critical in workflow performance level classification than other steps. Based on our frequent pattern mining and pairwise-option comparison of workflow contents, we derived optimal workflow recommendation rules and developed them into a freely accessible webserver (http://www.xxx.com). Additionally, we proposed integrating the top-ranked workflows through ensemble inference to achieve better performance than the individual best performance.
"))
)

tabWf_intro<-fluidRow(
  shinydashboard::box(
    title = "Workflow of a DEA process for label-free proteomics data and available tools for each step in the workflow",
    background = "navy",
    width = 12,
    img(src = "img/Figure 1_workflows_1.png",
        height = "500px", width = "1100px", align = "center"))
)

tabLOPOCV<-fluidRow(
  shinydashboard::box(
    title = "Good Generalizability confirmed by Leave-One-Project-Out Cross-Validation",
    background = "purple",
    width = 12,
    textOutput("LOPOCV performances"),
    br(),
    column(4,plotOutput("fgcv")),
    column(4,plotOutput("mqcv")),
    column(4,plotOutput("diacv")))
)

tabcls<-fluidRow(
  shinydashboard::box(
    title = "Workflow performance levels are predictable",
    background = "light-blue",
    width = 12,
    textOutput("10-fold cross-validation tests results of CatBoost classifiers"),
    br(),
    column(6,plotOutput("cvCls")),
    column(6,plotOutput("feaImp"))
  )
)

tabAC<-fluidRow(
  shinydashboard::box(
    title = "Acknowledgement and Citation",
    background = "orange",
    width = 12,
    h3("Acknowledgement"),
    p("This research/project is supported by the National Research Foundation, Singapore under its Industry Alignment Fund – Prepositioning (IAF-PP) Funding Initiative. Any opinions, findings and conclusions or recommendations expressed in this material are those of the author(s) and do not reflect the views of National Research Foundation, Singapore."),
    p("WWBG also acknowledges support from an MOE Tier 1 award (RS08/21)."),
    br(),
    h3('Publication'),
    p("Please cite the following paper:"),
    p("Hui Peng, He Wang, Weijia Kong, Jinyan Li*, Wilson Wen Bin Goh*. (2023). Optimizing Proteomics Data Differential Expression Analysis via Unveiling High-Performing Rules and Ensemble Inference")
  )
)
