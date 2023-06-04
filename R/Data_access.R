Datas1<-fluidRow(
  shinydashboard::box(
    title = "Raw data links",
    background = "black",
    width = 12,
    strong(h4("DDA data")),
    p("yeast2099:",a(href = "https://www.ebi.ac.uk/pride/archive/projects/PXD002099",
                     "PXD002099")),
    p("yeast1819:", a(href = "https://www.ebi.ac.uk/pride/archive/projects/PXD001819",
                       "PXD001819")),
    p("CPTAC:",a(href = "https://proteomic.datacommons.cancer.gov/pdc/TechnologyAdvancementStudies",
                "CPTAC-LTQ86; "), a(href = "https://proteomic.datacommons.cancer.gov/pdc/TechnologyAdvancementStudies",
                                    "CPTAC-LTQO65; "),
      a(href = "https://proteomic.datacommons.cancer.gov/pdc/TechnologyAdvancementStudies",
        "CPTAC-LTQP65; "),
      a(href = "https://proteomic.datacommons.cancer.gov/pdc/TechnologyAdvancementStudies",
        "CPTAC-LTQW56")),
    p("human_ecoli:", a(href = "https://www.ebi.ac.uk/pride/archive/projects/PXD018408",
                        "PXD018408")),
    strong(h4("DIA data")),
    h5("RD139:", a(href = "https://www.ebi.ac.uk/pride/archive/projects/PXD026600",
                   "Wide; "),
       a(href = "https://www.ebi.ac.uk/pride/archive/projects/PXD026600",
         "Narrow; "),
       a(href = "https://www.ebi.ac.uk/pride/archive/projects/PXD026600",
         "Overlap; "),
       a(href = "https://www.ebi.ac.uk/pride/archive/projects/PXD026600",
         "Variable")),
    p("human_ecoli:", a(href = "https://www.ebi.ac.uk/pride/archive/projects/PXD018408",
                        "PXD018408"))

  )
)

Datas2<-fluidRow(
  shinydashboard::box(title = "Quantification data",
      background = "black",
      width = 12,
      h4("Please download via the following link"),
      a(href = "https://drive.google.com/drive/folders/1qv-P_0Jhpe1ZevSVL84-ORLkG4B_vcuB?usp=sharing",
        "Quantification data under FragPipe, maxquant, DIA-NN")
  )
  )

Datas3<-fluidRow(
  shinydashboard::box(title = "benchmarking result data",
      background = "black",
      width = 12,
      h4("FragPipe workflow benchmarking"),
      h5("metrics:"),
      downloadLink("dD_fg1", "pAUC0.01s; "),
      downloadLink("dD_fg2", "pAUC0.05s; "),
      downloadLink("dD_fg3", "pAUC0.1s; "),
      downloadLink("dD_fg4", "nMCCs; "),
      downloadLink("dD_fg5", "G-means"),
      h5("ranks:"),
      downloadLink("dD_fg6", "ranks"),
      h4("maxquant workflow benchmarking"),
      h5("metrics:"),
      downloadLink("dD_mq1", "pAUC0.01s; "),
      downloadLink("dD_mq2", "pAUC0.05s; "),
      downloadLink("dD_mq3", "pAUC0.1s; "),
      downloadLink("dD_mq4", "nMCCs; "),
      downloadLink("dD_mq5", "G-means"),
      h5("ranks:"),
      downloadLink("dD_mq6", "ranks"),
      h4("DIA-NN workflow benchmarking"),
      h5("metrics:"),
      downloadLink("dD_dia1", "pAUC0.01s; "),
      downloadLink("dD_dia2", "pAUC0.05s; "),
      downloadLink("dD_dia3", "pAUC0.1s; "),
      downloadLink("dD_dia4", "nMCCs; "),
      downloadLink("dD_dia5", "G-means"),
      h5("ranks:"),
      downloadLink("dD_dia6", "ranks"),


  )
)
