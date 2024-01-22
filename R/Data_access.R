Datas1<-fluidRow(
  shinydashboard::box(
    title = "Raw data links",
    background = "black",
    width = 12,
    #strong(h3("Raw data links")),
    strong(h4("label-free DDA data")),
    column(2,p("HYE5600735_LFQ:",a(href = "https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD028735",
                                   "PXD028735"))),
    column(2,p("HYE6600735_LFQ:", a(href = "https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD028735",
                                    "PXD028735"))),
    column(2,p("HYEqe735_LFQ:",a(href = "https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD028735",
                                 "PXD028735"))),
    column(2,p("HYEtims735_LFQ:", a(href = "https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD028735",
                                    "PXD028735"))),
    column(2,p("HYtims134_LFQ:",a(href = "https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD036134",
                                  "PXD036134"))),
    column(2,p("HEtims425_LFQ:", a(href = "https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD021425",
                                   "PXD021425;"))),
    column(2,p("YUltq006_LFQ:",a(href = "https://proteomic.datacommons.cancer.gov/pdc/TechnologyAdvancementStudies/",
                                 "PDC000006"))),
    column(2,p("YUltq099_LFQ:", a(href = "https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD002099",
                                  "PXD002099"))),
    column(2,p("YUltq819_LFQ:",a(href = "https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD001819",
                                 "PXD001819"))),
    column(2,p("HEqe408_LFQ:", a(href = "https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD018408",
                                 "PXD018408"))),
    column(2,p("HYqfl683_LFQ:",a(href = "https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD007683",
                                 "PXD007683"))),
    column(2,p("HYEtims777_LFQ:", a(href = "https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD014777",
                                    "PXD014777;"))),
    column(12,br()),
    strong(h4("label-free DIA data")),
    #column(12,br()),
    column(3,p("HYEtims735_DIA:",a(href = "https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD028735",
                                   "PXD028735"))),
    column(3,p("MYtims709_DIA:",a(href = "http://www.iprox.org/page/project.html?id=IPX0004576000",
                                  "PXD034709"))),
    column(3,p("HEqe408_DIA:",a(href = "https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD018408",
                                "PXD018408"))),
    column(3,p("HEof_n600_DIA:",a(href = "https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD026600",
                                  "PXD026600"))),
    column(3,p("HEof_w600_DIA:",a(href = "https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD026600",
                                  "PXD026600"))),
    column(3,p("HYtims134_DIA:",a(href = "https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD036134",
                                  "PXD036134"))),
    column(3,p("HEqe777_DIA:",a(href = "https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD019777",
                                "PXD019777"))),
    column(3,p("HEqe777_DIA:",p(''))),
    column(12,br()),
    strong(h4("TMT data")),
    #column(12,br()),
    column(2,p("HEqe277_TMT10:",a(href = "https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD013277",
                                  "PXD013277"))),
    column(2,p("HYqfl683_TMT11:",a(href = "https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD007683",
                                   "PXD007683"))),
    column(2,p("HYms2faims815_TMT16:",a(href = "https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD020815",
                                        "PXD020815"))),
    column(2,p("HYsps2815_TMT16:",a(href = "https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD020815",
                                    "PXD020815"))),
    column(2,p("HYms2815_TMT16:",a(href = "https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD020815",
                                   "PXD020815")))
  )
)

Datas2<-fluidRow(
  shinydashboard::box(title = "Quantification data",
                      background = "black",
                      width = 12,
                      h4("Please download from zenodo via the following link"),
                      a(href = "https://zenodo.org/records/10482353",
                        "Zenodo: Raw quantification results from FragPipe, Maxquant, DIA-NN and Spectronaut"),
                      h4("Or download from google drive via the following link"),
                      a(href = "https://drive.google.com/file/d/1T3MHDMF5NokVCvFTN5IV3TABVAl-0qw1/view?usp=sharing",
                        "google Drive: Raw quantification results of DDA data; "),
                      a(href = "https://drive.google.com/file/d/1Q-u962AuCCAd8jQRVfiFqe0m03qAQI-8/view?usp=sharing",
                        "  DIA data; "),
                      a(href = "https://drive.google.com/file/d/1fYyHebPVW-t7iw67BWY-wvDziH00gEnB/view?usp=sharing",
                        "  TMT data")
  )
)

Datas3<-fluidRow(
  shinydashboard::box(title = "benchmarking result data",
      background = "black",
      width = 12,
      h4("expression matrices are available at: "),

      #column(4,downloadLink("DDA_exp", "DDA_expression_matrices; ")),
      column(4,h5("Download from zenodo via the following link"),
             a(href = "https://zenodo.org/records/10484253",
               "Zenodo: DDA_expression_matrices;")),
      #column(4,downloadLink("DIA_exp", "DIA_expression_matrices; ")),
      column(4,h5("Download from zenodo via the following link"),
             a(href = "https://zenodo.org/records/10484253",
               "Zenodo: DIA_expression_matrices;")),
      #column(4,downloadLink("TMT_exp", "TMT_expression_matrices; ")),
      column(4,h5("Download from zenodo via the following link"),
             a(href = "https://zenodo.org/records/10484253",
               "Zenodo: TMT_expression_matrices;")),

      h4("Performance metrices are available at: "),
      #column(4,downloadLink("DDA_met", "DDA_performance_metrics; ")),
      column(4,h5("Download from zenodo via the following link"),
             a(href = "https://zenodo.org/uploads/10484428",
               "Zenodo: DDA_performance_metrics;")),
      #column(4,downloadLink("DIA_met", "DIA_performance_metrics; ")),
      column(4,h5("Download from zenodo via the following link"),
             a(href = "https://zenodo.org/uploads/10484428",
               "Zenodo: DIA_performance_metrics;")),
      #column(4,downloadLink("TMT_met", "TMT_performance_metrics; ")),
      column(4,h5("Download from zenodo via the following link"),
             a(href = "https://zenodo.org/uploads/10484428",
               "Zenodo: TMT_performance_metrics;")),

      h4("Workflow ranks are available at: "),
      #column(4,downloadLink("DDA_rk", "DDA_workflow_ranks; ")),
      column(4,h5("Download from zenodo via the following link"),
             a(href = "https://zenodo.org/uploads/10484428",
               "Zenodo: DDA_workflow_ranks;")),
      #column(4,downloadLink("DIA_rk", "DIA_workflow_ranks; ")),
      column(4,h5("Download from zenodo via the following link"),
             a(href = "https://zenodo.org/uploads/10484428",
               "Zenodo: DIA_workflow_ranks;")),
      #column(4,downloadLink("TMT_rk", "TMT_workflow_ranks; ")),
      column(4,h5("Download from zenodo via the following link"),
             a(href = "https://zenodo.org/uploads/10484428",
               "Zenodo: TMT_workflow_ranks;")),


  )
)