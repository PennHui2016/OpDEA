root_fold<-'inst/app/www/'
b641 <- base64enc::dataURI(file=paste0(root_fold,"introduction_server.png"), mime="www/png")
b642 <- base64enc::dataURI(file=paste0(root_fold,"tutorial_viewing_benchmarking.png"), mime="www/png")
b643 <- base64enc::dataURI(file=paste0(root_fold,"toolkit_tutorial.png"), mime="www/png")
help_tab0<-fluidRow(
  shinydashboard::box(width = 12, height = 1200,
                      status = "primary", solidHeader = T,
                      strong("Introdcution of the Toolkit:"),
                      hr(),
                      p('This toolkit includes 5 main function panels including Introduction,  Benchmarking, Suggestion & DEA, Data and Help.
      In the [Introduction] panel, we first present the abstract of this work, then we show an overview of the proteomics data differential expression analysis workflow and the options available in each workflow steps. We later using the leave-one-dataset-out cross-validation results to show the good
      generalizability of our benchmarking which supports the realiability of our webserver for suggesting optimal workflows for newcoming data.
      We also usd the 10-fold cross-validation to prove the performance level of a workflow can be predict with a CatBoost classifier, which also
      lays the foundation of the ability in recommendation of optimal workflows by our webserver. At last, the Acknowledgement and Citation information are shown.'),
                      hr(),
                      strong('In the [Benchmarking] panel,'),
                      p('users can view our benchmarking results including checking the rank position and performance
      metric values of a workflow tested with our benchmaking datasets.'),
                      hr(),
                      strong('In the [Suggestion & DEA] pane, '),
                      p('we provide the tools for suggesting optimal workflows and conduct differential expression analysis with the
       suggested workflow directly.'),
                      hr(),
                      strong('In the [Data] pane, '),
                      p('the user can get the links where raw proteomics data are available. The raw quantification results of the raw data, the extracted expression matrices and our benchmarking results can be downloaded.'),
                      p('We also provide the link for downloading our offline toolkit with the same function as the webserver.',style = "color:red"),
                      hr(),
                      strong('In the [Help] pane, '),
                      p('We introduce the webserver and show what the users can do with our toolkit'),
                      # img(src = "img/introduction_server.png",
                      #          height = "500px", width = "1100px", align = "center")
                      img(src=b641, height = "500px", width = "1100px", align = "center")
  )
)

help_tab1<-fluidRow(
  shinydashboard::box(width = 12, height = 1200,
                      status = "primary", solidHeader = T,
                      strong("View workflow benchmarking results:"),
                      hr(),
                      strong('Step 1: Choosing an interested setting, e.g., label-free DDA data quantified with FragPipe, click the item of DDA_LFQ-FragPipe.'),
                      p('There are 7852, 7852, 6284, 6284, 4720 and 1568 workflows under settings of LFQ_DDA-FragPipe, LFQ_DDA-Maxquant, LFQ_DIA-DIANN,
      LFQ_DIA-Spectronaut, TMT-FragPipe and TMT-Maxquant.'),
                      hr(),
                      strong('Step 2: Filtering workflows, e.g., check the options in the 4 checkboxes that you are interested in, or by using a key word.'),
                      p('Drag the scroll bar below the table can check the ranking with a specific metric. The ranking of workflows is based-on the average
      ranking based on the five performance metrics (the column of avg_rank_mean).'),
                      hr(),
                      strong('Step 3: Click one row of the table to check the performance distributions testing on different datasets.'),
                      p('12 DDA datasets, 7 DIA datasets and 5 TMT datasets were used to evaluate the performances of workflows, the boxplot shows the performance
      distributions of the 5 metrics.'),
                      p('Only one row is permitted to be select each time, click the same row twice can cancel the selection.',style = "color:red" ),
                      hr(),
                      strong('Step 4: View the details of performance distrubutions of the five metrics in the right bottom figure.'),

                      # img(src = "img/tutorial_viewing_benchmarking.png",
                      #          height = "500px", width = "1100px", align = "center")
                      img(src=b642, height = "500px", width = "1100px", align = "center")
  )
)

help_tab2<-fluidRow(
  shinydashboard::box(width = 12, height = 1800,
                      status = "primary", solidHeader = T,
                      strong("Workflow suggestion and testing:"),
                      hr(),
                      strong('Step 1: Choose data type that the user has.'),
                      p('Lable-free DDA (DDA), label free DIA (DIA) and TMT data are supported. Choose the one correponding to the data you have.'),
                      hr(),
                      strong('Step 2: Choose a quantification platform.'),
                      p('For DDA data, FragPipe and Maxquant are supported. For DIA data, DIA-NN and Spectronaut are
      supported. For TMT data, again, FragPipe and Maxquant are supported. Just choose the one that you used to quantify your proteomics data.'),
                      hr(),
                      strong('Step 3: Choose to suggest single workflow or try the ensemble inference.'),
                      p('If single workflow is preferred. The user should choose whether the have some preferred options in the workflow step. If yes, then check Yes, and select the options in below checkboxes, otherwise select No.'),
                      p('For DDA and DIA data, our server also support to suggest ensemble inference where multiple workflow are integrated by a p-value integration method.
      The ens_multi-quant approach is used defaulty. '),
                      p('The imputation method missForest and MLE are quite time-consuming, if they are suggested, we will replace them with MinProb.',style = "color:red" ),
                      hr(),
                      strong('Step 4: Click the suggest worklfow button or Try ensemble inference button.'),
                      p('Our server will suggest the top 1st workflow after filtering the workflows with user specified option preference. The potential alternative
      options will also be shown according to our option comparsion results'),
                      hr(),
                      strong('Step 5: specify the path to the raw quantification result data, the path to python executable where drectlfq python package has been installed, and specify thresholds and submit the DEA task.'),
                      p('The user should specify the file paths of their raw quantification result file, e.g., the combined_protein.tsv file from FragPipe.
                      the directLFQ python package should be installed, see the tutorial at: https://github.com/MannLabs/directlfq.
                        The designed file showing the sample condition and group information must be uploaded at the same time.
                        The log2FC threshold and adj.pvalue (p-value is ajusted with BH method) threshold should be specified. At last, click DEA button to submit the task.'),
                      hr(),
                      strong('Step 6: Check DEA results.'),
                      p('After submitting the DEA task, the user should wait for a while till the task being completed. A volcano plot will be generated. The result file are stored locally, their file folder
                        will be shown above the volcano plot.'),
                      # img(src = "img/toolkit_tutorial.png",
                      #          height = "1000px", width = "1100px", align = "center")
                      img(src=b643, height = "1000px", width = "1100px", align = "center")
  )
)

help_tab3<-fluidRow(
  shinydashboard::box(width = 12, height = 400,
                      status = "primary", solidHeader = T,
                      h4("Contact:"),
                      p('Any problem please contact Hui Peng: '),
                      p('Email: hui.peng@ntu.edu.sg'),
                      br(),
                      p('Or, can email to the corresponding authors: '),
                      p('Jinyan Li: jinyan.li@siat.ac.cn; Wilson Wen Bin Goh: wilsongoh@ntu.edu.sg' )
  )
)
