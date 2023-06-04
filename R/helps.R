help_tab1<-fluidRow(
  shinydashboard::box(width = 12, height = 150,
      status = "primary", solidHeader = T,
      h4("1. About workflow benchmarking:"),
      p('Our benchmarking are based on three quantification platforms including FragPipe, maxquant and DIA-NN.
        There are 3004, 3004 and 4800 workflows for these three quantification platforms respectively, composed of choices of
        expression martrix types, normalization methods, missing value imputation algorithms and dfferential expression analysis
        tools (see more details of workflows in our Introduction page).'),
  ),
  shinydashboard::box(width = 12, height = 150,
      status = "primary", solidHeader = T,
      h4("2.  About contents shown in Benchmarking sheets:"),
      p("We have three sheets to present benchmarking results for workflows based on FragPipe, maxquant and DIA-NN. By clicking
         the corresponding platform name under the Benchmarking menuitem, users can browse different sheets. In each sheet,
         the top 4 boxes shows the filtering choices. The left-bottom box contains our benchmarking results showing in a table.
         the right-bottom box is used to show performance distributions of a selected workflow in the table of the left-bottom.
        ")
  ),
  shinydashboard::box(width = 12, height = 200,
      status = "primary", solidHeader = T,
      h4("3. How to view our benchmarking results"),
      p("If users hope to view benchmarking results of all the workflows, no item needs to be checked in the top 4 boxes. Just
      the tabke in the left-bottom box. In the table, each row presents detail information of a workflow, where the name of
      this workflow, the options chosed in each step of the workflow, the mean and
        median performances (indicated by metrics including pAUC0.01, pAUC0.05, pAUC0.1, nMCC and G-mean), ranks of it
        according to its performance scores, and its final ranks are shown."),
      p("Users also can view the workflows involve one or more specific expression matrix types (left-top box), one or more
        normalization method (right-top box), one or more missing value imputation algorithms (middle-left box), and one or
        more differential expression analysis tools (middle-right box), by checking correspond items in the boxes.")
  ),
  shinydashboard::box(width = 12, height = 150,
      status = "primary", solidHeader = T,
      h4("4. How to view a workflow's performance distributions"),
      p('When users click one row of the left-bottom table, then a boxplot will be shown in the right-bottom box. The
        distributions of five performance indicators including partial area under receiver operator characteristic curve
        at FPR threshold of 0.01 (pAUC0.01), pAUC0.05, pAUC0.1, normalizated Matthew’s correlation coefficient (nMCC) and
        geometric mean of specificity and recall (G-mean) are presented.'),
      p('We note that only one workflow can be chosen to present its performance distributions', style = "color:red"),
  ),
  shinydashboard::box(width = 12, height = 100,
      status = "primary", solidHeader = T,
      h4("5. How to obtain our benchmarking results"),
      p("Our benchmarking results can be downloaded via links in the page of Data.")

  )
)

help_tab2<-fluidRow(
  shinydashboard::box(width = 12, height = 150,
      status = "primary", solidHeader = T,
      h4("1. About workflow suggestion:"),
      p('In Introduction page, we show our leave-one-project-out cross-validation results to confirm that our benchmarking has
        good generalizability where optimal workflow could be suggested per our benchmarking. Then, we used the 10-fold
        cv to prove that a workflow performance level is predictable when we know the options in each step. We have found some
        frequent patterns that appear in high-performance workflows. In addition, we found that applying ensemble inference by
        integrating top-ranked workflows could improve differential expression analysis. Thus, our workflow suggestion is based
        on our benchmarking results the found frequent patterns and ensemble inference tests.'),
  ),
  shinydashboard::box(width = 12, height = 150,
      status = "primary", solidHeader = T,
      h4("2. About single workflow suggestion:"),
      p("We have three sheets to present benchmarking results for workflows based on FragPipe, maxquant and DIA-NN. By clicking
         the corresponding platform name under the Benchmarking menuitem, users can browse different sheets. In each sheet,
         the top 4 boxes shows the filtering choices. The left-bottom box contains our benchmarking results showing in a table.
         the right-bottom box is used to show performance distributions of a selected workflow in the table of the left-bottom.
        ")
  ),
  shinydashboard::box(width = 12, height = 200,
      status = "primary", solidHeader = T,
      h4("3. About ensemble inference suggestion"),
      p("If users hope to view benchmarking results of all the workflows, no item needs to be checked in the top 4 boxes. Just
      the tabke in the left-bottom box. In the table, each row presents detail information of a workflow, where the name of
      this workflow, the options chosed in each step of the workflow, the mean and
        median performances (indicated by metrics including pAUC0.01, pAUC0.05, pAUC0.1, nMCC and G-mean), ranks of it
        according to its performance scores, and its final ranks are shown."),
      p("Users also can view the workflows involve one or more specific expression matrix types (left-top box), one or more
        normalization method (right-top box), one or more missing value imputation algorithms (middle-left box), and one or
        more differential expression analysis tools (middle-right box), by checking correspond items in the boxes.")
  )
)

help_tab3<-fluidRow(
  shinydashboard::box(width = 12, height = 400,
      status = "primary", solidHeader = T,
      h4("About workflow benchmarking:"),
      p('Our benchmarking are based on three quantification platforms including FragPipe, maxquant and DIA-NN.
        There are 3004, 3004 and 4800 workflows for these three quantification platforms respectively, composed of choices of
        expression martrix types, normalization methods, missing value imputation algorithms and dfferential expression analysis
        tools (see more details of workflows in our Introduction page).')
  )
)
