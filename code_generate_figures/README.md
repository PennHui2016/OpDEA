The Figures in main text can be reproduced with the R scripts of Figure2.R, Figure3, Figure4.R and Figure5.R respectively.

To reproduce the workflow performance level classification, the following python packages are required to be installed (package version are behind package names):

   pandas-v2.0.3, numpy-v1.24.3, sklearn-v1.3.0, fpgrowth_py-v1.0.0, mlxtend-v0.23.0, catboost-v1.2.2
   
   then, run the codes with the following command with python(python 3.11 was used by the author), :
     
	  python Catboost_workflows_classification.py
	  python Catboost_workflows_classification_TMT_mq.py

To reproduce the figures, the following R packages (version numbers are behined package names) are required to be installed:
   
   readxl-v1.4.3, ggplot2-v3.4.3, reshape2-v1.4.4, ggpubr-v0.6.0, cowplot-v1.1.1, openxlsx-v4.2.5.2, ggalluvial-v0.12.5, ggsci-v3.0.0,
   ggrepel-v0.9.3, ggthemes-v5.0.0, ComplexHeatmap-v2.16.0, gridExtra-v2.3, UpSetR-v1.4.0, grid-v4.3.1, pROC-v1.18.4, ggradar-v0.2
   
   
 The codes were tested with following environment:
 R version 4.3.1 (2023-06-16 ucrt)
 Platform: x86_64-w64-mingw32/x64 (64-bit)
 Running under: Windows 11 x64 (build 22631)


