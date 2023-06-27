# OpDEA

## Installation
OpDEA is a R shiny application for:
1. presenting results of our benchmarking of proteomics data differential expression analysis workflow;
2. guiding users to select optimal workflows for analyzing their proteomics data;
3. providing links for downloading datasets used for benchmarking and for downloading our results.

It can be installed via downloading **"OpDEA_0.0.0.9000.tar.gz"** from this site, then installed with the following command:
    
    install.packages(pkgs = '~/OpDEA_0.0.0.9000.tar.gz', repos = NULL, type = "source")

Then, the shiny app can be launched via:

    library(shinydashboard)
    OpDEA::run_app()

## Requirements for installation of OpDEA
please install the packages with the following version numbers or higher:

R base: R-4.1.2;
R packages: shiny 1.7.4; shinydashboard 0.7.2; threejs 0.3.3; DT 0.28; ggplot2 3.3.6; reshape2 1.4.4; ggpubr 0.4.0; ggsci 2.9; readxl 1.4.0; ggalluvial 0.12.4; golem 0.4.0

## webserver
We prepared a freely-accessible webserver for helping users to use OpDEA without installation of the package, see http://www.ai4pro.tech:3838/OpDEA_test/.

## source codes for proteomics data differential expression analysis workflow benchmarking
The python and R source codes for benchmarking proteomics data differential expression analysis workflows are located in the folder **"source_codes_for_benchmarking"**.

## Contact
Any problems or requesting source codes for reproducing results in our paper please contact

    Hui Peng: hui.peng@ntu.edu.sg;  Wilson Wen Bin Goh: wilsongoh@ntu.edu.sg

