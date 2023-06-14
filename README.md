# OpDEA

## Installation
OpDEA is a R shiny application for:
1. presenting results of our benchmarking of proteomics data differential expression analysis workflow;
2. guiding users to select optimal workflows for analyzing their proteomics data;
3. providing links for downloading datasets used for benchmarking and for downloading our results.

It can be installed via downloading OpDEA_0.0.0.9000.tar.gz from this site, then install with the following command:
    
    install.packages(
  pkgs = './OpDEA_0.0.0.9000.tar.gz',
  lib = .libPaths()[length(.libPaths())],
  repos = NULL,
  dependencies = T
)

Then, the shiny app can be launched via:

OpDEA::run_app()

## Requirements for installation of OpDEA
R base: R-4.1.2 or higher;
R packages: shiny 1.7.4 or higher; shinydashboard 0.7.2 or higher; threejs; DT; ggplot2; reshape2; ggpubr; ggsci; readxl; ggalluvial

## webserver
We prepared a freely-accessible webserver for helping users using OpDEA without install the packages, see http://www.ai4pro.tech:3838/OpDEA_test/

## source codes for proteomics data differential expression analysis workflow benchmarking
The python and R source codes for benchmarking proteomics data differential expression analysis workflows are located in the folder ""
## Contact
Any problems or requesting source codes for reproducing results in our paper please contact

    Hui Peng: hui.peng@ntu.edu.sg;  Wilson Wen Bin Goh: wilsongoh@ntu.edu.sg

