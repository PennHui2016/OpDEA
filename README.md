# OpDEA

## installation
OpDEA is a R shiny application for:
1. presenting results of our benchmarking of proteomics data differential expression analysis workflow;
2. guiding users to select optimal workflows for analyzing their proteomics data;
3. providing links for downloading datasets used for benchmarking and for downloading our results.

It can be installed via two ways:
 1. download OpDEA_0.0.0.9000.tar.gz from this site, then install with the following command:
    
    install.packages(
  pkgs = './OpDEA_0.0.0.9000.tar.gz',
  lib = .libPaths()[length(.libPaths())],
  repos = NULL,
  dependencies = T
)

2. using devtools to install this package:

   install.packages("devtools")
   library(devtools)
   install_github("PennHui2016/OpDEA")

Then, the shiny app can be lanched via:

OpDEA::run_app()

## Required R package
