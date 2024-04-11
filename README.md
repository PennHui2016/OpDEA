# OpDEA
OpDEA is a R shiny application for:
1. presenting results of our benchmarking of proteomics data differential expression analysis workflow;
2. guiding users to select optimal workflows for analyzing their proteomics data;
3. For analyzing the user's proteomics data; 
4. providing links for downloading datasets used for benchmarking and for downloading our results.
   
## webserver and standalone toolkit (recommend)
We prepared a freely-accessible webserver for helping users to use OpDEA without installation of the package, 
see http://www.ai4pro.tech:3838/. The webserver requires uploading expression data, so we recommend using our
standalone toolkit (available at: https://doi.org/10.5281/zenodo.10958381, no R package needs to be installed and no need to upload data,
just decompress it and use it) instead. If you still hope to try our R package, please see  following installation instructions.

## Requirements for installation of OpDEA
You should install the following packages with the same versionor higher:

R base: R-4.2.0;
R packages: shiny 1.8.0; shinydashboard 0.7.2; threejs 0.3.3; DT 0.31; ggplot2 3.4.4; reshape2 1.4.4; ggpubr 0.6.0; ggsci 3.0.0; readxl 1.4.3; ggalluvial 0.12.5; golem 0.4.1;
iq 1.9.12; dplyr 1.1.4; aggregation 1.0.1; stringr 1.5.1; tidyverse 2.0.0; matrixStats 1.2.0; readr 2.1.4; rrcovNA 0.5.0; BiocManager 1.30.22; NormalyzerDE 1.16.0; limma 3.54.2;
ROTS 1.26.0; MSnbase 2.24.2; edgeR 3.40.2; proDA 1.12.0; DEqMS 1.16.0; plgem 1.70.0; DEP 1.20.0; MSstats 4.6.5; samr 3.0.0; mice 3.16.0; missForest 1.5; SeqKnn 1.0.1; 
GMSimpute 0.0.1.0 (see https://github.com/wangshisheng/NAguideR)

The whole R session of my R environment are as follows:

```
R version 4.2.0 (2022-04-22 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 22631)

Matrix products: default  

attached base packages:

[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:

[1] shinydashboard_0.7.2 shiny_1.8.0         

loaded via a namespace (and not attached):

 [1] circlize_0.4.15     shape_1.4.6         tidyselect_1.2.0    bslib_0.6.1        
 [5] purrr_1.0.2         colorspace_2.1-0    vctrs_0.6.5         generics_0.1.3     
 [9] htmltools_0.5.7     yaml_2.3.8          base64enc_0.1-3     utf8_1.2.4         
[13] rlang_1.1.2         jquerylib_0.1.4     later_1.3.2         pillar_1.9.0       
[17] glue_1.6.2          withr_2.5.2         lifecycle_1.0.4     munsell_0.5.0      
[21] gtable_0.3.4        ragg_1.2.7          fontawesome_0.5.2   ggalluvial_0.12.5  
[25] htmlwidgets_1.6.4   GlobalOptions_0.1.2 memoise_2.0.1       labeling_0.4.3     
[29] fastmap_1.1.1       golem_0.4.1         httpuv_1.6.13       fansi_1.0.6        
[33] Rcpp_1.0.11         xtable_1.8-4        promises_1.2.1      scales_1.3.0       
[37] DT_0.31             cachem_1.0.8        jsonlite_1.8.8      config_0.3.2       
[41] farver_2.1.1        mime_0.12           systemfonts_1.0.5   textshaping_0.3.7  
[45] ggplot2_3.4.4       digest_0.6.33       dplyr_1.1.4         grid_4.2.0         
[49] OpDEA_0.0.0.9000    cli_3.6.2           tools_4.2.0         magrittr_2.0.3     
[53] sass_0.4.8          tibble_3.2.1        tidyr_1.3.0         pkgconfig_2.0.3    
[57] ellipsis_0.3.2      rsconnect_1.2.0     attempt_0.3.1       rstudioapi_0.15.0  
[61] R6_2.5.1            compiler_4.2.0 
```
## Installation

It can be installed via two ways:
1. Install the package "devtools" if you have not installed it before

    ```
    if(!requireNamespace("devtools")){
       install.packages("devtools")
    }
    ```

Then, the package can be installed from github via the following code:

    library(devtools)
    install_github('PennHui2016/OpDEA')
   

2.Or via downloading **"OpDEA_0.0.0.9000.tar.gz"** from this site, then installed with the following command:
    
    install.packages(pkgs = '~/OpDEA_0.0.0.9000.tar.gz', repos = NULL, type = "source")

At last, the shiny app can be launched via:

    OpDEA::run_app()

If success, the page showing the introduction of our OpDEA will be presented. It can be used according to the contents in the help page.

## source codes for proteomics data differential expression analysis workflow benchmarking
The python and R source codes for benchmarking proteomics data differential expression analysis workflows are located in the folder **"codes_DEA_benchmarking"**.
Please read the **"README"** file inside the **"codes_DEA_benchmarking"** folder to find how to reproduce our benchmarking results.

## source codes for regenerating the figures in our paper
The R codes for regenerating our figures in our paper can be found in the folder **"code_generate_figures"**. Please read the 
**"README"** file inside the **"code_generate_figures"** folder to find how to regenerate our figures.

## Cite this article
Hui Peng, He Wang, Weijia Kong, Jinyan Li*, Wilson Wen Bin Goh*. (2024). Optimizing Differential Expression Analysis for Proteomics Data via High-Performing Rules and Ensemble Inference.

 
## Contact
Any problems or requesting source codes for reproducing results in our paper please contact

    Hui Peng: hui.peng@ntu.edu.sg;  Wilson Wen Bin Goh: wilsongoh@ntu.edu.sg

