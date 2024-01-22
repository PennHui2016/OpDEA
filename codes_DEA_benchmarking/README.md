The benchmarking codes are written with Python3, but the python codes will call R scripts to conduct differential expression analysis.

Please follow the below procedures to reproduce the benchmarking results:

1. download the expression matrix data from "https://zenodo.org/records/10484253", and
   decompress the files into the folder 'data/'
   
2. install required R and python packages including:

   python: 
     
	 pandas-v2.0.3, numpy-v1.24.3, sklearn-v1.3.0, scipy-v1.11.1
	 
   R:
     
	 NormalyzerDE-v1.18.1, limma-v3.56.2, ROTS-v1.28.0, MSnbase-v2.26.0, edgeR-v3.42.4,
	 proDA-v1.14.0, DEqMS-v1.18.0, plgem-v1.72.0, DEP-v1.22.0, MSstats-v4.8.7, samr-v3.0,
	 mice-v3.16.0, missForest-v1.5, SeqKnn-v1.0.1, reshape2-v1.4.4, dplyr-v1.1.3,
	 tidyverse-v2.0.0, matrixStats-v1.0.0, rrcovNA-v0.5.0, aggregation-v1.0.1, iq-v1.9.10
	 GMSimpute: see https://github.com/wangshisheng/NAguideR
	 
3. download the expression matrix data from: 

    https://zenodo.org/records/10484253

4. run the scripts using the folder containing the python scripts as working folder:
   
   for DDA data: python DEA_benchmarking_NC_DDA.py platform R_fold
   
   there are 2 required parameters including platform, R_fold
   platform should be "FragPipe" or "Maxquant" if your proteomics data are quantified
   with FragPipe or Maxquant. R_fold specifies the path to your RScript.exe which should
   be the installation path of your R.
   
   for DIA data: python DEA_benchmarking_NC_DIA.py platform R_fold
   
   the parameter platform should be "DIANN" or "spt" where "DIANN" corresponds to
   proteomics data quantified by DIA-NN and "spt" means the proteomics data quantified 
   by Spectronaut.
   
   for TMT data: python DEA_benchmarking_NC_TMT.py platform R_fold
   
   similarly, platform should be "FragPipe" or "Maxquant"
   
5. The ensemble inference should be conducted after finishing above step 4 as we need 
   to use the results from step3 to save time.
   
   firstly, the workflow benchmarking results should be cleared and formated to be acceptable
   by scripts for implementing ensemble inference. We provide the R scripts for clearing the
   benchmarking results and formating the results.
   
   The three R scripts "clear_res_DDA.R", "clear_res_DIA.R" and "clear_res_TMT.R" are for benchmarking
   data clearing.
   
   please ensure the benchmarking results should be saved in "benchmark_res/" folder, our python scripts
   for benchmarking already setting this folder as result storing folder, please don't change it.
   
   Then, the R scripts "combind_analysis_DDA.R", "combind_analysis_DIA.R" and "combind_analysis_TMT.R"
   are used to combine benchmarking results from multiple gold standard datasets.
   
   all these R scripts should be run using the current folder namely "codes_DEA_benchmarking/" as working folder.
   
   after the benchmarking results clearning and combining, now the ensemble inference can be implemented
   with the three python scripts
    
	"DEA_benchmarking_NC_DDA_ensemble.py" for DDA data
	"DEA_benchmarking_NC_DIA_ensemble.py" for DIA data
	"DEA_benchmarking_NC_TMT_ensemble.py" for TMT data
	
	commands:  
	  
	  python DEA_benchmarking_NC_DDA_ensemble.py platform ensemble_type Rscript_pth
	  
	  python DEA_benchmarking_NC_DIA_ensemble.py platform ensemble_type Rscript_pth
	  
	  python DEA_benchmarking_NC_TMT_ensemble.py platform ensemble_type Rscript_pth
	  
	  similarly, for DDA and TMT data, the parameter "platform" should be 'FragPipe' or "Maxquant", DIA data should
	  select "DIANN" or "spt";
	  the parameter "Rscript_pth" specifies the path to your RScript.exe which should
      be the installation path of your R. 
	  
	
	at last, the ensemble inference results can be cleared and combined by the R scripts of:
	
	"clear_ensemble_mean.R" and "combind_analysis_ensemble.R"
	
	the benchmarking results and ensemble inference results are stored under folder "benchmark_res/"
   
   
Any problem in running the scripts please contact:

  Hui Peng: hui.peng@ntu.edu.sg
   
  or:
  
  Wilson Goh: wilsongoh@ntu.edu.sg
   
   