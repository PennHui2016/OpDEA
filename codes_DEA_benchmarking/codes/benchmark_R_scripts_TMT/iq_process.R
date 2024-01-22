library(iq)

# process_long_format("D:/software/wechat/WeChat Files/wxid_ip7m0cn371qg22/FileStorage/File/2023-09/20230918_102252_50ng+15p625ng_hela_yeast_mix_Report.tsv",
#                     output_filename = "E:/MS_data/PXD034709/SearchResult-Benchmark-TIMS-Spectronaut-UniversalLib/iq-MaxLFQ.tsv", 
#                     sample_id  = "R.FileName",
#                     primary_id = "PG.ProteinGroups",
#                     secondary_id = c("EG.Library", "FG.Id", "FG.Charge", "F.FrgIon", "F.Charge", "F.FrgLossType"),
#                     intensity_col = "F.PeakArea",
#                     annotation_col = c("PG.Genes", "PG.ProteinNames", "PG.FastaFiles"),
#                     filter_string_equal = c("F.ExcludedFromQuantification" = "False"),
#                     filter_double_less = c("PG.Qvalue" = "0.01", "EG.Qvalue" = "0.01"),
#                     log2_intensity_cutoff = 0)

#modifiied
source('D:/data/benchmark/iq-master/R/iq-fast_MOD.R')
process_long_format("E:/MS_data/PXD036134/DIANN/report.tsv",
                    output_filename = "E:/MS_data/PXD036134/DIANN/test.tsv",
                    annotation_col = c("Protein.Names", "Genes"),
                    filter_double_less = c("Global.Q.Value" = "0.01", "Global.PG.Q.Value" = "0.01"),
                    method='topN', N=3)

#origianl
# library(iq)
# process_long_format("E:/MS_data/PXD036134/DIANN/report.tsv", 
#                     output_filename = "E:/MS_data/PXD036134/DIANN/test.tsv", 
#                     annotation_col = c("Protein.Names", "Genes"),
#                     filter_double_less = c("Global.Q.Value" = "0.01", "Global.PG.Q.Value" = "0.01"))

# raw <- read.delim("E:/MS_data/PXD036134/DIANN/report.tsv")
# 
# selected <-
#   !is.na(raw$Global.Q.Value) & (raw$Global.Q.Value < 0.01) &
#   !is.na(raw$Global.PG.Q.Value) & (raw$Global.PG.Q.Value < 0.01)
# 
# raw <- raw[selected,]
# norm_data <- iq::preprocess(raw)
data("spikeins")
norm_data <- iq::preprocess(spikeins, median_normalization = FALSE, pdf_out = NULL)
protein_list <- iq::create_protein_list(norm_data)
result <- iq::create_protein_table(protein_list, method = 'topN', N=3)


process_long_format("F:/PXD034709/PXD034709_SPN/20231013_092322_20210525_Report.tsv",
                    output_filename = "D:/data/benchmark/data/test.tsv", 
                    sample_id  = "R.FileName",
                    primary_id = "PG.ProteinGroups",
                    secondary_id = c("EG.Library", "FG.Id", "FG.Charge", "F.FrgIon", "F.Charge", "F.FrgLossType"),
                    intensity_col = "F.PeakArea",
                    annotation_col = c("PG.Genes", "PG.ProteinNames", "PG.FastaFiles"),
                    filter_string_equal = c("F.ExcludedFromQuantification" = "False"),
                    filter_double_less = c("PG.Qvalue" = "0.01", "EG.Qvalue" = "0.01"),
                    log2_intensity_cutoff = 0)
