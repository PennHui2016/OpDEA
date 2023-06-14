library(DEP)
library(dplyr)
library(tidyr)
library(stringr)
source('/home/hui/PycharmProjects/DL_proteomics/quant_funcs/preprocessing.R')


args <- commandArgs(trailingOnly = TRUE)
#print(args)
protein_file<-args[1]
peptide_file<-args[2]
design_file<-args[3]
msstats_file<-args[4]
in_type<-args[5]
inten_type<-args[6]
imput<-args[7]
normal<-args[8]
save_fold<-args[9]

if(imput=='blank'){
  imput=''
}

if(normal=='blank'){
  normal=''
}

if(inten_type=='basic'){
  inten_type=''
}

set.seed(123)

make_se2 <- function (proteins_unique, columns, expdesign) 
{
  assertthat::assert_that(is.data.frame(proteins_unique), is.integer(columns), 
                          is.data.frame(expdesign))
  if (any(!c("name", "ID") %in% colnames(proteins_unique))) {
    stop("'name' and/or 'ID' columns are not present in '", 
         deparse(substitute(proteins_unique)), "'.\nRun make_unique() to obtain the required columns", 
         call. = FALSE)
  }
  if (any(!c("label", "condition", "replicate") %in% colnames(expdesign))) {
    stop("'label', 'condition' and/or 'replicate' columns", 
         "are not present in the experimental design", call. = FALSE)
  }
  if (any(!apply(proteins_unique[, columns], 2, is.numeric))) {
    stop("specified 'columns' should be numeric", "\nRun make_se_parse() with the appropriate columns as argument", 
         call. = FALSE)
  }
  if (tibble::is_tibble(proteins_unique)) 
    proteins_unique <- as.data.frame(proteins_unique)
  if (tibble::is_tibble(expdesign)) 
    expdesign <- as.data.frame(expdesign)
  rownames(proteins_unique) <- proteins_unique$name
  raw <- proteins_unique[, columns]
  raw[raw == 0] <- NA
  #raw <- log2(raw)
  # expdesign <- mutate(expdesign, condition = make.names(condition)) %>% 
  #   tidyr::unite(ID, condition, replicate, remove = FALSE)
  expdesign <- mutate(expdesign, condition = make.names(condition), ID = label)
  rownames(expdesign) <- expdesign$ID
  matched <- match(make.names(DEP:::delete_prefix(expdesign$label)), 
                   make.names(DEP:::delete_prefix(colnames(raw))))
  if (any(is.na(matched))) {
    stop("None of the labels in the experimental design match ", 
         "with column names in 'proteins_unique'", "\nRun make_se() with the correct labels in the experimental design", 
         "and/or correct columns specification")
  }
  colnames(raw)[matched] <- expdesign$ID
  raw <- raw[, !is.na(colnames(raw))][rownames(expdesign)]
  row_data <- proteins_unique[, -columns]
  rownames(row_data) <- row_data$name
  se <- SummarizedExperiment:::SummarizedExperiment(assays = as.matrix(raw), colData = expdesign, 
                                                    rowData = row_data)
  return(se)
}

# design_file = '/home/hui/Documents/VM_share/DIA_data/yeast/design_yeast.txt'
#file_name = '/home/hui/Documents/VM_share/DIA_data/yeast/maxquant/txt/proteinGroups.txt'
file_name = protein_file # '/home/hui/Documents/VM_share/DIA_data/yeast/frager/combined_protein.tsv'
#save_fold =  '/home/hui/PycharmProjects/DL_proteomics/benchmark_res/yeast_'
designs = read.table(design_file, quote = "", sep = '\t', header = TRUE)

# in_type = 'maxquant'
# in_type = 'fragpipe'
# inten_type = ''
# imput = ''
# normal = ''
prepro_res = preprocessing_raw(file_name, design_file, in_type, inten_type, imput = imput, normal = normal)

#condition = as.factor(designs$condition)
#design = model.matrix(~0+condition) # fitting without intercept

all_conts<-function(design_data){
  conds<-design_data$condition
  uni_con<-unique(conds)
  conts<-c()
  groups<-vector()
  for (i in 1:(length(uni_con)-1)) {
    for(j in (i+1):length(uni_con)){
      conts<-c(conts, paste0('condition', uni_con[j], '-', 'condition', uni_con[i]))
      groups<-rbind(groups, c(uni_con[i], uni_con[j]))
    }
  }
  return(list(conts=conts, groups=groups))
}

consts<-all_conts(designs)

res_all<-vector()
for (i in 1:length(consts$conts)) {
  groups<-consts$groups
  g1<-groups[i,][1]
  g2<-groups[i,][2]
  
  intens<-prepro_res$normed[,c(grep(g1, designs$condition), grep(g2, designs$condition))]
  condition = as.factor(designs[c(grep(g1, designs$condition), grep(g2, designs$condition)),]$condition)
  design = model.matrix(~0+condition) # fitting without intercept
  intens<-as.data.frame(intens)
  
  intens$na_g1 = apply(intens[,c(grep(g1, colnames(intens)))],1,function(x) sum(is.na(x)))
  intens$na_g2 = apply(intens[,c(grep(g2, colnames(intens)))],1,function(x) sum(is.na(x)))
  # Filter protein table. DEqMS require minimum two values for each group.
  filter_idx = which(intens$na_g1<2 & intens$na_g2<2)
  intens.filter = as.matrix(intens[filter_idx,][,1:(length(intens)-2)])
  
  #row.names(intens)<-proteingroup$Protein#[unlist(test_idx)]
  dataset = cbind(prepro_res$filtered[filter_idx,][,1:7], intens.filter)
  
  experimental.design <- data.frame(designs[c(grep(g1, designs$condition), grep(g2, designs$condition)),])
  colnames(experimental.design)[1]<-'label'
  columns <- c((length(dataset[1,])-length(intens.filter[1,])+1):length(dataset[1,])) # get LFQ column numbers
  
  dataset$name<-dataset[,1]
  dataset$ID<-dataset[,1]
  data.se <- make_se2(dataset, columns, experimental.design)
  #data.norm <- normalize_vsn(data.se)
  #plot_normalization(data.se, data.norm)
  data_diff_all_contrasts <- DEP::test_diff(data.se, type = "all")
  
  dep <- add_rejections(data_diff_all_contrasts, alpha = 0.05, lfc = 0)
  res <- get_results(dep)
  
  res.DEP <-  res %>% 
    gather(statistic,value, contains('_vs_')) %>%
    mutate(statistic = str_replace(statistic, "_vs_",'vs')) %>%
    separate(statistic,c('contrast', 'statistic'),sep = '_') %>%
    spread(statistic,value) %>%
    transmute(protein = ID
              , feature = ID
              ## contrasts are differently specified
              , logFC = -ratio
              , pvalue = p.val
              , qvalue = p.adj
              ) %>%
    as_data_frame
  res.DEP$contrast<-consts$conts[i]
  res_all<-rbind(res_all, res.DEP)
}
#colnames(res_all)[length(res_all[1,])]<-'contrast'
colnames(res_all)[5]<-'adj.pvalue'
res_all<-res_all[order(res_all$pvalue),]
write.table(res_all, paste0(save_fold, 'DEP_', in_type, '_', inten_type, '_', imput, '_', normal, '.csv'), sep = ',', row.names = TRUE, col.names = TRUE)

