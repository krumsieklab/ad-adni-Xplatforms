# Catecholamine dataset - data cleaning and preprocessing
#
#
# by RB, ZW
# last update: 2022-09-03


# clear workspace
rm(list = setdiff(ls(),c("codes.makepath","data.makepath","results.makepath"))) 
# libraries
library(maplet) #
library(tidyverse) # %>%

# helper function
read_list_from_file <- function(file) {
  # check if file is Excel format
  if (!is.na(readxl::excel_format(path = file))) {
    # load first column from first sheet
    vec <- as.vector(read_excel(path=file, sheet = 1)[,1,drop=T])
  } else {
    # must be a text file
    vec <- readLines(file)
  }
  # remove comments
  vec <- gsub("#.*", "", vec)
  # trim whitespaces
  vec <- trimws(vec)
  # remove empty entries
  vec <- vec[nchar(vec)>0]
  # return
  vec
}

#### Read input files ----

# input
data_loc <- 'ADNI/Datasets/ADNI1_Catecholamines/'
anno_loc <- 'ADNI/Datasets/ADNI_metadata/'

# read input
raw_data <- read.csv(data.makepath(paste0(data_loc, 'raw_data/Catecholamines_unblinded_2019-02-25.csv')), header = T)
qc_data <- read.csv2(data.makepath(paste0(data_loc, 'raw_data/Catecholamine_QC_new.csv')), sep = ",",header = T)
fasting <- read.csv2(data.makepath(paste0(anno_loc, 'FastingStatusADNI1.txt')), sep = "\t", as.is = T)

# output files
output_pp_no_med <- 'catecholamines_preprocessing.xlsx'
output_pp_med <-  'catecholamines_preprocessing_medcor.xlsx'

#### Data cleaning ----

# define columns to be discarded
cols2discard <- c('X','RID','Pt_ID','Sample.Name','Batch_id')
qc_data1 <- qc_data[qc_data %>% select(Sample.Name) %>% as.matrix() %>% grep(pattern = 'NIST'),]
qc_data <- qc_data1 %>% mutate(RID = as.integer(rep('99999',nrow(qc_data1))))
raw_data1 <- raw_data %>% as.matrix() %>% gsub(' E','E',.) %>% gsub('^([1-9][A-Z])([1-9])$','\\10\\2',.) %>% 
  as.data.frame(stringsAsFactors = F) %>% arrange(Sample.Name) %>% mutate(RID = as.integer(RID)) %>%
  mutate(Batch_id = c(rep(1:10,each = 76),rep(11,32)))

raw_data <- bind_rows(raw_data1, qc_data %>% select(all_of(names(raw_data1 %>% select(-c(X, Pt_ID))))))

#sample_info
# colData
sample_info <- raw_data %>% select(all_of(cols2discard)) %>% 
  left_join(fast_stat, by='RID') %>% 
  left_join(meds, by='RID')

# rowData
met_info <- data.frame(col_ID=raw_data %>% select(-cols2discard) %>% names(), 
                       name=raw_data %>% select(-cols2discard) %>% names())

# assay
raw_data <- raw_data %>% select(-cols2discard) %>% data.frame()%>%
  mutate_all(as.matrix) %>% mutate_all(as.numeric) %>% t()


#### Write output files ----

# Summarazed Experiment
D <- SummarizedExperiment::SummarizedExperiment(assays = raw_data,
                                                colData = sample_info, rowData = met_info)

### QC ----

# MT pipeline (starting with %<>%, since D from se_qc file is already a SummarizedExperiment)
D %<>%
  # Nist based correction
  mt_pre_nist_based_correction(qc_samples= RID == 99999,
                               plate_col= 'Batch_id') %>%
  # Add CV from duplicated ids to rowData
  mt_modify_cv(qc_samples = (RID%in%RID[duplicated(RID)]), 
               col_lab = "Dup_cv", replicates=T, 
               id_col='RID') %>%
  # ICC for duplicates
  mt_modify_icc(qc_samples = (RID%in%RID[duplicated(RID)]), col_lab = "Dup_icc", id_col='RID') %>%
  mt_modify_mutate(anno_type = 'samples', col_name='Sample_ID', term=paste0('sample_', 1:length(RID))) %>%
  {.}

### Preprocess ----
D %<>%  
  # Headings
  mt_reporting_heading(strtitle = "Preprocessing", lvl = 1) %>%
  mt_reporting_heading(strtitle = "Filtering", lvl = 2) %>%
  
  # Filter out non-fasting samples 
  mt_modify_filter_samples(BIFAST == "Yes") %>%
  
  # Filter out missing data for both metabolites and samples using the wrapper function
  mtw_missingness(plot_options = list(met_max = 0.25,samp_max = 1),
                  filter_options = list(met_max = 0.25,sample_max = 1)) %>%
  
  # Filter out metabolites with cv > 25% or ICC < 65%
  #mt_modify_filter_metabolites(Dup_cv <0.25 || Dup_icc >0.65) %>%
  
  # Preprocess using the wrapper function
  mt_reporting_heading(strtitle = "Log transformation", lvl = 2) %>%
  
  # Log2 transformation
  mt_pre_trans_log(base = 2) %>%
  
  mt_reporting_heading(strtitle = "Imputation", lvl = 2) %>%
  
  # pre-imputation sample boxplot
  mt_plots_sampleboxplot(plottitle = "Before imputation") %>%
  
  mt_pre_impute_knn() %>% 
  
  # post-imputation sample boxplot
  mt_plots_sampleboxplot(plottitle = "After imputation") %>%
  
  # Average-combine duplicate samples
  mt_modify_averagesample(group_by = "RID") %>%
  
  # Add dataset info
  mt_reporting_heading(strtitle = "Dataset info", lvl = 2) %>%
  
  mt_logging_datasetinfo()
{.}

### Global statistics ----
D %<>%  
  
  # Global stats as data overview
  mt_reporting_heading(strtitle = "Global Statistics", lvl = 1) %>%
  
  mtw_global_stats(pca_options = list(scaledata = T, size = 2, ggadd=quote(scale_size_identity())),
                   umap_options = list(scaledata = T, size = 2, ggadd=quote(scale_size_identity())),
                   heatmap_options = list(scaledata = T, fontsize = 5)) %>%
  
  # Global outlier detection
  mt_pre_outlier(method='mahalanobis', pval=0.01) %>%
  
  mt_modify_filter_samples(outlier_mahalanobis != TRUE) %>%
  
  mtw_global_stats(pca_options = list(scaledata = T, size = 2, ggadd=quote(scale_size_identity())),
                   umap_options = list(scaledata = T, size = 2, ggadd=quote(scale_size_identity())),
                   do_pheatmap = F) %>%
  {.}