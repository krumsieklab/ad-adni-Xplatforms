# Catecholamine dataset - data cleaning and preprocessing
#
#
# by ZW and RB
# last update: 2022-09-04

# clear workspace
rm(list = setdiff(ls(),c("codes.makepath","data.makepath","results.makepath"))) 
# libraries
library(openxlsx) # for excel reading and writing
library(maplet) # maplet 
library(tidyverse) # %>%

# input
data_loc <- 'ADNI/Datasets/ADNI1_Catecholamines/'
anno_loc <- 'ADNI/Datasets/ADNI_metadata/'

# read input
raw_data <- read.csv(data.makepath(paste0(data_loc, 'raw_data/Catecholamines_unblinded_2019-02-25.csv')), header = T)
qc_data <- read.csv2(data.makepath(paste0(data_loc, 'raw_data/Catecholamine_QC_new.csv')), sep = ",",header = T)
fast_stat <- read.csv2(data.makepath(paste0(anno_loc, 'FastingStatusADNI1.txt')), sep = "\t", as.is = T)
meds <- read.csv2(data.makepath(paste0(anno_loc, 'MedicationsADNI1GO2.txt')), sep = "\t", as.is = T)
# output files
output_mt_ready <- '2022-09-03_catecholamines.xlsx'
output_pp_no_med <- '2022-09-03_catecholamines_preprocessing.xlsx'
output_pp_med <-  '2022-09-03_catecholamines_preprocessing_medcor.xlsx'

# sample measurements cleaning
  # remove space in samples names of raw data
raw_ids <- raw_data %>% mutate(Sample.Name=gsub('\\s*', '', Sample.Name)) %>%
  # add batch information 
  mutate(Batch_id = gsub('*.[A-Z]','',Sample.Name),
         Batch_id=as.numeric(as.matrix(Batch_id))) %>% 
  # remove uninformative samples
  select(-c('X', 'Pt_ID')) %>%
  # order columns
  select(RID, Sample.Name, Batch_id) %>%
  {.}
# function to remove less than sign from LOD values
rem_sp_char <- function(x){ gsub('<', '', x) } 
# measurements of samples, remove special characters and make numeric
raw_data %<>% select(-c(RID, X, Pt_ID, Sample.Name)) %>% 
 mutate_all(rem_sp_char)%>% 
  mutate_all(as.matrix) %>%
  mutate_all(as.numeric)
  # bind ids and measurements
raw_data <- bind_cols(raw_ids, raw_data)

# qc measurements cleaning
  # identify qc - nist samples
qc_data <-  qc_data[grep(pattern = 'NIST', qc_data$Sample.Name),]
  # separate ids of the qc samples 
qc_ids <- qc_data %>% mutate(RID = as.integer(rep('99999',nrow(qc_data)))) %>%
  # order columns
  select(RID, Sample.Name, Batch_id)
  # measurements of qc samples
qc_data %<>% select(-c(Sample.Name, Batch_id)) %>% mutate_all(as.matrix) %>%
  mutate_all(as.numeric)
  # bind ids and measurements
qc_data <- bind_cols(qc_ids, qc_data)

# add qc data to raw_data
raw_data <- bind_rows(raw_data, qc_data)

# define columns that are sample information and to be discarded from measurements sheet and metabolite information sheet
samp_info_cols <- c('RID','Sample.Name','Batch_id')

# colData
sample_info <- raw_data %>% select(all_of(samp_info_cols)) %>% 
  left_join(fast_stat, by='RID') %>% 
  left_join(meds, by='RID')
# rowData
met_info <- data.frame(col_ID=raw_data %>% select(-all_of(samp_info_cols)) %>% names(), 
                       name=raw_data %>% select(-all_of(samp_info_cols)) %>% names())

# assay
raw_data <- raw_data %>% select(-samp_info_cols) %>% data.frame()%>%
  mutate_all(as.matrix) %>% mutate_all(as.numeric) %>% t()

# create summarized Experiment
D <- SummarizedExperiment::SummarizedExperiment(assays = raw_data,
                                                colData = sample_info, rowData = met_info)

mt_write_se_xls(D, file=output_mt_ready)
# qc
D %<>%
  # mist based correction
  mt_pre_reference_correction(qc_samples= RID == 99999,
                               plate_col= 'Batch_id') %>%
  # add cv from duplicated ids to rowData
  mt_pre_cv(qc_samples = (RID%in%RID[duplicated(RID)]), 
               out_col = "Dup_cv", replicates=T, 
               id_col='RID') %>%
  # icc for duplicates
  mt_pre_icc(qc_samples = (RID%in%RID[duplicated(RID)]), out_col = "Dup_icc", id_col='RID') %>%
  mt_anno_mutate(anno_type ='samples', col_name='Sample_ID', term=paste0('sample_', 1:length(RID))) %>%
  {.}

# preprocessing 
D %<>%
  mt_reporting_heading(heading = "Preprocessing", lvl = 1) %>%
  # Filter out non-fasting samples 
  mt_modify_filter_samples(BIFAST == "Yes") %>%
  mt_reporting_heading(heading = "Missingness", lvl = 2) %>%
  mt_pre_zero_to_na() %>%
  mt_reporting_heading(heading = "Before filtering", lvl = 3) %>%
  # plot pre-filtering missingness distribution
  mt_plots_missingness(feat_max=0.4) %>%
  # Filter out metabolites with more than 40% missingness
  mt_pre_filter_missingness(feat_max = 0.4) %>%
  # Filter out samples with more than 40% missingness
  mt_pre_filter_missingness(samp_max = 0.4) %>%
  mt_reporting_heading(heading  = "After filtering", lvl = 3) %>%
  # plot pre-filtering missingness distribution
  mt_plots_missingness(feat_max=0.4) %>%
  # add missingness percentage as annotation to metabolites
  mt_anno_missingness(anno_type = "features", out_col = "missing") %>%
  # Log2 transformation
  mt_pre_trans_log() %>%
  mt_reporting_heading(heading ="Imputation", lvl = 2) %>%
  # plot pre-imputation sample boxplot
  mt_plots_sample_boxplot(title = "Before imputation") %>%
  # impute missing values
  mt_pre_impute_knn() %>%
  # plot post-imputation sample boxplot
  mt_plots_sample_boxplot(title = "After imputation") %>%
  # Average-combine duplicate samples
  mt_modify_avg_samples(group_col =  "RID") %>%
  
  {.}
# outlier detection and adjustments
# identify metabolite level outliers -> set to NA -> impute
D %<>%   # sample outlier detection
  #mt_pre_outlier_detection_mahalanobis(pval=0.01) %>%  The data matrix is not full-rank, Mahalanobis cannot be computed
  mt_pre_outlier_lof(seq_k = c(5, 10, 20, 30, 40, 50)) %>% # local outlier factor
  mt_pre_outlier_to_na(use_quant=TRUE, quant_thresh =0.025) %>% 
  mt_pre_impute_knn() %>%
  {.}

# data overview plots
D %<>% mt_reporting_heading(heading  = "Global Statistics", lvl = 2) %>%
  # plot PCA
  mt_plots_pca(scale_data = T, title = "PCA", size=2.5, ggadd=scale_size_identity()) %>%
  # plot UMAP
  mt_plots_umap(scale_data = T, title = "UMAP", size=2.5, ggadd=scale_size_identity()) %>%
  # plot heatmap
  mt_plots_heatmap(scale_data = T, fontsize = 5) %>%
  {.}

mt_write_se_xls(D, file=output_pp_no_med)

# medication correction [non maplet]
all_cols <- D %>% colData %>% as_tibble() %>% names()
med_cols <- all_cols[grep("Med", all_cols)] # column numbers of all meds
med_cols <- med_cols[which(med_cols%in%meds_to_exclude==F)] # meds to correct

# medication correction
Dmc <- D %>% mt_pre_confounding_correction_stepaic(cols_to_correct = med_cols, cols_to_exclude = meds_to_exclude)

mt_write_se_xls(Dmc, file=output_pp_med)