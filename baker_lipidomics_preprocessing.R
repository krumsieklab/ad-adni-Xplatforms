# baker lipidomics dataset - data cleaning and preprocessing
# 
# based on scripts written by: Matthias Arnold, Tyler Massaro
# Their comment: Data came batch normalized, so there is no batch removal included
# Adapted by RB, ZW 
# last update: 2022-09-04

# clear workspace
rm(list = setdiff(ls(),c("codes.makepath","data.makepath","results.makepath")))

# libraries
library(openxlsx) # for excel reading and writing
library(maplet) # MT 
library(tidyverse) # %>%

# helper function ----
combine_cvs <- function(D0){
  # select pqc_cv columns in the row data and create a vector
  pqc <- D0 %>% rowData() %>% data.frame() %>% 
    select(starts_with('PQC_cv')) %>% data.frame() %>%
    cbind(replicate(length(grep("PQC", D0$QCID)), .$'PQC_cv'))
  # select aqc_cv columns in the row data and create a vector
  aqc <- D0 %>% rowData() %>% data.frame() %>% 
    select(starts_with('AQC_cv')) %>% data.frame() %>%
    cbind(replicate(length(grep("A_QC", D0$QCID)), .$'AQC_cv'))
  # identify duplicates
  dups <- D0 %>% rowData() %>% data.frame() %>% 
    select(starts_with('Dup_cv')) %>% data.frame() %>%
    cbind(replicate(length((D0$QCID=='' & D0$RID%in%D0$RID[duplicated(D0$RID)])), .$'Dup_cv'))
  # bind together
  cvs <- bind_cols(pqc, aqc, dups) %>% rowMeans(na.rm = T)
  # add a mean
  rowData(D0)$Combine_cv <- cvs
  return(D0)
}

# input ----
data_loc <- 'ADNI/Datasets/ADNI1_Baker_Lipidomics/'
anno_loc <- '/ADNI/Datasets/ADNI_metadata/'
data_file <- data.makepath(paste0(data_loc,'raw_data/Baker_LIPIDOMICSDATABASE_03_05_19.csv'))
met_file <- data.makepath(paste0(data_loc,'raw_data/Baker_LIPIDOMICSDATABASE_DICT.csv'))
anno_file <- data.makepath(paste0(anno_loc,'/adni1go2_phenotypes_covariates.xlsx'))
fasting_file <- data.makepath(paste0(anno_loc,'FastingStatusADNI1.txt'))
med_file <- data.makepath(paste0(anno_loc,'MedicationsADNI1GO2.txt'))
meds_to_exclude <- c("Med.Anticholinesterases","Med.NMDAAntag") # meds to exclude

# output files ----
output_mt_ready <- '2022-09-03_baker_lipidomics.xlsx'
output_pp_no_med <- '2022-09-03_baker_lipidomics_preprocessing.xlsx'
output_pp_med <-  '2022-09-03_baker_lipidomics_preprocessing_medcor.xlsx'
output_html <- '2022-09-03_baker_lipidomics_preprocessing_medcor.html'

# read input ----
raw_data <- read.csv(data_file, header = T)
fast_stat <- read.csv2(fasting_file, sep = "\t")
meds <- read.csv2(med_file, sep = "\t")

# define columns that are sample information and to be discarded from measurements sheet and metabolite information sheet
samp_info_cols <- c("RID", "VIALID", "QCID", "ANALYSIS_START_DATE", "update_stamp")

# se ----
# colData
sample_info <- raw_data %>% select(all_of(samp_info_cols)) %>% 
  left_join(fast_stat, by='RID') %>% # add fasting info
  left_join(meds, by='RID') # add medication info
# assay
raw_data <- raw_data %>% select(-samp_info_cols) %>% t()
# rowData
met_info <- read.csv(met_file, header = T) %>%
  select(FLDNAME, NOTES) %>% dplyr::rename(col_ID=FLDNAME, name=NOTES) %>%
  filter(!col_ID%in%samp_info_cols)

D <- SummarizedExperiment::SummarizedExperiment(assays = raw_data,
                                                colData = sample_info, rowData = met_info)

# qc ----
D %<>%
  mt_pre_cv(qc_samples = QCID == "PQC", out_col = "PQC_cv")%>%
  mt_pre_cv(qc_samples = QCID == "A_QC", out_col = "AQC_cv") %>%
  mt_pre_cv(qc_samples = (QCID=='' & RID%in%RID[duplicated(RID)]), 
            out_col = "Dup_cv", replicates=T, id_col='RID') %>%
  # combine 3 CVs
  combine_cvs() %>% 
  # ICC for duplicates
  mt_pre_icc(qc_samples = (QCID=='' & RID%in%RID[duplicated(RID)]), out_col = "Dup_icc", 
             id_col='RID') %>%
  # select non qc samples
  mt_modify_filter_samples(QCID == '') %>% 
  mt_anno_mutate(anno_type = 'samples', col_name='Sample_ID', term=paste0('sample_', 1:length(RID)))

mt_write_se_xls(D, file=output_mt_ready)

# preprocessing ----
D %<>%
  mt_reporting_heading(heading= "Preprocessing", lvl = 1) %>%
  mt_reporting_heading(heading = "Filtering", lvl = 2) %>%
  # filter out non-fasting samples (as well as QC samples)
  mt_modify_filter_samples(BIFAST == "Yes") %>%
 # plot missingness distribution
  mt_plots_missingness(feat_max=0.4) %>%
  # filter metabolites with more than 40% missing values per group
  mt_pre_filter_missingness(feat_max = 0.4) %>%
  # plot missingness distribution after filtering
  mt_plots_missingness(feat_max=0.4) %>%
  # add missingness percentage as annotation to samples (remaining missing)
  mt_anno_missingness(anno_type = "samples", out_col = "missing") %>%
  # add missingness percentage as annotation to metabolites
  mt_anno_missingness(anno_type = "features", out_col = "missing") %>%
  # filter out metabolites with cv > 25% or ICC < 65%
  mt_modify_filter_features(Combine_cv < 0.25 || Dup_icc > 0.65) %>%
  # normalization
  mt_reporting_heading(heading = "Normalization", lvl = 2) %>%
  # pre-normalization sample boxplot
  mt_plots_sample_boxplot(title = "Before normalization") %>%
  # normalize abundances using probabilistic quotient
  mt_pre_norm_quot(feat_max = 0.4) %>%
  # post-normalization sample boxplot
  mt_plots_sample_boxplot(title = "After normalization") %>%
  # dilution plot showing dilution factors from quotient normalization
  mt_plots_dilution_factor(boxpl = T, in_col = 'BIFAST')  %>%
  # log2 transformation
  mt_pre_trans_log() %>%
  mt_reporting_heading(heading = "Imputation", lvl = 2) %>%
  # pre-imputation sample boxplot
  mt_plots_sample_boxplot(title = "Before imputation") %>%
  mt_pre_impute_knn() %>% 
  # post-imputation sample boxplot
  mt_plots_sample_boxplot(title = "After imputation") %>%
  # average-combine duplicate samples
  mt_modify_avg_samples(group_col = "RID") %>%
  {.}

# outlier detection and adjustments ----
D %<>%   # sample outlier detection
  # sample outlier detection
  # mt_pre_outlier_detection_leverage(thresh = 4) %>%
  #mt_pre_outlier(method='leverage', pval=0.01) %>%
  #mt_modify_filter_samples(outlier_leverage != TRUE) %>%
  # identify metabolite level outliers -> set to NA -> impute
  #mt_pre_outlier_detection_mahalanobis(pval=0.01) %>%  The data matrix is not full-rank, Mahalanobis cannot be computed
  #mt_pre_outlier_lof(seq_k = c(5, 10, 20, 30, 40, 50)) %>% # local outlier factor
  mt_pre_outlier_to_na(use_quant=TRUE, quant_thresh =0.025) %>% 
  mt_pre_impute_knn() %>%
  {.}

# data overview plots ----
D %<>% mt_reporting_heading(heading  = "Global Statistics", lvl = 2) %>%
  # plot PCA
  mt_plots_pca(scale_data = T, title = "PCA", size=2.5, ggadd=scale_size_identity()) %>%
  # plot UMAP
  mt_plots_umap(scale_data = T, title = "UMAP", size=2.5, ggadd=scale_size_identity()) %>%
  # plot heatmap
  mt_plots_heatmap(scale_data = T, fontsize = 5) %>%
  {.}

mt_write_se_xls(D, file=output_pp_no_med)

# medication correction ----
 # collect all column names of sample information
all_cols <- D %>% colData %>% as_tibble() %>% names()
 # grep columns with medication infor
med_cols <- all_cols[grep("Med", all_cols)] # column numbers of all meds
 # remove the AD samples from correction
med_cols <- med_cols[which(med_cols%in%meds_to_exclude==F)] # meds to correct
 # correction
Dmc <- D %>% mt_pre_confounding_correction_stepaic(cols_to_correct = med_cols, 
                                                   cols_to_exclude = meds_to_exclude, n_cores = 40)

mt_write_se_xls(Dmc, file=output_pp_med)
Dmc %>% mt_reporting_html(file=output_html)
# done ----