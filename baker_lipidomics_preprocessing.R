# Baker lipidomics dataset - data cleaning and preprocessing
#
#
# by RB, ZW
# last update: 2022-09-03


### Setup ----

# clear workspace
rm(list = setdiff(ls(),c("codes.makepath","data.makepath","results.makepath")))
# libraries
library(MetaboTools) # MT 
library(tidyverse) # %>%


### File paths ----

# Input files
data_loc <- 'ADNI/Datasets/ADNI1_Targeted_Lipidomics/'
anno_loc <- 'ADNI/Datasets/ADNI_metadata/'
se_qc <- data.makepath(paste0(data_loc, 'processed_data/ADNI1_Meikle_qc.rds'))
file_anno <- data.makepath(paste0(anno_loc, 'adni1go2_phenotypes_covariates.xlsx'))
file_outcome <- data.makepath(paste0(anno_loc, 'adni_outcomes.csv'))
file_meds <- data.makepath(paste0(anno_loc, 'adni_medications.csv'))

# Output files
se_output_pp_no_med <- 'Baker_Lipidomics_PP_no_med.xlsx'
se_output_pp_med <-  'Baker_Lipidomics_PP_med.xlsx'

### Load outcome and medication variable lists ----
med_cols <- read.csv(file_meds,header = T, as.is = T, stringsAsFactors = F) %>% as.matrix() %>% as.numeric()


### Start MT pipeline ----

# Load the SummarizedExperiment (=data)
load(se_qc)

### QC ----

# MT pipeline (starting with %<>%, since D from se_qc file is already a SummarizedExperiment)
D %<>%
  # CV for QC samples
  mt_modify_cv(qc_samples = QCID == "PQC", col_lab = "PQC_cv")%>%
  mt_modify_cv(qc_samples = QCID == "A_QC", col_lab = "AQC_cv") %>%
  mt_modify_cv(qc_samples = (QCID=='' & RID%in%RID[duplicated(RID)]), 
               col_lab = "Dup_cv", replicates=T, 
               id_col='RID') %>%
  # ICC for duplicates
  mt_modify_icc(qc_samples = (QCID=='' & RID%in%RID[duplicated(RID)]), col_lab = "Dup_icc", 
                id_col='RID') %>%
  # select non qc samples
  mt_modify_filter_samples(QCID == '') %>% 
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
  mtw_missingness(plot_options = list(met_max = 0.4,samp_max = 1),
                  filter_options = list(met_max = 0.4,sample_max = 1)) %>%
  
  # Filter out metabolites with cv greater than 25%
  mt_modify_filter_metabolites(PQC_cv <0.25 || AQC_cv <0.25) %>%
  
  # Normalization
  mt_reporting_heading(strtitle = "Normalization", lvl = 2) %>%
  
  # pre-normalization sample boxplot
  mt_plots_sampleboxplot(plottitle = "Before normalization") %>%
  
  # normalize abundances using probabilistic quotient
  mt_pre_norm_quot(met_max = 0.4) %>%
  
  # post-normalization sample boxplot
  mt_plots_sampleboxplot(plottitle = "After normalization") %>%
  
  # dilution plot showing dilution factors from quotient normalization
  mt_plots_qc_dilutionplot(boxpl = T,comp = 'BIFAST')  %>%
  
  # Log2 transformation
  mt_pre_trans_log() %>%
  
  mt_pre_trans_scale() %>%
  
  # WINSORIZING & RE-IMPUTATION : Local outlier correction
  mt_pre_outliercorrection(threshold = 4) %>%
  
  mt_reporting_heading(strtitle = "Imputation", lvl = 2) %>%
  
  # pre-imputation sample boxplot
  mt_plots_sample_boxplot(plottitle = "Before imputation") %>%
  
  mt_pre_impute_knn() %>% 
  
  # post-imputation sample boxplot
  mt_plots_sampleboxplot(plottitle = "After imputation") %>%
  
  # Average-combine duplicate samples
  mt_modify_averagesample(group_by = "RID") %>%
  {.}
  
### Global statistics ----
D %<>%  
  # Global stats as data overview
  mt_reporting_heading(heading  = "Global Statistics", lvl = 2) %>%
  # plot PCA
  mt_plots_pca(scale_data = T, title = "PCA", size=2.5, ggadd=scale_size_identity()) %>%
  # plot UMAP
  mt_plots_umap(scale_data = T, title = "UMAP", size=2.5, ggadd=scale_size_identity()) %>%
  # plot heatmap
  mt_plots_heatmap(scale_data = T, fontsize = 5)

  
  # Global outlier detection
  mt_pre_outlier(method='mahalanobis', pval=0.01) %>%
  
  {.}

mt_write_se_xls(D, file=se_output_pp_no_med)

# Medication correction ----
all_cols <- D %>% colData %>% as_tibble() %>% names()
meds_to_exclude <- c("Med.Anticholinesterases","Med.NMDAAntag") # meds to exclude
med_cols <- all_cols[grep("Med", all_cols)] # column numbers of all meds
med_cols <- med_cols[which(med_cols%in%meds_to_exclude==F)] # meds to correct

# medication correction
Dmc <- D %>% mt_pre_confounding_correction_stepaic(cols_to_correct = med_cols, cols_to_exclude = meds_to_exclude)

mt_write_se_xls(Dmc, file=se_output_pp_med)
