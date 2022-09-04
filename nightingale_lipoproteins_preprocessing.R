# nightingale lipoproteins dataset - data cleaning and preprocessing
#
#
# by RB
# last update: 2022-09-04

# clear workspace
rm(list = setdiff(ls(),c("codes.makepath","data.makepath","results.makepath")))
# libraries
library(openxlsx) # for excel reading and writing
library(maplet) # maplet 
library(tidyverse) # %>%

# input ----
data_loc <- '/ADNI/Datasets/ADNI1_GO2_Baseline_Nightingale/'
anno_loc <- '/ADNI/Datasets/ADNI_metadata/'
data_file <- data.makepath(paste0(data_loc, '/raw_data/16033-30-Sep-2020-Results.xlsx'))
metadata <- data.makepath(paste0(data_loc,'/meta_data/Nightingale_unblinding_sampleid_RID.xlsx'))
anno_file <- data.makepath(paste0(anno_loc,'/adni1go2_phenotypes_covariates.xlsx'))
batch_file <- data.makepath(paste0(anno_loc,'/ADNI_Batch.xlsx'))
meds_to_exclude <- c("Med.Anticholinesterases","Med.NMDAAntag") # meds to exclude

# output files ----
output_mt_ready <- '2022-09-03_ADNI_Nightingale_Baseline.xlsx'
output_pp_no_med <- '2022-09-03_ADNI_Nightingale_Baseline_preprocessed.xlsx'
output_pp_med <- '2022-09-03_ADNI_Nightingale_Baseline_preprocessed_medcor.xlsx'
output_html <- '2022-09-03_ADNI_Nightingale_Baseline_preprocessed_medcor.html'

# summarized experiment ----
D <- mt_load_nightingale (file=data_file, 
                          format_type = 'multiple_sheets_v1') %>%
  # print infos about dataset
  mt_reporting_data()%>%
  mt_anno_mutate(anno_type = 'samples', col_name = 'Sample_id', term=gsub('-', '_', Sample_id)) %>%
  # load sample annotations
  mt_anno_xls(file=metadata,sheet=1, 
              anno_type="samples", anno_id_col="sampleid", data_id_col="Sample_id") %>%
  # load batch annotations
  mt_anno_xls(file=batch_file,sheet=1, 
              anno_type="samples", anno_id_col="RID", data_id_col="RID") %>%
  {.}

mt_write_se_xls(D, file=output_mt_ready)

# qc ----
D <- D %>% mt_anno_mutate(anno_type = "samples", col_name = "Low_protein",
                          term = case_when(is.na(Low_protein) ~ '1', TRUE ~ Low_protein)) %>%
  mt_modify_filter_samples(Low_protein=='0') %>%
  {.}

# preprocessing ----
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

# outlier detection and adjustments ----
D %<>%  
 #sample outlier detection
 #mt_pre_outlier_detection_mahalanobis(pval=0.01) %>%  The data matrix is not full-rank, Mahalanobis cannot be computed
 #mt_pre_outlier_lof(seq_k = c(5, 10, 20, 30, 40, 50)) %>% # local outlier factor
 #identify metabolite level outliers -> set to NA -> impute
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