# Nightingale lipoproteins dataset - data cleaning and preprocessing
#
#
# by RB
# last update: 2022-09-03

# Set up ----

# clear workspace
rm(list = setdiff(ls(),c("codes.makepath","data.makepath","results.makepath")))
# libraries
library(openxlsx) # for excel reading and writing
library(maplet) # maplet 
library(tidyverse) # %>%

# Input files ----
data_loc <- '/ADNI/Datasets/ADNI1_GO2_Baseline_Nightingale/'
data_file <- 'raw_data/16033-30-Sep-2020-Results.xlsx'
metadata <- 'meta_data/Nightingale_unblinding_sampleid_RID.xlsx'
anno_file <- '/ADNI/Datasets/ADNI_metadata/adni1go2_phenotypes_covariates.xlsx'
outcome_file <- '/ADNI/Datasets/ADNI_metadata/traits.csv'
batch_file <- '/ADNI/Datasets/ADNI_metadata/ADNI_Batch.xlsx'

# Preparing fasting and medication annotaitons [to be done once]
# fasting_file_1 <- '/ADNI/Datasets/ADNI_metadata/FastingStatusADNI1.txt'
# fasting_file_2 <- '/ADNI/Datasets/ADNI_metadata/FastingStatusADNIGO2.txt'
# med_file <- '/ADNI/Datasets/ADNI_metadata/MedicationsADNI1GO2.txt'
# sample_info <- read.xlsx(data.makepath(paste0(data_loc, metadata)))
# fast_stat_1 <- read.csv2(data.makepath(fasting_file_1), sep = "\t", as.is = T)
# fast_stat_2 <- read.csv2(data.makepath(fasting_file_2), sep = "\t", as.is = T) %>%
#   group_by(RID) %>% summarise(BIFAST = BIFAST[which(!is.na(BIFAST))[1]]) %>% ungroup()
# fast_stat <- bind_rows(fast_stat_1, fast_stat_2) %>% mutate(RID=as.numeric(RID))
# meds <- read.csv2(data.makepath(med_file), sep = "\t", as.is = T) %>% mutate(RID=as.numeric(RID))
# sample_info  %<>% left_join(fast_stat, by='RID') %>% unique() %>% left_join(meds, by='RID') %>% unique()
# write.xlsx(sample_info, data.makepath(paste0(data_loc, metadata)))

# QC and PP Biomarker data ----
D <- mt_load_nightingale (file=data.makepath(paste0(data_loc, data_file)), 
                          format_type = 'multiple_sheets_v1') %>%
  # print infos about dataset
  mt_reporting_data()%>%
  mt_anno_mutate(anno_type = 'samples', col_name = 'Sample_id', term=gsub('-', '_', Sample_id)) %>%
  # load sample annotations
  mt_anno_xls(file=data.makepath(paste0(data_loc, metadata)),sheet=1, 
                    anno_type="samples", anno_id_col="sampleid", data_id_col="Sample_id") %>%
  # load batch annotations
  mt_anno_xls(file=data.makepath(batch_file),sheet=1, 
                    anno_type="samples", anno_id_col="RID", data_id_col="RID")
  
# QC 
D <- D %>% mt_anno_mutate(anno_type = "samples", col_name = "Low_protein",
                            term = case_when(is.na(Low_protein) ~ '1', TRUE ~ Low_protein)) %>%
  mt_modify_filter_samples(Low_protein=='0')
# Preprocessing 
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
  # Global stats as data overview
  mt_reporting_heading(heading  = "Global Statistics", lvl = 2) %>%
  # plot PCA
  mt_plots_pca(scale_data = T, title = "PCA", size=2.5, ggadd=scale_size_identity()) %>%
  # plot UMAP
  mt_plots_umap(scale_data = T, title = "UMAP", size=2.5, ggadd=scale_size_identity()) %>%
  # plot heatmap
  mt_plots_heatmap(scale_data = T, fontsize = 5)

# identify metabolite level outliers -> set to NA -> impute
D %<>% mt_pre_outlier_lof(seq_k = c(5, 10, 20, 30, 40, 50)) %>%
  mt_pre_outlier_to_na(use_quant=TRUE, quant_thresh =0.025) %>% 
  mt_pre_impute_knn()

mt_write_se_xls(Dmc, file='2022-09-03_ADNI_Nightingale_Baseline.xlsx')

# Medication correction ----
all_cols <- D %>% colData %>% as_tibble() %>% names()
meds_to_exclude <- c("Med.Anticholinesterases","Med.NMDAAntag") # meds to exclude
med_cols <- all_cols[grep("Med", all_cols)] # column numbers of all meds
med_cols <- med_cols[which(med_cols%in%meds_to_exclude==F)] # meds to correct

# medication correction
Dmc <- D %>% mt_pre_confounding_correction_stepaic(cols_to_correct = med_cols, cols_to_exclude = meds_to_exclude)

mt_write_se_xls(Dmc, file='2022-09-03_ADNI_Nightingale_Baseline_medcor.xlsx')
