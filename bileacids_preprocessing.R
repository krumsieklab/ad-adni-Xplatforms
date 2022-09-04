# Bileacids dataset - data cleaning and preprocessing
#
#
# by RB, ZW
# last update: 2022-09-04

# clear workspace
rm(list = setdiff(ls(),c("codes.makepath","data.makepath","results.makepath"))) 
# libraries
library(openxlsx) # for excel reading and writing
library(maplet) # MT 
library(tidyverse) # %>%

# input
data_loc <- 'ADNI/Datasets/ADNI1_GO2_BileAcids/'
anno_loc <- 'ADNI/Datasets/ADNI_metadata/'
se_qc <- data.makepath(paste0(data_loc, 'processed_data/ADNI1_BA_SE_qc.rds'))
file_anno <- data.makepath(paste0(anno_loc, 'adni1go2_phenotypes_covariates.xlsx'))
file_outcome <- data.makepath(paste0(anno_loc, 'adni_outcomes.csv'))
file_meds <- data.makepath(paste0(anno_loc, 'adni_medications.csv'))

# output file
output_mt_ready <- '2022-09-03_bileacids.xlsx'
output_pp_no_med <- '2022-09-03_bileacids_preprocessing.xlsx'
output_pp_med <-  '2022-09-03_bileacids_preprocessing_medcor.xlsx'

# read files
raw_data_1 <- read.csv(data.makepath(paste0(data_loc, 'raw_data/BA_ADNI1.csv')), header = T, as.is = T)
raw_data_2 <- read.csv(data.makepath(paste0(data_loc, 'raw_data/BA_ADNI2GO.csv')), header = T, as.is = T)
nist_data_1 <- read.csv2(data.makepath(paste0(data_loc, 'raw_data/BA_NIST_ADNI1.tsv')), header = T, sep='\t', as.is = T)
nist_data_2 <- read.csv2(data.makepath(paste0(data_loc, 'raw_data/BA_NIST_ADNI2GO.tsv')), header = T, sep='\t', as.is = T)
fast_stat_1 <- read.csv2(data.makepath(paste0(anno_loc, 'FastingStatusADNI1.txt')), sep = "\t", as.is = T)
fast_stat_2 <- read.csv2(data.makepath(paste0(anno_loc, 'FastingStatusADNIGO2.txt')), sep = "\t", as.is = T)
meds <- read.csv2(data.makepath(paste0(anno_loc, 'MedicationsADNI1GO2.txt')), sep = "\t", as.is = T)

#### Data cleaning ----

# define columns that are sample information and to be discarded from measurements sheet and metabolite information sheet
samp_info_cols <- c("RID","Plate.Bar.Code","SAMPLE.BAR.CODE","SAMPLE.TYPE","SPECIES","MATERIAL",
                 "WELL.POSITION","SAMPLE.VOLUME","RUN.NUMBER","INJECTION.NUMBER","Sample.Identification","Cohort")

# harmonize raw_datas
names(raw_data_1) <- toupper(names(raw_data_1))
com_cols <- intersect(names(raw_data_1), names(raw_data_2))
raw_data_1 %<>% select(all_of(com_cols)) %>% mutate(Cohort='ADNI1') %>% 
  dplyr::rename( "Plate.Bar.Code" ="PLATE.BAR.CODE") %>% mutate_all(as.character)
raw_data_2 %<>% select(all_of(com_cols)) %>% mutate(Cohort='ADNIGO2') %>% 
  dplyr::rename( "Plate.Bar.Code" ="PLATE.BAR.CODE") %>% mutate_all(as.character)
raw_data <- bind_rows(raw_data_1, raw_data_2) %>%
  filter(!RID==999999) %>% mutate('Sample.Identification'=as.character(0))

# harmonize nist datas
names(nist_data_1) <- toupper(names(nist_data_1))
names(nist_data_2) <- toupper(names(nist_data_2))

nist_data_1 %<>% mutate_all(as.character) %>%
  dplyr::rename("Plate.Bar.Code"="PLATE.BAR.CODE","Sample.Identification"="SAMPLE.IDENTIFICATION")
nist_data_2 %<>% mutate_all(as.character) %>%
  dplyr::rename("Plate.Bar.Code"="PLATE.BAR.CODE","Sample.Identification"="SAMPLE.IDENTIFICATION")

nist_data <- bind_rows(nist_data_1, nist_data_2)

# combine raw data and nist data
raw_data <- bind_rows(raw_data, nist_data %>% 
                        select(all_of(names(raw_data %>% select(-c(RID, Cohort))))))

# colData
fast_stat_2 %<>%
  group_by(RID) %>%
  summarise(BIFAST = BIFAST[which(!is.na(BIFAST))[1]]) %>% ungroup()

fast_stat <- bind_rows(fast_stat_1, fast_stat_2) %>% mutate(RID=as.character(RID))
meds %<>% mutate(RID=as.character(RID))

sample_info <- raw_data %>% select(all_of(samp_info_cols)) %>% 
  left_join(fast_stat, by='RID') %>% unique() %>%
  left_join(meds, by='RID') %>% unique()

# rowData
met_info <- data.frame(col_ID=raw_data %>% select(-samp_info_cols) %>% names(), 
                       name=raw_data %>% select(-samp_info_cols) %>% names())

# assay
raw_data <- raw_data %>% select(-samp_info_cols) %>% data.frame()%>%
  mutate_all(as.matrix) %>% mutate_all(as.numeric) %>% t()

# MT pipeline (starting with %<>%, since D from se_qc file is already a SummarizedExperiment)
D %<>%
  # Nist based correction
  mt_pre_nist_based_correction(qc_samples= Sample.Identification==102, 
                               plate_col='Plate.Bar.Code') %>%
  # Add CV from duplicated ids to rowData
  mt_pre_cv(qc_samples = (RID%in%RID[duplicated(RID)]), 
               out_col = "Dup_cv", replicates=T, 
               id_col='RID') %>%
  # ICC for duplicates
  mt_pre_icc(qc_samples = (RID%in%RID[duplicated(RID)]), out_col = "Dup_icc", id_col='RID') %>%
  mt_anno_mutate(anno_type = 'samples', col_name='Sample_ID', term=paste0('sample_', 1:length(RID))) %>%
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
  # filter out metabolites with cv > 25% or ICC < 65%
  mt_modify_filter_metabolites(Dup_cv < 0.25 || Dup_icc > 0.65) %>%
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
  
  # Global outlier detection
  #mt_pre_outlier(method='leverage', pval=0.01) %>%
  
  #mt_modify_filter_samples(outlier_leverage != TRUE) %>%
  #mt_pre_outlier_detection_mahalanobis(pval=0.01) %>%  The data matrix is not full-rank, Mahalanobis cannot be computed
  #mt_pre_outlier_lof(seq_k = c(5, 10, 20, 30, 40, 50)) %>% # local outlier factor
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