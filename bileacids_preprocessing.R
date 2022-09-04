# bileacids dataset - data cleaning and preprocessing
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

# input  ----
data_loc <- 'ADNI/Datasets/ADNI1_GO2_BileAcids/'
anno_loc <- 'ADNI/Datasets/ADNI_metadata/'
data_file1 <- data.makepath(paste0(data_loc, 'raw_data/BA_ADNI1.csv'))
data_file2 <- data.makepath(paste0(data_loc, 'raw_data/BA_ADNI2GO.csv'))
qc_file1 <- data.makepath(paste0(data_loc, 'raw_data/BA_NIST_ADNI1.tsv'))
qc_file2 <- data.makepath(paste0(data_loc, 'raw_data/BA_NIST_ADNI2GO.tsv'))
fasting_file1 <- data.makepath(paste0(anno_loc, 'FastingStatusADNI1.txt'))
fasting_file2 <- data.makepath(paste0(anno_loc, 'FastingStatusADNIGO2.txt'))
med_file <- data.makepath(paste0(anno_loc, 'MedicationsADNI1GO2.txt'))
meds_to_exclude <- c("Med.Anticholinesterases","Med.NMDAAntag") # meds to exclude

# output file ----
output_mt_ready <- '2022-09-03_bileacids.xlsx'
output_pp_no_med <- '2022-09-03_bileacids_preprocessing.xlsx'
output_pp_med <-  '2022-09-03_bileacids_preprocessing_medcor.xlsx'
output_html <- '2022-09-03_bileacids_preprocessing_medcor.html'

# read files  ----
raw_data_1 <- read.csv(data_file1, header = T, as.is = T)
raw_data_2 <- read.csv(data_file2, header = T, as.is = T)
nist_data_1 <- read.csv2(qc_file1, header = T, sep='\t', as.is = T)
nist_data_2 <- read.csv2(qc_file2, header = T, sep='\t', as.is = T)
fast_stat_1 <- read.csv2(fasting_file1, sep = "\t", as.is = T)
fast_stat_2 <- read.csv2(fasting_file2, sep = "\t", as.is = T)
meds <- read.csv2(med_file, sep = "\t", as.is = T)
# define columns that are sample information and to be discarded from measurements sheet and metabolite information sheet
samp_info_cols <- c("RID","Plate.Bar.Code","SAMPLE.BAR.CODE","SAMPLE.TYPE","SPECIES","MATERIAL",
                 "WELL.POSITION","SAMPLE.VOLUME","RUN.NUMBER","INJECTION.NUMBER","Sample.Identification","Cohort")

# harmonize raw_datas  ----
names(raw_data_1) <- toupper(names(raw_data_1)) # column names
com_cols <- intersect(names(raw_data_1), names(raw_data_2)) # file common columns
  # select common columns and harmonize names and datatype
raw_data_1 %<>% select(all_of(com_cols)) %>% mutate(Cohort='ADNI1') %>% 
  dplyr::rename( "Plate.Bar.Code" ="PLATE.BAR.CODE") %>% mutate_all(as.character)
  # same for second data file
raw_data_2 %<>% select(all_of(com_cols)) %>% mutate(Cohort='ADNIGO2') %>% 
  dplyr::rename( "Plate.Bar.Code" ="PLATE.BAR.CODE") %>% mutate_all(as.character)
  # combine the two data files
raw_data <- bind_rows(raw_data_1, raw_data_2) %>%
  filter(!RID==999999) %>% mutate('Sample.Identification'=as.character(0))

# harmonize nist datas  ----
names(nist_data_1) <- toupper(names(nist_data_1))
names(nist_data_2) <- toupper(names(nist_data_2))

nist_data_1 %<>% mutate_all(as.character) %>%
  dplyr::rename("Plate.Bar.Code"="PLATE.BAR.CODE","Sample.Identification"="SAMPLE.IDENTIFICATION")
nist_data_2 %<>% mutate_all(as.character) %>%
  dplyr::rename("Plate.Bar.Code"="PLATE.BAR.CODE","Sample.Identification"="SAMPLE.IDENTIFICATION")

nist_data <- bind_rows(nist_data_1, nist_data_2)

# combine raw data and nist data  ----
raw_data <- bind_rows(raw_data, nist_data %>% 
                        select(all_of(names(raw_data %>% select(-c(RID, Cohort))))))

# summarized experiment ----
  # colData
    # fasting data
fast_stat_2 %<>%
  group_by(RID) %>%
  summarise(BIFAST = BIFAST[which(!is.na(BIFAST))[1]]) %>% ungroup()
fast_stat <- bind_rows(fast_stat_1, fast_stat_2) %>% mutate(RID=as.character(RID))
    # medication data
meds %<>% mutate(RID=as.character(RID))
    # sample information
sample_info <- raw_data %>% select(all_of(samp_info_cols)) %>% 
  left_join(fast_stat, by='RID') %>% unique() %>%
  left_join(meds, by='RID') %>% unique()

  # rowData
met_info <- data.frame(col_ID=raw_data %>% select(-samp_info_cols) %>% names(), 
                       name=raw_data %>% select(-samp_info_cols) %>% names())

  # assay
raw_data <- raw_data %>% select(-samp_info_cols) %>% data.frame()%>%
  mutate_all(as.matrix) %>% mutate_all(as.numeric) %>% t()

D <- SummarizedExperiment::SummarizedExperiment(assays = raw_data,
                                                colData = sample_info, rowData = met_info)

# qc ----
D %<>%
  # Nist based correction
  mt_pre_reference_correction(qc_samples= Sample.Identification==102, 
                               plate_col='Plate.Bar.Code') %>%
  # Add CV from duplicated ids to rowData
  mt_pre_cv(qc_samples = (RID%in%RID[duplicated(RID)]), 
               out_col = "Dup_cv", replicates=T, 
               id_col='RID') %>%
  # ICC for duplicates
  mt_pre_icc(qc_samples = (RID%in%RID[duplicated(RID)]), out_col = "Dup_icc", id_col='RID') %>%
  mt_anno_mutate(anno_type = 'samples', col_name='Sample_ID', term=paste0('sample_', 1:length(RID))) %>%
  {.}

mt_write_se_xls(D, file=output_mt_ready)

# preprocessing  ----
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
  mt_modify_filter_features(Dup_cv < 0.25 || Dup_icc > 0.65) %>%
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