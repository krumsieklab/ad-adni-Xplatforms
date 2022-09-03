# Bileacids dataset - data cleaning and preprocessing
#
#
# by RB, ZW
# last update: 2022-09-03


#### Setup ----

# clear workspace
rm(list = setdiff(ls(),c("codes.makepath","data.makepath","results.makepath"))) 
# libraries
library(MetaboTools) # MT 
library(tidyverse) # %>%

# Helper function
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

# paths
data_loc <- 'ADNI/Datasets/ADNI1_GO2_BileAcids/'
anno_loc <- 'ADNI/Datasets/ADNI_metadata/'

# output file
output_se_qc <- 'processed_data/ADNI1_BA_SE_qc.rds'
se_qc <- data.makepath(paste0(data_loc, output_se_qc))

# Input files
data_loc <- 'ADNI/Datasets/ADNI1_GO2_BileAcids/'
anno_loc <- 'ADNI/Datasets/ADNI_metadata/'
se_qc <- data.makepath(paste0(data_loc, 'processed_data/ADNI1_BA_SE_qc.rds'))
file_anno <- data.makepath(paste0(anno_loc, 'adni1go2_phenotypes_covariates.xlsx'))
file_outcome <- data.makepath(paste0(anno_loc, 'adni_outcomes.csv'))
file_meds <- data.makepath(paste0(anno_loc, 'adni_medications.csv'))

# Output files
output_loc <- 'ADNI/Datasets/ADNI1_GO2_BileAcids/'
file_output_pp_no_med <- results.makepath(paste0(output_loc,  'ADNI1_BA_PP_no_med', '.html'))
file_output_pp_med <- results.makepath(paste0(output_loc, 'ADNI1_BA_PP_med', '.html'))
se_output_pp_no_med <- results.makepath(paste0(output_loc, 'ADNI1_BA_PP_no_med', '.rds'))
se_output_pp_med <- results.makepath(paste0(output_loc, 'ADNI1_BA_PP_med', '.rds'))
file_output_analysis_med <- results.makepath(paste0(output_loc, 'ADNI_BA_analysis' ,'.html'))
se_output_analysis_med <- results.makepath(paste0(output_loc, 'ADNI1_BA_analysis' ,'.rds'))


### Load outcome and medication variable lists ----
outcomes <- read.csv(file_outcome, header = T, as.is = T, stringsAsFactors = F)
med_cols <- read.csv(file_meds,header = T, as.is = T, stringsAsFactors = F) %>% as.matrix() %>% as.numeric()



# read files
raw_data_1 <- read.csv(data.makepath(paste0(data_loc, 'raw_data/BA_ADNI1.csv')), header = T, as.is = T)
raw_data_2 <- read.csv(data.makepath(paste0(data_loc, 'raw_data/BA_ADNI2GO.csv')), header = T, as.is = T)
nist_data_1 <- read.csv2(data.makepath(paste0(data_loc, 'raw_data/BA_NIST_ADNI1.tsv')), header = T, sep='\t', as.is = T)
nist_data_2 <- read.csv2(data.makepath(paste0(data_loc, 'raw_data/BA_NIST_ADNI2GO.tsv')), header = T, sep='\t', as.is = T)
fast_stat_1 <- read.csv2(data.makepath(paste0(anno_loc, 'FastingStatusADNI1.txt')), sep = "\t", as.is = T)
fast_stat_2 <- read.csv2(data.makepath(paste0(anno_loc, 'FastingStatusADNIGO2.txt')), sep = "\t", as.is = T)
meds <- read.csv2(data.makepath(paste0(anno_loc, 'MedicationsADNI1GO2.txt')), sep = "\t", as.is = T)
outcomes <- read.csv(data.makepath(paste0(anno_loc, 'traits.csv')), header = T, as.is = T,stringsAsFactors = F)
meds_to_exclude <- read_list_from_file(data.makepath(paste0(anno_loc, 'ADNI_meds_rm.txt'))) 


#### Data cleaning ----

# edit outcome name list
# Remove the one used in confounding correction
rows2rm <- c('ApoE4')
# Rename outcomes for analysis
outcomes1 <- outcomes %>% select(c('trait','var_type')) %>% rename(outcome=trait, type=var_type) %>% 
  subset(!(outcome %in% rows2rm)) %>% as.matrix() %>% gsub('<', '_', .) %>% gsub('\\(', '_', .) %>% 
  gsub('\\)', '_', .) %>% gsub('/', '_', .) %>% gsub('Diagnosis', 'Diag_num', .) 

outcomes1 %<>% data.frame(row.names = 1:nrow(outcomes1) ,stringsAsFactors = F) 
outcomes1 %>% write.table(data.makepath(paste0(anno_loc, 'adni_outcomes.csv')),sep = ",",row.names = F,quote = F)

# define colums to be discarded
cols2discard = c("RID","Plate.Bar.Code","SAMPLE.BAR.CODE","SAMPLE.TYPE","SPECIES","MATERIAL",
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
raw_data <- bind_rows(raw_data, nist_data %>% select(all_of(names(raw_data %>% select(-c(RID, Cohort))))))

# colData
fast_stat_2 %<>%
  group_by(RID) %>%
  summarise(BIFAST = BIFAST[which(!is.na(BIFAST))[1]]) %>% ungroup()

fast_stat <- bind_rows(fast_stat_1, fast_stat_2) %>% mutate(RID=as.character(RID))
meds %<>% mutate(RID=as.character(RID))

sample_info <- raw_data %>% select(all_of(cols2discard)) %>% 
  left_join(fast_stat, by='RID') %>% unique() %>%
  left_join(meds, by='RID') %>% unique()

# rowData
met_info <- data.frame(col_ID=raw_data %>% select(-cols2discard) %>% names(), 
                       name=raw_data %>% select(-cols2discard) %>% names())

# assay
raw_data <- raw_data %>% select(-cols2discard) %>% data.frame()%>%
  mutate_all(as.matrix) %>% mutate_all(as.numeric) %>% t()



#### Write output files ----

# Summarazed Experiment
D <- SummarizedExperiment::SummarizedExperiment(assays = raw_data,
                                                colData = sample_info, rowData = met_info) %>% 
  mt_files_write_SE(file = se_qc)

# List for Medication correction
meds_to_exclude <- D %>% colData %>% as_tibble() %>% names() %>% as.data.frame() %>% rownames_to_column() %>% 
  column_to_rownames(var = '.') %>% t() %>% as.data.frame() %>% select(all_of(meds_to_exclude)) %>% 
  as.matrix() %>% as.numeric()

med_cols <- D %>% colData %>% as_tibble() %>% names() %>% grep(pattern = 'Med') %>% 
  discard(~ .x %in% meds_to_exclude) %>% write.table(data.makepath(paste0(anno_loc, 'adni_medications.csv')),sep = ",",row.names = F,quote = F)

### QC ----
# MT pipeline (starting with %<>%, since D from se_qc file is already a SummarizedExperiment)
D %<>%
  # Nist based correction
  mt_pre_nist_based_correction(qc_samples= Sample.Identification==102, 
                               plate_col='Plate.Bar.Code') %>%
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
  mtw_missingness(plot_options = list(met_max = 0.4,samp_max = 1),
                  filter_options = list(met_max = 0.4,sample_max = 1)) %>%
  
  # Filter out metabolites with cv > 25% or ICC < 65%
  mt_modify_filter_metabolites(Dup_cv <0.25 || Dup_icc >0.65) %>%
  
  # Preprocess using the wrapper function
  mt_reporting_heading(strtitle = "Normalization", lvl = 2) %>%
  
  # Including quotient normalization, log transformation and kNN imputation
  mtw_preprocess(quot_options = list(met_max = 0.4),
                 log_base = 2,do_impute = T) %>%
  
  # Average-combine duplicate samples
  mt_modify_averagesample(group_by = "RID") %>%
  {.}

### Global statistics ----
D %<>%  
  
  # Global stats as data overview
  mt_reporting_heading(strtitle = "Global Statistics", lvl = 2) %>%
  
  mtw_global_stats(pca_options = list(scaledata = T, size = 2, ggadd=quote(scale_size_identity())),
                   umap_options = list(scaledata = T, size = 2, ggadd=quote(scale_size_identity())),
                   heatmap_options = list(scaledata = T, fontsize = 5)) %>%
  
  # Global outlier detection
  mt_pre_outlier(method='mahalanobis', pval=0.01) %>%
  
  # Report before medication correction
  mt_reporting_html(outfile=file_output_pp_no_med, output.calls = T) %>%
  # Save SE before medication correction
  mt_files_write_SE(file = se_output_pp_no_med) %>%
  
  # Medication correction
  mt_pre_confounding_correction_stepwise_aic(cols_to_cor = med_cols) %>%
  
  # Report after medication correction
  mt_reporting_html(outfile=file_output_pp_med, output.calls = T) %>%
  # Save SE after medication correction
  mt_files_write_SE(file = se_output_pp_med) %>%
  {.}
