# Wishart dataset - data cleaning and preprocessing
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
data_loc <- 'ADNI/Datasets/ADNI1_Wishart_HV/'
anno_loc <- 'ADNI/Datasets/ADNI_metadata/'

# output file
output_se_qc <- 'processed_data/ADNI1_Wishart_qc.rds'
se_qc <- data.makepath(paste0(data_loc, output_se_qc))

# read datas
raw_data <- read.csv(data.makepath(paste0(data_loc, 'raw_data/ADNI-1.hvmetabolites.csv')), header = T)
qc_data <- read.csv2(data.makepath(paste0(data_loc, 'raw_data/Catecholamine_QC_new.csv')), sep = ",",header = T)
fast_stat <- read.csv2(data.makepath(paste0(anno_loc, 'FastingStatusADNI1.txt')), sep = "\t", as.is = T)
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


# read datas
raw_data <- read.csv(file_data, header = T)
# read datas
fast_stat <- read.csv2(fast_data, sep = "\t")
meds <- read.csv2(med_data, sep = "\t")

# define columns to be discarded
cols2discard <- names(raw_data)[1:3]

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
                                                colData = sample_info, rowData = met_info) %>% 
  mt_files_write_SE(file = se_qc)

# List for Medication correction
meds_to_exclude <- D %>% colData %>% as_tibble() %>% names() %>% as.data.frame() %>% rownames_to_column() %>% 
  column_to_rownames(var = '.') %>% t() %>% as.data.frame() %>% select(all_of(meds_to_exclude)) %>% 
  as.matrix() %>% as.numeric()

med_cols <- D %>% colData %>% as_tibble() %>% names() %>% grep(pattern = 'Med') %>% 
  discard(~ .x %in% meds_to_exclude) %>% write.table(data.makepath(paste0(anno_loc, 'adni_medications.csv')),sep = ",",row.names = F,quote = F)


# QC 
# libraries
library(openxlsx) # for excel reading and writing
library(MetaboTools) # MT 
library(tidyverse) # %>%
library(magrittr)

# Input files ----
data_loc <- 'ADNI/Datasets/ADNI1_Wishart_HV/'
data_file <- 'raw_data/ADNI-1.hvmetabolites.csv'
anno_loc <- 'ADNI/Datasets/ADNI_metadata/'
fasting_file <- 'FastingStatusADNI1.txt'
med_file <- 'MedicationsADNI1GO2.txt'
output_loc <- 'ADNI/Datasets/ADNI1_Wishart_HV/'
output_xlsx <- 'raw_data/mt_ready_Wishart_HV.xlsx'
file_data <- data.makepath(paste0(data_loc, data_file))
xlsx_output <- data.makepath(paste0(data_loc, output_xlsx))
fast_data <- data.makepath(paste0(anno_loc, fasting_file))
med_data <- data.makepath(paste0(anno_loc, med_file))


# read datas
raw_data <- read.csv(file_data, header = T)
# read datas
fast_stat <- read.csv2(fast_data, sep = "\t")
meds <- read.csv2(med_data, sep = "\t")

# Data load ----

cols2discard <- names(raw_data)[1:3]

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

# SE ----

D0 <- SummarizedExperiment::SummarizedExperiment(assays = raw_data,
                                                 colData = sample_info, rowData = met_info)

# Data preprocessing ----

D0 %<>%
  
  # Add CV from duplicated ids to rowData
  mt_modify_cv(qc_samples = (RID%in%RID[duplicated(RID)]),
               col_lab = "Dup_cv", replicates=T,
               id_col='RID') %>%
  # ICC for duplicates
  mt_modify_icc(qc_samples = (RID%in%RID[duplicated(RID)]), col_lab = "Dup_icc", id_col='RID') %>%
  mt_modify_mutate(anno_type = 'samples', col_name='Sample_ID', term=paste0('sample_', 1:length(RID)))


# Write out the files ----
# take assay data, add metabolite names column and add RID to column names
assay_df <- D0 %>% assay() %>% t() %>% data.frame()
names(assay_df) <- D0 %>% rowData() %>% data.frame() %>% pull(name)

assay_df %<>% dplyr::bind_cols( D0 %>% colData() %>% data.frame() %>% select(Sample_ID), .)

# work book
wb <- createWorkbook()

addWorksheet(wb, "data")
writeData(wb, "data", assay_df, colNames = T, rowNames = F)

addWorksheet(wb, "sampleinfo")
writeData(wb, "sampleinfo", colData(D0))

addWorksheet(wb, "metinfo")
writeData(wb, "metinfo", rowData(D0))

saveWorkbook(wb, file = xlsx_output, overwrite = T)


### File paths ----

# Input files
data_loc <- 'ADNI/Datasets/ADNI1_Wishart_HV/'
anno_loc <- 'ADNI/Datasets/ADNI_metadata/'
se_qc <- data.makepath(paste0(data_loc, 'processed_data/ADNI1_Wishart_qc.rds'))
file_anno <- data.makepath(paste0(anno_loc, 'adni1go2_phenotypes_covariates.xlsx'))
file_outcome <- data.makepath(paste0(anno_loc, 'adni_outcomes.csv'))
file_meds <- data.makepath(paste0(anno_loc, 'adni_medications.csv'))

# Output files
output_loc <- 'ADNI/Datasets/ADNI1_Wishart_HV/'
file_output_pp_no_med <- results.makepath(paste0(output_loc,  'ADNI1_Wishart_PP_no_med', '.html'))
file_output_pp_med <- results.makepath(paste0(output_loc, 'ADNI1_Wishart_PP_med', '.html'))
se_output_pp_no_med <- results.makepath(paste0(output_loc, 'ADNI1_Wishart_PP_no_med', '.rds'))
se_output_pp_med <- results.makepath(paste0(output_loc, 'ADNI1_Wishart_PP_med', '.rds'))
file_output_analysis_med <- results.makepath(paste0(output_loc, 'ADNI_Wishart_analysis' ,'.html'))
se_output_analysis_med <- results.makepath(paste0(output_loc, 'ADNI1_Wishart_analysis' ,'.rds'))


### Load outcome and medication variable lists ----
outcomes <- read.csv(file_outcome, header = T, as.is = T, stringsAsFactors = F)
med_cols <- read.csv(file_meds,header = T, as.is = T, stringsAsFactors = F) %>% as.matrix() %>% as.numeric()


### Start MT pipeline ----

# Load the SummarizedExperiment (=data)
load(se_qc)

### QC ----

# MT pipeline (starting with %<>%, since D from se_qc file is already a SummarizedExperiment)
D %<>%
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
  
  mt_pre_zeroToNA() %>%
  # Filter out non-fasting samples 
  mt_modify_filter_samples(BIFAST == "Yes") %>%
  
  mtw_missingness(plot_options = list(met_max = 0.25,samp_max = 1),
                  filter_options = list(met_max = 0.25,sample_max = 1)) %>%
  
  # Filter out metabolites with cv > 25% or ICC < 65%
  mt_modify_filter_metabolites(Dup_cv < 0.25 || Dup_icc > 0.65) %>%
  
  # Preprocess
  mt_reporting_heading(strtitle = "Normalization", lvl = 2) %>%
  
  mtw_preprocess(quot_options = list(met_max = 0),
                 log_base = 2,do_impute = T) %>%
  
  # Average-combine duplicate samples
  mt_modify_averagesample(group_by = "RID") %>%
  
  # Global outlier detection
  mt_pre_outlier(method='leverage', pval=0.01) %>%
  
  mt_modify_filter_samples(outlier_leverage != TRUE) %>%
  {.}

### Global statistics ----
D %<>%  
  
  # Global stats as data overview
  mt_reporting_heading(strtitle = "Global Statistics", lvl = 1) %>%
  
  mtw_global_stats(annos_pca_umap = list(quote(outlier_leverage)),
                   pca_options = list(scaledata = T, ggadd=quote(scale_size_identity()), size = 2),
                   umap_options = list(scaledata = T, ggadd=quote(scale_size_identity()), size = 2),
                   heatmap_options = list(scaledata = T, fontsize = 5)) %>%
  
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

# libraries
library(openxlsx) # for excel reading and writing
library(MetaboTools) # MT 
library(tidyverse) # %>%
library(magrittr)

# Input files ----
data_loc <- 'ADNI/Datasets/ADNI1_Wishart_HV/'
data_file <- 'raw_data/mt_ready_Wishart_HV.xlsx'
output_loc <- 'ADNI/Datasets/ADNI1_Wishart_HV/'
output_file <- 'ADNI1_Wishart_PP_no_med.html'
output_file_1 <- 'ADNI1_Wishart_PP_med.html'
output_se <- 'processed_data/ADNI1_Wishart_SE_PP_no_med.rds'
output_se_1 <- 'processed_data/ADNI1_Wishart_SE_PP_med.rds'
file_data <- data.makepath(paste0(data_loc, data_file))
file_output <- results.makepath(paste0(output_loc,output_file))
file_output <- results.makepath(paste0(output_loc,output_file_1))
se_output <- data.makepath(paste0(data_loc, output_se))
se_output_1 <- data.makepath(paste0(data_loc, output_se_1))


# Data load ----
D <- 
  # load data
  mt_files_data_xls(file=file_data, sheet="data", samples_in_row=T, ID_col="Sample_ID", zero_to_NA = T) %>% 
  # validate checksum
  # mt_files_checksum(file=file_data, checksum = "80afcd72481c6cf3dcf83342e3513699") %>%
  # load metabolite annotations
  mt_files_anno_xls(file=file_data, sheet="metinfo",anno_type="metabolites", anno_ID="name", data_ID="name") %>%
  # load clinical annotations
  mt_files_anno_xls(file=file_data, sheet="sampleinfo", anno_type="samples", anno_ID="Sample_ID", data_ID="Sample_ID") %>%
  # print infos about dataset
  mt_logging_datasetinfo()


# Data preprocessing ----

D %<>%
  
  mt_pre_zeroToNA() %>%
  
  mt_reporting_heading(strtitle = "Preprocessing", lvl = 1) %>%
  
  mt_reporting_heading(strtitle = "Filtering", lvl = 2) %>%
  
  # Filter out non-fasting samples 
  mt_modify_filter_samples(BIFAST == "Yes") %>%
  
  mtw_missingness(plot_options = list(met_max = 0.25,samp_max = 1),
                  filter_options = list(met_max = 0.25,sample_max = 1)) %>%
  
  # Filter out metabolites with cv > 25% or ICC < 65%
  mt_modify_filter_metabolites(Dup_cv < 0.25 || Dup_icc > 0.65) %>%
  
  # Preprocess
  mt_reporting_heading(strtitle = "Normalization", lvl = 2) %>%
  
  mtw_preprocess(quot_options = list(met_max = 0),
                 log_base = 2,do_impute = T) %>%
  
  # Average-combine duplicate samples
  mt_modify_averagesample(group_by = "RID") %>%
  
  # Global outlier detection
  mt_pre_outlier(method='leverage', pval=0.01) %>%
  
  mt_modify_filter_samples(outlier_leverage != TRUE) %>%
  
  # Global stats as data overview
  mt_reporting_heading(strtitle = "Global Statistics", lvl = 2) %>%
  
  mtw_global_stats(annos_pca_umap = list(quote(outlier_leverage)),
                   pca_options = list(scaledata = T, ggadd=quote(scale_size_identity()), size = 2),
                   umap_options = list(scaledata = T, ggadd=quote(scale_size_identity()), size = 2),
                   heatmap_options = list(scaledata = T, fontsize = 5))

# Medication correction ----

# meds for correction
meds_to_cor <- c("Med.Fish.Oils","Med.Nicotinic.acid","Med.COQ10","Med.AcetylCarn",
                 "Med.ThyroidHormInhib", "Med.ThyroidHormRecAg",
                 "Med.NSAIDs","Med.Acetaminaphen",
                 "Med.LoopDiur","Med.Thiazide","Med.K.SparAg","Med.K.AldostAntag","Med.Biguanides","Med.Sulfonylureas","Med.Insulin",
                 "Med.Thiazolidinediones","Med.Statin","Med.Resins","Med.Fibrates","Med.OtherLipids","Med.ACE.inhib",
                 "Med.AngiotensinIIant","Med.Alphaantagonis","Med.BetablockSelec","Med.BetablockNonSelec","Med.Alphaandbetablock","Med.Imidazoline.RAgon",
                 "Med.Verapamil","Med.Amiodarone","Med.Organic.nitrates","Med.Digitalis","Med.AntiEpileptic","Med.AntipsycAgents",
                 "Med.Bupropion","Med.Benzodiazepine","Med.OtherHypnotic","Med.SerotoninReupInh","Med.SerotNorepReupInh","Med.TrycicAntidep", 
                 "Med.MelatoninAntag") 

# all column numbers
med_cols <- D %>% colData %>% as_tibble() %>% names() %>% as.data.frame() %>% rownames_to_column() %>% 
  column_to_rownames(var = '.') %>% t() %>% as.data.frame() %>% select(all_of(meds_to_cor)) %>% as.matrix() %>% as.numeric()

D1 <- D %>% mt_pre_confounding_correction_stepwise_aic(cols_to_cor = med_cols)

# Data output ----

# Report of preprocessing
D %>% mt_reporting_html(outfile=file_output, output.calls = T)
D1 %>% mt_reporting_html(outfile=file_output_1, output.calls = T)

# Save the SE
save(D, file = se_output)
save(D1, file = se_output_1)

