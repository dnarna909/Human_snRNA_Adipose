#' change "D:/" to "/media/jianie/Extreme SSD1/" # if in LUNIX computer

rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.
memory.limit(size = 80000000)

# set parameters
source(paste0("/media/jianie/Extreme SSD1/2022-09-01 STARR_SGLT2 Combine/Project Parameters.R"), local = knitr::knit_global())
# source(paste0("D:/2022-09-01 STARR_SGLT2 Combine/Project Parameters.R"), local = knitr::knit_global())

library(dplyr)

# add patient information ----------------------------------------------------------------
# combine
combine.df <- readRDS(paste0(sample.info.disk, sample.info.path, sample.file.name)) %>% 
  mutate(Group = ifelse(`Group Assignment:` == "Nutritional counseling", "C",
                        ifelse(`Group Assignment:` == "Drug - dapagliflozin", "D", NA)))

library(readxl)
# insulin result -1
Insulin_Results_1 <- rio::import(paste0(sample.info.disk, sample.info.path, "STARR_Insulin Measurement.xlsx"),  
                                sheet = "Insulin_1", col_types = c("text", "text", "numeric", "numeric"))

# insulin result - 2
Insulin_Results_2 <- rio::import(paste0(sample.info.disk, sample.info.path, "STARR_Insulin Measurement.xlsx"),  
                                sheet = "Insulin_2", col_types = c("text", "numeric"))

combine.df <- combine.df %>% left_join(Insulin_Results_1, by = c("sample_id" = "id") ) %>% 
  left_join(Insulin_Results_2, by = c("sample_id" = "id") ) 
colnames(combine.df)
combine.df$`Insulin(µIU/mL)_MSD`

# combine
# patient <- sglt2 %>% full_join(starr, by = intersect(colnames(sglt2), colnames(starr))) 
patient <- combine.df %>% arrange(Agg.List) %>% 
  mutate(Dataset = paste0("Sample", stringr::str_pad( 1: n(), 2, pad = "0") ),
         `HOMA-IR` = Glucose * `Insulin(µIU/mL)_MSD`/405) 
#, HOMA-IR = Glucose_mM * Insulin / 22.5, Glucose_mg/dL * Insulin / 405
# HOMA-IR was calculated according to the formula: fasting insulin (microU/L) x fasting glucose (nmol/L)/22.5.
unique(patient$Dataset)
table(patient$Dataset, patient$Agg.List)
saveRDS(patient, file = paste0(sample.info.disk, sample.info.path, "patients.rds"))

library(gt)
library(tidyverse)
library(glue)
patient %>%
  dplyr::filter() %>% rename(`Fasting Glucose` = Glucose, `HbA1c Value` = `A1c Value`,  `Diabetes Diagostic` = dbc4_diabetes) %>%
  select(all_of(c("Dataset", "Gender" ,"Age", "BMI", "Fasting Glucose", "HbA1c Value"))) %>%
  gt() %>%
  tab_header(
    title = "Patient Information",
    # subtitle = glue("{start_date} to {end_date}")
  ) %>%
  # fmt_date(
  #   columns = date,
  #   date_style = 3
  # ) %>%
  # fmt_currency(
  #   columns = c(open, high, low, close),
  #   currency = "USD"
  # ) %>%
  fmt_number(
    columns = c(Age, BMI, `Fasting Glucose`, `HbA1c Value`),
    suffixing = TRUE
  ) %>% gtsave(paste0(Disk, Project.folder, "/", figures.folder, "/", "Figure 1_patient information.png"), expand = 10)

