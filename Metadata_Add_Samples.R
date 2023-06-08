###' DISCLAIMER
#' Some of the algorithms are non-deterministic making the results slightly different from run to run.
#' Many of the algorithms are frequently updated with new releases making the results slightly different from version to version.
#' Depending on your system, this code below may not produce an exact 1:1 copy of the results. 
#' This adversely affects clustering, filtering etc. if the code is naively copy-pasted.
#' For better reproducibility, please download quality-filtered Rds objects and work on downstream analyses. 
#' change "D:/" to "/media/jianie/Extreme SSD1/" # if in LUNIX computer

##' Load libraries-------------------------------------------------------------------------------------------------------------
# rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memory and report the memory usage.
memory.limit(size = 1e+13) # 'memory.limit()' is Windows-specific
source(paste0(Disk, Project.folder, "/","Project Parameters.R"), local = knitr::knit_global())

# import packages
# library("velocyto.R") # linux -specific


# Import file name -----------------------------------------------------------------------------------------------
basic.meta.data <- readRDS(paste0(Disk, Project.folder, "/",  Rds.folder, "/","SAT_basic.meta.data.Rds"))
table(basic.meta.data$Dataset)

# import sample.file.name1: Aggregation Information
sample.data1 <- data.frame()
for (ss in sample.file.name1) {
  sample <-  gsub('\\.', '-', stringr::str_split(ss, "_")[[1]][1])
  sample1 <- rio::import(paste0(sample.info.disk1, sample.info.path1, ss)) %>% 
    mutate(
      Project.Aggr = sample,
      # List = 1:n(),
      Dataset = paste0(sample, "_", "Sample", stringr::str_pad( 1:n(), 2, pad = "0") )
    )%>% 
    rename_with(make.names) %>% 
    as.data.frame() %>%
    dplyr::select(!one_of(c("...16"))) %>%
    rename(List = Agg.List,
           Seq.Status = Status)
  sample.data1 <- sample.data1 %>% bind_rows(sample1)
  #  colnames(sample.data1 %>% dplyr::select_if(is.numeric) )
}
rm(sample1, ss)


# import sample.file.name2: All Sample Information
sample.data2 <- data.frame()
for (ss in sample.file.name2) {
  sample2 <- rio::import(paste0(sample.info.disk2, sample.info.path2, ss)) %>% 
    rename_with(make.names) %>% 
    as.data.frame()
  sample.data2 <- sample.data2 %>% bind_rows(sample2)
  colnames(sample.data2 %>% select_if(is.numeric) )
  # sglt2$Samples
}
rm(sample2, ss)

# merge to sample.data
intersect(colnames(sample.data1), colnames(sample.data2))
sample.data <- left_join(sample.data1, sample.data2, 
                         by = intersect(colnames(sample.data1), colnames(sample.data2)) )
# sample.data <- full_join(sample.data1, sample.data2, 
#                          by = intersect(colnames(sample.data1), colnames(sample.data2)))

# further modify
colnames(sample.data)[sapply(sample.data, is.character)]
table(sample.data$Dataset)
table(sample.data$Gender)
table(sample.data$Diabetes)
table(sample.data$Pre_Post)
table(sample.data$Race)
table(sample.data$Group.Assignment.)
table(sample.data$Group)
table(sample.data$Age_Group)
table(sample.data$Status)
table(sample.data$Age_Group, sample.data$Status)
table(sample.data$Project)
table(sample.data$Project.Aggr)
table(sample.data$sample_id, sample.data$Group.Assignment.)
table(sample.data$Ethnicity)
table(sample.data$Ethnicity, sample.data$Project)
table(sample.data$Prediabetes1, sample.data$Prediabetes2)
table(sample.data$Prediabetes1, sample.data$Prediabetes2, sample.data$Project)

sample.data$Status <- factor(sample.data$Status, 
                             levels = c("Obese", "Overweight", "Lean", "Underweight"))
sample.data$Age_Group <- factor(sample.data$Age_Group, 
                                levels = c( "Older","Middle", "Young"))
sample.data$Group <- factor(sample.data$Group, 
                            levels = c("Notreatment", "NUT", "SGLT2i"))
sample.data$Diabetes <- factor(sample.data$Diabetes, 
                               levels = c("Yes", "No"))
sample.data$Pre_Post <- factor(sample.data$Pre_Post, 
                               levels = c("Notreatment", "Pre", "Post"))
sample.data$Prediabetes1 <- factor(sample.data$Prediabetes1, 
                                   levels = c("No", "PreD", "Diabetes"))
sample.data$Prediabetes2 <- factor(sample.data$Prediabetes2, 
                                   levels = c("No", "PreD", "Diabetes"))

sample.data <- sample.data %>% mutate(Treatment_Group1 = factor(paste0(Group, "_", Pre_Post), 
                                                               levels = c("Notreatment_Notreatment", 
                                                                          "NUT_Pre", "NUT_Post",
                                                                          "SGLT2i_Pre", "SGLT2i_Post")),
                                      BA_Group = ifelse(Status == "Lean" & Prediabetes1 == "No" & Age_Group == "Older" & Diabetes == "No", "Older_Lean",
                                                        ifelse(Status == "Overweight" & Prediabetes1 == "No" & Age_Group == "Older" & Diabetes == "No", "Older_Overweight",
                                                               ifelse(Status == "Obese" & Prediabetes1 == "PreD" & Age_Group == "Older" & Diabetes == "No", "Older_PreD_Obese",
                                                                      ifelse(Status == "Underweight" & Prediabetes1 == "No" & Age_Group == "Older" & Diabetes == "No", "Older_Underweight",
                                                                             ifelse(Status == "Lean" & Prediabetes1 == "PreD" & Age_Group == "Older" & Diabetes == "No", "Older_PreD_Lean",
                                                                                    ifelse(Status == "Overweight" & Prediabetes1 == "PreD" & Age_Group == "Older" & Diabetes == "No", "Older_PreD_Overweight",
                                                                                           
                                                                             ifelse(Status == "Lean" & Prediabetes1 == "No" & Age_Group == "Middle" & Diabetes == "No", "Middle_Lean",
                                                                                    ifelse(Status == "Overweight" & Prediabetes1 == "No" & Age_Group == "Middle" & Diabetes == "No", "Middle_Overweight",
                                                                                           ifelse(Status == "Obese" & Prediabetes1 == "No" & Age_Group == "Middle" & Diabetes == "No", "Middle_Obese",
                                                                                                  
                                                                                                  ifelse(Prediabetes1 == "Diabetes" | Diabetes == "Yes", "Diabetes", "No"
                                                                                           )))))))))))
table(sample.data$BA_Group); unique(sample.data$BA_Group); is.na(sample.data$BA_Group)
sample.data$BA_Group <- factor(sample.data$BA_Group, levels = c("Older_Underweight", 
                                                                "Older_Lean", "Older_Overweight", "Older_PreD_Obese", "Diabetes",
                                                                "Older_PreD_Lean", "Older_PreD_Overweight", 
                                                                "Middle_Lean", "Middle_Overweight", "Middle_Obese"))# "Older_PreD_Lean", "Older_PreD_Overweight",
table(sample.data$BA_Group); unique(sample.data$BA_Group); is.na(sample.data$BA_Group)

pd.ids <- sample.data %>%
  group_by(sample_id) %>% dplyr::filter(row_number(sample_id) == 1)%>%
  ungroup() %>%
  as.data.frame()  %>% 
  dplyr::filter(
  Pre_Post %in% c("Pre", "Post") 
) %>% 
  group_by(patient_id) %>%
  mutate(n = n()) %>%
  dplyr::filter(n == 2) %>%
  ungroup() %>%
  dplyr::select( -n) %>% arrange(sample_id) %>% pull(sample_id) 
sample.data <- sample.data %>% 
  mutate(Project_Pd = ifelse( (sample_id %in% pd.ids) & Project == "SGLT2", "SGLT2_Pd", "Non_Pd")) %>%
  mutate(Treatment_Group = ifelse(Project_Pd == "SGLT2_Pd" & Pre_Post == "Pre" & Group == "NUT", "NUT_Pre", 
                                  ifelse(Project_Pd == "SGLT2_Pd" & Pre_Post == "Post" & Group == "NUT", "NUT_Post",
                                         ifelse(Project_Pd == "SGLT2_Pd" & Pre_Post == "Pre" & Group == "SGLT2i", "SGLT2i_Pre", 
                                                ifelse(Project_Pd == "SGLT2_Pd" & Pre_Post == "Post" & Group == "SGLT2i", "SGLT2i_Post", "Notreatment") ) )))

table(sample.data$Project_Pd)
table(sample.data$BA_Group, sample.data$Project_Pd)
table(sample.data$BA_Group, sample.data$Treatment_Group)

df_sg <- sample.data %>% dplyr::filter(
  Pre_Post %in% c("Pre", "Post") 
) %>% 
  group_by(patient_id) %>%
  mutate(n = n()) %>%
  dplyr::filter(n >= 2) %>%
  ungroup() %>%
  dplyr::select( -n)%>% 
  distinct(sample_id, .keep_all = TRUE)
df_Pre <- sample.data %>% dplyr::filter( 
  Age_Group == "Older" & 
    Diabetes == "No"&
    Pre_Post %in% c("Pre", "Notreatment")  &
    Prediabetes1 %in% c("PreD", "No")   &
    # sample_id != c("STARR_081")
    BMI >= 18.5
) 
sample.data <- sample.data %>% 
 mutate(Compare_Group1 = ifelse( (BA_Group %in% c("Older_Lean", "Older_Overweight", "Older_PreD_Obese", "Diabetes")) &
                                  (Treatment_Group %in% c("Notreatment", "NUT_Pre", "SGLT2i_Pre")) , TRUE, FALSE),
        Compare_Group2 = ifelse( (sample_id %in% df_sg$sample_id) , TRUE, FALSE),
        Compare_Group3 = ifelse( (BA_Group %in% c("Middle_Lean", "Middle_Overweight", "Older_Lean", "Older_Overweight")) &
                                   (Treatment_Group %in% c("Notreatment", "NUT_Pre", "SGLT2i_Pre")) , TRUE, FALSE)
        )
table(sample.data$BA_Group, sample.data$Compare_Group1)
table(sample.data$Treatment_Group, sample.data$Compare_Group2)
sample.data <- sample.data %>% 
  mutate(Select_Group = ifelse( (Compare_Group1 == TRUE) | (Compare_Group2 == TRUE) | (Compare_Group3 == TRUE), TRUE, FALSE))
table(sample.data$Treatment_Group, sample.data$Select_Group)
table(sample.data$BA_Group, sample.data$Select_Group)

table(sample.data$BA_Group, sample.data$Medicine.Status);
table(sample.data$BA_Group, sample.data$medication_1_name);
table(sample.data$BA_Group, sample.data$medication_2_name);
table(sample.data$BA_Group, sample.data$medication_3_name);
table(sample.data$BA_Group, sample.data$medication_4_name);
table(sample.data$BA_Group, sample.data$medication_5_name);
table(sample.data$BA_Group, sample.data$medication_6_name);
table(sample.data$BA_Group, sample.data$medication_7_name);
table(sample.data$BA_Group, sample.data$medication_8_name);


sample.data <- sample.data %>% # Remove columns from dataframe where ALL values are NA
  mutate(
    `HOMA_IR` = Glucose * `Insulin`/405,
    Weight_kg = Weight_pounds/2.20462
    # ,BMI = 703*Weight_pounds/(Height_inches^2)
    , Fat_Lean_Ratio = Fat..g./Lean..g.
  )
colnames(sample.data)
saveRDS(sample.data, paste0(Disk, Project.folder, "/",  Rds.folder, "/", "Samples.df.Rds")) ## This object is downloadable from Open Science Framework, with further annotations as subsequent scripts add to the file
xlsx::write.xlsx(sample.data, file = paste0(Disk, Project.folder, "/",  Rds.folder, "/", "Samples.df.xlsx"), sheetName = "SGLT2_STARR", row.names=FALSE)

# sample.data %>% dplyr::select(one_of(c("patient_id", "sample_id", "Pre_Post", "Group", "Project", "Project.Aggr", 
#                                 "Ethnicity", "Race", "A1c")))

# generate sample meta data ------------------------------------------------------
metadata.new <- basic.meta.data %>% dplyr::select(one_of(c("Dataset"))) %>% 
  tibble::rownames_to_column(var = "rowname") %>% 
  left_join(sample.data, by= "Dataset")%>% 
  tibble::column_to_rownames(var = "rowname") %>% 
  rename_with(make.names) 
# metadata.new$Treatment_Group <- factor(metadata.new$Treatment_Group, 
#                                        levels = c("Control", "Rapamycin"))
# print(table(metadata.new$Treatment_Group, metadata.new$Dataset))

saveRDS(metadata.new, paste0(Disk, Project.folder, "/",  Rds.folder, "/", "SAT_Sample.meta.data.Rds")) ## This object is downloadable from Open Science Framework, with further annotations as subsequent scripts add to the file

rm(sample.data, sample.data1, sample.data2, basic.meta.data, df_Pre, df_sg)
gc() #free up memory and report the memory usage.


rm(metadata.new)



