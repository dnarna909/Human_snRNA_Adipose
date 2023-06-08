df <- readRDS("/media/jianie/Extreme SSD/2022-09-01 STARR_SGLT2 Combine/Rds files_Aggr_all/Samples.df.Rds")
library(dplyr)

df1 <- df %>%
dplyr::filter(Treatment_Group %in% c("Notreatment", "NUT_Pre", "SGLT2i_Pre"),
BA_Group %in% c("Older_Lean", "Older_Overweight", "Older_PreD_Obese", "Diabetes")#
) %>%
mutate(BA_Group = factor(BA_Group, levels = c("Older_Lean", "Older_Overweight", "Older_PreD_Obese", "Diabetes"))) %>%
  select(one_of(c("sample_id" , "Group.Assignment.", "BA_Group" ,
                  "Project_Pd" , "Treatment_Group"  ,"Age" , "Gender"  ,"recruited_from" ,"A1c_Value" 
                  )))

df2 <- df %>% dplyr::filter(
  Pre_Post %in% c("Pre", "Post") 
) %>% 
  group_by(patient_id) %>%
  mutate(n = n()) %>%
  dplyr::filter(n >= 2) %>%
  ungroup() %>%
  dplyr::select( -n)  %>% 
  distinct(sample_id, .keep_all = TRUE)%>%
  select(one_of(c("sample_id" , "Group.Assignment.", "BA_Group" ,"patient_id",
                  "Project_Pd" , "Treatment_Group"  ,"Age" , "Gender"  ,"A1c_Value" 
  ))) %>% as.data.frame()

xlsx::write.xlsx(df1, file = paste0("/media/jianie/Extreme SSD/2022-09-01 STARR_SGLT2 Combine/Rds files_Aggr_all/", "Samples_NoTreatment.df.xlsx"), sheetName = "NoTreatment", row.names=FALSE)
xlsx::write.xlsx(df2, file = paste0("/media/jianie/Extreme SSD/2022-09-01 STARR_SGLT2 Combine/Rds files_Aggr_all/", "Samples_SGLT2.df.xlsx"), sheetName = "SGLT2", row.names=FALSE)
