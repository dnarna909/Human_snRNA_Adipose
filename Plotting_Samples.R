gc() #free up memory and report the memory usage.
memory.limit(size = 1e+13) # 'memory.limit()' is Windows-specific
source(paste0(Disk, Project.folder, "/","Project Parameters.R"), local = knitr::knit_global())

# import packages
# library("velocyto.R") # linux -specific


# Import file name -----------------------------------------------------------------------------------------------
# df <- readRDS(paste0(Disk, Project.folder, "/",  Rds.folder, "/", "Samples.df.Rds")) %>% 
#   dplyr::filter(!patient_id %in% c("76615")) 

df <- readRDS(paste0(Disk, Project.folder, "/",  Rds.folder, "/", "Samples.df.Rds")) %>% 
  mutate(
    `HOMA_IR` = Glucose * `Insulin`/405,
    Weight_kg = Weight_pounds/2.20462
    # ,BMI = 703*Weight_pounds/(Height_inches^2)
  ) %>% 
  group_by(sample_id) %>% dplyr::filter(row_number(sample_id) == 1)%>%
  ungroup() %>%
  as.data.frame() 

# df$Group <- factor(df$Group, levels = c("C", "D"))
# df$Pre_Post <- factor(df$Pre_Post, levels = c("Visit 5", "Visit 11"))
# df$Treatment_Group <- factor(df$Treatment_Group, levels = c("C_Visit 5", "C_Visit 11",
#                                                             "D_Visit 5", "D_Visit 11"))

dir.create(file.path(paste0(Disk, Project.folder, "/", figures.folder, "/", export.folder, ".plotting")), 
           showWarnings = FALSE)
dir1 <- paste0(Disk, Project.folder, "/", figures.folder, "/", export.folder, ".plotting", "/")

# Plot All by project --------------------------
library("RColorBrewer")
library("ggpubr")
library("ggplot2")
## Age_BMI  ------------
table(is.na(df$Age))
table(is.na(df$BMI))
table(is.na(df$SingleNuclei))
table(is.na(df$Gender))
table(is.na(substr(df$sample_id, 4,9)))

p <- ggplot(df, aes(x = Age, y = BMI)) + 
  geom_point(aes(color=Project, shape = Gender), size = ifelse(df$Diabetes == "Yes", 4, 2)) +
  labs(title="Patients distribution", y = "Body mass index", x = "Age (years)") +  
  ylim(min(df$BMI), max(df$BMI[which(df$BMI < Inf)])+1) + xlim(30, max(df$Age)+1) + 
  theme_classic() +
  theme(panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
        axis.title.x = element_text(color = "black", size = 13, face = "plain"),
        axis.title.y = element_text(color = "black", size = 13, face = "plain"),
        axis.text.x = element_text(color = "black", size = 11.5, face = "plain", angle = 0, hjust = 1),
        plot.title = element_text(color = "black", size = 14, face = "bold", hjust = 0.5),
        legend.text = element_text(color = "black", size = 10.5, face = "plain"),
        #legend.position = c(0.9, 0.8),
        legend.position = c(.11, .81),
        #legend.position = "none",
        #legend.justification = c("left", "bottom"),
        #legend.box.just = "left",
        legend.spacing.y = unit(0.15, 'cm') ,
        legend.key.size = unit(0.45, "cm"),
        legend.background = element_rect( fill = "grey98", color = "grey98", linewidth = 1)
  ) + 
  geom_vline(xintercept = 60, color = "black", linetype="dashed", size=0.2) +
  geom_vline(xintercept = 45, color = "black", linetype="dashed", size=0.2) +
  geom_hline(yintercept = 29.9, color = "black", linetype="dashed", size=0.2) +
  geom_hline(yintercept = 24.9, color = "black", linetype="dashed", size=0.2) + 
  geom_text(aes(color=Project),
            label= substr(df$sample_id, 4,9),
            nudge_x=0.5, nudge_y=0.5,
            check_overlap=F
  )
p
png(file=paste0(dir1, "All_PatientInfor_Age_BMI_distribution.png"), 
    width=8, height=8, units = "in", res = 300)
p
dev.off()

## A1c_Glucose -----
p2 <- ggplot(df , aes(x = `A1c_Value`, y = Glucose)) + 
  geom_point(aes(color=Project, shape = Gender), size = ifelse(df$Diabetes == "Yes", 4, 2)) +
  labs(title="Patients distribution", y = "Fasting Glucose (mg/dL)", x= "Hemoglobin A1C (%)") + # , y = "BMI", y = "BMI" 
  ylim(min(df$Glucose), max(df$Glucose)*1.02) + xlim(min(df$`A1c_Value`), max(df$`A1c_Value`)*1.02) + 
  theme_classic() +
  theme(panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
        axis.title.x = element_text(color = "black", size = 13, face = "plain"),
        axis.title.y = element_text(color = "black", size = 13, face = "plain"),
        axis.text.x = element_text(color = "black", size = 11.5, face = "plain", angle = 0, hjust = 1),
        plot.title = element_text(color = "black", size = 14, face = "bold", hjust = 0.5),
        legend.text = element_text(color = "black", size = 10.5, face = "plain"),
        #legend.position = c(0.9, 0.8),
        legend.position = c(.11, .812),
        #legend.position = "none",
        #legend.justification = c("left", "bottom"),
        #legend.box.just = "left",
        legend.spacing.y = unit(0.15, 'cm') ,
        legend.key.size = unit(0.45, "cm"),
        legend.background = element_rect( fill = "grey98", color = "grey98", size = 1)
  ) + 
  geom_vline(xintercept = 5.7, color = "black", linetype="dashed", size=0.2) +
  geom_vline(xintercept = 6.4, color = "black", linetype="dashed", size=0.2) +
  geom_hline(yintercept = 100, color = "black", linetype="dashed", size=0.2) +
  geom_hline(yintercept = 125, color = "black", linetype="dashed", size=0.2) + 
  geom_text(aes(color=Project), 
            label= substr(df$sample_id, 4,9),
            nudge_x=0.1, nudge_y=0.2, 
            check_overlap=F
  )
p2
png(file=paste0(dir1,  "All_PatientInfor_A1c_Glucose_distribution.png"), 
    width=8, height=8, units = "in", res = 300)
p2
dev.off()

## A1c and BMI ------------
p <- ggplot(df, aes(x = A1c_Value, y = BMI)) + 
  geom_point(aes(color=Project, shape = Gender), size = ifelse(df$Diabetes == "Yes", 4, 2)) +
  labs(title="Patients distribution", y = "Body mass index", x = "Hemoglobin A1C (%)") +  
  ylim(min(df$BMI), max(df$BMI[which(df$BMI < Inf)])*1.1) + xlim(min(df$A1c_Value), max(df$A1c_Value)*1.1) + 
  theme_classic() +
  theme(panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
        axis.title.x = element_text(color = "black", size = 13, face = "plain"),
        axis.title.y = element_text(color = "black", size = 13, face = "plain"),
        axis.text.x = element_text(color = "black", size = 11.5, face = "plain", angle = 0, hjust = 1),
        plot.title = element_text(color = "black", size = 14, face = "bold", hjust = 0.5),
        legend.text = element_text(color = "black", size = 10.5, face = "plain"),
        #legend.position = c(0.9, 0.8),
        legend.position = c(.11, .81),
        #legend.position = "none",
        #legend.justification = c("left", "bottom"),
        #legend.box.just = "left",
        legend.spacing.y = unit(0.15, 'cm') ,
        legend.key.size = unit(0.45, "cm"),
        legend.background = element_rect( fill = "grey98", color = "grey98", linewidth = 1)
  ) + 
  geom_vline(xintercept = 5.7, color = "black", linetype="dashed", size=0.2) +
  geom_vline(xintercept = 6.4, color = "black", linetype="dashed", size=0.2) +
  geom_hline(yintercept = 29.9, color = "black", linetype="dashed", size=0.2) +
  geom_hline(yintercept = 24.9, color = "black", linetype="dashed", size=0.2) + 
  geom_text(aes(color=Project),
            label= substr(df$sample_id, 4,9),
            nudge_x=0.1, nudge_y=0.1,
            check_overlap=F
  )
p
png(file=paste0(dir1, "All_PatientInfor_A1c_Value_BMI_distribution.png"), 
    width=8, height=8, units = "in", res = 300)
p
dev.off()


## Insulin, FFA -------------
p2 <- ggplot(df , aes(x = Insulin, y = FFA)) + 
  geom_point(aes(color=Project, shape = Gender), size = ifelse(df$Diabetes == "Yes", 4, 2)) +
  labs(title="Patients distribution", y = "Fasting FFA (mM)", x= "Fasting Inuslin (µIU/mL)") + # , y = "BMI", y = "BMI" 
  ylim(min(df$FFA[which(!is.na(df$FFA))])*0.8, max(df$FFA[which(!is.na(df$FFA))])*1.3) + 
  xlim(min(df$Insulin[which(!is.na(df$Insulin))])*0.8, max(df$Insulin[which(!is.na(df$Insulin))])*1.3) + 
  theme_classic() +
  theme(panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
        axis.title.x = element_text(color = "black", size = 13, face = "plain"),
        axis.title.y = element_text(color = "black", size = 13, face = "plain"),
        axis.text.x = element_text(color = "black", size = 11.5, face = "plain", angle = 0, hjust = 1),
        plot.title = element_text(color = "black", size = 14, face = "bold", hjust = 0.5),
        legend.text = element_text(color = "black", size = 10.5, face = "plain"),
        #legend.position = c(0.9, 0.8),
        legend.position = c(0.95, .812),
        #legend.position = "none",
        #legend.justification = c("left", "bottom"),
        #legend.box.just = "left",
        legend.spacing.y = unit(0.15, 'cm') ,
        legend.key.size = unit(0.45, "cm"),
        legend.background = element_rect( fill = "grey98", color = "grey98", size = 1)
  ) + 
  geom_vline(xintercept = 25, color = "black", linetype="dashed", size=0.2) +
  # geom_vline(xintercept = 6.4, color = "black", linetype="dashed", size=0.2) +
  # geom_hline(yintercept = 0.45, color = "black", linetype="dashed", size=0.2) +
  # geom_hline(yintercept = 0.6, color = "black", linetype="dashed", size=0.2) + 
  geom_text(aes(color=Project), 
            label= substr(df$sample_id, 4,9),
            nudge_x=0.01, nudge_y=0.02, 
            check_overlap=F
  )
p2
png(file=paste0(dir1,  "All_PatientInfor_Insulin_FFA_distribution.png"), 
    width=8, height=8, units = "in", res = 300)
p2
dev.off()

## `HOMA_IR` and Age ----------------
p <- ggplot(df, aes(x = Age, y = `HOMA_IR`)) + 
  geom_point(aes(color=Project, shape = Gender), size = ifelse(df$Diabetes == "Yes", 4, 2)) +
  labs(title="Patients distribution", y = "HOMA_IR", x = "Age (years)") +  
  ylim(min(df$`HOMA_IR`), max(df$`HOMA_IR`[which(df$`HOMA_IR` < Inf)])+1) + xlim(30, max(df$Age)+1) + 
  theme_classic() +
  theme(panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
        axis.title.x = element_text(color = "black", size = 13, face = "plain"),
        axis.title.y = element_text(color = "black", size = 13, face = "plain"),
        axis.text.x = element_text(color = "black", size = 11.5, face = "plain", angle = 0, hjust = 1),
        plot.title = element_text(color = "black", size = 14, face = "bold", hjust = 0.5),
        legend.text = element_text(color = "black", size = 10.5, face = "plain"),
        #legend.position = c(0.9, 0.8),
        legend.position = c(.11, .81),
        #legend.position = "none",
        #legend.justification = c("left", "bottom"),
        #legend.box.just = "left",
        legend.spacing.y = unit(0.15, 'cm') ,
        legend.key.size = unit(0.45, "cm"),
        legend.background = element_rect( fill = "grey98", color = "grey98", linewidth = 1)
  ) + 
  geom_vline(xintercept = 60, color = "black", linetype="dashed", size=0.2) +
  geom_vline(xintercept = 45, color = "black", linetype="dashed", size=0.2) +
  geom_hline(yintercept = quantile(df$`HOMA_IR`[which(df$`HOMA_IR` < Inf)],prob=1-(20/100) )[[1]], color = "black", linetype="dashed", size=0.2) +
  geom_hline(yintercept = quantile(df$`HOMA_IR`[which(df$`HOMA_IR` < Inf)],prob=1-(90/100) )[[1]], color = "black", linetype="dashed", size=0.2) + 
  geom_text(aes(color=Project),
            label= substr(df$sample_id, 4,9),
            nudge_x=0.5, nudge_y=0.5,
            check_overlap=F
  )
p
png(file=paste0(dir1, "All_PatientInfor_Age_HOMA_IR_distribution.png"), 
    width=8, height=8, units = "in", res = 300)
p
dev.off()

# Publication_Patients in Compare.Group1 and Compare.Group3 ------------
df1 <- df %>%
  dplyr::filter(Compare_Group3 == TRUE | 
                  Compare_Group1 == TRUE )

label.id = unique(df1%>%dplyr::filter(Diabetes == "Yes" | Prediabetes1 == c("Diabetes") )%>% 
                    pull(sample_id))

label.id2 = unique(df%>%dplyr::filter(Pre_Post %in% c("Post"))%>% pull(patient_id))

## Age and BMI ---------------------------
p <- ggplot(df1, aes(x = Age, y = BMI)) + 
  geom_point(aes(color=Gender, shape = Gender), size = ifelse(df1$Diabetes == "Yes" | df1$Prediabetes1 == c("Diabetes") , 4, 2) ) +
  labs(title="Patients distribution", y = "Body mass index", x = "Age (years)") +  
  # ylim(min(df1$BMI), max(df1$BMI[which(df1$BMI < Inf)])*1.1) + 
  ylim(18, 45) +
  xlim(min(df1$Age), max(df1$Age)) + 
  theme_classic() +
  theme(panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
        axis.title.x = element_text(color = "black", size = 13, face = "plain"),
        axis.title.y = element_text(color = "black", size = 13, face = "plain"),
        axis.text.x = element_text(color = "black", size = 11.5, face = "plain", angle = 0, hjust = 1),
        plot.title = element_text(color = "black", size = 14, face = "bold", hjust = 0.5),
        legend.text = element_text(color = "black", size = 10.5, face = "plain"),
        #legend.position = c(0.9, 0.8),
        legend.position = c(.14, .88),
        #legend.position = "none",
        #legend.justification = c("left", "bottom"),
        #legend.box.just = "left",
        legend.spacing.y = unit(0.15, 'cm') ,
        legend.key.size = unit(0.45, "cm"),
        legend.background = element_rect( fill = "grey98", color = "grey98", linewidth = 1)
  ) + 
  # geom_vline(xintercept = 60, color = "black", linetype="dashed", size=0.2) +
  # geom_vline(xintercept = 45, color = "black", linetype="dashed", size=0.2) +
  geom_hline(yintercept = 29.99, color = "black", linetype="dashed", size=0.2) +
  geom_hline(yintercept = 24.99, color = "black", linetype="dashed", size=0.2) 
p
png(file=paste0(dir1, "Publication_PatientInfor_Age_BMI_distribution.png"), 
    width=4, height=4, units = "in", res = 300)
p
dev.off()  

p <- ggplot(df1, aes(x = Age, y = BMI)) + 
  geom_point(aes(color=Gender, shape = Gender), size = ifelse(df1$Diabetes == "Yes" | df1$Prediabetes1 == c("Diabetes") , 4, 2) ) +
  geom_point(data = df1 %>%  dplyr::filter(patient_id %in% label.id2),
             pch=21, 
             size=3.8,
             colour="black") +
  labs(title="Patients distribution", y = "Body mass index", x = "Age (years)") +  
  # ylim(min(df1$BMI), max(df1$BMI[which(df1$BMI < Inf)])*1.1) + 
  ylim(18, 45) +
  xlim(min(df1$Age), max(df1$Age)) + 
  theme_classic() +
  theme(panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
        axis.title.x = element_text(color = "black", size = 13, face = "plain"),
        axis.title.y = element_text(color = "black", size = 13, face = "plain"),
        axis.text.x = element_text(color = "black", size = 11.5, face = "plain", angle = 0, hjust = 1),
        plot.title = element_text(color = "black", size = 14, face = "bold", hjust = 0.5),
        legend.text = element_text(color = "black", size = 10.5, face = "plain"),
        #legend.position = c(0.9, 0.8),
        legend.position = c(.14, .88),
        #legend.position = "none",
        #legend.justification = c("left", "bottom"),
        #legend.box.just = "left",
        legend.spacing.y = unit(0.15, 'cm') ,
        legend.key.size = unit(0.45, "cm"),
        legend.background = element_rect( fill = "grey98", color = "grey98", linewidth = 1)
  ) + 
  # geom_vline(xintercept = 60, color = "black", linetype="dashed", size=0.2) +
  # geom_vline(xintercept = 45, color = "black", linetype="dashed", size=0.2) +
  geom_hline(yintercept = 29.99, color = "black", linetype="dashed", size=0.2) +
  geom_hline(yintercept = 24.99, color = "black", linetype="dashed", size=0.2) 
p
png(file=paste0(dir1, "Publication_PatientInfor_Age_BMI_distribution_Circle.png"), 
    width=4, height=4, units = "in", res = 300)
p
dev.off()  

## A1c and BMI --------------------------------------
p <- ggplot(df1, aes(x = A1c_Value, y = BMI)) + 
  geom_point(aes(color=Gender, shape = Gender), size = ifelse(df1$Diabetes == "Yes" | df1$Prediabetes1 == c("Diabetes") , 4, 2)) +
  geom_point(data = df1 %>%  dplyr::filter(patient_id %in% label.id2),
             pch=21, 
             size=3.8,
             colour="black") +
  labs(title="Patients distribution", y = "Body mass index", x = "Hemoglobin A1C (%)") +  
  # ylim(min(df1$BMI), max(df1$BMI[which(df1$BMI < Inf)])*1.1) + 
  ylim(18,45) + 
  xlim(min(df1$A1c_Value), max(df1$A1c_Value)) + 
  theme_classic() +
  theme(panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
        axis.title.x = element_text(color = "black", size = 13, face = "plain"),
        axis.title.y = element_text(color = "black", size = 13, face = "plain"),
        axis.text.x = element_text(color = "black", size = 11.5, face = "plain", angle = 0, hjust = 1),
        plot.title = element_text(color = "black", size = 14, face = "bold", hjust = 0.5),
        legend.text = element_text(color = "black", size = 10.5, face = "plain"),
        #legend.position = c(0.9, 0.8),
        legend.position = c(.14, .88),
        #legend.position = "none",
        #legend.justification = c("left", "bottom"),
        #legend.box.just = "left",
        legend.spacing.y = unit(0.15, 'cm') ,
        legend.key.size = unit(0.45, "cm"),
        legend.background = element_rect( fill = "grey98", color = "grey98", linewidth = 1)
  ) + 
  geom_vline(xintercept = 5.69, color = "black", linetype="dashed", size=0.2) +
  geom_vline(xintercept = 6.4, color = "black", linetype="dashed", size=0.2) +
  geom_hline(yintercept = 29.99, color = "black", linetype="dashed", size=0.2) +
  geom_hline(yintercept = 24.99, color = "black", linetype="dashed", size=0.2) 
p
png(file=paste0(dir1, "Publication_PatientInfor_A1c_BMI_distribution_4Groups_Circle.png"), 
    width=4.2, height=4, units = "in", res = 300)
p
dev.off()

p <- ggplot(df1, aes(x = A1c_Value, y = BMI)) + 
  geom_point(aes(color=Gender, shape = Gender), size = ifelse(df1$Diabetes == "Yes" | df1$Prediabetes1 == c("Diabetes") , 4, 2)) +
  labs(title="Patients distribution", y = "Body mass index", x = "Hemoglobin A1C (%)") +  
  # ylim(min(df1$BMI), max(df1$BMI[which(df1$BMI < Inf)])*1.1) + 
  ylim(18,45) + 
  xlim(min(df1$A1c_Value), max(df1$A1c_Value)) + 
  theme_classic() +
  theme(panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
        axis.title.x = element_text(color = "black", size = 13, face = "plain"),
        axis.title.y = element_text(color = "black", size = 13, face = "plain"),
        axis.text.x = element_text(color = "black", size = 11.5, face = "plain", angle = 0, hjust = 1),
        plot.title = element_text(color = "black", size = 14, face = "bold", hjust = 0.5),
        legend.text = element_text(color = "black", size = 10.5, face = "plain"),
        #legend.position = c(0.9, 0.8),
        legend.position = c(.14, .88),
        #legend.position = "none",
        #legend.justification = c("left", "bottom"),
        #legend.box.just = "left",
        legend.spacing.y = unit(0.15, 'cm') ,
        legend.key.size = unit(0.45, "cm"),
        legend.background = element_rect( fill = "grey98", color = "grey98", linewidth = 1)
  ) + 
  geom_vline(xintercept = 5.69, color = "black", linetype="dashed", size=0.2) +
  geom_vline(xintercept = 6.4, color = "black", linetype="dashed", size=0.2) +
  geom_hline(yintercept = 29.99, color = "black", linetype="dashed", size=0.2) +
  geom_hline(yintercept = 24.99, color = "black", linetype="dashed", size=0.2) 
p
png(file=paste0(dir1, "Publication_PatientInfor_A1c_BMI_distribution_4Groups.png"), 
    width=4.2, height=4, units = "in", res = 300)
p
dev.off()

# Healthy_Patient_Plot samples for Compare_Group3 ------------------------
df1 <- df %>%
  dplyr::filter(Compare_Group3 == TRUE)
table(is.na(df1$Age))
table(is.na(df1$BMI))
table(is.na(df1$SingleNuclei))
table(is.na(df1$Gender))
table(is.na(substr(df1$sample_id, 4,9)))


## Age_BMI  ------------
p <- ggplot(df1, aes(x = Age, y = BMI)) + 
  geom_point(aes(color=Gender, shape = Gender), size = 2) +
  labs(title="Patients distribution", y = "Body mass index", x = "Age (years)") +  
  # ylim(min(df1$BMI), max(df1$BMI[which(df1$BMI < Inf)])*1.1) + 
  ylim(18, 45) +
  xlim(min(df1$Age), max(df1$Age)) + 
  theme_classic() +
  theme(panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
        axis.title.x = element_text(color = "black", size = 13, face = "plain"),
        axis.title.y = element_text(color = "black", size = 13, face = "plain"),
        axis.text.x = element_text(color = "black", size = 11.5, face = "plain", angle = 0, hjust = 1),
        plot.title = element_text(color = "black", size = 14, face = "bold", hjust = 0.5),
        legend.text = element_text(color = "black", size = 10.5, face = "plain"),
        #legend.position = c(0.9, 0.8),
        legend.position = c(.14, .88),
        #legend.position = "none",
        #legend.justification = c("left", "bottom"),
        #legend.box.just = "left",
        legend.spacing.y = unit(0.15, 'cm') ,
        legend.key.size = unit(0.45, "cm"),
        legend.background = element_rect( fill = "grey98", color = "grey98", linewidth = 1)
  ) + 
  # geom_vline(xintercept = 60, color = "black", linetype="dashed", size=0.2) +
  # geom_vline(xintercept = 45, color = "black", linetype="dashed", size=0.2) +
  geom_hline(yintercept = 29.99, color = "black", linetype="dashed", size=0.2) +
  geom_hline(yintercept = 24.99, color = "black", linetype="dashed", size=0.2) 
p
png(file=paste0(dir1, "Healthy_PatientInfor_Age_BMI_distribution.png"), 
    width=4, height=4, units = "in", res = 300)
p
dev.off()  

## A1c_Glucose -----
p2 <- ggplot(df1 , aes(x = `A1c_Value`, y = Glucose)) + 
  geom_point(aes(color=Gender, shape = Gender), size = ifelse(df1$Diabetes == "Yes", 4, 2)) +
  labs(title="Patients distribution", y = "Fasting Glucose (mg/dL)", x= "Hemoglobin A1C (%)") + # , y = "BMI", y = "BMI" 
  ylim(min(df1$Glucose), max(df1$Glucose)*1.02) + xlim(min(df1$`A1c_Value`), max(df1$`A1c_Value`)*1.02) + 
  theme_classic() +
  theme(panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
        axis.title.x = element_text(color = "black", size = 13, face = "plain"),
        axis.title.y = element_text(color = "black", size = 13, face = "plain"),
        axis.text.x = element_text(color = "black", size = 11.5, face = "plain", angle = 0, hjust = 1),
        plot.title = element_text(color = "black", size = 14, face = "bold", hjust = 0.5),
        legend.text = element_text(color = "black", size = 10.5, face = "plain"),
        #legend.position = c(0.9, 0.8),
        legend.position = c(.11, .812),
        #legend.position = "none",
        #legend.justification = c("left", "bottom"),
        #legend.box.just = "left",
        legend.spacing.y = unit(0.15, 'cm') ,
        legend.key.size = unit(0.45, "cm"),
        legend.background = element_rect( fill = "grey98", color = "grey98", size = 1)
  ) + 
  geom_vline(xintercept = 5.7, color = "black", linetype="dashed", size=0.2) +
  geom_vline(xintercept = 6.4, color = "black", linetype="dashed", size=0.2) +
  geom_hline(yintercept = 100, color = "black", linetype="dashed", size=0.2) +
  geom_hline(yintercept = 125, color = "black", linetype="dashed", size=0.2) 
# + 
#   geom_text(aes(color=Gender), 
#             label= substr(df1$sample_id, 4,9),
#             nudge_x=0.1, nudge_y=0.2, 
#             check_overlap=F
#   )
p2
png(file=paste0(dir1,  "Healthy_PatientInfor_A1c_Glucose_distribution.png"), 
    width=4, height=4, units = "in", res = 300)
p2
dev.off()

## A1c and BMI ------------
p <- ggplot(df1, aes(x = A1c_Value, y = BMI)) + 
  geom_point(aes(color=Gender, shape = Gender), size = ifelse(df1$Diabetes == "Yes", 4, 2)) +
  labs(title="Patients distribution", y = "Body mass index", x = "Hemoglobin A1C (%)") +  
  ylim(min(df1$BMI), max(df1$BMI[which(df1$BMI < Inf)])*1.1) + xlim(min(df1$A1c_Value), max(df1$A1c_Value)*1.1) + 
  theme_classic() +
  theme(panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
        axis.title.x = element_text(color = "black", size = 13, face = "plain"),
        axis.title.y = element_text(color = "black", size = 13, face = "plain"),
        axis.text.x = element_text(color = "black", size = 11.5, face = "plain", angle = 0, hjust = 1),
        plot.title = element_text(color = "black", size = 14, face = "bold", hjust = 0.5),
        legend.text = element_text(color = "black", size = 10.5, face = "plain"),
        #legend.position = c(0.9, 0.8),
        legend.position = c(.86, .88),
        #legend.position = "none",
        #legend.justification = c("left", "bottom"),
        #legend.box.just = "left",
        legend.spacing.y = unit(0.15, 'cm') ,
        legend.key.size = unit(0.45, "cm"),
        legend.background = element_rect( fill = "grey98", color = "grey98", linewidth = 1)
  ) + 
  geom_vline(xintercept = 5.7, color = "black", linetype="dashed", size=0.2) +
  geom_vline(xintercept = 6.4, color = "black", linetype="dashed", size=0.2) +
  geom_hline(yintercept = 29.9, color = "black", linetype="dashed", size=0.2) +
  geom_hline(yintercept = 24.9, color = "black", linetype="dashed", size=0.2) 
# + 
#   geom_text(aes(color=Gender),
#             label= substr(df1$sample_id, 4,9),
#             nudge_x=0.1, nudge_y=0.1,
#             check_overlap=F
#   )
p
png(file=paste0(dir1, "Healthy_PatientInfor_A1c_Value_BMI_distribution.png"), 
    width=4, height=4, units = "in", res = 300)
p
dev.off()


## Insulin, FFA -------------
p2 <- ggplot(df1 , aes(x = Insulin, y = FFA)) + 
  geom_point(aes(color=Gender, shape = Gender), size = ifelse(df1$Diabetes == "Yes", 4, 2)) +
  labs(title="Patients distribution", y = "Fasting FFA (mM)", x= "Fasting Inuslin (µIU/mL)") + # , y = "BMI", y = "BMI" 
  ylim(min(df1$FFA[which(!is.na(df1$FFA))])*0.8, max(df1$FFA[which(!is.na(df1$FFA))])*1.3) + 
  xlim(min(df1$Insulin[which(!is.na(df1$Insulin))])*0.8, max(df1$Insulin[which(!is.na(df1$Insulin))])*1.3) + 
  theme_classic() +
  theme(panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
        axis.title.x = element_text(color = "black", size = 13, face = "plain"),
        axis.title.y = element_text(color = "black", size = 13, face = "plain"),
        axis.text.x = element_text(color = "black", size = 11.5, face = "plain", angle = 0, hjust = 1),
        plot.title = element_text(color = "black", size = 14, face = "bold", hjust = 0.5),
        legend.text = element_text(color = "black", size = 10.5, face = "plain"),
        #legend.position = c(0.9, 0.8),
        legend.position = c(0.95, .812),
        #legend.position = "none",
        #legend.justification = c("left", "bottom"),
        #legend.box.just = "left",
        legend.spacing.y = unit(0.15, 'cm') ,
        legend.key.size = unit(0.45, "cm"),
        legend.background = element_rect( fill = "grey98", color = "grey98", size = 1)
  ) + 
  geom_vline(xintercept = 25, color = "black", linetype="dashed", size=0.2) # +
# geom_vline(xintercept = 6.4, color = "black", linetype="dashed", size=0.2) +
# geom_hline(yintercept = 0.45, color = "black", linetype="dashed", size=0.2) +
# geom_hline(yintercept = 0.6, color = "black", linetype="dashed", size=0.2) 
# + 
# geom_text(aes(color=Gender), 
#           label= substr(df1$sample_id, 4,9),
#           nudge_x=0.01, nudge_y=0.02, 
#           check_overlap=F
# )
p2
png(file=paste0(dir1,  "Healthy_PatientInfor_Insulin_FFA_distribution.png"), 
    width=4, height=4, units = "in", res = 300)
p2
dev.off()

## HOMA_IR and Age ----------------
p <- ggplot(df1, aes(x = Age, y = `HOMA_IR`)) + 
  geom_point(aes(color=Gender, shape = Gender), size = ifelse(df1$Diabetes == "Yes", 4, 2)) +
  labs(title="Patients distribution", y = "HOMA_IR", x = "Age (years)") +  
  ylim(min(df1$`HOMA_IR`), max(df1$`HOMA_IR`[which(df1$`HOMA_IR` < Inf)])+1) + xlim(45, max(df1$Age)+1) + 
  theme_classic() +
  theme(panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
        axis.title.x = element_text(color = "black", size = 13, face = "plain"),
        axis.title.y = element_text(color = "black", size = 13, face = "plain"),
        axis.text.x = element_text(color = "black", size = 11.5, face = "plain", angle = 0, hjust = 1),
        plot.title = element_text(color = "black", size = 14, face = "bold", hjust = 0.5),
        legend.text = element_text(color = "black", size = 10.5, face = "plain"),
        #legend.position = c(0.9, 0.8),
        legend.position = c(.11, .81),
        #legend.position = "none",
        #legend.justification = c("left", "bottom"),
        #legend.box.just = "left",
        legend.spacing.y = unit(0.15, 'cm') ,
        legend.key.size = unit(0.45, "cm"),
        legend.background = element_rect( fill = "grey98", color = "grey98", linewidth = 1)
  ) + 
  geom_vline(xintercept = 60, color = "black", linetype="dashed", size=0.2) +
  geom_vline(xintercept = 45, color = "black", linetype="dashed", size=0.2) +
  geom_hline(yintercept = quantile(df1$`HOMA_IR`[which(df1$`HOMA_IR` < Inf)],prob=1-(20/100) )[[1]], color = "black", linetype="dashed", size=0.2) +
  geom_hline(yintercept = quantile(df1$`HOMA_IR`[which(df1$`HOMA_IR` < Inf)],prob=1-(90/100) )[[1]], color = "black", linetype="dashed", size=0.2) 
# + 
#   geom_text(aes(color=Gender),
#             label= substr(df1$sample_id, 4,9),
#             nudge_x=0.5, nudge_y=0.5,
#             check_overlap=F
#   )
p
png(file=paste0(dir1, "Healthy_PatientInfor_Age_HOMA_IR_distribution.png"), 
    width=4, height=4, units = "in", res = 300)
p
dev.off()

# 4 Groups_Plot Prediabetes and diabetes: STARR (OLD, non Diabetes, diabetes) and SGLT2_V5 --------------------------
df_Pre <- df %>% dplyr::filter( 
  Age_Group == "Older" & 
    Diabetes == "No"&
    Pre_Post %in% c("Pre", "Notreatment")  &
    Prediabetes1 %in% c("PreD", "No")   &
    # sample_id != c("STARR_081")
    BMI >= 18.5
) 
df_Pre <- df %>% dplyr::filter( 
  Compare_Group1 == TRUE  &
    Age_Group == "Older"
) 
label.id = unique(df_Pre%>%dplyr::filter(Diabetes == "Yes" | Prediabetes1 == c("Diabetes") )%>% 
                    pull(sample_id))

label.id2 = unique(df%>%dplyr::filter(Pre_Post %in% c("Post"))%>% pull(patient_id))

## A1C_BMI  ------------
p <- ggplot(df_Pre, aes(x = A1c_Value, y = BMI)) + 
  geom_point(aes(color=Gender, shape = Gender), size = ifelse(df_Pre$Diabetes == "Yes" | df_Pre$Prediabetes1 == c("Diabetes") , 4, 2)) +
  labs(title="Patients distribution", y = "Body mass index", x = "Hemoglobin A1C (%)") +  
  ylim(min(df_Pre$BMI), max(df_Pre$BMI[which(df_Pre$BMI < Inf)])*1.1) + 
  xlim(min(df_Pre$A1c_Value), max(df_Pre$A1c_Value)) + 
  theme_classic() +
  theme(panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
        axis.title.x = element_text(color = "black", size = 13, face = "plain"),
        axis.title.y = element_text(color = "black", size = 13, face = "plain"),
        axis.text.x = element_text(color = "black", size = 11.5, face = "plain", angle = 0, hjust = 1),
        plot.title = element_text(color = "black", size = 14, face = "bold", hjust = 0.5),
        legend.text = element_text(color = "black", size = 10.5, face = "plain"),
        #legend.position = c(0.9, 0.8),
        legend.position = c(.14, .88),
        #legend.position = "none",
        #legend.justification = c("left", "bottom"),
        #legend.box.just = "left",
        legend.spacing.y = unit(0.15, 'cm') ,
        legend.key.size = unit(0.45, "cm"),
        legend.background = element_rect( fill = "grey98", color = "grey98", linewidth = 1)
  ) + 
  geom_vline(xintercept = 5.69, color = "black", linetype="dashed", size=0.2) +
  geom_vline(xintercept = 6.4, color = "black", linetype="dashed", size=0.2) +
  geom_hline(yintercept = 29.99, color = "black", linetype="dashed", size=0.2) +
  geom_hline(yintercept = 24.99, color = "black", linetype="dashed", size=0.2) 
p
png(file=paste0(dir1, "Basal_PatientInfor_A1c_BMI_distribution_4Groups.png"), 
    width=4.2, height=4, units = "in", res = 300)
p
dev.off()

p <- ggplot(df_Pre, aes(x = A1c_Value, y = BMI)) + 
  geom_point(aes(color=Gender, shape = Gender), size = ifelse(df_Pre$Diabetes == "Yes" | df_Pre$Prediabetes1 == c("Diabetes") , 4, 2)) +
  geom_point(data = df_Pre %>%  dplyr::filter(patient_id %in% label.id2),
             pch=21, 
             size=3.8,
             colour="black") +
  labs(title="Patients distribution", y = "Body mass index", x = "Hemoglobin A1C (%)") +  
  ylim(min(df_Pre$BMI), max(df_Pre$BMI[which(df_Pre$BMI < Inf)])*1.1) + 
  xlim(min(df_Pre$A1c_Value), max(df_Pre$A1c_Value)) + 
  theme_classic() +
  theme(panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
        axis.title.x = element_text(color = "black", size = 13, face = "plain"),
        axis.title.y = element_text(color = "black", size = 13, face = "plain"),
        axis.text.x = element_text(color = "black", size = 11.5, face = "plain", angle = 0, hjust = 1),
        plot.title = element_text(color = "black", size = 14, face = "bold", hjust = 0.5),
        legend.text = element_text(color = "black", size = 10.5, face = "plain"),
        #legend.position = c(0.9, 0.8),
        legend.position = c(.14, .88),
        #legend.position = "none",
        #legend.justification = c("left", "bottom"),
        #legend.box.just = "left",
        legend.spacing.y = unit(0.15, 'cm') ,
        legend.key.size = unit(0.45, "cm"),
        legend.background = element_rect( fill = "grey98", color = "grey98", linewidth = 1)
  ) + 
  geom_vline(xintercept = 5.69, color = "black", linetype="dashed", size=0.2) +
  geom_vline(xintercept = 6.4, color = "black", linetype="dashed", size=0.2) +
  geom_hline(yintercept = 29.99, color = "black", linetype="dashed", size=0.2) +
  geom_hline(yintercept = 24.99, color = "black", linetype="dashed", size=0.2) 
p
png(file=paste0(dir1, "Basal_PatientInfor_A1c_BMI_distribution_4Groups_Circle.png"), 
    width=4.2, height=4, units = "in", res = 300)
p
dev.off()

## Age_BMI  ------------
p <- ggplot(df_Pre, aes(x = Age, y = BMI)) + 
  geom_point(aes(color=Gender, shape = Gender), size = ifelse(df_Pre$Diabetes == "Yes" | df_Pre$Prediabetes1 == c("Diabetes") , 4, 2)) +
  labs(title="Patients distribution", y = "Body mass index", x = "Age (years)") +  
  # ylim(min(df_Pre$BMI), max(df_Pre$BMI[which(df_Pre$BMI < Inf)])*1.1) + 
  ylim(18, 45) +
  xlim(min(df_Pre$Age), max(df_Pre$Age)) + 
  theme_classic() +
  theme(panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
        axis.title.x = element_text(color = "black", size = 13, face = "plain"),
        axis.title.y = element_text(color = "black", size = 13, face = "plain"),
        axis.text.x = element_text(color = "black", size = 11.5, face = "plain", angle = 0, hjust = 1),
        plot.title = element_text(color = "black", size = 14, face = "bold", hjust = 0.5),
        legend.text = element_text(color = "black", size = 10.5, face = "plain"),
        #legend.position = c(0.9, 0.8),
        legend.position = c(.14, .88),
        #legend.position = "none",
        #legend.justification = c("left", "bottom"),
        #legend.box.just = "left",
        legend.spacing.y = unit(0.15, 'cm') ,
        legend.key.size = unit(0.45, "cm"),
        legend.background = element_rect( fill = "grey98", color = "grey98", linewidth = 1)
  ) + 
  geom_vline(xintercept = 60, color = "black", linetype="dashed", size=0.2) +
  geom_vline(xintercept = 45, color = "black", linetype="dashed", size=0.2) +
  geom_hline(yintercept = 29.99, color = "black", linetype="dashed", size=0.2) +
  geom_hline(yintercept = 24.99, color = "black", linetype="dashed", size=0.2) 
p
png(file=paste0(dir1, "Basal_PatientInfor_Age_BMI_distribution_4Groups.png"), 
    width=4, height=4, units = "in", res = 300)
p
dev.off()

p <- ggplot(df_Pre, aes(x = Age, y = BMI)) + 
  geom_point(aes(color=Gender, shape = Gender), size = ifelse(df_Pre$Diabetes == "Yes" | df_Pre$Prediabetes1 == c("Diabetes") , 4, 2)) +
  geom_point(data = df_Pre %>%  dplyr::filter(patient_id %in% label.id2),
             pch=21, 
             size=3.8,
             colour="black") +
  labs(title="Patients distribution", y = "Body mass index", x = "Age (years)") +  
  # ylim(min(df_Pre$BMI), max(df_Pre$BMI[which(df_Pre$BMI < Inf)])*1.1) + 
  ylim(18, 45) +
  xlim(min(df_Pre$Age), max(df_Pre$Age)) + 
  theme_classic() +
  theme(panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
        axis.title.x = element_text(color = "black", size = 13, face = "plain"),
        axis.title.y = element_text(color = "black", size = 13, face = "plain"),
        axis.text.x = element_text(color = "black", size = 11.5, face = "plain", angle = 0, hjust = 1),
        plot.title = element_text(color = "black", size = 14, face = "bold", hjust = 0.5),
        legend.text = element_text(color = "black", size = 10.5, face = "plain"),
        #legend.position = c(0.9, 0.8),
        legend.position = c(.14, .88),
        #legend.position = "none",
        #legend.justification = c("left", "bottom"),
        #legend.box.just = "left",
        legend.spacing.y = unit(0.15, 'cm') ,
        legend.key.size = unit(0.45, "cm"),
        legend.background = element_rect( fill = "grey98", color = "grey98", linewidth = 1)
  ) + 
  geom_vline(xintercept = 60, color = "black", linetype="dashed", size=0.2) +
  geom_vline(xintercept = 45, color = "black", linetype="dashed", size=0.2) +
  geom_hline(yintercept = 29.99, color = "black", linetype="dashed", size=0.2) +
  geom_hline(yintercept = 24.99, color = "black", linetype="dashed", size=0.2) 
p
png(file=paste0(dir1, "Basal_PatientInfor_Age_BMI_distribution_4Groups_Circle.png"), 
    width=4, height=4, units = "in", res = 300)
p
dev.off()

# 4 Groups_plot tables in four group ------------------ 
cols.list <- list(
  Others = c( "Blood.Urea.Nitrogen..BUN.", "Creatinine", "Sodium", "Potassium", 
              "Calcium", "Protein..Total", "Albumin", "Alkaline.Phosphatase", "AST", "ALT",  
              "Red.Blood.Cell...RBC", "Hemoglobin", "Hematocrit", "MCV", "MCH", 
              "Platelets", "White.Blood.Cell...WBC", "HDL", "LDL", "Weight_pounds"
  ), #  
  # miss = c(),#"Basophils",  "Chloride", "RDW", "MCHC","Z.Score", "Granulocytes", "Manual.eGFR.calculation","PT", "INR", "PTT", 
  # "Temperature..F.", "CO2",
  # "Monocyte", "Neutrophils", "Eosinophils" ,"TSH", "Urine.pH",  "MoCA.Score" , "VLDL",
  # "Tissue...Fat.", "Region...Fat.", "Fat..g.","Lean..g.", "Tissue..g.",  "BMC..g.", "Fat.Free..g.", "T.Mass..kg.",
  #           "Heart.Rate..bpm.", "Respirations.per.minute",
  Basal = c(
    "Age","A1c_Value", "Glucose", "Weight_kg","Height_inches", "BMI", 
    "Blood.Pressure.Systolic", "Blood.Pressure.Diastolic", 
    "HOMA_IR",  "Insulin", "FFA",
    "Cholesterol..Total",  "Triglycerides",  "Bilirubin..Total"),
  Plasma = c( "TNF.alpha", "IL.6", "IL.8.CXCL8", "VEGF", "IL.7", "CCL3.MIP.1alpha", 
              "CCL4.MIP.1", "RAGE.AGER", "SOST.Sclerostin", "ADAMTS13", "Osteropontin.OPN", "ICAM.1.CD54", 
              "IL.15", "GDF.15", "TNF.RI.TNFRSF1A", "Fas.TNFRSF6.CD95", "GDNF", "ActivinA" )
)

compare.groups <- list(
  compare1 = c("Older_Lean",  "Older_Overweight", "Older_PreD_Obese", "Diabetes" )
)
for (cc in names(cols.list)){
  for (ll in names(compare.groups)) {
    rm(demtable)
    demtable <- df %>%
      dplyr::filter( Compare_Group1 == TRUE & 
                       Pre_Post %in% c("Pre", "Notreatment") )  %>%
      dplyr::filter(BA_Group %in% compare.groups[[ll]]) %>%
      mutate(BA_Group = factor(BA_Group, levels = c(compare.groups[[ll]][1],  compare.groups[[ll]][2],  
                                                    compare.groups[[ll]][3],  compare.groups[[ll]][4]) )) %>%
      dplyr::select(one_of(c("sample_id", "BA_Group", cols.list[[cc]]))) %>%
      dplyr::filter(complete.cases(.)) %>% # delete missing values
      tbl_summary(by = "BA_Group", include = -c("sample_id"), 
                  type = everything() ~ "continuous",
                  statistic = list(all_continuous() ~ "{mean} ({sd})",
                                   all_categorical() ~ "{n} / {N} ({p}%)"),
                  digits = all_continuous() ~ 2,
                  # label = grade ~ "Tumor Grade", # change names
                  missing_text = "(Missing)"
      ) %>%
      add_p(all_continuous() ~ "aov")%>% # "oneway.test" # "ancova"
      add_n() %>%
      # add_overall(last = FALSE) %>%
      add_significance_stars() %>%
      modify_header(label ~ "**Variable**") %>%
      modify_spanning_header(c("stat_1", "stat_2", "stat_3", "stat_4") ~ "**Group**") %>%
      modify_footnote(
        all_stat_cols() ~ "Mean (SD)" #" or Frequency (%)"
      ) %>%
      modify_caption("**Table 2. Type of Assessments in Four Groups**") %>%
      bold_labels()
    demtable
    gt::gtsave(as_gt(demtable),file=paste0(dir1, "STARR_", cc, "_", "4Groups", "_Group_Event2.png"), vwidth = 3000, vheight = 1500)
  }
}

for (cc in names(cols.list)){
  for (ll in names(compare.groups)) {
    rm(demtable)
    demtable <- df %>%
      dplyr::filter( Compare_Group1 == TRUE &
                       Pre_Post %in% c("Pre", "Notreatment") )  %>%
      dplyr::filter(BA_Group %in% compare.groups[[ll]]) %>%
      mutate(BA_Group = factor(BA_Group, levels = c(compare.groups[[ll]][1],  compare.groups[[ll]][2],  
                                                    compare.groups[[ll]][3],  compare.groups[[ll]][4]))) %>%
      dplyr::select(one_of(c("sample_id", "BA_Group","Gender", cols.list[[cc]]))) %>%
      dplyr::filter(complete.cases(.)) %>% # delete missing values
      tbl_strata(
        strata = Gender,
        .tbl_fun =
          ~ .x %>%
          tbl_summary(by = BA_Group, include = -c(sample_id),
                      type = everything() ~ "continuous",
                      statistic = list(all_continuous() ~ "{mean} ({sd})",
                                       all_categorical() ~ "{n} / {N} ({p}%)"),
                      digits = all_continuous() ~ 2,
                      # label = grade ~ "Tumor Grade", # change names
                      missing_text = "(Missing)") %>%
          add_p(all_continuous() ~ "aov")%>% # "oneway.test" # "ancova"
          # add_n() %>%
          # add_overall(last = FALSE) %>%
          add_significance_stars() %>%
          modify_header(label ~ "**Variable**") %>%
          modify_spanning_header(c("stat_1", "stat_2", "stat_3", "stat_4") ~ "**Timeline**") %>%
          modify_footnote(
            all_stat_cols() ~ "Mean (SD)" #" or Frequency (%)"
          ) %>%
          modify_caption("**Table 2. Type of Assessments in Three Groups by Gender**") %>%
          bold_labels(),
        .header = "**{strata}**, N = {n}"
      )
    demtable
    gt::gtsave(as_gt(demtable),file=paste0(dir1, "STARR_Gender_", cc, "_", "4Groups", "_Group_Event2.png"), vwidth = 3500, vheight = 1500)
  }
}


for (cc in names(cols.list)){
  for (ll in names(compare.groups)) {
    rm(demtable)
    demtable <- df %>%
      dplyr::filter( Compare_Group1 == TRUE & 
                       Pre_Post %in% c("Pre", "Notreatment") &
                       Gender == "Female")  %>%
      dplyr::filter(BA_Group %in% compare.groups[[ll]]) %>%
      mutate(BA_Group = factor(BA_Group, levels = c(compare.groups[[ll]][1],  compare.groups[[ll]][2],  
                                                    compare.groups[[ll]][3],  compare.groups[[ll]][4]) )) %>%
      dplyr::select(one_of(c("sample_id", "BA_Group", cols.list[[cc]]))) %>%
      dplyr::filter(complete.cases(.)) %>% # delete missing values
      tbl_summary(by = "BA_Group", include = -c("sample_id"), 
                  type = everything() ~ "continuous",
                  statistic = list(all_continuous() ~ "{mean} ({sd})",
                                   all_categorical() ~ "{n} / {N} ({p}%)"),
                  digits = all_continuous() ~ 2,
                  # label = grade ~ "Tumor Grade", # change names
                  missing_text = "(Missing)"
      ) %>%
      add_p(all_continuous() ~ "aov")%>% # "oneway.test" # "ancova"
      add_n() %>%
      # add_overall(last = FALSE) %>%
      add_significance_stars() %>%
      modify_header(label ~ "**Variable**") %>%
      modify_spanning_header(c("stat_1", "stat_2", "stat_3", "stat_4") ~ "**Group**") %>%
      modify_footnote(
        all_stat_cols() ~ "Mean (SD)" #" or Frequency (%)"
      ) %>%
      modify_caption("**Table 2. Type of Assessments in Four Groups in Male**") %>%
      bold_labels()
    demtable
    gt::gtsave(as_gt(demtable),file=paste0(dir1, "STARR_Gender_Female_", cc, "_", "4Groups", "_Group_Event2.png"), vwidth = 3000, vheight = 1500)
  }
}

for (cc in names(cols.list)){
  for (ll in names(compare.groups)) {
    rm(demtable)
    demtable <- df %>%
      dplyr::filter( Compare_Group1 == TRUE & 
                       Pre_Post %in% c("Pre", "Notreatment") &
                       Gender == "Male")  %>%
      dplyr::filter(BA_Group %in% compare.groups[[ll]]) %>%
      mutate(BA_Group = factor(BA_Group, levels = c(compare.groups[[ll]][1],  compare.groups[[ll]][2],  
                                                    compare.groups[[ll]][3],  compare.groups[[ll]][4]) )) %>%
      dplyr::select(one_of(c("sample_id", "BA_Group", cols.list[[cc]]))) %>%
      dplyr::filter(complete.cases(.)) %>% # delete missing values
      tbl_summary(by = "BA_Group", include = -c("sample_id"), 
                  type = everything() ~ "continuous",
                  statistic = list(all_continuous() ~ "{mean} ({sd})",
                                   all_categorical() ~ "{n} / {N} ({p}%)"),
                  digits = all_continuous() ~ 2,
                  # label = grade ~ "Tumor Grade", # change names
                  missing_text = "(Missing)"
      ) %>%
      add_p(all_continuous() ~ "aov")%>% # "oneway.test" # "ancova"
      add_n() %>%
      # add_overall(last = FALSE) %>%
      add_significance_stars() %>%
      modify_header(label ~ "**Variable**") %>%
      modify_spanning_header(c("stat_1", "stat_2", "stat_3", "stat_4") ~ "**Group**") %>%
      modify_footnote(
        all_stat_cols() ~ "Mean (SD)" #" or Frequency (%)"
      ) %>%
      modify_caption("**Table 2. Type of Assessments in Four Groups in Male**") %>%
      bold_labels()
    demtable
    gt::gtsave(as_gt(demtable),file=paste0(dir1, "STARR_Gender_Male_", cc, "_", "4Groups", "_Group_Event2.png"), vwidth = 3000, vheight = 1500)
  }
}

# 4 Groups_plot numeric.variables by barplot -------------------
cat(paste(shQuote(colnames(df)[sapply(df, is.numeric)], type = "cmd"), collapse = ", "))
numeric.variables = c("White.Blood.Cell...WBC", "Red.Blood.Cell...RBC", "Hemoglobin", "Hematocrit", "MCV", "MCH", 
                      "Platelets", "Glucose", "Blood.Urea.Nitrogen..BUN.", "Creatinine", "Sodium", "Potassium", 
                      "Calcium", "Protein..Total", "Albumin", "Alkaline.Phosphatase", "Bilirubin..Total", "AST", "ALT", 
                      "Cholesterol..Total", "Triglycerides", "HDL", "LDL", "VLDL", "A1c_Value", 
                      "Weight_pounds", "Height_inches", "BMI", "Blood.Pressure.Systolic", "Blood.Pressure.Diastolic", 
                      "Insulin", "FFA", "TNF.alpha", "IL.6", "IL.8.CXCL8", "VEGF", "IL.7", "CCL3.MIP.1alpha", 
                      "CCL4.MIP.1", "RAGE.AGER", "SOST.Sclerostin", "ADAMTS13", "Osteropontin.OPN", "ICAM.1.CD54", 
                      "IL.15", "GDF.15", "TNF.RI.TNFRSF1A", "Fas.TNFRSF6.CD95", "GDNF", "ActivinA", 
                      "Heart.Rate..bpm.", "Respirations.per.minute", "Temperature..F.", 
                      "Monocyte",  "Neutrophils", "Eosinophils", "Basophils", 
                      "MCHC", "RDW", "Chloride", "CO2", "Manual.eGFR.calculation",
                      "MoCA.Score",  
                      "Tissue...Fat.", "T.Mass..kg.", "Region...Fat.", "Tissue..g.", "Fat..g.", 
                      "Lean..g.", "BMC..g.", "Fat.Free..g.", "HOMA_IR", "Weight_kg")
compare.groups <- list(
  # compare1 = c("Older_Lean",  "Older_Overweight", "Older_PreD_Obese" ),
  compare2= c("Older_Lean",  "Older_Overweight", "Older_PreD_Obese", "Diabetes" )
  
)

result.t.test <- data.frame()
result.anova <- data.frame()
result.summary <- data.frame()

result.anova.Gender <- data.frame()
result.t.test.Gender <- data.frame()
result.summary.Gender <- data.frame()

result.t.test.sg <- data.frame()
result.anova.sg <- data.frame()
result.summary.sg <- data.frame()

library(rstatix)

for (ll in names(compare.groups)) {
  # ll = names(compare.groups)[1]
  # comparisons = expand.grid(a = compare.groups[[ll]], b = compare.groups[[ll]]) %>% 
  #   dplyr::filter(a != b, ) %>% 
  #   arrange(a)
  comparisons = combn(compare.groups[[ll]],2)
  my_comparisons <- list()
  for (nn in 1:ncol(comparisons)) {
    my_comparisons[[nn]] <- as.vector(c(comparisons[1, nn], comparisons[2, nn]))
  }
  df.test <- df %>%
    dplyr::filter( Compare_Group1 == TRUE & 
                     Pre_Post %in% c("Pre", "Notreatment") )  %>%
    dplyr::filter(BA_Group %in% compare.groups[[ll]]) %>%
    mutate(BA_Group = factor(BA_Group, levels = c(compare.groups[[ll]][1],  compare.groups[[ll]][2],  compare.groups[[ll]][3],  compare.groups[[ll]][4]) )) 
  
  for (vv in numeric.variables) {
    # vv = numeric.variables[5]
    df.test2 <- df.test[!is.na(df.test[[vv]]), ] %>% group_by(BA_Group, Gender) %>% dplyr::filter(n()>1) %>%
      ungroup() %>%
      dplyr::filter(n_distinct(BA_Group) >1) %>%
      ungroup() %>%
      as.data.frame()  
    
    if(nrow(df.test2) > 2){
      # summary
      summary.df.Gender <-df.test2 %>% 
        group_by(Gender, BA_Group) %>%
        get_summary_stats(as.name(vv), type = "common")
      result.summary.Gender <- rbind(result.summary.Gender, summary.df.Gender)
      
      # anova
      # stat.anova.Gender <- compare_means(
      #   formula(paste(vv, "~", "BA_Group")) , data = df.test,
      #   method = "anova", group.by = "Gender", p.adjust.method = "bonferroni"
      # )
      stat.anova.Gender <- df.test2 %>%
        group_by(Gender) %>%
        anova_test(formula(paste(vv, "~", "BA_Group"))) %>%
        adjust_pvalue(method = "bonferroni") %>%
        add_significance("p.adj")%>%
        add_significance("p")%>%
        as.data.frame() 
      stat.anova.Gender <- stat.anova.Gender %>%
        mutate(y.position = (max(df.test2[[vv]], na.rm = T))* 1.09) 
      result.anova.Gender <- rbind(result.anova.Gender, stat.anova.Gender)
      
      # t.test
      # stat.test.Gender <- compare_means(
      #   formula(paste(vv, "~", "BA_Group")) , data = df.test %>% group_by(BA_Group, Gender) %>% dplyr::filter(n()>1),
      #   method = "t.test", p.adjust.method = "bonferroni", group.by = "Gender"
      # )%>%
      #   add_significance("p.adj") 
      stat.test.Gender <- df.test2 %>%
        group_by(Gender) %>%
        t_test(formula(paste(vv, "~", "BA_Group"))) %>%
        adjust_pvalue(method = "bonferroni") %>%
        add_significance("p.adj")%>%
        add_significance("p")
      result.t.test.Gender <- rbind(result.t.test.Gender, stat.test.Gender)
      
      if (TRUE %in% c(stat.anova.Gender$p <= 0.06)) {
        p1 <- ggbarplot(df.test2, x = "BA_Group", y = vv,
                        add = c("mean_se", "jitter"),
                        position = position_dodge(),
                        color = "BA_Group", 
                        # id = "patient_id", line.color = "gray", line.size = 0.4,
                        palette = "npg",
                        facet.by = c( "Gender"), short.panel.labs = T,
                        ncol = length(unique(df.test[["Gender"]])),
                        xlab = "", ylab = "Values",
                        ylim = c(0, (max(df.test[[vv]], na.rm = T)*1.1))
        )+ labs(color=NULL, title = vv)+
          #stat_compare_means(method = "anova", label.y = 90 ) +
          # stat_pvalue_manual(stat.anova.Gender, label = "p.adj")+ 
          stat_compare_means(aes_string(group= "BA_Group"),
                             method = "anova",
                             # paired = T,
                             label = "p.format", size = 4,
                             #hide.ns = TRUE, show.legend = F,
                             label.y = (max(df.test[[vv]], na.rm = T)*1.09), )+
          theme(
            legend.position="right",
            axis.text.x=element_blank(),
            #axis.text.x = element_text(angle = 0, hjust = 0.5),
            axis.ticks.x=element_blank()
          ) + guides(color=guide_legend(nrow=length(compare.groups[[ll]]),byrow=TRUE))
        print(p1)
        png(filename = paste0(dir1, "STARR", "_Gender_", length(compare.groups[[ll]]), "Groups_", vv, ".barplot.png"), 
            width = (length(unique(df.test[["Gender"]])) * length(compare.groups[[ll]])/2) + (max(nchar(compare.groups[[ll]]))/8), 
            height = 4, res = 300, units = "in")
        print(p1)
        dev.off()
      }
    }
    
    # no gender
    df.test2 <- df.test[!is.na(df.test[[vv]]), ] %>% group_by(BA_Group) %>% dplyr::filter(n()>1) %>%
      ungroup() %>%
      as.data.frame()  %>% 
      dplyr::filter(n_distinct(BA_Group) >1) 
    if(nrow(df.test2) > 2){
      # summary
      summary.df <-df.test2 %>% 
        group_by( BA_Group) %>%
        get_summary_stats(as.name(vv), type = "common")
      result.summary <- rbind(result.summary, summary.df)
      
      # t.test
      # stat.test <- compare_means(
      #   formula(paste(vv, "~", "BA_Group")) , data = df.test,
      #   method = "t.test", p.adjust.method = "bonferroni"
      # )
      library(rstatix)
      stat.test <- df.test2 %>% 
        t_test(formula(paste(vv, "~", "BA_Group"))) %>%
        adjust_pvalue(method = "bonferroni") %>%
        add_significance("p.adj")%>%
        add_significance("p")%>%
        as.data.frame() 
      #         add_xy_position(x = "BA_Group") 
      stat.test <- stat.test %>%
        mutate(y.position = (max(df.test2[[vv]], na.rm = T))* seq(1,1.5, 0.5/(length(my_comparisons)-1)) ) 
      result.t.test <- rbind(result.t.test, stat.test)
      
      # anova
      # stat.anova <- compare_means(
      #   formula(paste(vv, "~", "BA_Group")) , data = df.test,
      #   method = "anova"
      # )
      stat.anova <- df.test2 %>%
        anova_test(formula(paste(vv, "~", "BA_Group"))) %>%
        adjust_pvalue(method = "bonferroni") %>%
        add_significance("p.adj")%>%
        add_significance("p")%>%
        as_tibble() 
      stat.anova <- stat.anova %>%
        mutate(y.position = (max(df.test2[[vv]], na.rm = T))* 1.09)
      result.anova <- rbind(result.anova, stat.anova)
      
      if (TRUE %in% c(stat.test$p <= 0.06)) {
        p1 <- ggbarplot(df.test2, x = "BA_Group", y = vv,
                        add = c("mean_se", "jitter"),
                        position = position_dodge(),
                        color = "BA_Group", 
                        # id = "patient_id", line.color = "gray", line.size = 0.4,
                        palette = "npg",
                        # facet.by = c( "Gender"), short.panel.labs = T,
                        ncol = length(unique(df.test[["Gender"]])),
                        xlab = "", ylab = "Values",
                        ylim = c(0, (max(df.test[[vv]], na.rm = T)*1.5))
        )+labs(color=NULL, title = vv)+
          # stat_compare_means(method = "anova", label.y = 90 ) +
          # stat_pvalue_manual(stat.anova, label = "p.adj" , size = 3)+ 
          stat_compare_means(aes_string(group= "BA_Group"),
                             method = "anova",
                             # paired = T,
                             # label = "p.format", size = 4,
                             hide.ns = TRUE,  show.legend = F,
                             label.y = (max(df.test[[vv]], na.rm = T)*1.5))+
          # stat_compare_means(comparisons = my_comparisons, 
          #                    mapping=aes(label = format.pval(..p.adj.., digits = 1)),
          #                    method = "t.test", 
          #                    # paired = T,
          #                    # label = "p.signif", size = 4, 
          #                    label = "p.format", size = 2, 
          #                    hide.ns = TRUE,  show.legend = F, 
          #                    step.increase = ((max(df.test[[vv]], na.rm = T))*0.5)/500,
          #                    # label.y = (max(df.test[[vv]], na.rm = T)*1) 
          #                    ) +
          stat_pvalue_manual(stat.test #%>% dplyr::filter(p.adj <= 0.05)
                             # , label = "p = {p.adj}"
                             , label = "{p.adj}{p.signif}"
                             # , label = "p.adj"
                             , hide.ns = TRUE
                             , size = 3
          )+ 
          theme(
            legend.position="right", # legend.justification="middle",
            axis.text.x=element_blank(),
            #axis.text.x = element_text(angle = 0, hjust = 0.5),
            axis.ticks.x=element_blank()
          ) + guides(color=guide_legend(nrow=length(compare.groups[[ll]]),byrow=TRUE))
        print(p1)
        png(filename = paste0(dir1, "STARR", "_", length(compare.groups[[ll]]), "Groups_", vv, ".barplot.png"), 
            width = (length(compare.groups[[ll]])/1.5)+ (max(nchar(compare.groups[[ll]]))/8), 
            height = 5, res = 300, units = "in")
        print(p1)
        dev.off()
      }
    }
    
    
    print(vv)
  }
}
rm(vv, p1, numeric.variables, df.test,comparisons, my_comparisons)


# 3 Groups_Plot Prediabetes: STARR (OLD, non Diabetes) and SGLT2_V5 --------------------------
# df_Pre <- df %>% dplyr::filter( 
#   Age_Group == "Older" & 
#     Diabetes == "No"&
#     Pre_Post %in% c("Pre", "Notreatment")  &
#     Prediabetes1 %in% c("PreD", "No")   &
#     # sample_id != c("STARR_081")
#     BMI >= 18.5
# ) 
# df_Pre <- df %>% dplyr::filter( 
#   Compare_Group1 == TRUE & 
#     Diabetes == "No" &
#     Prediabetes1 %in% c("PreD", "No")  
# ) 
# label.id = unique(df%>%dplyr::filter(Pre_Post %in% c("Post"))%>% pull(patient_id))
# 
## point plot for A1c and BMI ------------
# p <- ggplot(df_Pre, aes(x = A1c_Value, y = BMI)) + 
#   geom_point(aes(color=Gender, shape = Gender), size = 2) +
#   geom_point(data = df_Pre %>%  dplyr::filter(patient_id %in% label.id),
#              pch=21, 
#              size=5,
#              colour="purple") +
#   labs(title="Patients distribution", y = "Body mass index", x = "Hemoglobin A1C (%)") +  
#   ylim(min(df_Pre$BMI), max(df_Pre$BMI[which(df_Pre$BMI < Inf)])*1.1) + xlim(min(df_Pre$A1c_Value), 6.4) + 
#   theme_classic() +
#   theme(panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
#         axis.title.x = element_text(color = "black", size = 13, face = "plain"),
#         axis.title.y = element_text(color = "black", size = 13, face = "plain"),
#         axis.text.x = element_text(color = "black", size = 11.5, face = "plain", angle = 0, hjust = 1),
#         plot.title = element_text(color = "black", size = 14, face = "bold", hjust = 0.5),
#         legend.text = element_text(color = "black", size = 10.5, face = "plain"),
#         #legend.position = c(0.9, 0.8),
#         legend.position = c(.14, .88),
#         #legend.position = "none",
#         #legend.justification = c("left", "bottom"),
#         #legend.box.just = "left",
#         legend.spacing.y = unit(0.15, 'cm') ,
#         legend.key.size = unit(0.45, "cm"),
#         legend.background = element_rect( fill = "grey98", color = "grey98", linewidth = 1)
#   ) + 
#   geom_vline(xintercept = 5.69, color = "black", linetype="dashed", size=0.2) +
#   geom_vline(xintercept = 6.4, color = "black", linetype="dashed", size=0.2) +
#   geom_hline(yintercept = 29.99, color = "black", linetype="dashed", size=0.2) +
#   geom_hline(yintercept = 24.99, color = "black", linetype="dashed", size=0.2) 
# p
# png(file=paste0(dir1, "3Groups_PatientInfor_A1c_BMI_distribution_Circle.png"), 
#     width=4, height=4, units = "in", res = 300)
# p
# dev.off()
# 
# p <- ggplot(df_Pre, aes(x = A1c_Value, y = BMI)) + 
#   geom_point(aes(color=Gender, shape = Gender), size = 2) +
#   labs(title="Patients distribution", y = "Body mass index", x = "Hemoglobin A1C (%)") +  
#   ylim(min(df_Pre$BMI), max(df_Pre$BMI[which(df_Pre$BMI < Inf)])*1.1) + 
#   xlim(min(df_Pre$A1c_Value), 6.4) + 
#   theme_classic() +
#   theme(panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
#         axis.title.x = element_text(color = "black", size = 13, face = "plain"),
#         axis.title.y = element_text(color = "black", size = 13, face = "plain"),
#         axis.text.x = element_text(color = "black", size = 11.5, face = "plain", angle = 0, hjust = 1),
#         plot.title = element_text(color = "black", size = 14, face = "bold", hjust = 0.5),
#         legend.text = element_text(color = "black", size = 10.5, face = "plain"),
#         #legend.position = c(0.9, 0.8),
#         legend.position = c(.14, .88),
#         #legend.position = "none",
#         #legend.justification = c("left", "bottom"),
#         #legend.box.just = "left",
#         legend.spacing.y = unit(0.15, 'cm') ,
#         legend.key.size = unit(0.45, "cm"),
#         legend.background = element_rect( fill = "grey98", color = "grey98", linewidth = 1)
#   ) + 
#   geom_vline(xintercept = 5.69, color = "black", linetype="dashed", size=0.2) +
#   geom_vline(xintercept = 6.4, color = "black", linetype="dashed", size=0.2) +
#   geom_hline(yintercept = 29.99, color = "black", linetype="dashed", size=0.2) +
#   geom_hline(yintercept = 24.99, color = "black", linetype="dashed", size=0.2) 
# p
# png(file=paste0(dir1, "3Groups_PatientInfor_A1c_BMI_distribution.png"), 
#     width=4, height=4, units = "in", res = 300)
# p
# dev.off()
# 
# p <- ggplot(df_Pre, aes(x = Age, y = BMI)) + 
#   geom_point(aes(color=Gender, shape = Gender), size = 2) +
#   labs(title="Patients distribution", y = "Body mass index", x = "Age (years)") +  
#   ylim(min(df_Pre$BMI), max(df_Pre$BMI[which(df_Pre$BMI < Inf)])*1.1) + 
#   xlim(60, max(df_Pre$Age)) + 
#   theme_classic() +
#   theme(panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
#         axis.title.x = element_text(color = "black", size = 13, face = "plain"),
#         axis.title.y = element_text(color = "black", size = 13, face = "plain"),
#         axis.text.x = element_text(color = "black", size = 11.5, face = "plain", angle = 0, hjust = 1),
#         plot.title = element_text(color = "black", size = 14, face = "bold", hjust = 0.5),
#         legend.text = element_text(color = "black", size = 10.5, face = "plain"),
#         #legend.position = c(0.9, 0.8),
#         legend.position = c(.14, .88),
#         #legend.position = "none",
#         #legend.justification = c("left", "bottom"),
#         #legend.box.just = "left",
#         legend.spacing.y = unit(0.15, 'cm') ,
#         legend.key.size = unit(0.45, "cm"),
#         legend.background = element_rect( fill = "grey98", color = "grey98", linewidth = 1)
#   ) + 
#   geom_vline(xintercept = 60, color = "black", linetype="dashed", size=0.2) +
#   # geom_vline(xintercept = 45, color = "black", linetype="dashed", size=0.2) +
#   geom_hline(yintercept = 29.99, color = "black", linetype="dashed", size=0.2) +
#   geom_hline(yintercept = 24.99, color = "black", linetype="dashed", size=0.2) 
# p
# png(file=paste0(dir1, "3Groups_PatientInfor_Age_BMI_distribution.png"), 
#     width=4, height=4, units = "in", res = 300)
# p
# dev.off()  


## plot tables in two group ------------------ 
# compare.groups <- list(
#   compare1 = c("Older_Lean",  "Older_Overweight" ),
#   compare2 = c("Older_Overweight", "Older_PreD_Obese" ),
#   compare3 = c("Older_PreD_Obese", "Diabetes")
# )
# for (cc in names(cols.list)){
#   for (ll in names(compare.groups)) {
#     rm(demtable)
#     demtable <- df %>%
#       dplyr::filter( Compare_Group1 == TRUE & 
#                        Pre_Post %in% c("Pre", "Notreatment") )  %>%
#       dplyr::filter(BA_Group %in% compare.groups[[ll]]) %>%
#       mutate(BA_Group = factor(BA_Group, levels = c(compare.groups[[ll]][1],  compare.groups[[ll]][2]))) %>%
#       dplyr::select(one_of(c("sample_id", "BA_Group",cols.list[[cc]]))) %>%
#       dplyr::filter(complete.cases(.)) %>% # delete missing values
#       tbl_summary(by = "BA_Group", include = -c("sample_id"), 
#                   type = everything() ~ "continuous",
#                   statistic = list(all_continuous() ~ "{mean} ({sd})",
#                                    all_categorical() ~ "{n} / {N} ({p}%)"),
#                   digits = all_continuous() ~ 2,
#                   # label = grade ~ "Tumor Grade", # change names
#                   missing_text = "(Missing)"
#       ) %>%
#       add_p(all_continuous() ~ "t.test")%>% 
#       add_n() %>%
#       add_overall(last = FALSE) %>%
#       add_significance_stars() %>%
#       modify_header(label ~ "**Variable**") %>%
#       modify_spanning_header(c("stat_1", "stat_2") ~ "**Group**") %>%
#       modify_footnote(
#         all_stat_cols() ~ "Mean (SD)" #" or Frequency (%)"
#       ) %>%
#       # modify_caption("**Drug - dapagliflozin**") %>%
#       bold_labels()
#     demtable
#     gt::gtsave(as_gt(demtable),file=paste0(dir1, "STARR_", cc, "_", compare.groups[[ll]][1], ".",  compare.groups[[ll]][2], "_Group_Event2.png"), vwidth = 1500, vheight = 1500)
#   }
# }
# 
# compare.groups <- list(
#   compare1 = c("Older_Lean",  "Older_Overweight" ),
#   compare2 = c("Older_Overweight", "Older_PreD_Obese" )
# )
# for (cc in names(cols.list)){
#   for (ll in names(compare.groups)) {
#     rm(demtable)
#     demtable <- df %>%
#       dplyr::filter( Compare_Group1 == TRUE & 
#                        Pre_Post %in% c("Pre", "Notreatment") )  %>%
#       dplyr::filter(BA_Group %in% compare.groups[[ll]]) %>%
#       mutate(BA_Group = factor(BA_Group, levels = c(compare.groups[[ll]][1],  compare.groups[[ll]][2]))) %>%
#       dplyr::select(one_of(c("sample_id", "BA_Group","Gender", cols.list[[cc]]))) %>%
#       dplyr::filter(complete.cases(.)) %>% # delete missing values
#       tbl_strata(
#         strata = Gender,
#         .tbl_fun =
#           ~ .x %>%
#           tbl_summary(by = BA_Group, include = -c(sample_id), 
#                       type = everything() ~ "continuous",
#                       statistic = list(all_continuous() ~ "{mean} ({sd})",
#                                        all_categorical() ~ "{n} / {N} ({p}%)"),
#                       digits = all_continuous() ~ 2,
#                       # label = grade ~ "Tumor Grade", # change names
#                       missing_text = "(Missing)") %>%
#           add_p(all_continuous() ~ "t.test"
#           ) %>% 
#           # add_n() %>%
#           add_overall(last = FALSE) %>%
#           add_significance_stars() %>%
#           modify_header(label ~ "**Variable**") %>%
#           modify_spanning_header(c("stat_1", "stat_2") ~ "**Timeline**") %>%
#           modify_footnote(
#             all_stat_cols() ~ "Mean (SD)" #" or Frequency (%)"
#           ) %>%
#           modify_caption("**Table 2. Type of Assessments**") %>%
#           bold_labels(),
#         .header = "**{strata}**, N = {n/2}"
#       ) 
#     demtable
#     gt::gtsave(as_gt(demtable),file=paste0(dir1, "STARR_Gender_", cc, "_", compare.groups[[ll]][1], ".",  compare.groups[[ll]][2], "_Group_Event2.png"), vwidth = 2000, vheight = 1500)
#     
#   }
# }

## plot tables in three group ------------------ 
# compare.groups <- list(
#   compare1 = c("Older_Lean",  "Older_Overweight", "Older_PreD_Obese" )
# )
# for (cc in names(cols.list)){
#   for (ll in names(compare.groups)) {
#     rm(demtable)
#     demtable <- df %>%
#       dplyr::filter( Compare_Group1 == TRUE & 
#                        Pre_Post %in% c("Pre", "Notreatment") )  %>%
#       dplyr::filter(BA_Group %in% compare.groups[[ll]]) %>%
#       mutate(BA_Group = factor(BA_Group, levels = c(compare.groups[[ll]][1],  compare.groups[[ll]][2],  compare.groups[[ll]][3]) )) %>%
#       dplyr::select(one_of(c("sample_id", "BA_Group", cols.list[[cc]]))) %>%
#       dplyr::filter(complete.cases(.)) %>% # delete missing values
#       tbl_summary(by = "BA_Group", include = -c("sample_id"), 
#                   type = everything() ~ "continuous",
#                   statistic = list(all_continuous() ~ "{mean} ({sd})",
#                                    all_categorical() ~ "{n} / {N} ({p}%)"),
#                   digits = all_continuous() ~ 2,
#                   # label = grade ~ "Tumor Grade", # change names
#                   missing_text = "(Missing)"
#       ) %>%
#       add_p(all_continuous() ~ "aov")%>% # "oneway.test" # "ancova"
#       add_n() %>%
#       add_overall(last = FALSE) %>%
#       add_significance_stars() %>%
#       modify_header(label ~ "**Variable**") %>%
#       modify_spanning_header(c("stat_1", "stat_2") ~ "**Group**") %>%
#       modify_footnote(
#         all_stat_cols() ~ "Mean (SD)" #" or Frequency (%)"
#       ) %>%
#       modify_caption("**Table 2. Type of Assessments in Three Groups**") %>%
#       bold_labels()
#     demtable
#     gt::gtsave(as_gt(demtable),file=paste0(dir1, "STARR_", cc, "_", "3Groups", "_Group_Event2.png"), vwidth = 3000, vheight = 1500)
#   }
# }

# for (cc in names(cols.list)){
#   for (ll in names(compare.groups)) {
#     rm(demtable)
#     demtable <- df %>%
#       dplyr::filter( Compare_Group1 == TRUE & 
#                        Pre_Post %in% c("Pre", "Notreatment") )  %>%
#       dplyr::filter(BA_Group %in% compare.groups[[ll]]) %>%
#       mutate(BA_Group = factor(BA_Group, levels = c(compare.groups[[ll]][1],  compare.groups[[ll]][2],  compare.groups[[ll]][3]))) %>%
#       dplyr::select(one_of(c("sample_id", "BA_Group","Gender", cols.list[[cc]]))) %>%
#       dplyr::filter(complete.cases(.)) %>% # delete missing values
#       tbl_strata(
#         strata = Gender,
#         .tbl_fun =
#           ~ .x %>%
#           tbl_summary(by = BA_Group, include = -c(sample_id), 
#                       type = everything() ~ "continuous",
#                       statistic = list(all_continuous() ~ "{mean} ({sd})",
#                                        all_categorical() ~ "{n} / {N} ({p}%)"),
#                       digits = all_continuous() ~ 2,
#                       # label = grade ~ "Tumor Grade", # change names
#                       missing_text = "(Missing)") %>%
#           add_p(all_continuous() ~ "aov")%>% # "oneway.test" # "ancova"
#           # add_n() %>%
#           add_overall(last = FALSE) %>%
#           add_significance_stars() %>%
#           modify_header(label ~ "**Variable**") %>%
#           modify_spanning_header(c("stat_1", "stat_2") ~ "**Timeline**") %>%
#           modify_footnote(
#             all_stat_cols() ~ "Mean (SD)" #" or Frequency (%)"
#           ) %>%
#           modify_caption("**Table 2. Type of Assessments in Three Groups by Gender**") %>%
#           bold_labels(),
#         .header = "**{strata}**, N = {n/2}"
#       ) 
#     demtable
#     gt::gtsave(as_gt(demtable),file=paste0(dir1, "STARR_Gender_", cc, "_", "3Groups", "_Group_Event2.png"), vwidth = 3500, vheight = 1500)
#     
#   }
# }


# Plot sglt2 only ----------------------
library(tidyverse)
df_sg <- df %>% dplyr::filter(
  Pre_Post %in% c("Pre", "Post") 
) %>% 
  group_by(patient_id) %>%
  mutate(n = n()) %>%
  dplyr::filter(n >= 2) %>%
  ungroup() %>%
  dplyr::select( -n)%>% 
  distinct(sample_id, .keep_all = TRUE)

## point plot for A1c and BMI ------------
p <- ggplot(df_sg, aes(x = A1c_Value, y = BMI)) + 
  geom_point(aes(color=Group, shape = Group) , size = 2)  +# shape = Gender, 
  # scale_colour_discrete(guide = "none")+ 
  # guides(color = "none") +
  geom_segment(data = reshape(df_sg %>% dplyr::select(one_of(c("A1c_Value", "BMI", "Pre_Post","patient_id", "Group" )))%>% as.data.frame(), 
                              idvar = c("patient_id", "Group"), timevar = "Pre_Post", direction = "wide"),
               aes(x=A1c_Value.Pre, xend=A1c_Value.Post, y=BMI.Pre, yend=BMI.Post, color=Group), size = 0.25,
               arrow = arrow(length = unit(0.20, "cm")), show.legend=FALSE) +
  labs(title="Patients distribution", y = "Body mass index", x = "Hemoglobin A1C (%)") +  
  ylim(min(df_sg$BMI), max(df_sg$BMI[which(df_sg$BMI < Inf)])*1) + xlim(min(df_sg$A1c_Value), 6.2) + 
  theme_classic() +
  theme(panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
        axis.title.x = element_text(color = "black", size = 13, face = "plain"),
        axis.title.y = element_text(color = "black", size = 13, face = "plain"),
        axis.text.x = element_text(color = "black", size = 11.5, face = "plain", angle = 0, hjust = 1),
        plot.title = element_text(color = "black", size = 14, face = "bold", hjust = 0.5),
        legend.text = element_text(color = "black", size = 10.5, face = "plain"),
        #legend.position = c(0.9, 0.8),
        legend.position = c(.14, .88),
        #legend.position = "none",
        #legend.justification = c("left", "bottom"),
        #legend.box.just = "left",
        legend.spacing.y = unit(0.15, 'cm') ,
        legend.key.size = unit(0.45, "cm"),
        legend.background = element_rect( fill = "grey98", color = "grey98", linewidth = 1)
  ) 
p
png(file=paste0(dir1, "SGLT2i_PatientInfor_A1c_BMI_distribution.png"), 
    width=4, height=4, units = "in", res = 300)
p
dev.off()


## plot numeric.variables by ggpaired and point plot-----------------------------------------------------------------------------------------
cat(paste(shQuote(colnames(df)[sapply(df, is.numeric)], type = "cmd"), collapse = ", "))
numeric.variables = c("White.Blood.Cell...WBC", "Red.Blood.Cell...RBC", "Hemoglobin", "Hematocrit", "MCV", "MCH", 
                      "Platelets", "Glucose", "Blood.Urea.Nitrogen..BUN.", "Creatinine", "Sodium", "Potassium", 
                      "Calcium", "Protein..Total", "Albumin", "Alkaline.Phosphatase", "Bilirubin..Total", "AST", "ALT", 
                      "Cholesterol..Total", "Triglycerides", "HDL", "LDL", "VLDL", "A1c_Value", 
                      "Weight_pounds", "Height_inches", "BMI", "Blood.Pressure.Systolic", "Blood.Pressure.Diastolic", 
                      "Insulin", "FFA", "TNF.alpha", "IL.6", "IL.8.CXCL8", "VEGF", "IL.7", "CCL3.MIP.1alpha", 
                      "CCL4.MIP.1", "RAGE.AGER", "SOST.Sclerostin", "ADAMTS13", "Osteropontin.OPN", "ICAM.1.CD54", 
                      "IL.15", "GDF.15", "TNF.RI.TNFRSF1A", "Fas.TNFRSF6.CD95", "GDNF", "ActivinA", 
                      "Heart.Rate..bpm.", "Respirations.per.minute", "Temperature..F.", 
                      "Monocyte",  "Neutrophils", "Eosinophils", "Basophils", 
                      "MCHC", "RDW", "Chloride", "CO2", "Manual.eGFR.calculation",
                      "MoCA.Score",  
                      "Tissue...Fat.", "T.Mass..kg.", "Region...Fat.", "Tissue..g.", "Fat..g.", 
                      "Lean..g.", "BMC..g.", "Fat.Free..g.", "HOMA_IR", "Weight_kg")# "Granulocytes","TSH", "Urine.pH", "PTT",  "PT", "INR", "Z.Score", 
# numeric.variables = c("Weight_kg", "BMI", "Glucose", "A1c_Value")

for (vv in numeric.variables) {
  # vv = numeric.variables[5]
  
  df.test2 <- df_sg[!is.na(df_sg[[vv]]), ] %>% group_by(Group, patient_id) %>% dplyr::filter(n()>1) %>%
    ungroup() %>%
    as.data.frame()  %>% 
    dplyr::filter(n_distinct(Group) >1) %>% 
    arrange(sample_id, patient_id)
  
  if(nrow(df.test2) > 2){
    # summary
    summary.df <-df.test2 %>% 
      group_by( Group,  (!!sym("Pre_Post")) ) %>%
      get_summary_stats(as.name(vv), type = "common")
    result.summary.sg <- rbind(result.summary.sg, summary.df)
    
    # t.test
    library(rstatix)
    stat.test.sg <- df.test2 %>% 
      group_by( Group) %>%
      t_test(formula(paste(vv, "~", "Pre_Post")), paired = T) %>%
      adjust_pvalue(method = "bonferroni") %>%
      add_significance("p.adj")%>%
      add_significance("p")%>%
      as.data.frame() 
    stat.test.sg <- stat.test.sg %>%
      mutate(y.position = (max(df.test2[[vv]], na.rm = T))* 1.05 ) 
    result.t.test.sg <- rbind(result.t.test.sg, stat.test.sg)
    
    # anova
    stat.anova.sg <- df.test2 %>%
      anova_test(formula(paste(vv, "~", "Pre_Post", "*", "Group"))) %>%
      adjust_pvalue(method = "bonferroni") %>%
      add_significance("p.adj")%>%
      add_significance("p")%>%
      as_tibble() 
    stat.anova.sg <- stat.anova.sg %>%
      mutate(y.position = (max(df.test2[[vv]], na.rm = T))* 1.0,
             variable = vv) 
    result.anova.sg <- rbind(result.anova.sg, stat.anova.sg)
    
    
    if (TRUE %in% c(stat.anova.sg$p <= 0.06)) {
      p1 <- ggpaired(df.test2 ,
                     x = "Pre_Post", y = vv,
                     color = "Pre_Post", palette = "npg",
                     id = "patient_id", line.color = "gray", line.size = 0.4,
                     facet.by = c("Group"), short.panel.labs = T,  scales="fixed",
                     xlab = "", ylab = vv,
                     ylim = c(range(df.test2[[vv]])[1]*0.9, range(df.test2[[vv]])[2]*1.1)
      )+
        # stat_compare_means(
        #   aes_string(group="Pre_Post"),
        #   method = "t.test", paired = T,
        #   label = "p.format",  hide.ns = TRUE, show.legend = F, size = 5,
        #   label.y = range(df_sg[[vv]])[1]*0.91 )+
        stat_pvalue_manual(stat.test.sg #%>% dplyr::filter(p.adj <= 0.05)
                           # , label = "p = {p.adj}"
                           , label = "{p.adj}{p.signif}"
                           # , label = "p.adj"
                           , hide.ns = TRUE
                           , size = 3
        )+ 
        theme(axis.text.x=element_blank(),
              axis.ticks.x=element_blank())
      print(p1)
      png(filename = paste0(dir1, "SGLT2i_",vv, ".png"),
          width = 4, height = 5, res = 300, units = "in")
      print(p1)
      dev.off()
    }
  }
}
# for (vv in numeric.variables) {
#   # vv = numeric.variables[5]  
#   p <- ggplot(data = reshape(df_sg %>% dplyr::select(one_of(c(vv, "Pre_Post","patient_id", "Group" )))%>% as.data.frame(), 
#                              idvar = c("patient_id", "Group"), timevar = "Pre_Post", direction = "wide"), 
#               aes(x = .data[[paste0(vv, ".Pre")]], y = .data[[paste0(vv, ".Post") ]])) + 
#     geom_point(aes(color=Group, shape = Group) , size = 2)  +# shape = Gender, 
#     # scale_colour_discrete(guide = "none")+ 
#     # guides(color = "none") +
#     labs(title= vv, y = "Post Treatment", x = "Pre Treatment") +  
#     # ylim(min(Composition_3$BMI), max(Composition_3$BMI[which(Composition_3$BMI < Inf)])*1) + 
#     ylim(min(df_sg[[vv]]), max(df_sg[[vv]])) + 
#     xlim(min(df_sg[[vv]]), max(df_sg[[vv]])) + 
#     theme_classic() +
#     theme(panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
#           axis.title.x = element_text(color = "black", size = 13, face = "plain"),
#           axis.title.y = element_text(color = "black", size = 13, face = "plain"),
#           axis.text.x = element_text(color = "black", size = 11.5, face = "plain", angle = 0, hjust = 1),
#           plot.title = element_text(color = "black", size = 14, face = "bold", hjust = 0.5),
#           legend.text = element_text(color = "black", size = 10.5, face = "plain"),
#           #legend.position = c(0.9, 0.8),
#           legend.position = c(.14, .88),
#           #legend.position = "none",
#           #legend.justification = c("left", "bottom"),
#           #legend.box.just = "left",
#           legend.spacing.y = unit(0.15, 'cm') ,
#           legend.key.size = unit(0.45, "cm"),
#           legend.background = element_rect( fill = "grey98", color = "grey98", linewidth = 1)
#     ) +
#     geom_abline(intercept = 0, slope = 1, colour = "black", linetype = "dashed")
#   print(p)
#   png(file=paste0(dir1, "SGLT2i_",vv, "_Pre_Post.png"), 
#       width=4, height=4, units = "in", res = 300)
#   print(p)
#   dev.off()
# }
rm(vv, p1, numeric.variables)



## plot table -------------------
cols.list <- list(
  Others = c("White.Blood.Cell...WBC", "Red.Blood.Cell...RBC", "Hemoglobin", "Hematocrit", "MCV", "MCH", 
             "Platelets",  "Blood.Urea.Nitrogen..BUN.", "Creatinine", "Sodium", "Potassium", 
             "Calcium", "Protein..Total", "Albumin", "Alkaline.Phosphatase", "AST", "ALT", 
             "Temperature..F.", 
             "Monocyte", "Neutrophils", "Eosinophils", 
             "MCHC", "RDW", "Chloride", "CO2"
  ), #  
  miss = c("Triglycerides","Z.Score", "Granulocytes", "Manual.eGFR.calculation","PT", "INR", "PTT" ,"TSH", "Urine.pH",  "MoCA.Score" , "Age", 
           "Basophils", "Bilirubin..Total"),
  Basal = c("Cholesterol..Total",  "HDL", "LDL", "VLDL", 
            "A1c_Value", "Glucose", 
            "Tissue...Fat.", "Region...Fat.", "Fat..g.","Lean..g.", "Tissue..g.",  "BMC..g.", "Fat.Free..g.", "T.Mass..kg.", "Weight_kg",
            "Weight_pounds", "Height_inches", "BMI", "Blood.Pressure.Systolic", "Blood.Pressure.Diastolic",
            "Heart.Rate..bpm.", "Respirations.per.minute"),
  Plasma = c("HOMA_IR",  "Insulin", "FFA", "TNF.alpha", "IL.6", "IL.8.CXCL8", "VEGF", "IL.7", "CCL3.MIP.1alpha", 
             "CCL4.MIP.1", "RAGE.AGER", "SOST.Sclerostin", "ADAMTS13", "Osteropontin.OPN", "ICAM.1.CD54", 
             "IL.15", "GDF.15", "TNF.RI.TNFRSF1A", "Fas.TNFRSF6.CD95", "GDNF", "ActivinA" ))

### table for all -------------------
for (cc in names(cols.list)){
  demtable <- df_sg %>%
    dplyr::filter(Pre_Post %in% c("Pre", "Post")) %>%
    mutate(Pre_Post = factor(Pre_Post, levels = c("Pre", "Post"))) %>%
    # dplyr::filter(Group == "SGLT2i") %>%
    dplyr::select(one_of(c("patient_id", "Pre_Post", "Group", cols.list[[cc]]))) %>%
    # dplyr::filter(complete.cases(.)) %>% # delete missing values
    group_by(patient_id) %>% # keep IDs with both measurements
    dplyr::filter(n() == 2) %>%
    ungroup() %>% as.data.frame() %>%
    tbl_strata(
      strata = Group,
      .tbl_fun =
        ~ .x %>%
        tbl_summary(by = Pre_Post, include = -c(patient_id), 
                    type = everything() ~ "continuous",
                    statistic = list(all_continuous() ~ "{mean} ({sd})",
                                     all_categorical() ~ "{n} / {N} ({p}%)"),
                    digits = all_continuous() ~ 2,
                    # label = grade ~ "Tumor Grade", # change names
                    missing_text = "(Missing)") %>%
        add_p(all_continuous() ~ "paired.t.test", 
              group = "patient_id"
        ) %>% 
        # add_n() %>%
        # add_overall(last = FALSE) %>%
        add_significance_stars() %>%
        modify_header(label ~ "**Variable**") %>%
        modify_spanning_header(c("stat_1", "stat_2") ~ "**Timeline**") %>%
        modify_footnote(
          all_stat_cols() ~ "Mean (SD)" #" or Frequency (%)"
        ) %>%
        modify_caption("**Table 2. Type of Assessments**") %>%
        bold_labels(),
      .header = "**{strata}**, N = {n/2}"
    ) 
  demtable
  gt::gtsave(as_gt(demtable),file=paste0(dir1, "SGLT2i_",cc, "_Group_Event2.png"), vwidth = 1500, vheight = 1500)
  
}

### table by Gender ------------------
# for (cc in names(cols.list)){
#   demtable <- df_sg %>%
#     dplyr::filter(Pre_Post %in% c("Pre", "Post")) %>%
#     mutate(Pre_Post = factor(Pre_Post, levels = c("Pre", "Post"))) %>%
#     # dplyr::filter(Group == "SGLT2i") %>%
#     dplyr::select(one_of(c("patient_id", "Pre_Post", "Group","Gender", cols.list[[cc]]))) %>%
#     # dplyr::filter(complete.cases(.)) %>% # delete missing values
#     group_by(patient_id) %>% # keep IDs with both measurements
#     dplyr::filter(n() == 2) %>%
#     ungroup() %>% as.data.frame() %>%
#     tbl_strata(
#       strata = c( Gender, Group),
#       .tbl_fun =  ~ .x %>%
#         tbl_summary(by = Pre_Post, include = -c(patient_id), 
#                     type = everything() ~ "continuous",
#                     statistic = list(all_continuous() ~ "{mean} ({sd})",
#                                      all_categorical() ~ "{n} / {N} ({p}%)"),
#                     digits = all_continuous() ~ 2,
#                     # label = grade ~ "Tumor Grade", # change names
#                     missing_text = "(Missing)") %>%
#         # add_p(all_continuous() ~ "paired.t.test", 
#         #       group = "patient_id") %>% 
#         # add_n() %>%
#         add_overall(last = FALSE) %>%
#         add_significance_stars() %>%
#         modify_header(label ~ "**Variable**") %>%
#         modify_spanning_header(c("stat_1", "stat_2") ~ "**Timeline**") %>%
#         modify_footnote(
#           all_stat_cols() ~ "Mean (SD)" #" or Frequency (%)"
#           ) %>%
#         modify_caption("**Table 2. Type of Assessments**") %>%
#         bold_labels(),
#       .header = "**{strata}**, N = {n/2}"
#     ) 
#   demtable
#   gt::gtsave(as_gt(demtable),file=paste0(dir1, "SGLT2i_Gender_",cc, "_Group_Event2.png"), vwidth = 3000, vheight = 1500)
# }

for (cc in names(cols.list)){
  demtable <- df_sg %>%
    dplyr::filter(Pre_Post %in% c("Pre", "Post")) %>%
    dplyr::filter(Gender %in% c("Male")) %>%
    mutate(Pre_Post = factor(Pre_Post, levels = c("Pre", "Post"))) %>%
    # dplyr::filter(Group == "SGLT2i") %>%
    dplyr::select(one_of(c("patient_id", "Pre_Post", "Group", cols.list[[cc]]))) %>%
    # dplyr::filter(complete.cases(.)) %>% # delete missing values
    group_by(patient_id) %>% # keep IDs with both measurements
    dplyr::filter(n() == 2) %>%
    ungroup() %>% as.data.frame() %>%
    tbl_strata(
      strata = c( Group),
      .tbl_fun =
        ~ .x %>%
        tbl_summary(by = Pre_Post, include = -c(patient_id), 
                    type = everything() ~ "continuous",
                    statistic = list(all_continuous() ~ "{mean} ({sd})",
                                     all_categorical() ~ "{n} / {N} ({p}%)"),
                    digits = all_continuous() ~ 2,
                    # label = grade ~ "Tumor Grade", # change names
                    missing_text = "(Missing)") %>%
        add_p(all_continuous() ~ "paired.t.test", 
              group = "patient_id"
        ) %>% 
        # add_n() %>%
        # add_overall(last = FALSE) %>%
        add_significance_stars() %>%
        modify_header(label ~ "**Variable**") %>%
        modify_spanning_header(c("stat_1", "stat_2") ~ "**Timeline**") %>%
        modify_footnote(
          all_stat_cols() ~ "Mean (SD)" #" or Frequency (%)"
        ) %>%
        modify_caption("**Table 2. Type of Assessments in Male**") %>%
        bold_labels(),
      .header = "**{strata}**, N = {n/2}"
    ) 
  demtable
  gt::gtsave(as_gt(demtable),file=paste0(dir1, "SGLT2i_Gender_Male_",cc, "_Group_Event2.png"), vwidth = 3000, vheight = 1500)
  
  
  demtable <- df_sg %>%
    dplyr::filter(Pre_Post %in% c("Pre", "Post")) %>%
    dplyr::filter(Gender %in% c("Female")) %>%
    mutate(Pre_Post = factor(Pre_Post, levels = c("Pre", "Post"))) %>%
    # dplyr::filter(Group == "SGLT2i") %>%
    dplyr::select(one_of(c("patient_id", "Pre_Post", "Group", cols.list[[cc]]))) %>%
    # dplyr::filter(complete.cases(.)) %>% # delete missing values
    group_by(patient_id) %>% # keep IDs with both measurements
    dplyr::filter(n() == 2) %>%
    ungroup() %>% as.data.frame() %>% group_by(Group, Pre_Post) %>% dplyr::filter(n()>1) 
  
  if (nrow(demtable) > 1) {
    demtable <- demtable %>%
      tbl_strata(
        strata = c( Group),
        .tbl_fun =
          ~ .x %>%
          tbl_summary(by = Pre_Post, include = -c(patient_id), 
                      type = everything() ~ "continuous",
                      statistic = list(all_continuous() ~ "{mean} ({sd})",
                                       all_categorical() ~ "{n} / {N} ({p}%)"),
                      digits = all_continuous() ~ 2,
                      # label = grade ~ "Tumor Grade", # change names
                      missing_text = "(Missing)") %>%
          add_p(all_continuous() ~ "paired.t.test", 
                group = "patient_id"
          ) %>% 
          # add_n() %>%
          # add_overall(last = FALSE) %>%
          add_significance_stars() %>%
          modify_header(label ~ "**Variable**") %>%
          modify_spanning_header(c("stat_1", "stat_2") ~ "**Timeline**") %>%
          modify_footnote(
            all_stat_cols() ~ "Mean (SD)" #" or Frequency (%)"
          ) %>%
          modify_caption("**Table 2. Type of Assessments in Female**") %>%
          bold_labels(),
        .header = "**{strata}**, N = {n/2}"
      ) 
    demtable
    gt::gtsave(as_gt(demtable),file=paste0(dir1, "SGLT2i_Gender_Female_",cc, "_Group_Event2.png"), vwidth = 3000, vheight = 1500)
  }
  
}

save(result.anova, result.anova.Gender,result.anova.sg,
     result.t.test, result.t.test.Gender,result.t.test.sg,
     result.summary, result.summary.Gender, result.summary.sg, 
     file = paste0(dir1, "All.test", "_", "Result.RData"))

# End -------------------------------
rm(df, df_Pre, df_sg, df1, p ,p2 , cols.list, compare.groups, demtable, cc, dir1, label.id, ll)
rm(result.anova, result.anova.Gender,result.anova.sg,
    result.t.test, result.t.test.Gender,result.t.test.sg,
    result.summary, result.summary.Gender,result.summary.sg, 
   summary.df, summary.df.Gender,
   label.id2, nn, df.test2,
   stat.anova, stat.anova.Gender,stat.anova.sg,
   stat.test, stat.test.Gender, stat.test.sg
   )
