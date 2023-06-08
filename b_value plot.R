
library(readxl)
library(ggplot2)
library(dplyr)
b_gal <- read_excel("C:/Users/nieji/Desktop/2023_SGLT2_X_gal_b_value.xlsx") %>% 
  mutate(b_value2 = 100- b_value)

# No Treatment
df <- metadata.all %>% 
  dplyr::filter(Treatment_Group %in% c("Notreatment", "NUT_Pre", "SGLT2i_Pre"),
                BA_Group %in% c("Older_Lean", "Older_Overweight", "Older_PreD_Obese", "Diabetes")#
  ) %>% 
  mutate(BA_Group = factor(BA_Group, levels = c("Older_Lean", "Older_Overweight", "Older_PreD_Obese", "Diabetes"))) %>% 
  arrange(vars(one_of(c("BA_Group", col.group, "sample_id")))) %>% 
  group_by_at(vars(one_of(c("BA_Group", col.group, "sample_id")))) %>%
  dplyr::filter(get(col.group1) > 0) %>% # , !is.na(get(col.group1))
  summarise(Avg = mean(get(col.group1), na.rm = TRUE)) 

# combine
df.merge <- df %>% 
  group_by_at(vars(one_of(c("BA_Group", "sample_id")))) %>%
  summarise(Avg = mean(get('Avg'), na.rm = TRUE)) %>% 
  inner_join(b_gal, by = "sample_id" ) 


# SGLT2i
df <- metadata.all %>% 
  dplyr::filter(Pre_Post %in% c("Pre", "Post") ) %>% 
  dplyr::filter(!patient_id %in% c("76615")) %>% 
  arrange(vars(one_of(c("Pre_Post", "Group", col.group, "patient_id", "sample_id")))) %>% 
  group_by_at(vars(one_of(c("Pre_Post", "Group", col.group, "patient_id", "sample_id")))) %>%
  dplyr::filter(get(col.group1) > 0) %>% # , !is.na(get(col.group1))
  summarise(Avg = mean(get(col.group1), na.rm = TRUE)) %>% 
  mutate(Group_Event = factor(paste0(Group, "_", Pre_Post), levels = c("NUT_Pre", "NUT_Post",
                                                                       "SGLT2i_Pre", "SGLT2i_Post")))  %>%
  dplyr::filter(
    !sample_id %in% c("76658_V5", "76658_V11", "76656_V11")
    )#, "76664_V5", "76660_V11"

# combine
df.merge <- df %>% 
  group_by_at(vars(one_of(c("sample_id")))) %>%
  summarise(Avg = mean(get('Avg'), na.rm = TRUE)) %>% 
  inner_join(b_gal, by = "sample_id" ) 



# plot compare
p <- ggplot(df.merge, aes( x = Avg, y = b_value2)) + 
  geom_point(aes(color=sample_id)) +
  labs(title="", y = "Beta Gal Staining", x = "Senecense Score") +  
  ylim(min(df.merge$b_value2), max(df.merge$b_value2)*1.1) + 
  xlim(min(df.merge$Avg), max(df.merge$Avg)*1.1) + 
  theme_classic() +
  theme(panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
        axis.title.x = element_text(color = "black", size = 13, face = "plain"),
        axis.title.y = element_text(color = "black", size = 13, face = "plain"),
        axis.text.x = element_text(color = "black", size = 11.5, face = "plain", angle = 0, hjust = 1),
        plot.title = element_text(color = "black", size = 14, face = "bold", hjust = 0.5),
        legend.text = element_text(color = "black", size = 10.5, face = "plain"),
        #legend.position = c(0.9, 0.8),
        #legend.position = c(.11, .81),
        legend.position = "none",
        #legend.justification = c("left", "bottom"),
        #legend.box.just = "left",
        legend.spacing.y = unit(0.15, 'cm') ,
        legend.key.size = unit(0.45, "cm"),
        legend.background = element_rect( fill = "grey98", color = "grey98", linewidth = 1)
  )  +
  stat_cor(label.x = 1, label.y =90) +
  stat_regline_equation(label.x = 1, label.y = 90-2) +
   geom_smooth(method=lm,  linetype="dashed",se=FALSE,
               color="black", fill="blue") +
  geom_text(
            label= substr(df.merge$sample_id, 4,9),
            nudge_x=0.01, nudge_y=0.01,
            check_overlap=F)
p
png(file=paste0( "Senecense Score.png"), 
    width=4, height=4, units = "in", res = 300)
p
dev.off()




# PLOT whole group  

df <- metadata.all %>% 
  dplyr::filter(Treatment_Group %in% c("Notreatment", "NUT_Pre", "SGLT2i_Pre"),
                BA_Group %in% c("Older_Lean", "Older_Overweight", "Older_PreD_Obese", "Diabetes")#
  ) %>% 
  mutate(BA_Group = factor(BA_Group, levels = c("Older_Lean", "Older_Overweight", "Older_PreD_Obese", "Diabetes"))) %>% 
  arrange(vars(one_of(c("BA_Group", col.group, "sample_id")))) %>% 
  group_by_at(vars(one_of(c("BA_Group", col.group, "sample_id")))) %>%
  dplyr::filter(get(col.group1) > 0) %>% # , !is.na(get(col.group1))
  summarise(Avg = mean(get(col.group1), na.rm = TRUE)) %>% 
  group_by_at(vars(one_of(c("BA_Group", "sample_id")))) %>%
  summarise(Avg = mean(get('Avg'), na.rm = TRUE)) 

library(rstatix)
ano.result <- df  %>%
  as.data.frame() %>%
  rstatix::anova_test(Avg ~ BA_Group, type = 3) %>% 
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj") %>% 
  ungroup() 

test.result <- df%>%
  as.data.frame() %>%
  rstatix::t_test(Avg ~ BA_Group) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj") %>% 
  mutate(method = "T-test")%>% 
  add_xy_position(x = "BA_Group") 

p <-ggline(df,
           x ="BA_Group" , y = "Avg", color = "BA_Group",
           add = c("mean_se", "dotplot"),
           palette = "BA_Group",
           xlab = "", ylab = col.group1
           # ,ylim = c(0, max(df_test1[["Avg"]])*1)
) + theme_bw() +
  theme(
    axis.text.x=element_blank()
    , legend.position="right"
   , legend.title = element_blank()
  ) + guides(color=guide_legend(nrow=4,byrow=TRUE)) +
  stat_pvalue_manual(
    test.result, label = "p.adj.signif", # "p = {p.adj}",
    hide.ns = TRUE,
    tip.length = 0.01
  )+
  labs(
    title = "GG",
    subtitle = get_test_label(ano.result, detailed = TRUE),
    caption = get_pwc_label(test.result)
  )
print(p) 
png(
#  filename = paste0(dir1, file.name, "_",col.group,"_", col.group1, "_",gg,".4Groups.postive.ggline.avg.png"), 
  filename = paste0("4Groups.postive.ggline.avg.png"), 
    width = 6, 
    height = 4.5, res = 300, units = "in")
print(p)
dev.off()

