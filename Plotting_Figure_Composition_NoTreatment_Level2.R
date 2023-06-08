
# rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.
memory.limit(size = 1e+13)

# import parameters and packages --------------------------------------------------------------------------------------------
# source(paste0("/media/jianie/Extreme SSD1/2022-09-01 STARR_SGLT2 Combine/Project Parameters.R"), local = knitr::knit_global())
source(paste0(paste0(Disk, Project.folder, "/", "Project Parameters.R")), local = knitr::knit_global())
library("RColorBrewer")
library("ggpubr")
library(rstatix)
# Import data -----------------------------------------------------------------------------------------------------------------------------------------------------
# export.folder <- "Compositions"
# dir.create(file.path(paste0(Disk, Project.folder, "/", Rds.folder, "/"), export.folder), showWarnings = FALSE)
dir0 <- paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder, "/")
files <- list.files(dir0)

dir.create(file.path(paste0(Disk, Project.folder, "/", figures.folder, "/", export.folder, ".plotting")), showWarnings = FALSE)
dir1 <- paste0(Disk, Project.folder, "/", figures.folder, "/", export.folder, ".plotting", "/")

for (file in files) {
  # file = files[1]
  load(paste0(dir0, file))
  
  # prepare file name
  file.name <- stringr::str_split(file, ".metadata_Compositions.RData")[[1]][1]
  Annotation.file.name <- paste(stringr::str_split(file.name, "_")[[1]][c(1,2,3)], collapse = "_")
  Subtype.file.name <- stringr::str_split(file.name, paste0(Annotation.file.name, "_"))[[1]][2]
  if( length(stringr::str_split(Subtype.file.name, paste0("_"))[[1]]) > 1){
    celltype.name <- stringr::str_split(Subtype.file.name, paste0("_"))[[1]][1]
    #label.name <- stringr::str_split(Subtype.file.name, paste0("_"))[[1]][2]
  } else {
    celltype.name <- Annotation.file.name
    # label.name <- stringr::str_split(Subtype.file.name, paste0("_"))[[1]][1]
  }
  
  # import data
  for (list.name in names(Composition.lists)) {
    # list.name = names(Composition.lists)[1]
    label.name <- stringr::str_split(list.name, paste0("_"))[[1]][length(stringr::str_split(list.name, paste0("_"))[[1]])]
    x.name <- stringr::str_split(list.name, paste0("_"))[[1]][2]
    Composition_3 <- Composition.lists[[list.name]][["Composition_3All"]] %>% 
      dplyr::filter(Treatment_Group %in% c("Notreatment", "NUT_Pre", "SGLT2i_Pre"),
                    BA_Group %in% c("Older_Lean", "Older_Overweight", "Older_PreD_Obese", "Diabetes")#
      ) %>% 
      mutate(BA_Group = factor(BA_Group, levels = c("Older_Lean", "Older_Overweight", "Older_PreD_Obese", "Diabetes")))# 
    Composition <- Composition.lists[[list.name]][["Composition_All"]] 
    
    clusters = colnames(Composition)
    # display.brewer.pal(length(clusters), "Set1")
    # Composition_3$x.label <- factor(Composition_3[[label.name]], levels=clusters)
    
    if (x.name == label.name){
      p1 <- ggplot(Composition_3, aes(x = .data[[x.name]] , y = .data[["composition"]], group =  .data[[x.name]])) + 
        geom_boxplot( size=0.5) +
        geom_point(aes(color= .data[["sample_id"]]), size=2) +
        labs(title="Cell type distribution", y = "Average percentage of all nuclei") +
        ylim(0, max(Composition_3[["composition"]])) +
        theme_classic()+
        theme(panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
              axis.title.x = element_blank(),
              axis.title.y = element_text(color = "black", size = 13, face = "plain"),
              axis.text.x = element_text(color = "black", size = 11.5, face = "plain", angle = 45, hjust = 1),
              plot.title = element_text(color = "black", size = 14, face = "bold", hjust = 0.5),
              #legend.text = element_text(color = "black", size = 10.5, face = "plain"),
              #legend.position = c(0.9, 0.8),
              #legend.position = c(.99, .99),
              legend.position = "none",
              #legend.justification = c("right", "top"),
              #legend.box.just = "right",
              legend.spacing.y = unit(0.15, 'cm') ,
              legend.key.size = unit(0.45, "cm"),
              legend.background = element_rect( fill = "grey98", color = "grey98", linewidth  = 1)
        ) 
      print(p1)
      png(file=paste0(dir1, file.name, "_", x.name, "_", "NoTreatment_4Groups", ".Boxplot.png"), 
          width=3.7, height=4, res = 300, units = "in")
      print(p1)
      dev.off()
      
      # plotting
      p1 <- ggbarplot(Composition_3, x = "BA_Group", y = "composition",
                      add = c("mean_se", "jitter"),
                      position = position_dodge(),
                      color = "BA_Group", 
                      # id = "patient_id", line.color = "gray", line.size = 0.4,
                      palette = "npg",
                      facet.by = c( x.name), short.panel.labs = T,
                      ncol = length(unique(Composition_3[[x.name]])),
                      xlab = "", ylab = "Percentage of Cells")+
        #stat_compare_means(method = "anova", label.y = 90 ) +
        stat_compare_means(aes_string(group= "BA_Group"), 
                           method = "anova", 
                           # paired = T,
                           label = "p.format", size = 4, 
                           #hide.ns = TRUE, show.legend = F, 
                           label.y = (max(Composition_3[["composition"]])*1.1), )+
        theme(
          axis.text.x=element_blank(),
          #axis.text.x = element_text(angle = 0, hjust = 0.5),
          axis.ticks.x=element_blank()
        )
      print(p1)
      png(filename = paste0(dir1, file.name, "_", x.name, "_", "NoTreatment_4Groups", ".barplot.png"), 
          width = length(unique(Composition_3[[x.name]]))*length(unique(Composition_3[["BA_Group"]]))/2, 
          height = 4, res = 300, units = "in")
      print(p1)
      dev.off()
      
      # plot each one
      library(car)
      df1 <- Composition_3 %>% data.frame() %>% 
        group_by_at(vars(one_of(c(x.name, "BA_Group"))))   %>%
        # nest() %>%
        dplyr::filter(n() > 2 ) %>%
        ungroup() %>% as.data.frame() %>% 
        group_by_at(vars(one_of(c(x.name))))   %>%
        dplyr::filter(n_distinct(BA_Group) >1)%>% 
        dplyr::filter(BA_Group %in% c("Older_Lean", "Older_Overweight", "Older_PreD_Obese", "Diabetes")) %>%
        mutate(BA_Group = factor(BA_Group, levels = c("Older_Lean", "Older_Overweight", "Older_PreD_Obese", "Diabetes")))
      
      for (gg in unique(df1[[x.name]]) ) {
        df_test1 <- df1 %>% 
          dplyr::filter(get(x.name) == gg)
        
        if (sum(df_test1$composition) > 0){
          
          
          ano.result <- df_test1  %>%
            as.data.frame() %>%
            rstatix::anova_test(composition ~ BA_Group, type = 3) %>% 
            adjust_pvalue(method = "bonferroni") %>%
            add_significance("p.adj") %>% 
            ungroup() 
          
          if (ano.result$p <= 0.06) {
            test.result <- df_test1%>%
              as.data.frame() %>% # nest() %>% pluck("data", 3) %>% 
              rstatix::t_test(composition ~ BA_Group) %>%
              adjust_pvalue(method = "bonferroni") %>%
              add_significance("p.adj") %>% 
              mutate(method = "T-test")%>% 
              add_xy_position(x = "BA_Group")
            
            # Plot 1 way ANOVA
            p <-ggline(df_test1,
                       x ="BA_Group" , y = "composition", color = "BA_Group",
                       add = c("mean_se", "dotplot"),
                       palette = "BA_Group",
                       xlab = "", ylab = "Percentage of Cells"
                       # ,ylim = c(0, max(df_test1[["Avg"]])*1)
            ) + theme_bw() +
              theme(
                axis.text.x=element_blank()
                , legend.position="right"
                , legend.title = element_blank()
              ) +
              stat_pvalue_manual(
                test.result, label = "p.adj.signif", # "p = {p.adj}",
                hide.ns = TRUE,
                tip.length = 0.01
              )+
              labs(
                title = gg,
                subtitle = get_test_label(ano.result, detailed = TRUE),
                caption = get_pwc_label(test.result)
              )
            print(p)
            png(
              filename = paste0(dir1, file.name, "_", x.name, "_", gg, "_NoTreatment_4Groups", ".ggline.png"), 
              width = 6, 
              height = 4.5, res = 300, units = "in")
            print(p)
            dev.off()
          }}
      }
      
      # df1 <- Composition_3%>% 
      #   dplyr::filter(BA_Group %in% c("Older_Lean", "Older_Overweight"))
      # p1 <- ggbarplot(df1, 
      #                 x = "BA_Group", y = "composition",
      #                 add = c("mean_se", "jitter"),
      #                 position = position_dodge(),
      #                 color = "BA_Group", 
      #                 # id = "patient_id", line.color = "gray", line.size = 0.4,
      #                 palette = "npg",
      #                 facet.by = c( x.name), short.panel.labs = T,
      #                 ncol = length(unique(df1[[x.name]])),
      #                 xlab = "", ylab = "Percentage of Cells")+
      #   #stat_compare_means(method = "anova", label.y = 90 ) +
      #   stat_compare_means(aes_string(group= "BA_Group"), 
      #                      method = "t.test", 
      #                      # paired = T,
      #                      label = "p.format", size = 4, 
      #                      #hide.ns = TRUE, show.legend = F, 
      #                      label.y = (max(df1[["composition"]])*1.1), )+
      #   theme(
      #     axis.text.x=element_blank(),
      #     #axis.text.x = element_text(angle = 0, hjust = 0.5),
      #     axis.ticks.x=element_blank()
      #   )
      # print(p1)
      # png(filename = paste0(dir1, file.name, "_", x.name, "_", "NoTreatment_LeanvsOver", ".barplot.png"), 
      #     width = length(unique(df1[[x.name]]))*length(unique(df1[["BA_Group"]]))/2, 
      #     height = 4, res = 300, units = "in")
      # print(p1)
      # dev.off()
      # 
      # 
      # df1 <- Composition_3%>% 
      #   dplyr::filter(BA_Group %in% c("Older_Overweight", "Older_PreD_Obese"))
      # p1 <- ggbarplot(df1, 
      #                 x = "BA_Group", y = "composition",
      #                 add = c("mean_se", "jitter"),
      #                 position = position_dodge(),
      #                 color = "BA_Group", 
      #                 # id = "patient_id", line.color = "gray", line.size = 0.4,
      #                 palette = "npg",
      #                 facet.by = c( x.name), short.panel.labs = T,
      #                 ncol = length(unique(df1[[x.name]])),
      #                 xlab = "", ylab = "Percentage of Cells")+
      #   #stat_compare_means(method = "anova", label.y = 90 ) +
      #   stat_compare_means(aes_string(group= "BA_Group"), 
      #                      method = "t.test", 
      #                      # paired = T,
      #                      label = "p.format", size = 4, 
      #                      #hide.ns = TRUE, show.legend = F, 
      #                      label.y = (max(df1[["composition"]])*1.1), )+
      #   theme(
      #     axis.text.x=element_blank(),
      #     #axis.text.x = element_text(angle = 0, hjust = 0.5),
      #     axis.ticks.x=element_blank()
      #   )
      # print(p1)
      # png(filename = paste0(dir1, file.name, "_", x.name, "_", "NoTreatment_OvervsPreD", ".barplot.png"), 
      #     width = length(unique(df1[[x.name]]))*length(unique(df1[["BA_Group"]]))/2, 
      #     height = 4, res = 300, units = "in")
      # print(p1)
      # dev.off()
      
      df1 <- Composition_3%>% 
        dplyr::filter(BA_Group %in% c("Older_Lean", "Older_Overweight", "Older_PreD_Obese"))
      p1 <- ggplot(df1, aes(x = .data[[x.name]] , y = .data[["composition"]], group =  .data[[x.name]])) + 
        geom_boxplot( size=0.5) +
        geom_point(aes(color= .data[["sample_id"]]), size=2) +
        labs(title="Cell type distribution", y = "Average percentage of all nuclei") +
        ylim(0, max(df1[["composition"]])) +
        theme_classic()+
        theme(panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
              axis.title.x = element_blank(),
              axis.title.y = element_text(color = "black", size = 13, face = "plain"),
              axis.text.x = element_text(color = "black", size = 11.5, face = "plain", angle = 45, hjust = 1),
              plot.title = element_text(color = "black", size = 14, face = "bold", hjust = 0.5),
              #legend.text = element_text(color = "black", size = 10.5, face = "plain"),
              #legend.position = c(0.9, 0.8),
              #legend.position = c(.99, .99),
              legend.position = "none",
              #legend.justification = c("right", "top"),
              #legend.box.just = "right",
              legend.spacing.y = unit(0.15, 'cm') ,
              legend.key.size = unit(0.45, "cm"),
              legend.background = element_rect( fill = "grey98", color = "grey98", linewidth  = 1)
        ) 
      print(p1)
      png(file=paste0(dir1, file.name, "_", x.name, "_", "NoTreatment_3Groups", ".Boxplot.png"), 
          width=3.7, height=4, res = 300, units = "in")
      print(p1)
      dev.off()
      
      df1 <- Composition_3%>% 
        dplyr::filter(BA_Group %in% c("Older_Lean", "Older_Overweight", "Older_PreD_Obese"))
      p1 <- ggbarplot(df1, 
                      x = "BA_Group", y = "composition",
                      add = c("mean_se", "jitter"),
                      position = position_dodge(),
                      color = "BA_Group", 
                      # id = "patient_id", line.color = "gray", line.size = 0.4,
                      palette = "npg",
                      facet.by = c( x.name), short.panel.labs = T, # scale = "free",
                      ncol = length(unique(df1[[x.name]])),
                      xlab = "", ylab = "Percentage of Cells")+
        #stat_compare_means(method = "anova", label.y = 90 ) +
        stat_compare_means(aes_string(group= "BA_Group"), 
                           method = "anova" , # "t.test", 
                           # paired = T,
                           label = "p.format", size = 4, 
                           #hide.ns = TRUE, show.legend = F, 
                           label.y = (max(df1[["composition"]])*1.1), )+
        theme(
          axis.text.x=element_blank(),
          #axis.text.x = element_text(angle = 0, hjust = 0.5),
          axis.ticks.x=element_blank()
        )
      print(p1)
      png(
        filename = paste0(dir1, file.name, "_", x.name, "_", "NoTreatment_3Groups", ".barplot.png"), 
          width = length(unique(df1[[x.name]]))*length(unique(df1[["BA_Group"]]))/2, 
          height = 4, res = 300, units = "in")
      print(p1)
      dev.off()
      
      
      # df1 <- Composition_3%>% 
      #   dplyr::filter(BA_Group %in% c("Older_PreD_Obese", "Diabetes"))
      # p1 <- ggbarplot(df1, 
      #                 x = "BA_Group", y = "composition",
      #                 add = c("mean_se", "jitter"),
      #                 position = position_dodge(),
      #                 color = "BA_Group", 
      #                 # id = "patient_id", line.color = "gray", line.size = 0.4,
      #                 palette = "npg",
      #                 facet.by = c( x.name), short.panel.labs = T,
      #                 ncol = length(unique(df1[[x.name]])),
      #                 xlab = "", ylab = "Percentage of Cells")+
      #   #stat_compare_means(method = "anova", label.y = 90 ) +
      #   stat_compare_means(aes_string(group= "BA_Group"), 
      #                      method = "t.test", 
      #                      # paired = T,
      #                      label = "p.format", size = 4, 
      #                      #hide.ns = TRUE, show.legend = F, 
      #                      label.y = (max(df1[["composition"]])*1.1), )+
      #   theme(
      #     #axis.text.x=element_blank(),
      #     axis.text.x = element_text(angle = 0, hjust = 0.5),
      #     axis.ticks.x=element_blank()
      #   )
      # print(p1)
      # png(filename = paste0(dir1, file.name, "_", x.name, "_", "NoTreatment_PreDvsD", ".barplot.png"), 
      #     width = length(unique(df1[[x.name]]))*length(unique(df1[["BA_Group"]]))/2, 
      #     height = 4, res = 300, units = "in")
      # print(p1)
      # dev.off()
      
      # plot each one
      library(car)
      df1 <- Composition_3 %>% data.frame() %>% 
        group_by_at(vars(one_of(c(x.name, "BA_Group"))))   %>%
        # nest() %>%
        dplyr::filter(n() > 2 ) %>%
        ungroup() %>% as.data.frame() %>% 
        group_by_at(vars(one_of(c(x.name))))   %>%
        dplyr::filter(n_distinct(BA_Group) >1)%>% 
        dplyr::filter(BA_Group %in% c("Older_Lean", "Older_Overweight", "Older_PreD_Obese")) %>%
        mutate(BA_Group = factor(BA_Group, levels = c("Older_Lean", "Older_Overweight", "Older_PreD_Obese")))
      
      for (gg in unique(df1[[x.name]]) ) {
        df_test1 <- df1 %>% 
          dplyr::filter(get(x.name) == gg)
        
        if (sum(df_test1$composition) > 0){
          
        
        ano.result <- df_test1  %>%
          as.data.frame() %>%
          rstatix::anova_test(composition ~ BA_Group, type = 3) %>% 
          adjust_pvalue(method = "bonferroni") %>%
          add_significance("p.adj") %>% 
          ungroup() 
        
        if (ano.result$p <= 0.06) {
          test.result <- df_test1%>%
            as.data.frame() %>% # nest() %>% pluck("data", 3) %>% 
            rstatix::t_test(composition ~ BA_Group) %>%
            adjust_pvalue(method = "bonferroni") %>%
            add_significance("p.adj") %>% 
            mutate(method = "T-test")%>% 
            add_xy_position(x = "BA_Group")

          # Plot 1 way ANOVA
          p <-ggline(df_test1,
                     x ="BA_Group" , y = "composition", color = "BA_Group",
                     add = c("mean_se", "dotplot"),
                     palette = "BA_Group",
                     xlab = "", ylab = "Percentage of Cells"
                     # ,ylim = c(0, max(df_test1[["Avg"]])*1)
          ) + theme_bw() +
            theme(
              axis.text.x=element_blank()
              , legend.position="right"
              , legend.title = element_blank()
            ) +
            stat_pvalue_manual(
              test.result, label = "p.adj.signif", # "p = {p.adj}",
              hide.ns = TRUE,
              tip.length = 0.01
            )+
            labs(
              title = gg,
              subtitle = get_test_label(ano.result, detailed = TRUE),
              caption = get_pwc_label(test.result)
            )
          print(p)
          png(
            filename = paste0(dir1, file.name, "_", x.name, "_", gg, "_NoTreatment_3Groups", ".ggline.png"), 
              width = 6, 
              height = 4.5, res = 300, units = "in")
          print(p)
          dev.off()
        }}
      }
      rm(gg)
      
      
      
    }

    
    if (x.name != label.name){
      for (ll in unique(Composition_3[[label.name]]) ) {
        df <- Composition_3 %>% dplyr::filter((!!sym(label.name)) == ll)
        # plotting
        p1 <- ggbarplot(df, x = "BA_Group", y = "composition",
                        add = c("mean_se", "jitter"),
                        position = position_dodge(),
                        color = "BA_Group", 
                        # id = "patient_id", line.color = "gray", line.size = 0.4,
                        palette = "npg",
                        facet.by = c(x.name), short.panel.labs = T,
                        ncol = length(unique(df[[x.name]]))*length(unique(df[["BA_Group"]])),
                        xlab = "", ylab = "Percentage of Cells")+
          #stat_compare_means(method = "anova", label.y = 90 ) +
          stat_compare_means(aes_string(group= "BA_Group"), 
                             method = "anova", 
                             # paired = T,
                             label = "p.format", size = 4, 
                             #hide.ns = TRUE, show.legend = F, 
                             label.y = (max(df[["composition"]])*1.1 ), )+
          theme(
            axis.text.x=element_blank(),
            #axis.text.x = element_text(angle = 0, hjust = 0.5),
            axis.ticks.x=element_blank()
          )
        print(p1)
        png(filename = paste0(dir1, file.name, "_", x.name, "_", ll, "_", "NoTreatment_4Groups", ".barplot.png"), 
            width = length(unique(df[[x.name]]))*length(unique(df[["BA_Group"]]))/2, 
            height = 4, res = 300, units = "in")
        print(p1)
        dev.off()
        
        # df1 <- Composition_3%>% 
        #   dplyr::filter(BA_Group %in% c("Older_Lean", "Older_Overweight"))
        # p1 <- ggbarplot(df1, 
        #                 x = "BA_Group", y = "composition",
        #                 add = c("mean_se", "jitter"),
        #                 position = position_dodge(),
        #                 color = "BA_Group", 
        #                 # id = "patient_id", line.color = "gray", line.size = 0.4,
        #                 palette = "npg",
        #                 facet.by = c( x.name), short.panel.labs = T,
        #                 ncol = length(unique(df1[[x.name]])),
        #                 xlab = "", ylab = "Percentage of Cells")+
        #   #stat_compare_means(method = "anova", label.y = 90 ) +
        #   stat_compare_means(aes_string(group= "BA_Group"), 
        #                      method = "t.test", 
        #                      # paired = T,
        #                      label = "p.format", size = 4, 
        #                      #hide.ns = TRUE, show.legend = F, 
        #                      label.y = (max(df1[["composition"]])*1.1), )+
        #   theme(
        #     #axis.text.x=element_blank(),
        #     axis.text.x = element_text(angle = 0, hjust = 0.5),
        #     axis.ticks.x=element_blank()
        #   )
        # print(p1)
        # png(filename = paste0(dir1, file.name, "_", x.name, "_", "NoTreatment_LeanvsOver", ".barplot.png"), 
        #     width = length(unique(df1[[x.name]]))*length(unique(df1[["BA_Group"]]))/2, 
        #     height = 4, res = 300, units = "in")
        # print(p1)
        # dev.off()
        
        
        # df1 <- Composition_3%>% 
        #   dplyr::filter(BA_Group %in% c("Older_Overweight", "Older_PreD_Obese"))
        # p1 <- ggbarplot(df1, 
        #                 x = "BA_Group", y = "composition",
        #                 add = c("mean_se", "jitter"),
        #                 position = position_dodge(),
        #                 color = "BA_Group", 
        #                 # id = "patient_id", line.color = "gray", line.size = 0.4,
        #                 palette = "npg",
        #                 facet.by = c( x.name), short.panel.labs = T,
        #                 ncol = length(unique(df1[[x.name]])),
        #                 xlab = "", ylab = "Percentage of Cells")+
        #   #stat_compare_means(method = "anova", label.y = 90 ) +
        #   stat_compare_means(aes_string(group= "BA_Group"), 
        #                      method = "t.test", 
        #                      # paired = T,
        #                      label = "p.format", size = 4, 
        #                      #hide.ns = TRUE, show.legend = F, 
        #                      label.y = (max(df1[["composition"]])*1.1), )+
        #   theme(
        #     #axis.text.x=element_blank(),
        #     axis.text.x = element_text(angle = 0, hjust = 0.5),
        #     axis.ticks.x=element_blank()
        #   )
        # print(p1)
        # png(filename = paste0(dir1, file.name, "_", x.name, "_", "NoTreatment_OvervsPreD", ".barplot.png"), 
        #     width = length(unique(df1[[x.name]]))*length(unique(df1[["BA_Group"]]))/2, 
        #     height = 4, res = 300, units = "in")
        # print(p1)
        # dev.off()
        
        # df1 <- Composition_3%>% 
        #   dplyr::filter(BA_Group %in% c("Older_PreD_Obese", "Diabetes"))
        # p1 <- ggbarplot(df1, 
        #                 x = "BA_Group", y = "composition",
        #                 add = c("mean_se", "jitter"),
        #                 position = position_dodge(),
        #                 color = "BA_Group", 
        #                 # id = "patient_id", line.color = "gray", line.size = 0.4,
        #                 palette = "npg",
        #                 facet.by = c( x.name), short.panel.labs = T,
        #                 ncol = length(unique(df1[[x.name]])),
        #                 xlab = "", ylab = "Percentage of Cells")+
        #   #stat_compare_means(method = "anova", label.y = 90 ) +
        #   stat_compare_means(aes_string(group= "BA_Group"), 
        #                      method = "t.test", 
        #                      # paired = T,
        #                      label = "p.format", size = 4, 
        #                      #hide.ns = TRUE, show.legend = F, 
        #                      label.y = (max(df1[["composition"]])*1.1), )+
        #   theme(
        #     #axis.text.x=element_blank(),
        #     axis.text.x = element_text(angle = 0, hjust = 0.5),
        #     axis.ticks.x=element_blank()
        #   )
        # print(p1)
        # png(filename = paste0(dir1, file.name, "_", x.name, "_", "NoTreatment_PreDvsD", ".barplot.png"), 
        #     width = length(unique(df1[[x.name]]))*length(unique(df1[["BA_Group"]]))/2, 
        #     height = 4, res = 300, units = "in")
        # print(p1)
        # dev.off()
        
      }
    }
  }
}

rm(dir1, dir0, p1, files, file,
   file.name, Annotation.file.name, Subtype.file.name,
   celltype.name, list.name, label.name, x.name,
   Composition_3, Composition,
   clusters.lists, Composition.lists, clusters
)








