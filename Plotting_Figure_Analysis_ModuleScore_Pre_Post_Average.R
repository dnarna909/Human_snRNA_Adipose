
# rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.
memory.limit(size = 1e+13)

# import parameters and packages --------------------------------------------------------------------------------------------
# source(paste0("/media/jianie/Extreme SSD1/2022-09-01 STARR_SGLT2 Combine/Project Parameters.R"), local = knitr::knit_global())
source(paste0(paste0(Disk, Project.folder, "/", "Project Parameters.R")), local = knitr::knit_global())
library("RColorBrewer")
library("ggpubr")

# Import data -----------------------------------------------------------------------------------------------------------------------------------------------------
# export.folder <- "Compositions"
# dir.create(file.path(paste0(Disk, Project.folder, "/", Rds.folder, "/"), export.folder), showWarnings = FALSE)

dir.create(file.path(paste0(Disk, Project.folder, "/", figures.folder, "/", export.folder, ".plotting")), showWarnings = FALSE)
dir1 <- paste0(Disk, Project.folder, "/", figures.folder, "/", export.folder, ".plotting", "/")


ano.results <- data.frame()
t.test.results <- data.frame()
Plot2W.list <- list()

for (seurat.file in seurat.files) {
  # seurat.file = seurat.files[1]
  
  file.name <- stringr::str_split(seurat.file, ".Rds")[[1]][1]
  Annotation.file.name <- paste(stringr::str_split(file.name, "_")[[1]][c(1,2,3)], collapse = "_")
  Subtype.file.name <- stringr::str_split(file.name, paste0(Annotation.file.name, "_"))[[1]][2]
  if (is.na(Subtype.file.name)) {Subtype.file.name = ""}
  
  # Import sample data
  sample.df <- readRDS(paste0(Disk, Project.folder, "/",  Rds.folder, "/", "Samples.df.Rds")) 
  
  # Import metadata
  sample.meta <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/", 
                                stringr::str_split(file.name, "_")[[1]][1], "_Sample.meta.data.Rds")) %>% 
    tibble::rownames_to_column(var = "rowname")
  
  Subtype.meta <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/", 
                                 file.name, "_", file.type)) %>% 
    tibble::rownames_to_column(var = "rowname")
  
  Subtype.meta1 <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/",  
                                  file.type1)) %>% 
    tibble::rownames_to_column(var = "rowname")
  
  metadata.all <- left_join(Subtype.meta, sample.meta, by = "rowname") %>% 
    left_join(Subtype.meta1, by = "rowname")
  
  row.groups.new <- paste(row.groups, col.groups, sep = "_")
  for (row.group.new in row.groups.new) {
    metadata.all[[row.group.new]] <-
      paste(metadata.all[[stringr::str_split(row.group.new, "_")[[1]][1]]], 
            metadata.all[[stringr::str_split(row.group.new, "_")[[1]][2]]], sep = "_")
  }
  
  for (col.group in col.groups) {
    # col.group = col.groups[2]
    for (col.group1 in col.groups1) { 
      # col.group1 = col.groups1[1]
      
      options(dplyr.summarise.inform = FALSE)
      df <- metadata.all %>% 
        dplyr::filter(Pre_Post %in% c("Pre", "Post") ) %>% 
        dplyr::filter(!patient_id %in% c("76615")) %>% 
        arrange(vars(one_of(c("Pre_Post", "Group", col.group, "patient_id")))) %>% 
        group_by_at(vars(one_of(c("Pre_Post", "Group", col.group, "patient_id")))) %>%
        dplyr::filter(get(col.group1) > 0) %>% # , !is.na(get(col.group1))
        summarise(Avg = mean(get(col.group1), na.rm = TRUE)) %>% 
        mutate(Group_Event = factor(paste0(Group, "_", Pre_Post), levels = c("NUT_Pre", "NUT_Post",
                                                                             "SGLT2i_Pre", "SGLT2i_Post")),
               Group = factor(Group, levels = c("NUT", "SGLT2i")),
               Pre_Post = factor(Pre_Post, levels = c("Pre", "Post")))
      my_comparisons <- list(c("NUT_Pre", "NUT_Post"),
                             c("SGLT2i_Pre", "SGLT2i_Post"),
                             c("NUT_Pre", "SGLT2i_Pre"),
                             c("NUT_Post", "SGLT2i_Post"))
      
      # Comparisons by two way ANOVA
      library(rstatix)
      
      # # simple t.test
      test.result <- df %>% data.frame() %>% 
        group_by_at(vars(one_of(c(col.group,  "Group_Event"))))   %>%
        # nest() %>%
        dplyr::filter(n() > 2 ) %>%
        ungroup() %>% as.data.frame() %>% 
        group_by_at(vars(one_of(c(col.group))))   %>%
        dplyr::filter(n_distinct(Group_Event) >1)  %>%
        rstatix::t_test(Avg ~ Group_Event) %>%
        add_xy_position(x = "Group_Event") %>%
        adjust_pvalue(method = "bonferroni") %>%
        rename(Cells = !!sym(col.group)) %>%
        mutate(method = "T-test", Group = "Group_Event")%>% 
        add_significance("p.adj")
      t.test.results <- rbind(t.test.results, test.result )
      # 
      # 
      # # paired t.test by patient ID, where is subject id?
      test.result <- df %>% data.frame() %>% 
        group_by_at(vars(one_of(c(col.group,  "Group", "Pre_Post"))))   %>%
        dplyr::filter(n() > 2 ) %>%
        ungroup() %>% as.data.frame() %>% 
        group_by_at(vars(one_of(c(col.group,  "patient_id"))))   %>%
        dplyr::filter(n() > 1 ) %>%
        ungroup() %>% as.data.frame() %>% 
        group_by_at(vars(one_of(c(col.group, "Group"))))   %>%
        dplyr::filter(n_distinct(Pre_Post) >1)  %>%
        rename(id = patient_id) %>%
        arrange(!!sym(col.group), Group, id, Pre_Post) %>% #nest() %>% pluck("data", 59) %>% 
        rstatix::t_test(Avg ~ Pre_Post, paired = TRUE) %>%
        add_xy_position(x = "Pre_Post")%>%
        adjust_pvalue(method = "bonferroni") %>%
        mutate(method = "paired.T-test")%>% 
        add_significance("p.adj")  %>% 
        rename(Cells = !!sym(col.group))
      t.test.results <- rbind(t.test.results, test.result )
      
      # # 2 ANOVA
      ano.result <- df %>% data.frame() %>% 
        group_by_at(vars(one_of(c(col.group, "Pre_Post"))))   %>%
        dplyr::filter(n_distinct(Group) >1 ) %>%
        ungroup() %>% as.data.frame() %>% 
        group_by_at(vars(one_of(c(col.group, "Group"))))   %>%
        dplyr::filter(n_distinct(Pre_Post) >1)  %>%
        ungroup() %>% as.data.frame() %>% 
        group_by_at(vars(one_of(c(col.group))))   %>%
        dplyr::filter(n_distinct(Group_Event) >1)  %>%# nest()  %>% pluck("data", 59) %>%
        rstatix::anova_test(Avg ~ Group*Pre_Post , type = 3) %>%
        adjust_pvalue(method = "bonferroni") %>%
        add_significance("p.adj") %>% 
        ungroup() %>% as.data.frame() %>% 
        rename(Cells = !!sym(col.group)) %>%
        mutate(stats.type = "ANOVA_III",
               file.type.1 = file.name,
               file.type.2 = file.type,
               file.type.3 = stringr::str_split(file.type1, "_")[[1]][1],
               Y.value = col.group1,
               Factor1 = "Treatment.Group",
               Factor2 = "Pre_Post",
               cell.group.type = col.group, 
               # cell.group.name = gg,
               sig = ifelse(.data[["p"]] < 0.05 & !is.na(.data[["p"]]) , "Yes", "ns"),
               Test.type = "2ANOVA"
        )
      ano.results <- rbind(ano.results, ano.result)
      
      # Mixed  ANOVA by patient ID
      # ano.result <- df %>% data.frame() %>% 
      #   group_by_at(vars(one_of(c(col.group,  "Group", "Pre_Post"))))   %>%
      #   dplyr::filter(n() > 2 ) %>%
      #   ungroup() %>% as.data.frame() %>% 
      #   group_by_at(vars(one_of(c(col.group))))   %>%
      #   dplyr::filter(n_distinct(Group_Event) >1)  %>%
      #   arrange(!!sym(col.group), Group, patient_id, Pre_Post) %>%
      #   rename(id = patient_id, Event = Pre_Post) %>%
      #   group_by_at(vars(one_of(c(col.group))))  %>%
      #   # nest() %>%
      #   # pluck("data", 1) %>% 
      #  rstatix::anova_test(Avg ~ Group*Event + Error(id/(Event)) ) %>% # error
      #    # rstatix::anova_test(dv = Avg, wid = id, within = c(Event)) %>% #
      #   #add_xy_position(x = "Pre_Post") %>%
      #   adjust_pvalue(method = "bonferroni") %>%
      #   add_significance("p.adj")
      
      ano.result <- df %>% data.frame() %>% 
        group_by_at(vars(one_of(c(col.group, "patient_id"))))   %>%
        dplyr::filter(n() >1 ) %>%
        ungroup() %>% as.data.frame() %>% 
        group_by_at(vars(one_of(c(col.group, "Pre_Post"))))   %>%
        dplyr::filter(n_distinct(Group) >1 ) %>%
        ungroup() %>% as.data.frame() %>% 
        group_by_at(vars(one_of(c(col.group, "Group"))))   %>%
        dplyr::filter(n_distinct(Pre_Post) >1)  %>%
        ungroup() %>% as.data.frame() %>% 
        group_by_at(vars(one_of(c(col.group, "Group", "Pre_Post"))))   %>%
        dplyr::filter(n() >1 ) %>%
        ungroup() %>% as.data.frame() %>% 
        group_by_at(vars(one_of(c(col.group))))   %>%
        dplyr::filter(n_distinct(Group_Event) >1)  %>%  # nest()  %>% pluck("data", 59) %>%
        arrange(!!sym(col.group), Group, patient_id, Pre_Post) %>%
        rstatix::anova_test(dv = Avg, wid = patient_id,
                            between = c(Group), within = c(Pre_Post), type = 3) %>% #
        #add_xy_position(x = "Pre_Post") %>%
        adjust_pvalue(method = "bonferroni") %>%
        add_significance("p.adj")%>%
        ungroup() %>% as.data.frame() %>%
        rename(Cells = !!sym(col.group)) %>%
        mutate(stats.type = "ANOVA_III",
               file.type.1 = file.name,
               file.type.2 = file.type,
               file.type.3 = stringr::str_split(file.type1, "_")[[1]][1],
               Y.value = col.group1,
               Factor1 = "Treatment.Group",
               Factor2 = "Pre_Post",
               cell.group.type = col.group,
               # cell.group.name = gg,
               sig = ifelse(.data[["p.adj"]] < 0.05 & !is.na(.data[["p.adj"]]) , "Yes", "ns"),
               Test.type = "Mixed ANOVA"
        )
      ano.results <- rbind(ano.results, ano.result)
      
      library(car)
      df_test <- df %>% data.frame() %>% 
        group_by_at(vars(one_of(c(col.group, "Pre_Post"))))   %>%
        dplyr::filter(n_distinct(Group) >1 ) %>%
        ungroup() %>% as.data.frame() %>% 
        group_by_at(vars(one_of(c(col.group, "Group"))))   %>%
        dplyr::filter(n_distinct(Pre_Post) >1)  %>%
        ungroup() %>% as.data.frame() %>% 
        group_by_at(vars(one_of(c(col.group))))   %>%
        dplyr::filter(n_distinct(Group_Event) >1) 
      
      for (gg in unique(df_test[[col.group]]) ) {
        df_test1 <- df_test %>%
          dplyr::filter(get(col.group) == gg)
        
        ano.result <- df_test1  %>%
          as.data.frame() %>%
          rstatix::anova_test(Avg ~ Group*Pre_Post, type = 3) %>% 
          adjust_pvalue(method = "bonferroni") %>%
          add_significance("p.adj") %>% 
          ungroup() 
        
        if ( ano.result$p[1] <= 0.05 | ano.result$p[2] <= 0.05 ) {
          ano.result2 <- ano.result [c((ano.result$p[1:2] < 0.05), FALSE),]
          print(gg)
          test.result <- df_test1%>%
            as.data.frame() %>% 
            group_by_at(vars(one_of(c(col.group,  "Group"))))   %>%
            rstatix::t_test(Avg ~ Pre_Post) %>%
            adjust_pvalue(method = "bonferroni") %>%
            add_significance("p.adj") %>% 
            mutate(method = "T-test")%>% 
            add_xy_position(x = "Pre_Post")
          
          # Plot 2 way ANOVA
          library(CGPfunctions)
          library(pwr)
          # Plot2WayANOVA(Avg ~ Group*Pre_Post, df_test1, plottype = "line")
          # Plot2WayANOVA(Avg ~ Group*Pre_Post, df_test1, plottype = "line",
          #               overlay.type = "box",
          #               mean.label = TRUE)
          graphics.off()
          Plot2W.list[[paste0(col.group,"_", col.group1, "_",gg)]] <- Plot2WayANOVA(Avg ~ Group*Pre_Post, df_test1, plottype = "line", 
                                                                                    title = gg, ylab = col.group1,
                                                                                    confidence = .99,
                                                                                    ggplot.component = theme(axis.text.x = element_text(size=13, color="darkred")))
          ggsave(filename = paste0(dir1, file.name, "_",col.group,"_", col.group1, "_",gg,".SGLT2i.postive.avg_2ANOVA.Plot2Way.png"), 
                 width = 6, 
                 height = 4.5, dpi = 300, units = "in")
          graphics.off()
          
          # Plot 2 way ANOVA
          p <-ggline(df_test1,
                     x ="Pre_Post" , y = "Avg", color = "Group",
                     add = c("mean_se", "dotplot"),
                     palette = "Pre_Post",
                     xlab = "", ylab = col.group1
                     # ,ylim = c(0, max(df_test1[["Avg"]])*1)
          ) + theme_bw() +
            theme(
              # axis.text.x=element_blank(),
              legend.position="right"
              , legend.title = element_blank()
            ) +
            stat_pvalue_manual(
              test.result, label = "p.adj.signif", # "p = {p.adj}",
              hide.ns = TRUE,
              tip.length = 0.01
            )+
            labs(
              title = paste0(gg, ", ",ano.result2$Effect),
              subtitle = get_test_label(ano.result2, detailed = TRUE),
              caption = get_pwc_label(test.result)
            )
          print(p)
          png(filename = paste0(dir1, file.name, "_",col.group,"_", col.group1, "_",gg,".SGLT2i.postive.ggline.avg_2ANOVA.png"),
              width = 6,
              height = 4.5, res = 300, units = "in")
          print(p)
          dev.off()
        }
      }
      rm(gg)
      
      p0 <- ggline(df, 
                   x ="Pre_Post" , y = "Avg", color = "Group",
                   add = c("mean_se", "dotplot"),
                   facet.by = c(col.group), short.panel.labs = T,
                   palette = "Pre_Post",
                   ncol = length(unique(metadata.all[[col.group]])),
                   xlab = "", ylab = col.group1,
                   ylim = c(0, max(df[["Avg"]])*1)
      ) +theme_bw()
      print(p0)
      png(filename = paste0(dir1, file.name, "_",col.group,"_", col.group1, ".SGLT2i.postive.bar.avg_2ANOVA.png"), 
          width = (length(unique(df[[col.group]]))*1.3 + 0.5), 
          height = 4, res = 300, units = "in")
      print(p0)
      dev.off()
      
      # Comparisons by one way ANOVA
      p0 <- ggbarplot(df, x = "Group_Event", y = "Avg",
                      #label = TRUE, label.pos = "out", # lab.vjust = -1.6,
                      add = c("mean_se", "jitter"),
                      position = position_dodge(),
                      color ="Pre_Post", shape = "Group", 
                      palette = "Group",
                      facet.by = c(col.group), short.panel.labs = T,
                      ncol = length(unique(metadata.all[[col.group]])),
                      xlab = "", ylab = col.group1,
                      ylim = c(0, max(df[["Avg"]])*1.24)
      ) +
        stat_compare_means(method = "anova", size = 2, label.y = max(df[["Avg"]]) *1.237 ) + 
        theme_bw() +
        stat_compare_means(mapping=aes(label = format.pval(..p.adj.., digits = 1)), 
                           #aes_string(group="Group"), 
                           comparisons=my_comparisons,
                           method = "t.test", 
                           #label = "p.format",                            
                           size = 2, 
                           label.y = c(max(df[["Avg"]]) *1,
                                       max(df[["Avg"]]) *1,
                                       max(df[["Avg"]]) *1.07,
                                       max(df[["Avg"]]) *1.155), 
                           hide.ns = TRUE, # show.legend = F 
        )+theme_bw() +
        theme(
          #axis.text.x=element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.ticks.x=element_blank()
        )
      print(p0)
      png(filename = paste0(dir1, file.name, "_",col.group,"_", col.group1, ".SGLT2i.postive.bar.avg_1ANOVA.png"), 
          width = (length(unique(df[[col.group]]))*1.3 + 0.8), 
          height = 4, res = 300, units = "in")
      print(p0)
      dev.off()
      
      
      p1 <- ggbarplot(df, x = "Group", y = "Avg",
                      #label = TRUE, label.pos = "out", # lab.vjust = -1.6,
                      add = c("mean_se", "jitter"),
                      position = position_dodge(),
                      color = "Group", 
                      palette = "Group",
                      facet.by = c("Pre_Post", col.group), short.panel.labs = T,
                      ncol = length(unique(metadata.all[[col.group]])),
                      xlab = "", ylab = col.group1,
                      ylim = c(0, max(df[["Avg"]])*1.2)
      ) + 
        stat_compare_means(aes_string(group="Group"), 
                           method = "t.test", 
                           label = "p.format", 
                           size = 4, 
                           label.y = (max(df[["Avg"]]) *1), 
                           hide.ns = TRUE, # show.legend = F 
        ) + theme_bw() +
        theme(
          #axis.text.x=element_blank(),
          axis.text.x = element_text(angle = 0, hjust = 0.5),
          axis.ticks.x=element_blank()
        )
      print(p1)
      png(filename = paste0(dir1, file.name, "_",col.group,"_", col.group1, ".SGLT2i.postive.bar.avg_Group.png"), 
          width = (length(unique(df[[col.group]]))*1.3 + 0.8), 
          height = (length(unique(metadata.all[["Pre_Post"]]))*3.5 + 2), res = 300, units = "in")
      print(p1)
      dev.off()
      
      
      p2 <- ggbarplot(df, x = "Pre_Post", y = "Avg",
                      #label = TRUE, label.pos = "out", # lab.vjust = -1.6,
                      add = c("mean_se", "jitter"),
                      position = position_dodge(),
                      color = "Pre_Post", 
                      palette = "Pre_Post",
                      facet.by = c("Group", col.group), short.panel.labs = T,
                      ncol = length(unique(metadata.all[[col.group]])),
                      xlab = "", ylab = col.group1,
                      ylim = c(0, max(df[["Avg"]])*1.1)
      ) +
        stat_compare_means(aes_string(group="Pre_Post"), 
                           method = "t.test", 
                           label = "p.format", size = 4, 
                           label.y = (max(df[["Avg"]]) *1), 
                           hide.ns = TRUE, # show.legend = F 
        )+ theme_bw() +
        theme(
          #axis.text.x=element_blank(),
          axis.text.x = element_text(angle = 0, hjust = 0.5),
          axis.ticks.x=element_blank()
        )
      print(p2)
      png(filename = paste0(dir1, file.name, "_",col.group,"_", col.group1, ".SGLT2i.postive.bar.avg.png"), 
          width = (length(unique(df[[col.group]]))*1.3 + 0.8), 
          height = (length(unique(metadata.all[["Group"]]))*3.5 + 2), res = 300, units = "in")
      print(p2)
      dev.off()
      
      # options(digits=2)
      p3 <- ggpaired(df, x = "Pre_Post", y = "Avg",
                     # label = TRUE, label.pos = "out", # lab.vjust = -1.6,
                     add = c("mean_se", "jitter"),
                     position = position_dodge(),
                     color = "Pre_Post", 
                     id = "patient_id", line.color = "gray", line.size = 0.4,
                     palette = "Pre_Post",
                     facet.by = c("Group", col.group), short.panel.labs = T,
                     ncol = length(unique(metadata.all[[col.group]])),
                     xlab = "", ylab = col.group1,
                     ylim = c(0, max(df[["Avg"]])*1.1)
      ) + 
        stat_compare_means(aes_string(group=c("Pre_Post")), 
                           method = "t.test", paired = T,
                           label = "p.format", size = 4, 
                           hide.ns = TRUE, # show.legend = F, 
                           label.y = (max(df[["Avg"]]) *1) 
        )+ # theme_bw() +
        theme(
          #axis.text.x=element_blank(),
          axis.text.x = element_text(angle = 0, hjust = 0.5),
          axis.ticks.x=element_blank()
        )
      print(p3)
      png(filename = paste0(dir1, file.name, "_",col.group,"_", col.group1, ".SGLT2i.postive.bar.paired.avg.png"), 
          width = (length(unique(df[[col.group]]))*1.3 + 0.8), 
          height = (length(unique(metadata.all[["Group"]]))*3.5 + 2), res = 300, units = "in")
      print(p3)
      dev.off()
      
      print(paste(col.group, col.group1, sep = ":"))
    }
  }
}

save(ano.results, t.test.results, Plot2W.list, 
        file = paste0(Disk, Project.folder, "/", Rds.folder, "/", 
                      stringr::str_split(file.type1, "/")[[1]][1], "_SGLT2i.ModuleScoreAve_2ano_Results.RData"))
rm(dir1, p1, p0, p2, p3, seurat.file , df, my_comparisons,
   file.name, Annotation.file.name, Subtype.file.name,
   sample.df,  sample.meta, Subtype.meta, Subtype.meta1, 
   metadata.all, row.groups.new, row.group.new, 
   col.group, col.group1,
   ano.result, ano.results,
   df_test, df_test1, ano.result2, p
   
   
)








