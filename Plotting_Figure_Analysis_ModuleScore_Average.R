
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
    # col.group = col.groups[1]
    for (col.group1 in col.groups1) { 
      # col.group1 = col.groups1[1]
      
      options(dplyr.summarise.inform = FALSE)
      df <- metadata.all %>% arrange(vars(one_of(c("Event.Name", "Group", col.group, "patient_id")))) %>% 
        group_by_at(vars(one_of(c("Event.Name", "Group", col.group, "patient_id")))) %>%
        dplyr::filter(get(col.group1) > 0) %>% # , !is.na(get(col.group1))
        summarise(Avg = mean(get(col.group1), na.rm = TRUE)) %>% 
        dplyr::filter(!patient_id %in% c("76615")) %>% 
        mutate(Group_Event = factor(paste0(Group, "_", Event.Name), levels = c("C_Visit 5", "C_Visit 11",
                                                                                 "D_Visit 5", "D_Visit 11")))
      my_comparisons <- list(c("C_Visit 5", "C_Visit 11"),
                             c("D_Visit 5", "D_Visit 11"),
                             c("C_Visit 5", "D_Visit 5"),
                             c("C_Visit 11", "D_Visit 11"))
      
      # Comparisons by two way ANOVA
      
      library(rstatix)
      # simple t.test
      test.result <- df %>% data.frame() %>% 
        group_by_at(vars(one_of(c(col.group))))  %>%
        rstatix::t_test(Avg ~ Group_Event) %>% 
        add_xy_position(x = "Event.Name") %>%
        adjust_pvalue(method = "bonferroni") %>%
        add_significance("p.adj")
      
      # t.test by group
      test.result <- df %>% data.frame() %>% 
        group_by_at(vars(one_of(c(col.group, "Event.Name"))))  %>%
        rstatix::t_test(Avg ~ Group) %>% 
        add_xy_position(x = "Event.Name")%>%
        adjust_pvalue(method = "bonferroni") %>%
        add_significance("p.adj")
      
      # paired t.test by patient ID, where is subject id?
      test.result <- df %>% data.frame() %>% rename(id = patient_id) %>%
        group_by_at(vars(one_of(c(col.group, "Group"))))  %>%
        rstatix::t_test(Avg ~ Event.Name, paired = TRUE) %>% 
        add_xy_position(x = "Event.Name")%>%
        adjust_pvalue(method = "bonferroni") %>%
        add_significance("p.adj")
      
      # 2 ANOVA
      ano.result <- df %>% data.frame() %>% 
        group_by_at(vars(one_of(c(col.group))))  %>% 
        rstatix::anova_test(Avg ~ Group*Event.Name, type = 3) %>% 
        #add_xy_position(x = "Event.Name") %>%
        adjust_pvalue(method = "bonferroni") %>%
        add_significance("p.adj")
      
      # paired 2 ANOVA by patient ID
      ano.result <- df %>% data.frame() %>% rename(id = patient_id, Event = Event.Name) %>% 
        group_by_at(vars(one_of(c(col.group))))  %>% 
        rstatix::anova_test(Avg ~ Group*Event+Error(id/(Group*Event))) %>% # 
        #add_xy_position(x = "Event.Name") %>%
        adjust_pvalue(method = "bonferroni") %>%
        add_significance("p.adj")
      
      ano.result <- df %>% data.frame() %>% 
        group_by_at(vars(one_of(c(col.group))))  %>% 
        rstatix::anova_test(dv = Avg, wid = patient_id, 
                            between = c(Group), within = c(Event.Name), type = 3) %>% # 
        #add_xy_position(x = "Event.Name") %>%
        adjust_pvalue(method = "bonferroni") %>%
        add_significance("p.adj")%>% 
        mutate(stats.type = "ANOVA_III",
               file.type.1 = file.name,
               file.type.2 = file.type,
               file.type.3 = stringr::str_split(file.type1, "_")[[1]][1],
               Y.value = col.group1,
               Factor1 = "Treatment.Group",
               Factor2 = "Event.Name",
               cell.group.type = col.group, 
               cell.group.name = gg,
               sig = ifelse(.data[["Pr(>F)"]] < 0.05 & !is.na(.data[["Pr(>F)"]]) , "Yes", "ns")
        )
      
      res.aov <- anova_test(data = df%>% data.frame() %>% 
                              group_by_at(vars(one_of(c(col.group)))) , dv = Avg, wid = patient_id, 
                            between = c(Group), within = c(Event.Name), type = 3)
      get_anova_table(res.aov)
      ggline(df, 
             x ="Event.Name" , y = "Avg", color = "Group",
             add = c("mean_se", "dotplot"),
             facet.by = c(col.group), short.panel.labs = T,
             palette = "Event.Name",
             ncol = length(unique(metadata.all[[col.group]])),
             xlab = "", ylab = col.group1,
             ylim = c(0, max(df[["Avg"]])*1)
      ) +
        stat_pvalue_manual(
          test.result, label = "p.adj",
          #y.position = c(1, 1.2, 1.5)
        )+
        labs(
          subtitle = get_test_label(ano.result, detailed = TRUE),
          caption = get_pwc_label(test.result)
        )
      
      library(car)
      for (gg in unique(df[[col.group]]) ) {
        if(df %>% dplyr::filter(get(col.group) == gg) %>% nrow() > 12){
          my_anova <- aov(Avg ~ Group * Event.Name, 
                          data = df %>% dplyr::filter(get(col.group) == gg) %>% data.frame())
          ano.result <- Anova(my_anova, type = "III") %>% 
            tibble::rownames_to_column(var = "rowname") %>% 
            mutate(stats.type = "ANOVA_III",
                   file.type.1 = file.name,
                   file.type.2 = file.type,
                   file.type.3 = stringr::str_split(file.type1, "_")[[1]][1],
                   Y.value = col.group1,
                   Factor1 = "Treatment.Group",
                   Factor2 = "Event.Name",
                   cell.group.type = col.group, 
                   cell.group.name = gg,
                   sig = ifelse(.data[["Pr(>F)"]] < 0.05 & !is.na(.data[["Pr(>F)"]]) , "Yes", "ns")
            ) 
          ano.results <- rbind(ano.results, ano.result)
          rm(my_anova)
        }
      }
      rm(gg)
      
      p0 <- ggline(df, 
                   x ="Event.Name" , y = "Avg", color = "Group",
                   add = c("mean_se", "dotplot"),
                   facet.by = c(col.group), short.panel.labs = T,
                   palette = "Event.Name",
                   ncol = length(unique(metadata.all[[col.group]])),
                   xlab = "", ylab = col.group1,
                   ylim = c(0, max(df[["Avg"]])*1)
      )
      png(filename = paste0(dir1, file.name, "_",col.group,"_", col.group1, ".postive.bar.avg_2ANOVA.png"), 
          width = (length(unique(df[[col.group]]))*1.3 + 0.5), 
          height = 4, res = 300, units = "in")
      print(p0)
      dev.off()
      
      
      
      # interaction.plot(x.factor = df[["Group"]], trace.factor = df[["Event.Name"]],
      #                  response = df[["Avg"]], fun = mean,
      #                  type = "b", legend = TRUE,
      #                  xlab = "Group", ylab=col.group1,
      #                  pch=c(1,19), col = c("#00AFBB", "#E7B800"))
      # 
      # CGPfunctions::Plot2WayANOVA(Avg ~ Group * Event.Name, df, plottype = "line",
      #               overlay.type = "box",show.dots =T,
      #               mean.label = TRUE)
      # CGPfunctions::Plot2WayANOVA(Avg ~ Group * Event.Name, 
      #                             df %>% dplyr::filter(get(col.group) == "FAP") %>% data.frame()
      #                             , confidence = .99,show.dots =T,
      #                             ggplot.component = theme(axis.text.x = element_text(size=13, color="darkred")))
      
      
      # Comparisons by one way ANOVA
      p0 <- ggbarplot(df, x = "Group_Event", y = "Avg",
                      #label = TRUE, label.pos = "out", # lab.vjust = -1.6,
                      add = c("mean_se", "jitter"),
                      position = position_dodge(),
                      color ="Event.Name", shape = "Group", 
                      palette = "Group",
                      facet.by = c(col.group), short.panel.labs = T,
                      ncol = length(unique(metadata.all[[col.group]])),
                      xlab = "", ylab = col.group1,
                      ylim = c(0, max(df[["Avg"]])*1.24)
      ) +
        stat_compare_means(method = "anova", size = 2, label.y = max(df[["Avg"]]) *1.237 ) + 
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
        )+
        theme(
          #axis.text.x=element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.ticks.x=element_blank()
        )
      png(filename = paste0(dir1, file.name, "_",col.group,"_", col.group1, ".postive.bar.avg_1ANOVA.png"), 
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
                      facet.by = c("Event.Name", col.group), short.panel.labs = T,
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
        )+
        theme(
          #axis.text.x=element_blank(),
          axis.text.x = element_text(angle = 0, hjust = 0.5),
          axis.ticks.x=element_blank()
        )
      png(filename = paste0(dir1, file.name, "_",col.group,"_", col.group1, ".postive.bar.avg_Group.png"), 
          width = (length(unique(df[[col.group]]))*1.3 + 0.8), 
          height = (length(unique(metadata.all[["Event.Name"]]))*3.5 + 2), res = 300, units = "in")
      print(p1)
      dev.off()
      
      
      p2 <- ggbarplot(df, x = "Event.Name", y = "Avg",
                      #label = TRUE, label.pos = "out", # lab.vjust = -1.6,
                      add = c("mean_se", "jitter"),
                      position = position_dodge(),
                      color = "Event.Name", 
                      palette = "Event.Name",
                      facet.by = c("Group", col.group), short.panel.labs = T,
                      ncol = length(unique(metadata.all[[col.group]])),
                      xlab = "", ylab = col.group1,
                      ylim = c(0, max(df[["Avg"]])*1.1)
                      ) +
        stat_compare_means(aes_string(group="Event.Name"), 
                           method = "t.test", 
                           label = "p.format", size = 4, 
                           label.y = (max(df[["Avg"]]) *1), 
                           hide.ns = TRUE, # show.legend = F 
        )+
        theme(
          #axis.text.x=element_blank(),
          axis.text.x = element_text(angle = 0, hjust = 0.5),
          axis.ticks.x=element_blank()
        )
      png(filename = paste0(dir1, file.name, "_",col.group,"_", col.group1, ".postive.bar.avg.png"), 
          width = (length(unique(df[[col.group]]))*1.3 + 0.8), 
          height = (length(unique(metadata.all[["Group"]]))*3.5 + 2), res = 300, units = "in")
      print(p2)
      dev.off()
      
      # options(digits=2)
      p3 <- ggpaired(df, x = "Event.Name", y = "Avg",
                     # label = TRUE, label.pos = "out", # lab.vjust = -1.6,
                     add = c("mean_se", "jitter"),
                     position = position_dodge(),
                     color = "Event.Name", 
                     id = "patient_id", line.color = "gray", line.size = 0.4,
                     palette = "Event.Name",
                     facet.by = c("Group", col.group), short.panel.labs = T,
                     ncol = length(unique(metadata.all[[col.group]])),
                     xlab = "", ylab = col.group1,
                     ylim = c(0, max(df[["Avg"]])*1.1)
                     ) + 
        stat_compare_means(aes_string(group=c("Event.Name")), 
                           method = "t.test", paired = T,
                           label = "p.format", size = 4, 
                           hide.ns = TRUE, # show.legend = F, 
                           label.y = (max(df[["Avg"]]) *1) 
        )+
        theme(
          #axis.text.x=element_blank(),
          axis.text.x = element_text(angle = 0, hjust = 0.5),
          axis.ticks.x=element_blank()
        )
      #print(p1)
      png(filename = paste0(dir1, file.name, "_",col.group,"_", col.group1, ".postive.bar.paired.avg.png"), 
          width = (length(unique(df[[col.group]]))*1.3 + 0.8), 
          height = (length(unique(metadata.all[["Group"]]))*3.5 + 2), res = 300, units = "in")
      print(p3)
      dev.off()
      
      # gghistogram(metadata.all, x = col.group1,
      #             add = "mean", rug = TRUE,
      #             color = "Event.Name", fill = "Event.Name",
      #             facet.by = c(col.group), short.panel.labs = T, 
      #             ncol = length(unique(metadata.all[[col.group]])),
      #             ggtheme = theme_light(), legend = "top")
      print(paste(col.group, col.group1, sep = ":"))
    }
  }
}
saveRDS(ano.results, file = paste0(Disk, Project.folder, "/", Rds.folder, "/", 
                                  stringr::str_split(file.type1, "/")[[1]][1], "_ModuleScoreAve_2ano_Results.Rds"))
rm(dir1, p1, p0, p2, p3, seurat.file , df, my_comparisons,
   file.name, Annotation.file.name, Subtype.file.name,
   sample.df,  sample.meta, Subtype.meta, Subtype.meta1, 
   metadata.all, row.groups.new, row.group.new, 
   col.group, col.group1
   
)








