
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
      # col.group1 = col.groups1[2]
      
      options(dplyr.summarise.inform = FALSE)
      df <- metadata.all %>% 
        dplyr::filter(Treatment_Group %in% c("Notreatment", "NUT_Pre", "SGLT2i_Pre"),
                      BA_Group %in% c("Older_Lean", "Older_Overweight", "Older_PreD_Obese", "Diabetes")#
        ) %>% 
        mutate(BA_Group = factor(BA_Group, levels = c("Older_Lean", "Older_Overweight", "Older_PreD_Obese", "Diabetes"))) %>% 
        arrange(vars(one_of(c("BA_Group", col.group, "sample_id")))) %>% 
        group_by_at(vars(one_of(c("BA_Group", col.group, "sample_id")))) %>%
        dplyr::filter(get(col.group1) > 0) %>% # , !is.na(get(col.group1))
        summarise(Avg = mean(get(col.group1), na.rm = TRUE)) 
      my_comparisons <- list(c("Older_Lean", "Older_Overweight"),
                             c("Older_Lean", "Older_PreD_Obese"),
                             c("Older_Lean", "Diabetes"),
                             c("Older_Overweight", "Older_PreD_Obese"),
                             c("Older_Overweight", "Diabetes"),
                             c("Older_PreD_Obese", "Diabetes")
      )
      
      # Comparisons by one way ANOVA
      
      library(rstatix)
      df_test <- df %>% data.frame() %>% 
        group_by_at(vars(one_of(c(col.group, "BA_Group"))))   %>%
        # nest() %>%
        dplyr::filter(n() > 2 ) %>%
        ungroup() %>% as.data.frame() %>% 
        group_by_at(vars(one_of(c(col.group))))   %>%
        dplyr::filter(n_distinct(BA_Group) >1)
      
      # simple t.test
      test.result <- df_test%>%
        # nest() %>%
        # pluck("data", 46) %>% 
        rstatix::t_test(Avg ~ BA_Group) %>%
        adjust_pvalue(method = "bonferroni") %>%
        add_significance("p.adj") %>% 
        rename(Cells = !!sym(col.group)) %>%
        mutate(method = "T-test")%>% 
        add_xy_position(x = "BA_Group") 
      t.test.results <- rbind(t.test.results, test.result )
      # rm(test.result)
      
      # 2 ANOVA
      ano.result <- df_test %>%
        # as.data.frame() %>%
        rstatix::anova_test(Avg ~ BA_Group, type = 3) %>%
        adjust_pvalue(method = "bonferroni") %>%
        add_significance("p.adj") %>% 
        ungroup() %>% as.data.frame() %>% 
        rename(Cells = !!sym(col.group)) %>%
        mutate(stats.type = "ANOVA_III",
               file.type.1 = file.name,
               file.type.2 = file.type,
               file.type.3 = stringr::str_split(file.type1, "_")[[1]][1],
               Y.value = col.group1,
               Factor1 =  "BA_Group",
               cell.group.type = col.group, 
               # cell.group.name = gg,
               sig = ifelse(.data[["p"]] < 0.05 & !is.na(.data[["p"]]) , "Yes", "ns"),
               Test.type = "2ANOVA"
        )
      ano.results <- rbind(ano.results, ano.result)
      
      library(car)
      for (gg in unique(df_test[[col.group]]) ) {
        df_test1 <- df_test %>% dplyr::filter(get(col.group) == gg)
        ano.result <- df_test1  %>%
          as.data.frame() %>%
          rstatix::anova_test(Avg ~ BA_Group, type = 3) %>% 
          adjust_pvalue(method = "bonferroni") %>%
          add_significance("p.adj") %>% 
          ungroup() 
        
        if (ano.result$p < 0.05) {
          test.result <- df_test1%>%
            rstatix::t_test(Avg ~ BA_Group) %>%
            adjust_pvalue(method = "bonferroni") %>%
            add_significance("p.adj") %>% 
            mutate(method = "T-test")%>% 
            add_xy_position(x = "BA_Group") 
          
          p <-ggline(df_test1,
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
          png(filename = paste0(dir1, file.name, "_",col.group,"_", col.group1, "_",gg,".4Groups.postive.ggline.avg.png"), 
              width = 6, 
              height = 4.5, res = 300, units = "in")
          print(p)
          dev.off()
          
        }

      }
      rm(gg)
      
      # Comparisons by one way ANOVA
      p0 <- ggbarplot(df, x = "BA_Group", y = "Avg",
                      #label = TRUE, label.pos = "out", # lab.vjust = -1.6,
                      add = c("mean_se", "jitter"),
                      position = position_dodge(),
                      color ="BA_Group", shape = "BA_Group", 
                      palette = "BA_Group",
                      facet.by = c(col.group), short.panel.labs = T, # scales="free",
                      ncol = length(unique(metadata.all[[col.group]])),
                      xlab = "", ylab = col.group1,
                      ylim = c(0, max(df[["Avg"]])*1.24)
      ) +
        stat_compare_means(method = "anova", size = 3, label.y = max(df[["Avg"]]) *1.237 ) + 
        theme_bw() +
        theme(
          axis.text.x=element_blank(),
          #axis.text.x = element_text(angle = 45, hjust = 1),
          axis.ticks.x=element_blank()
          , legend.title = element_blank()
        )
      print(p0)
      png(filename = paste0(dir1, file.name, "_",col.group,"_", col.group1, ".4Groups.postive.bar.avg_1ANOVA.png"), 
          width = (length(unique(df[[col.group]]))*1.3 + 0.8), 
          height = 4, res = 300, units = "in")
      print(p0)
      dev.off()
      
      print(paste(col.group, col.group1, sep = ":"))
    }
  }
}
save(ano.results,  t.test.results, df, df_test,
     file = paste0(Disk, Project.folder, "/", Rds.folder, "/", 
                                  stringr::str_split(file.type1, "/")[[1]][1], "_ModuleScoreAve_Results.RData"))
rm(dir1, p0, seurat.file , df, my_comparisons,
   file.name, Annotation.file.name, Subtype.file.name,
   sample.df,  sample.meta, Subtype.meta, Subtype.meta1, 
   metadata.all, row.groups.new, row.group.new, 
   col.group, col.group1
   
)








