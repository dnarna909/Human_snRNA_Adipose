
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
      DF <- metadata.all %>% 
        dplyr::filter(Pre_Post %in% c("Pre", "Post") ) %>% 
        dplyr::filter(!patient_id %in% c("76615"))
      p1 <- ggdensity(DF, x = col.group1,
            add = "mean", rug = TRUE,
            add.params = list(size = 2, linetype = "longdash"),
            color = "Pre_Post", fill = "Pre_Post",
            facet.by = c("Group", col.group), short.panel.labs = T, 
            ncol = length(unique(DF[[col.group]])),
            ggtheme = theme_light(), legend = "top"
            )+   
        theme(strip.text = element_text(size = 20, face = "bold")# Change font size
              , axis.text.x = element_text(size=15)
              , legend.text = element_text(size=20)
              , axis.title = element_text(size=20, face = "bold")
              )  
      print(p1)
      png(filename = paste0(dir1, file.name, "_",col.group,"_", col.group1, "_SGLT2i.png"), 
          width = length(unique(DF[[col.group]]))*2, 
          height = length(unique(DF[["Group"]]))*4, res = 300, units = "in")
      print(p1)
      dev.off()
      
      df <- DF[ DF[[col.group1]] > 0, ] 
     #  %>%  mutate(across(where(is.numeric), ~ round(., 2)))
      
      p1 <- ggdensity(df, x = col.group1,
                      add = "mean", rug = TRUE, 
                      add.params = list(size = 2, linetype = "longdash"),
                      color = "Pre_Post", fill = "Pre_Post",
                      facet.by = c("Group", col.group), short.panel.labs = T, 
                      ncol = length(unique(metadata.all[[col.group]])),
                      ggtheme = theme_light(), legend = "top"
      ) + 
        theme(strip.text = element_text(size = 20, face = "bold")# Change font size
              , axis.text.x = element_text(size=15)
              , legend.text = element_text(size=20)
              , axis.title = element_text(size=20, face = "bold")
              )
      #print(p1)
      png(filename = paste0(dir1, file.name, "_",col.group,"_", col.group1, "_SGLT2i.postive.png"), 
          width = length(unique(metadata.all[[col.group]]))*2, 
          height = length(unique(metadata.all[["Group"]]))*4, res = 300, units = "in")
      print(p1)
      dev.off()
      
      p2 <- ggviolin(df, x = "Pre_Post", y = col.group1,
               fill = "Pre_Post", add = "boxplot", add.params = list(fill = "white"),
               # id = "patient_id", line.color = "gray", line.size = 0.4,
               palette = "Pre_Post",
               facet.by = c("Group", col.group), short.panel.labs = T,
               ncol = length(unique(metadata.all[[col.group]])),
               xlab = "", ylab = col.group1)+
        #stat_compare_means(method = "anova", label.y = 90 ) +
        stat_compare_means(aes_string(group="Pre_Post"), 
                           method = "t.test", 
                           # paired = T,
                           label = "p.format", size = 4, 
                           label.y = (max(df[[col.group1]]) *1.1), 
                           hide.ns = TRUE, # show.legend = F 
                           )+
        theme(
          #axis.text.x=element_blank(),
          axis.text.x = element_text(angle = 0, hjust = 0.5),
          axis.ticks.x=element_blank()
        )
      png(filename = paste0(dir1, file.name, "_",col.group,"_", col.group1, "_SGLT2i.postive.violin.png"), 
          width = length(unique(metadata.all[[col.group]]))*1.5, 
          height = length(unique(metadata.all[["Group"]]))*4, res = 300, units = "in")
      print(p2)
      dev.off()
      
      # options(digits=2)
      p3 <- ggbarplot(df, x = "Pre_Post", y = col.group1,
                      # label = TRUE, label.pos = "out", # lab.vjust = -1.6,
                      add = c("mean_se", "jitter"),
                      position = position_dodge(),
                      color = "Pre_Post", 
                      # id = "patient_id", line.color = "gray", line.size = 0.4,
                      palette = "Pre_Post",
                      facet.by = c("Group", col.group), short.panel.labs = T,
                      ncol = length(unique(metadata.all[[col.group]])),
                      xlab = "", ylab = col.group1)+
        # label(sprintf("%0.2f", round(a, digits = 2)))+
        #stat_compare_means(method = "anova", label.y = 90 ) +
        stat_compare_means(aes_string(group="Pre_Post"), 
                           method = "t.test", 
                           # paired = T,
                           label = "p.format", size = 4, 
                           hide.ns = TRUE, # show.legend = F, 
                           label.y = (max(df[[col.group1]])*1.1) 
                           )+
        theme(
          #axis.text.x=element_blank(),
          axis.text.x = element_text(angle = 0, hjust = 0.5),
          axis.ticks.x=element_blank()
        )
      #print(p1)
      png(filename = paste0(dir1, file.name, "_",col.group,"_", col.group1, "_SGLT2i.postive.bar.png"), 
          width = length(unique(metadata.all[[col.group]]))*1.5, 
          height = length(unique(metadata.all[["Group"]]))*4, res = 300, units = "in")
      print(p3)
      dev.off()
      
      library(ggpubr)
      tt <- split(df , df[, c("Group", col.group)])
      tt <- tt[sapply(tt, nrow)>0]
             num.position <- c()
        for (tx in 1:length(tt)){
          num.position <- c(num.position, mean(tt[[tx]][[col.group1]]))
        }
        max(num.position)
        
        p3 <- ggbarplot(df, x = "Pre_Post", y = col.group1,
                      # label = TRUE, label.pos = "out", # lab.vjust = -1.6,
                      # add = c("mean_se", "jitter"),
                      add = "mean_se", 
                      # error.plot = "pointrange", 
                      # label = TRUE, lab.vjust = -1.6,
                      # position = position_dodge(),
                      color = "Pre_Post", 
                      # id = "patient_id", line.color = "gray", line.size = 0.4,
                      palette = "Pre_Post",
                      facet.by = c("Group", col.group), short.panel.labs = T,
                      ncol = length(unique(metadata.all[[col.group]])),
                      xlab = "", ylab = col.group1)+ # ylim(c(0, mean(df[[col.group1]])*1.5)) +
        # label(sprintf("%0.2f", round(a, digits = 2)))+
        #stat_compare_means(method = "anova", label.y = 90 ) +
        stat_compare_means(aes_string(group="Pre_Post"), 
                           method = "t.test", 
                           # paired = T,
                           label = "p.format", size = 4, 
                           hide.ns = TRUE, # show.legend = F, 
                           label.y = max(num.position)*1.1 
        )+
        theme(
          #axis.text.x=element_blank(),
          axis.text.x = element_text(angle = 0, hjust = 0.5),
          axis.ticks.x=element_blank()
        )
      print(p3)
      png(filename = paste0(dir1, file.name, "_",col.group,"_", col.group1, "_SGLT2i.postive.bar.nodot.png"), 
          width = length(unique(metadata.all[[col.group]]))*1.5, 
          height = length(unique(metadata.all[["Group"]]))*3, res = 300, units = "in")
      print(p3)
      dev.off()
      
      # gghistogram(metadata.all, x = col.group1,
      #             add = "mean", rug = TRUE,
      #             color = "Pre_Post", fill = "Pre_Post",
      #             facet.by = c(col.group), short.panel.labs = T, 
      #             ncol = length(unique(metadata.all[[col.group]])),
      #             ggtheme = theme_light(), legend = "top")
      print(paste(col.group, col.group1, sep = ":"))
    }
  }
}
rm(dir1, p1, p2, p3, seurat.file , df,
   file.name, Annotation.file.name, Subtype.file.name,
   sample.df,  sample.meta, Subtype.meta, Subtype.meta1, 
   metadata.all, row.groups.new, row.group.new, 
   col.group, col.group1
)








