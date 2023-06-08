
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
    # list.name = names(Composition.lists)[2]
    label.name <- stringr::str_split(list.name, paste0("_"))[[1]][length(stringr::str_split(list.name, paste0("_"))[[1]])]
    x.name <- stringr::str_split(list.name, paste0("_"))[[1]][2]
    Composition_3 <- Composition.lists[[list.name]][["Composition_3All"]] %>% 
      dplyr::filter(Pre_Post %in% c("Pre", "Post") ) %>% 
      dplyr::filter(!patient_id %in% c("76615"))
    Composition <- Composition.lists[[list.name]][["Composition_All"]] 
    
    clusters = colnames(Composition)
    # display.brewer.pal(length(clusters), "Set1")
    # Composition_3$x.label <- factor(Composition_3[[label.name]], levels=clusters)
    
    if (x.name == label.name){
      # plotting
      p1 <- ggpaired(Composition_3, x = "Pre_Post", y = "composition",
                     color = "Pre_Post", id = "patient_id", line.color = "gray", line.size = 0.4,
                     palette = "npg",
                     facet.by = c("Group", x.name), short.panel.labs = T,
                     ncol = length(unique(Composition_3[[x.name]])),
                     xlab = "", ylab = "Percentage of Cells")+
        stat_compare_means(aes_string(group="Pre_Post"), 
                           method = "t.test", paired = T,
                           # label = "p.signif", size = 8, 
                           label = "p.format", size = 4, 
                           show.legend = F, hide.ns = TRUE, 
                           # label.y = 0,
                           label.y = (max(Composition_3[["composition"]])*1.1)
                           )+
        theme(#axis.text.x=element_blank(),
          axis.text.x = element_text(angle = 0, hjust = 0.5),
              axis.ticks.x=element_blank())
      print(p1)
      png(filename = paste0(dir1, file.name, "_", x.name,"_", "SGLT2i", ".ggpaired.png"), 
          width = length(unique(Composition_3[[x.name]]))*length(unique(Composition_3[["Pre_Post"]])), 
          height = length(unique(Composition_3[["Group"]]))*6, res = 300, units = "in")
      print(p1)
      dev.off()

    }
    
    if (x.name != label.name){
      for (ll in unique(Composition_3[[label.name]]) ) {
        df <- Composition_3 %>% dplyr::filter((!!sym(label.name)) == ll)
        
        p1 <- ggpaired(df, x = "Pre_Post", y = "composition",
                       color = "Pre_Post", id = "patient_id", line.color = "gray", line.size = 0.4,
                       palette = "npg",
                       facet.by = c("Group", x.name), short.panel.labs = T,
                       ncol = length(unique(df[[x.name]]))*length(unique(df[["Group"]])),
                       xlab = "", ylab = "Percentage of Cells")+
          stat_compare_means(aes_string(group="Pre_Post"), 
                             method = "t.test", paired = T,
                             # label = "p.signif", size = 8, 
                             label = "p.format", size = 4, 
                             show.legend = F, hide.ns = TRUE, 
                             # label.y = 0,
                             label.y = (max(df[["composition"]])*1.1)
                             )+
          theme(#axis.text.x=element_blank(),
            axis.text.x = element_text(angle = 0, hjust = 0.5),
            axis.ticks.x=element_blank())
        
        # ggsave(filename = paste0(Disk, Project.folder, "/", figures.folder, "/", file.name, "_",label.name, ".png"), plot = p1,
        #        width = 12, height = 8)
        print(p1)
        
        png(filename = paste0(dir1, file.name, "_", x.name, "_", ll, "_", "SGLT2i", ".ggpaired.png"), 
            width = length(unique(df[[x.name]]))*length(unique(df[["Pre_Post"]])), 
            height = length(unique(df[["Group"]]))*6, res = 300, units = "in")
        print(p1)
        dev.off()
      }
    }
  }
}

rm(dir1, dir0, p1, files, file,df,ll,
   file.name, Annotation.file.name, Subtype.file.name,
   celltype.name, list.name, label.name, x.name,
   Composition_3, Composition,
   clusters.lists, Composition.lists, clusters
)








