
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
    label.name <- stringr::str_split(list.name, paste0("_"))[[1]][2]
    Composition_3 <- Composition.lists[[list.name]][["Composition_3All"]] %>% dplyr::filter(!patient_id %in% c("76615"))
    Composition <- Composition.lists[[list.name]][["Composition_All"]] 
    
    clusters = colnames(Composition)
    # display.brewer.pal(length(clusters), "Set1")
    # Composition_3$x.label <- factor(Composition_3[[label.name]], levels=clusters)
    
    # plotting
    p1 <- ggpaired(Composition_3, x = "Event.Name", y = "composition",
                   color = "Event.Name", id = "patient_id", line.color = "gray", line.size = 0.4,
                   palette = "npg",
                   facet.by = c("Group", label.name), short.panel.labs = T,
                   xlab = "", ylab = "Percentage of Cells")+
      stat_compare_means(aes_string(group="Event.Name"), 
                         method = "t.test", paired = T,
                         label = "p.signif",  hide.ns = TRUE, show.legend = F, size = 8,
                         label.y = 0 )+
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank())
    print(p1)
    png(filename = paste0(dir1, file.name, "_",label.name, ".ggpaired.png"), 
        width = length(clusters), height = 8, res = 300, units = "in")
    print(p1)
    dev.off()
    
    select.variables <- c("Changes_Weight..pounds._Group", "Changes_BMI_Group" , 
                          "Changes_Glucose_Group",  "Changes_A1c.Value_Group")
    for (ss in select.variables) {
      p2 <- ggpaired(Composition_3, x = "Event.Name", y = "composition",
                     color = "Event.Name", id = "patient_id", line.color = "gray", line.size = 0.4,
                     palette = "npg",
                     facet.by = c(ss, label.name), short.panel.labs = T,
                     xlab = "", ylab = "Percentage of Cells")+
        stat_compare_means(aes_string(group="Event.Name"), 
                           method = "t.test", paired = T,
                           label = "p.signif",  hide.ns = TRUE, show.legend = F, size = 8,
                           label.y = 0 )+
        theme(axis.text.x=element_blank(),
              axis.ticks.x=element_blank())
      print(p2)
      
      png(filename = paste0(dir1, file.name, "_",label.name, "_", ss, ".ggpaired.png"), 
          width = length(clusters), height = 12, res = 300, units = "in")
      print(p2)
      dev.off()
    }
    
    
    ggboxplot(Composition_3, x = "Event.Name", y = "composition",
              color = "Event.Name", palette = "npg",
              add = "jitter",
              facet.by = c("Group", label.name), short.panel.labs = T)+
      stat_compare_means(aes_string(group="Event.Name"), 
                         method = "t.test", paired = F,
                         label = "p.signif",  hide.ns = TRUE, show.legend = F, size = 8,
                         label.y = 0 )+
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank())
    
  }
}

rm(dir1, dir0, p1, files, file,
   file.name, Annotation.file.name, Subtype.file.name,
   celltype.name, list.name, label.name, Composition_3, Composition,
   clusters.lists, Composition.lists, clusters
)








