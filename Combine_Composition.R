
# rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.
memory.limit(size = 1e+13)

# import parameters and packages --------------------------------------------------------------------------------------------
# source(paste0("/media/jianie/Extreme SSD1/2022-09-01 STARR_SGLT2 Combine/Project Parameters.R"), local = knitr::knit_global())
source(paste0(paste0(Disk, Project.folder, "/", codes.folder, "/", "Project Parameters.R")), local = knitr::knit_global())
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


combine.ls <- list()
for (file in files) {
  # file = files[1]
  # file = files[5]
  load(paste0(dir0, file))
  
  # prepare file name
  file.name <- stringr::str_split(file, ".metadata_Compositions.RData")[[1]][1]
  file.name
  Annotation.file.name <- paste(stringr::str_split(file.name, "_")[[1]][c(1,2,3)], collapse = "_")
  Annotation.file.name
  Subtype.file.name <- stringr::str_split(file.name, paste0(Annotation.file.name, "_"))[[1]][2]
  Subtype.file.name
  if( length(stringr::str_split(Subtype.file.name, paste0("_"))[[1]]) > 1){
    celltype.name <- stringr::str_split(Subtype.file.name, paste0("_"))[[1]][1]
    #label.name <- stringr::str_split(Subtype.file.name, paste0("_"))[[1]][2]
  } else {
    celltype.name <- Annotation.file.name
    # label.name <- stringr::str_split(Subtype.file.name, paste0("_"))[[1]][1]
  }
  celltype.name
  
  # import data
  for (list.name in names(Composition.lists)) {
    # list.name = names(Composition.lists)[1]
    list.name
    label.name <- stringr::str_split(list.name, paste0("_"))[[1]][length(stringr::str_split(list.name, paste0("_"))[[1]])] # the last one in list.name
    label.name
    x.name <- stringr::str_split(list.name, paste0("_"))[[1]][2] # the second one in list name
    x.name
    Composition <- Composition.lists[[list.name]][["Composition_All"]] 
    colnames(Composition) <- paste0(celltype.name, "_",colnames(Composition))
    Composition <- Composition %>% tibble::rownames_to_column(var = "Dataset")
    clusters = colnames(Composition)
    # display.brewer.pal(length(clusters), "Set1")
    # Composition_3$x.label <- factor(Composition_3[[label.name]], levels=clusters)
    
    if (x.name == label.name){
      combine.ls[[paste0(celltype.name, "_",label.name)]] <- Composition
    }
  }
}

Composition.combine <- Reduce(full_join, combine.ls)

saveRDS(Composition.combine, file = paste0(Disk, Project.folder, "/", Rds.folder, "/", file.name, "_", "Composition.combine.Rds") )

rm(dir1, dir0, files, file,
   file.name, Annotation.file.name, Subtype.file.name,
   celltype.name, list.name, label.name, x.name,
   Composition, combine.ls, Composition.combine,
   clusters.lists, Composition.lists, clusters
)








