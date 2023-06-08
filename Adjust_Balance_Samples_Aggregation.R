###' DISCLAIMER
#' Some of the algorithms are non-deterministic making the results slightly different from run to run.
#' Many of the algorithms are frequently updated with new releases making the results slightly different from version to version.
#' Depending on your system, this code below may not produce an exact 1:1 copy of the results. 
#' This adversely affects clustering, filtering etc. if the code is naively copy-pasted.
#' For better reproducibility, please download quality-filtered Rds objects and work on downstream analyses. 
#' change "D:/" to "/media/jianie/Extreme SSD1/" # if in LUNIX computer

##' Load libraries-------------------------------------------------------------------------------------------------------------
# rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memory and report the memory usage.
memory.limit(size = 1e+13) # 'memory.limit()' is Windows-specific
if (serve == "GCCRI"){
  source(paste0(Disk, Project.folder, "/", "Project Parameters_GCCRI.R"), local = knitr::knit_global())
} else {
  source(paste0(Disk, Project.folder, "/", "Project Parameters.R"), local = knitr::knit_global())
}

# import packages
# library("velocyto.R") # linux -specific
library(Seurat)

# Import file name -----------------------------------------------------------------------------------------------
metadata.new <- readRDS(paste0(Disk, Project.folder, "/",  Rds.folder, "/", "SAT_Sample.meta.data.Rds")) ## This object is downloadable from Open Science Framework, with further annotations as subsequent scripts add to the file

# Import data and subset ---------------------------------------------------------------------------------------------------------------------------------
for (seurat.file in seurat.files) {
  file.name <- stringr::str_split(seurat.file, ".Rds")[[1]][1]
  seurat.data <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/", seurat.file))
  gc() #free up memory and report the memory usage.
  seurat.data <- AddMetaData( object = seurat.data, metadata = metadata.new)
  
  # put into seurat.data
  print(seurat.file)
  g <-DimPlot(seurat.data, group.by = "Label", split.by = "Treatment_Group", raster = F)
  print(g)
  # png(file=paste0(Disk, Project.folder, "/", figures.folder, "/", file.name, ".umap.png"), 
  #     width= length(unique(seurat.data@meta.data[["Treatment_Group"]])) * 2, height=3, res = 300, units = "in")
  # print(g)
  # dev.off()
  
  seurat.data <- SetIdent(seurat.data, value = "Dataset")
  unique(seurat.data@meta.data[[ "Dataset"]])
  seurat.data <- subset(seurat.data, idents = c("STARR-SGLT2-1_Sample01", "STARR-SGLT2-2_Sample01"), invert = TRUE)
  file.name2 <- paste(stringr::str_split(file.name, "_")[[1]][c(1,2,3)], collapse = "_")
  saveRDS(seurat.data, paste0(Disk, Project.folder, "/", Rds.folder, "/", file.name2, ".Rds") )
  print(seurat.file)
  g <-DimPlot(seurat.data, group.by = "Label", split.by = "Treatment_Group", raster = F)
  print(g)
  
  # DimPlot(seurat.data, group.by = "Label", split.by = "sample_id", raster = F)
  rm(seurat.data)
  gc() #free up memory and report the memory usage.
}
rm(metadata.new, g, file.name2)

