# DISCLAIMER
# Some of the algorithms are non-deterministic making the results slightly different from run to run.
# Many of the algorithms are frequently updated with new releases making the results slightly different from version to version.
# Depending on your system, this code below may not produce an exact 1:1 copy of the results. 
# This adversely affects clustering, filtering etc. if the code is naively copy-pasted.
# For better reproducibility, please download quality-filtered Rds objects and work on downstream analyses. 
# change "D:/" to "/media/jianie/Extreme SSD1/" # if in LUNIX computer

# rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.
memory.limit(size = 1e+13)
source(paste0(paste0(Disk, Project.folder, "/", "Project Parameters.R")), local = knitr::knit_global())
# Load libraries 
# library("velocyto.R") # linux -specific

# Import file name -----------------------------------------------------------------------------------------------
file.name <- stringr::str_split(seurat.file, ".Rds")[[1]][1]


# combine Subtype metadata --------------------------------------------
metadata <- bind_rows(readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/",
                                     file.name, "_", "Immune", "_", "Subtype.metadata.Rds")),
                      readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/",
                                     file.name, "_", "FAP", "_", "Subtype.metadata.Rds")))
metadata <- bind_rows(metadata,
                      readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/",
                                     file.name, "_", "Adipocytes", "_", "Subtype.metadata.Rds")))
metadata <- bind_rows(metadata, readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/",
                                     file.name, "_", "Endothelial", "_", "Subtype.metadata.Rds"))
                      )

table(metadata$Subtype, metadata$Annotation)
table(metadata$Subtype, metadata$Annotation_Label)

## save meta.data ----------------------------------------------
saveRDS(metadata, paste0(Disk, Project.folder, "/", Rds.folder, "/",file.name, "_", "Subtype.metadata.Rds"))

# add metadata ------------------------------------------
seurat.data <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/",file.name, ".Rds"))
seurat.data <- AddMetaData( object = seurat.data, metadata = metadata)
gc() #free up memory and report the memory usage.

DimPlot(seurat.data, group.by = "Subtype", label = F, raster=FALSE)
DimPlot(seurat.data, group.by = "Annotation_Label", label = T, raster=FALSE)
DimPlot(seurat.data, group.by = "Annotation", label = T, raster=FALSE)

# DimPlot(seurat.data, group.by = "Annotation", label = F, raster=FALSE)

rm(file.name, 
   metadata, seurat.data)

