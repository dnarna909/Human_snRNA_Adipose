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

# in GCCRI terminal:
# source /opt/conda/etc/profile.d/conda.sh
# conda activate biocore
# R
# library(reticulate)

# import parameters and packages --------------------------------------------------------------------------------------------
if (serve == "GCCRI"){
  source(paste0(Disk, Project.folder, "/", "Project Parameters_GCCRI.R"), local = knitr::knit_global())
  library(Seurat)
  # First install phate in Python by running the following code from a terminal:
  # pip install --user phate
  # Then install phateR from CRAN by running the following code in R:
  # install.packages("phateR")
  library(dplyr)
  library(harmony)
  reticulate::use_python("/opt/conda/envs/biocore/bin/python")
  # scipy <- reticulate::import("scipy");phate <- reticulate::import("phate") # not need
  reticulate::py_discover_config(required_module="phate")
  library(phateR)
  tree.phate <- phate(tree.data$data)# test
  summary(tree.phate)
} else {
  source(paste0(Disk, Project.folder, "/", "Project Parameters.R"), local = knitr::knit_global())
}
source(paste0(Disk, "00_Functions_Refs/BasicFunction.R"), local = knitr::knit_global())

# source(paste0(Disk, "00_Functions_Refs/Functions_Human Fat snRNA.R"), local = knitr::knit_global())
# or sys.source("your-script.R", envir = knitr::knit_global())

# Load libraries ------------------------------------------------------------------------
# reticulate::py_discover_config(required_module = "phate")
# reticulate::import("phate")
#library("velocyto.R") # linux -specific


# Import data -----------------------------------------------------------------------------------------------------------------------------------------------------
#' ### Combine data
file.name <- stringr::str_split(seurat.file, ".Rds")[[1]][1]
Annotation.file.name <- paste(stringr::str_split(file.name, "_")[[1]][c(1,2,3)], collapse = "_")
Subtype.file.name <- stringr::str_split(file.name, paste0(Annotation.file.name, "_"))[[1]][2]
if (is.na(Subtype.file.name)) {Subtype.file.name = ""}
Subtype.meta <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/", Annotation.file.name, "_", file.type)) %>% tibble::rownames_to_column(var = "rowname")
metadata.all <- Subtype.meta %>% tibble::column_to_rownames(var = "rowname")

if(!file.exists(paste0(Disk, Project.folder, "/", Rds.folder, "/","Adipogenesis_PHATE.CPA1.Rds"))) {
  seurat.data <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/", seurat.file ))
  # n <- ncol(seurat.data)
  # seurat.data2 <- seurat.data[, 1:trunc(n / 10)]
  # Adipogenesis <- seurat.data2
  
  # Import seurat data
  if (exists("sub.folder")){
    seurat.data <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/", sub.folder, "/", seurat.file)) 
  }
  if (!exists("sub.folder")){
    seurat.data <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/", seurat.file)) 
  }
  seurat.data <- AddMetaData(seurat.data, metadata = metadata.all)
  if (!exists("Select.Group")){
    seurat.data <- SetIdent(seurat.data, value = Select.Group)
    seurat.data <- subset(seurat.data, idents = TRUE)
  }
  gc() #free up memory and report the memory usage.
  
  
  print(DimPlot(seurat.data, group.by = "Subtype",raster=FALSE, label = TRUE))
  # Amb <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/","Ambient.Rds"))
  seurat.data <- subset(seurat.data, subset = Annotation %in% c("Adipocytes", "FAP"))
  gc() #free up memrory and report the memory usage.
  gc() #free up memrory and report the memory usage.
  # seurat.data <- subset(seurat.data, subset = Annotation %in% c("Adipocytes", "FAP") |  Subtype %in% c("CEM1", "CEM2"))
  Idents(object = seurat.data) <- "Subtype"
  seurat.data <- subset(seurat.data, subset = Subtype %in% c("NewP2","NewP3","NewP4","NewP5","NewP6","CPA2","FIP1","FIP2","APC","NewP1"), invert = TRUE)# ,"CPA1"
  dim(seurat.data)
  print(DimPlot(seurat.data, group.by = "Annotation", raster=FALSE))
  print(DimPlot(seurat.data, group.by = "Subtype", raster=FALSE))
  
  
  # Trajectory analysis of FAP and adipocytes --------------------------------------------------------------------------------------------------------------------
  #' ## Trajectory
  #' ### PHATE embedding  of adipocytes and FAPs.
  #' NOTE: THIS STEP IS NON-DETERMINISTIC (Harmony) - RESULTS VARY FROM RUN TO RUN. 
  #' Thus, in order to reproduce downstream analyses using the same Harmony results, please download Adipogenesis.Rds
  seurat.data <- FindVariableFeatures(seurat.data, nfeatures= dim(seurat.data)[1]-1 ) # (dim(seurat.data)[1]-1 ) # 20000
  # VariableFeatures(seurat.data) <- VariableFeatures(seurat.data)[!(VariableFeatures(seurat.data) %in% Amb)]
  seurat.data <- ScaleData(seurat.data)
  gc() #free up memrory and report the memory usage.
  seurat.data <- RunPCA(seurat.data)
  gc() #free up memrory and report the memory usage.
  seurat.data <- RunHarmony(seurat.data, group.by.vars="Dataset", assay.use = "originalexp")
  gc() #free up memrory and report the memory usage.
  saveRDS(seurat.data, paste0(Disk, Project.folder, "/", Rds.folder, "/","Adipogenesis_PHATE.CPA1.Rds"))
  gc() #free up memrory and report the memory usage.
} else {
  seurat.data <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/","Adipogenesis_PHATE.CPA1.Rds"))
}
gc() #free up memrory and report the memory usage.


# phate analysis

seurat.data <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/","Adipogenesis_PHATE.CPA1.Rds"))
# PHATE calculation
PHATE <- phateR::phate(Embeddings(seurat.data, "harmony"), knn = 3, gamma = 0, t = 8) # updated # Seurat..ProjectDim.originalexp.harmony
seurat.data[["phate"]] <- CreateDimReducObject(PHATE$embedding, assay = "originalexp", key = "phate_") # updated
DimPlot(seurat.data, reduction = "phate", raster=FALSE)
DimPlot(seurat.data, reduction = "phate", group.by = "Subtype", raster=FALSE)
seurat.data@project.name <- "Adipogenesis"
DimPlot(seurat.data, reduction = "phate", group.by = "Subtype", label = TRUE, raster=FALSE)
DimPlot(seurat.data, reduction = "phate", group.by = "Annotation", label = TRUE, raster=FALSE)

# save files ----------------------------------------------------------------------------------------------------------------------
# saveRDS(PHATE, paste0(Disk, Project.folder, "/", Rds.folder, "/","Adipogenesis.PHATE.Rds"))
saveRDS(PHATE, paste0(Disk, Project.folder, "/", Rds.folder, "/","PHATE.CPA1.Rds"))
saveRDS(seurat.data, paste0(Disk, Project.folder, "/", Rds.folder, "/","Adipogenesis_PHATE.CPA1.Rds"))
gc() #free up memrory and report the memory usage.

