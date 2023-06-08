# Load libraries -------------------------------------------------------------------------------------------------------------
rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.
memory.limit(size = 1e+13) # 'memory.limit()' is Windows-specific

# set parameters
# source(paste0("/media/jianie/Extreme SSD1/2022-09-01 STARR Combine/Project Parameters.R"), local = knitr::knit_global())
source(paste0("D:/2022-09-01 STARR Combine/Project Parameters.R"), local = knitr::knit_global())

# import packages
setwd(paste0(Disk, Project.folder, "/"))
packageList <- scan( paste0(Disk, "00_Functions_Refs/SingleCellPreprocessPackageList.txt") , character(), quote = "")
for(package in packageList){
  if(!require(package,character.only = TRUE)){
    install.packages(package);require(package,character.only = TRUE);}
}


eWAT <- readRDS("D:/2022-09-01 STARR Combine/Rds files_Aggr_all/SAT.Rds")
##' 4. Final embedding and clustering using ambient gene removal and Harmony (best performance in above test). NOTE: THIS STEP IS NON-DETERMINISTIC (Harmony) - RESULTS VARY FROM RUN TO RUN
eWAT <- FindVariableFeatures(eWAT, nfeature=5000)
# VariableFeatures(eWAT) <- VariableFeatures(eWAT)[!(VariableFeatures(eWAT) %in% Amb)]
eWAT <- ScaleData(eWAT)
eWAT <- RunPCA(eWAT)
eWAT <- RunHarmony(eWAT, group.by.vars="Dataset", dims.use=1:20, assay.use = "originalexp")
eWAT <- RunUMAP(eWAT, dims=1:15, reduction="harmony")
FeaturePlot(eWAT, features = "DoubletScore", raster=FALSE)
FeaturePlot(eWAT, features = "sum", raster=FALSE)
FeaturePlot(eWAT, features = "sizeFactor", raster=FALSE)
table(eWAT$outlier)
saveRDS(eWAT, paste0(Disk, Project.folder, "/",  Rds.folder, "/SAT_5000.Rds")) ## This object is downloadable from Open Science Framework, with further annotations as subsequent scripts add to the file
saveRDS(eWAT@meta.data, paste0(Disk, Project.folder, "/",  Rds.folder, "/SAT.meta.data.Rds"))


eWAT <- FindNeighbors(object = eWAT, dims = 1:2, reduction = "umap")
eWAT <- FindClusters(eWAT, resolution = 0.05, algorithm = 1) # 17
DimPlot(eWAT,  raster=FALSE, label = T) # final
eWAT <- FindClusters(eWAT, resolution = 0.02, algorithm = 1) # 11
DimPlot(eWAT,  raster=FALSE, label = T) # final
eWAT <- FindClusters(eWAT, resolution = 0.015, algorithm = 1) # 9
DimPlot(eWAT,  raster=FALSE, label = T) # final
eWAT <- FindClusters(eWAT, resolution = 0.012, algorithm = 1) # 8
DimPlot(eWAT,  raster=FALSE, label = T) # final
eWAT <- FindClusters(eWAT, resolution = 0.01, algorithm = 1) # 8
DimPlot(eWAT,  raster=FALSE, label = T) # final
eWAT <- FindClusters(eWAT, resolution = 0.007, algorithm = 1) # 8
eWAT <- FindClusters(eWAT, resolution = 0.006, algorithm = 1) # 7
eWAT <- FindClusters(eWAT, resolution = 0.005, algorithm = 1) # 7

DimPlot(eWAT,  group.by = "originalexp_snn_res.0.007", raster=FALSE, label = T) # final


DimPlot(eWAT, group.by = "Dataset") # final
# DimPlot(eWAT)
# DimPlot(eWAT, group.by = "Gender")
# DimPlot(eWAT, group.by = "Age")
# DimPlot(eWAT, group.by = "BMI")
# DimPlot(eWAT, split.by = "Dataset", raster=FALSE) # final
VlnPlot(eWAT, features = "DoubletScore", raster=FALSE)
# VlnPlot(eWAT, features = "sizeFactor", raster=FALSE)

saveRDS(eWAT, paste0(Disk, Project.folder, "/",  Rds.folder, "/SAT_5000.Rds")) ## This object is downloadable from Open Science Framework, with further annotations as subsequent scripts add to the file
