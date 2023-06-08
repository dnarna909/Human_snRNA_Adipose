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
# source(paste0("/media/jianie/Extreme SSD1/2022-09-01 STARR Combine/Project Parameters.R"), local = knitr::knit_global())
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
  library("ElPiGraph.R")
  # tree.phate <- phate(tree.data$data)# test
  # summary(tree.phate)
} else {
  source(paste0(Disk, Project.folder, "/", "Project Parameters.R"), local = knitr::knit_global())
}
source(paste0(Disk, "00_Functions_Refs/BasicFunction.R"), local = knitr::knit_global())
# source(paste0(Disk, "00_Functions_Refs/Functions_Human Fat snRNA.R"), local = knitr::knit_global())
# or sys.source("your-script.R", envir = knitr::knit_global())
PlotPG <- LoadFunction(file=paste0(Disk, "00_Functions_Refs/Functions_Human Fat snRNA.R"), PlotPG)# PlotPG function has an error in orignal package, updated here

# Load libraries ------------------------------------------------------------------------
# reticulate::py_discover_config(required_module = "phate")
# reticulate::import("phate")
#library("velocyto.R") # linux -specific

dir.create(file.path(paste0(Disk, Project.folder, "/", figures.folder, "/",  ".Trajectory.phate", "/")), showWarnings = FALSE)
dir0 <- paste0(Disk, Project.folder, "/", figures.folder, "/", ".Trajectory.phate",  "/")
dir.create(file.path(paste0(Disk, Project.folder, "/", Rds.folder, "/",  ".Trajectory.phate", "/")), showWarnings = FALSE)
dir1 <- paste0(Disk, Project.folder, "/", Rds.folder, "/",  ".Trajectory.phate",  "/")

# Import data -----------------------------------------------------------------------------------------------------------------------------------------------------
Differentiation <- readRDS( paste0(Disk, Project.folder, "/", Rds.folder, "/","Differentiation_PHATE_", groups, ".Rds") )
meta.data<- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/","Differentiation_PHATE_", groups, ".Adipogesis.meta.Rds") )

# DimPlots
DimPlot(Differentiation, reduction = "phate", group.by = "Subtype", label = TRUE, split.by="Subjects")
DimPlot(Differentiation, reduction = "phate", group.by = "Subtype", label = TRUE)
DimPlot(Differentiation, reduction = "phate", group.by = "Adipogenesis.Type")
DimPlot(Differentiation, reduction = "phate", group.by = "Adipogenesis.Type", split.by="Age_Group")
# ggsave(filename = paste0("D:/2021-01-11 Figures for Grant resubmission/Rds files/Figure 6_Differentiation_human_Type_Age_Group.pdf"), height = 3.5, width = 7)

# plotting -----------
p <-DimPlot(Differentiation, reduction = "phate", group.by = "Adipogenesis.Type", raster=FALSE)
p
png(file=paste0("Figure 6_Differentiation_Adipogenesis.Type_DimPlot.png"), 
    width= 7, height=6, res = 300, units = "in")
print(p)
dev.off()

tt <- Differentiation@meta.data %>% dplyr::filter(Subtype == "CA3") %>% tibble::rownames_to_column() %>% pull(rowname)
p <-DimPlot(Differentiation, reduction = "phate", cells.highlight = tt, raster=FALSE) + NoLegend()
p
png(file=paste0("Figure 6_Differentiation_CA3_DimPlot.png"), 
    width= 7, height=6, res = 300, units = "in")
print(p)
dev.off()

tt <- Differentiation@meta.data %>% dplyr::filter(Subtype == "CPA1") %>% tibble::rownames_to_column() %>% pull(rowname)
p <-DimPlot(Differentiation, reduction = "phate", cells.highlight = tt, raster=FALSE) + NoLegend()
p
png(file=paste0("Figure 6_Differentiation_CPA1_DimPlot.png"), 
    width= 7, height=6, res = 300, units = "in")
print(p)
dev.off()

tt <- Differentiation@meta.data %>% dplyr::filter(Subtype == "FIP1") %>% tibble::rownames_to_column() %>% pull(rowname)
p <-DimPlot(Differentiation, reduction = "phate", cells.highlight = tt, raster=FALSE) + NoLegend()
p
png(file=paste0("Figure 6_Differentiation_FIP1_DimPlot.png"), 
    width= 7, height=6, res = 300, units = "in")
print(p)
dev.off()


tt <- Differentiation@meta.data %>% dplyr::filter(Compare_Group1) %>% tibble::rownames_to_column() %>% pull(rowname)
p <-DimPlot(Differentiation[, tt], reduction = "phate", group.by = "Adipogenesis.Type", split.by="BA_Group", raster=FALSE)
p
png(file=paste0("Figure 6_Differentiation_Type_Compare_Group1_BA_Group_DimPlot.png"), 
    width= 18, height=6, res = 300, units = "in")
print(p)
dev.off()
dt <- as.data.frame.matrix(table(Differentiation[, tt]@meta.data$BA_Group, Differentiation[, tt]@meta.data$Adipogenesis.Type)) %>% as.matrix()
dt <- as.data.frame.matrix(table(Differentiation[, tt]@meta.data$BA_Group, Differentiation[, tt]@meta.data$EdgeID)) %>% filter_all(any_vars(. != 0))
prop.table(dt%>% filter_all(any_vars(. != 0)) %>% as.matrix(), margin = 1)

plist  <-FeaturePlot(Differentiation[, tt], reduction = "phate", features = "Senescence.score", raster=FALSE, min.cutoff = 0, split.by="BA_Group") 
plist[[1]] + xlim(c(-0.01,0.025)) + ylim(c(-0.04,0.04))
png(file=paste0("Figure 6_Differentiation_Senescence.score_BA_Group_1_FeaturePlot.png"), 
    width= 7, height=6, res = 300, units = "in")
print(plist[[1]] + xlim(c(-0.01,0.025)) + ylim(c(-0.04,0.04)))
dev.off()
png(file=paste0("Figure 6_Differentiation_Senescence.score_BA_Group_2_FeaturePlot.png"), 
    width= 7, height=6, res = 300, units = "in")
print(plist[[2]] + xlim(c(-0.01,0.025)) + ylim(c(-0.04,0.04)))
dev.off()
png(file=paste0("Figure 6_Differentiation_Senescence.score_BA_Group_3_FeaturePlot.png"), 
    width= 7, height=6, res = 300, units = "in")
print(plist[[3]] + xlim(c(-0.01,0.025)) + ylim(c(-0.04,0.04)))
dev.off()
png(file=paste0("Figure 6_Differentiation_Senescence.score_BA_Group_4_FeaturePlot.png"), 
    width= 7, height=6, res = 300, units = "in")
print(plist[[4]] + xlim(c(-0.01,0.025)) + ylim(c(-0.04,0.04)))
dev.off()



table(Differentiation@meta.data$Treatment_Group)
tt <- Differentiation@meta.data %>% dplyr::filter(Compare_Group2) %>% tibble::rownames_to_column() %>% pull(rowname)
p <-DimPlot(Differentiation[, tt], reduction = "phate", group.by = "Adipogenesis.Type", split.by="Treatment_Group", raster=FALSE)
p
png(file=paste0("Figure 6_Differentiation_Type_Compare_Group2_Treatment_Group_DimPlot.png"), 
    width= 18, height=6, res = 300, units = "in")
print(p)
dev.off()
dt <- as.data.frame.matrix(table(Differentiation[, tt]@meta.data$Treatment_Group, Differentiation[, tt]@meta.data$Adipogenesis.Type)) %>% as.matrix()
prop.table(dt%>% as.matrix(), margin = 1)


p <-FeaturePlot(Differentiation, reduction = "phate", features = "Senescence.score", raster=FALSE, min.cutoff = 0) + xlim(c(-0.01,0.025)) + ylim(c(-0.04,0.04))
p
png(file=paste0("Figure 6_Differentiation_Senescence.score_FeaturePlot.png"), 
    width= 7, height=6, res = 300, units = "in")
print(p)
dev.off()

plist  <-FeaturePlot(Differentiation[, tt], reduction = "phate", features = "Senescence.score", raster=FALSE, min.cutoff = 0, split.by="Treatment_Group") 
plist[[1]] + xlim(c(-0.01,0.025)) + ylim(c(-0.04,0.04))
png(file=paste0("Figure 6_Differentiation_Senescence.score_Treatment_Group_1_FeaturePlot.png"), 
    width= 7, height=6, res = 300, units = "in")
print(plist[[1]] + xlim(c(-0.01,0.025)) + ylim(c(-0.04,0.04)))
dev.off()
png(file=paste0("Figure 6_Differentiation_Senescence.score_Treatment_Group_2_FeaturePlot.png"), 
    width= 7, height=6, res = 300, units = "in")
print(plist[[2]] + xlim(c(-0.01,0.025)) + ylim(c(-0.04,0.04)))
dev.off()
png(file=paste0("Figure 6_Differentiation_Senescence.score_Treatment_Group_3_FeaturePlot.png"), 
    width= 7, height=6, res = 300, units = "in")
print(plist[[3]] + xlim(c(-0.01,0.025)) + ylim(c(-0.04,0.04)))
dev.off()
png(file=paste0("Figure 6_Differentiation_Senescence.score_Treatment_Group_4_FeaturePlot.png"), 
    width= 7, height=6, res = 300, units = "in")
print(plist[[4]] + xlim(c(-0.01,0.025)) + ylim(c(-0.04,0.04)))
dev.off()


DimPlot(Differentiation, reduction = "phate", group.by = "Subtype", label = FALSE)

p <- DimPlot(Differentiation, reduction = "phate", group.by = "Annotation", label = FALSE)
p
png(file=paste0("Figure 6_Differentiation_Annotation_DimPlot.png"), 
    width= 7, height=6, res = 300, units = "in")
print(p)
dev.off()

p <- DimPlot(Differentiation, reduction = "phate", group.by = "Subtype", label = FALSE)
p
png(file=paste0("Figure 6_Differentiation_Subtype_DimPlot.png"), 
    width= 7, height=6, res = 300, units = "in")
print(p)
dev.off()

FeaturePlot(Differentiation, features = paste0(path, "_Pt"), raster=FALSE)
p <- VlnPlot(Differentiation, features = paste0(path, "_Pt"), group.by = "Subtype", raster=FALSE)

RidgePlot(Differentiation, features = paste0(path, "_Pt"), group.by = "Subtype")

p <- RidgePlot(Differentiation, features = paste0(path, "_Pt"), group.by = "Annotation")
p
png(file=paste0("Figure 6_Differentiation_Annotation_RidgePlot.png"), 
    width= 7, height=6, res = 300, units = "in")
print(p)
dev.off()

tt <- Differentiation@meta.data %>% dplyr::filter(Annotation == "Adipocytes") %>% tibble::rownames_to_column() %>% pull(rowname)
p <- RidgePlot(Differentiation[, tt], features = paste0(path, "_Pt"), group.by = "Subtype")
p
png(file=paste0("Figure 6_Differentiation_Adipocytes_Subtype_Pt_RidgePlot.png"), 
    width= 7, height=6, res = 300, units = "in")
print(p)
dev.off()
p <- VlnPlot(Differentiation[, tt], features = paste0(path, "_Pt"), group.by = "Subtype", raster=FALSE, sort = "increasing", pt.size = 0.01)
p
png(file=paste0("Figure 6_Differentiation_Adipocytes_Subtype_Pt_VlnPlot.png"), 
    width= 7, height=6, res = 300, units = "in")
print(p)
dev.off()


tt <- Differentiation@meta.data %>% dplyr::filter(Annotation == "FAP") %>% tibble::rownames_to_column() %>% pull(rowname)
p <- RidgePlot(Differentiation[, tt], features = paste0(path, "_Pt"), group.by = "Subtype")
p
png(file=paste0("Figure 6_Differentiation_FAP_Subtype_Pt_RidgePlot.png"), 
    width= 7, height=6, res = 300, units = "in")
print(p)
dev.off()
p <- VlnPlot(Differentiation[, tt], features = paste0(path, "_Pt"), group.by = "Subtype", raster=FALSE, sort = "decreasing", pt.size = 0.01)
p
png(file=paste0("Figure 6_Differentiation_FAP_Subtype_Pt_VlnPlot.png"), 
    width= 7, height=6, res = 300, units = "in")
print(p)
dev.off()

