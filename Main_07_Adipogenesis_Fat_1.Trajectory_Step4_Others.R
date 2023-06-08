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
meta_data_merge <- readRDS(file=paste0(Disk, Project.folder, "/", Rds.folder, "/","Adipogenesis.meta.data.PHATE_", groups, "Pt..Rds"))
colnames(meta_data_merge)
load( file = paste0(Disk, Project.folder, "/", Rds.folder, "/", "Adipogenesis.PHATE_", groups,"Branch..RData"))
Adipogenesis <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/","Adipogenesis_PHATE_", groups, ".Rds"))
Adipogenesis <- AddMetaData(Adipogenesis, metadata = meta_data_merge)

# Meta-stable states -------------------------------------------------------------------------------------------------------------------
#' ## Subset to only the differentiation branch
Branch
path = Branch
Idents(object = Adipogenesis) <- path
Differentiation <- subset(Adipogenesis, idents = path)
hist(Differentiation@meta.data[[paste0(path, "_Pt")]])
rm(Adipogenesis)
gc() #free up memrory and report the memory usage.

# Define the intermediary state
PHATE <- Embeddings(Differentiation, "phate")
EdgeID <- Differentiation@meta.data[, c("Subtype", "EdgeID", paste0(path, "_Pt") )] %>% arrange((!!sym(paste0(path, "_Pt"))))
DimPlot(Differentiation, reduction = "phate", group.by = "Subtype", label = TRUE, raster=FALSE)
DimPlot(Differentiation, reduction = "phate", group.by = "EdgeID", label = TRUE, raster=FALSE)
EdgeID.df <- EdgeID[duplicated(EdgeID$EdgeID) ==F, c("EdgeID", paste0(path, "_Pt"))] %>% arrange((!!sym(paste0(path, "_Pt")))) %>%
  mutate(BranchID = 1:n())
EdgeID.df$EdgeID
Differentiation@meta.data$Adipogenesis.Type = "NA"
Differentiation@meta.data[ rownames(Differentiation@meta.data) %in% rownames(EdgeID[EdgeID$EdgeID %in% c(EdgeID.df$EdgeID[1: floor(length(EdgeID.df$EdgeID)/4)]),]), "Adipogenesis.Type"] <- "Early_Pre"
Differentiation@meta.data[ rownames(Differentiation@meta.data) %in% rownames(EdgeID[EdgeID$EdgeID %in% c(EdgeID.df$EdgeID[(floor(length(EdgeID.df$EdgeID)/4)+1): floor(length(EdgeID.df$EdgeID)*2/4)]),]), "Adipogenesis.Type"] <- "Late_Pre"
Differentiation@meta.data[ rownames(Differentiation@meta.data) %in% rownames(EdgeID[EdgeID$EdgeID %in% c(EdgeID.df$EdgeID[(floor(length(EdgeID.df$EdgeID)*2/4)+1): floor(length(EdgeID.df$EdgeID)*3/4)]),]), "Adipogenesis.Type"] <- "Transitioning"
Differentiation@meta.data[ rownames(Differentiation@meta.data) %in% rownames(EdgeID[EdgeID$EdgeID %in% c(EdgeID.df$EdgeID[(floor(length(EdgeID.df$EdgeID)*3/4)+1): length(EdgeID.df$EdgeID)]),]), "Adipogenesis.Type"] <- "Mature"
unique(Differentiation@meta.data[["Adipogenesis.Type"]])
Differentiation@meta.data[["Adipogenesis.Type"]] <- factor(Differentiation@meta.data[["Adipogenesis.Type"]], levels = c( "Early_Pre", "Late_Pre", "Transitioning", "Mature"))
DimPlot(Differentiation, reduction = "phate", group.by = "Adipogenesis.Type", raster=FALSE)
DimPlot(Differentiation, reduction = "phate", group.by = "Adipogenesis.Type", raster=FALSE, split.by = "Subtype")

EdgeID <- EdgeID %>% left_join(EdgeID.df[, c("EdgeID", "BranchID")], by = "EdgeID" )
Differentiation <- AddMetaData(Adipogenesis, metadata = EdgeID[, c("EdgeID", "BranchID")])

# Combine metadata 
file.type = "SAT_Label_2000_Subtype.metadata.Rds"; # with Annotation.file.name # subset seurat.data
sample.meta <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/", file.type)) 
colnames(sample.meta)
Differentiation <- AddMetaData(Differentiation, metadata = sample.meta)

sample.file.type = "SAT_Sample.meta.data.Rds"; # sample files
sample.meta <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/", sample.file.type)) 
colnames(sample.meta)
Differentiation <- AddMetaData(Differentiation, metadata = sample.meta)

file.type1 = "Senescence_SAT_Label_2000/Senescence.Marker_ModuleScore_Label_metadata.rds";
sen.meta <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/", file.type1)) 
colnames(sen.meta)
Differentiation <- AddMetaData(Differentiation, metadata = sen.meta)

table(Differentiation@meta.data[["Annotation"]] , Differentiation@meta.data[["Adipogenesis.Type"]] )
Differentiation@meta.data[["Annotation.Adipogenesis"]] <- paste0(Differentiation@meta.data[["Annotation"]] , ".", Differentiation@meta.data[["Adipogenesis.Type"]] )
table(Differentiation@meta.data[["Annotation.Adipogenesis"]])

table(Differentiation@meta.data[["Subtype"]] , Differentiation@meta.data[["Adipogenesis.Type"]] )
Differentiation@meta.data[["Subtype.Adipogenesis"]] <- paste0(Differentiation@meta.data[["Subtype"]] , ".", Differentiation@meta.data[["Adipogenesis.Type"]] )
table(Differentiation@meta.data[["Subtype.Adipogenesis"]])


table(Differentiation@meta.data[["Annotation"]] , Differentiation@meta.data[["BranchID"]] )
Differentiation@meta.data[["Annotation.Adipogenesis.BranchID"]] <- paste0(Differentiation@meta.data[["Annotation"]] , ".", Differentiation@meta.data[["BranchID"]] )
table(Differentiation@meta.data[["Annotation.Adipogenesis.BranchID"]])

table(Differentiation@meta.data[["Subtype"]] , Differentiation@meta.data[["BranchID"]] )
Differentiation@meta.data[["Subtype.Adipogenesis.BranchID"]] <- paste0(Differentiation@meta.data[["Subtype"]] , ".", Differentiation@meta.data[["BranchID"]] )
table(Differentiation@meta.data[["Subtype.Adipogenesis.BranchID"]])


FeaturePlot(Differentiation,  reduction = "phate", features = "DoubletScore", raster=FALSE, order = T, min.cutoff = 0)+ xlim(c(-0.01,0.025)) + ylim(c(-0.04,0.04))

# Save the object -----------------------------------------------------------------------------------------------------
saveRDS(Differentiation, paste0(Disk, Project.folder, "/", Rds.folder, "/","Differentiation_PHATE_", groups, ".Rds") )
saveRDS(Differentiation@meta.data[, c("EdgeID", path, paste0(path, "_Pt"), "Adipogenesis.Type",
                                      "Annotation.Adipogenesis", "Subtype.Adipogenesis")], 
        paste0(Disk, Project.folder, "/", Rds.folder, "/","Differentiation_PHATE_", groups, ".Adipogesis.meta.Rds") )
gc()

# DimPlots
DimPlot(Differentiation, reduction = "phate", group.by = "Subtype", label = TRUE, split.by="Dataset")
DimPlot(Differentiation, reduction = "phate", group.by = "Subtype", label = TRUE)
DimPlot(Differentiation, reduction = "phate", group.by = "Adipogenesis.Type")
DimPlot(Differentiation, reduction = "phate", group.by = "Adipogenesis.Type", split.by="Age_Group")
# ggsave(filename = paste0("D:/2021-01-11 Figures for Grant resubmission/Rds files/Figure 6_Differentiation_human_Type_Age_Group.pdf"), height = 3.5, width = 7)


# plot FAP and Adipocyte markers --------------------------
## Gene modules
# Overall_Markers <- readRDS(paste0("eWAT_Overall_Markers.Rds")) # Generated during analysis of the main datasets
Overall_Markers <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/", "FindMarkers_SAT_Label_2000/Annotation.Level2.Markers/SAT_Label_2000_Annotation.Select_Group_Combined.Mean.Rds"))

# Extract marker gene results
Clust1 <- Overall_Markers[["FAP"]][["All"]] # FAPs
Clust2 <- Overall_Markers[["Adipocytes"]][["All"]]  # Adipocytes

# Subset to genes that are enriched or exclusive (!= NS)
# Clust1 <- Clust1[Clust1$Marker != "NS", ]
# Clust2 <- Clust2[Clust2$Marker != "NS", ]

# Extracting the minimum logFC across replicates in pairwise cluster comparisons
# RATIONALE: The genes with the highest minimum fold change are the most specific ones for any given cluster
# Clust1$logFC_OvO <- apply(Clust1[, grep("logFC_Cluster", colnames(Clust1))], 1, FUN = "min")
# Clust2$logFC_OvO <- apply(Clust2[, grep("logFC_Cluster", colnames(Clust2))], 1, FUN = "min")

# Computing gene module scores using the top 50 most specific marker genes
Differentiation <- AddModuleScore(Differentiation, features = list(Cluster1 = Clust1[order(Clust1$Rank.ID), ][1:50, ]$rowname,
                                                                   Cluster2 = Clust2[order(Clust2$Rank.ID), ][1:50, ]$rowname))

# Scale the module scores (across all nuclei)
Differentiation$Cluster1 <- scale(Differentiation$Cluster1)
Differentiation$Cluster2 <- scale(Differentiation$Cluster2)

# Boxplots
par(mfcol=c(1,2))
boxplot(
  Differentiation@meta.data[ Differentiation@meta.data$Adipogenesis.Type == "Early_Pre", "Cluster1"][,1],
  Differentiation@meta.data[ Differentiation@meta.data$Adipogenesis.Type == "Late_Pre", "Cluster1"][,1],
  Differentiation@meta.data[ Differentiation@meta.data$Adipogenesis.Type == "Transitioning", "Cluster1"][,1],
  Differentiation@meta.data[ Differentiation@meta.data$Adipogenesis.Type == "Mature", "Cluster1"][,1],
  outline=F, las=1, main = "FAP signature", names = c("Early_Pre", "Late_Pre","Transitioning","Mature"), ylab="Scaled score")
boxplot(
  Differentiation@meta.data[ Differentiation@meta.data$Adipogenesis.Type == "Early_Pre", "Cluster2"][,1],
  Differentiation@meta.data[ Differentiation@meta.data$Adipogenesis.Type == "Late_Pre", "Cluster2"][,1],
  Differentiation@meta.data[ Differentiation@meta.data$Adipogenesis.Type == "Transitioning", "Cluster2"][,1],
  Differentiation@meta.data[ Differentiation@meta.data$Adipogenesis.Type == "Mature", "Cluster2"][,1],
  outline=F, las=1, main = "Adipocyte signature", names = c("Early_Pre", "Late_Pre","Transitioning","Mature"), ylab="Scaled score")


DimPlot(Differentiation, reduction = "phate", group.by = "EdgeID", label = TRUE, raster=FALSE)
EdgeID.df$EdgeID
DimPlot(Differentiation, reduction = "phate", group.by = "Annotation", label = TRUE, raster=FALSE)
par(mfcol=c(1,2))
boxplot(
  Differentiation@meta.data[ Differentiation@meta.data$EdgeID == EdgeID.df$EdgeID[1], "Cluster1"][,1],
  Differentiation@meta.data[ Differentiation@meta.data$EdgeID == EdgeID.df$EdgeID[2], "Cluster1"][,1],
  Differentiation@meta.data[ Differentiation@meta.data$EdgeID == EdgeID.df$EdgeID[3], "Cluster1"][,1],
  Differentiation@meta.data[ Differentiation@meta.data$EdgeID == EdgeID.df$EdgeID[4], "Cluster1"][,1],
  Differentiation@meta.data[ Differentiation@meta.data$EdgeID == EdgeID.df$EdgeID[5], "Cluster1"][,1],
  Differentiation@meta.data[ Differentiation@meta.data$EdgeID == EdgeID.df$EdgeID[6], "Cluster1"][,1],
  Differentiation@meta.data[ Differentiation@meta.data$EdgeID == EdgeID.df$EdgeID[7], "Cluster1"][,1],
  Differentiation@meta.data[ Differentiation@meta.data$EdgeID == EdgeID.df$EdgeID[8], "Cluster1"][,1],
  Differentiation@meta.data[ Differentiation@meta.data$EdgeID == EdgeID.df$EdgeID[9], "Cluster1"][,1],
  Differentiation@meta.data[ Differentiation@meta.data$EdgeID == EdgeID.df$EdgeID[10], "Cluster1"][,1],
  Differentiation@meta.data[ Differentiation@meta.data$EdgeID == EdgeID.df$EdgeID[11], "Cluster1"][,1],
  Differentiation@meta.data[ Differentiation@meta.data$EdgeID == EdgeID.df$EdgeID[12], "Cluster1"][,1],
  Differentiation@meta.data[ Differentiation@meta.data$EdgeID == EdgeID.df$EdgeID[13], "Cluster1"][,1],
  Differentiation@meta.data[ Differentiation@meta.data$EdgeID == EdgeID.df$EdgeID[14], "Cluster1"][,1],
  Differentiation@meta.data[ Differentiation@meta.data$EdgeID == EdgeID.df$EdgeID[15], "Cluster1"][,1],
  Differentiation@meta.data[ Differentiation@meta.data$EdgeID == EdgeID.df$EdgeID[16], "Cluster1"][,1],
  outline=F, las=1, main = "FAP signature", names = as.character(EdgeID.df$EdgeID), ylab="Scaled score")
boxplot(
  Differentiation@meta.data[ Differentiation@meta.data$EdgeID == EdgeID.df$EdgeID[1], "Cluster2"][,1],
  Differentiation@meta.data[ Differentiation@meta.data$EdgeID == EdgeID.df$EdgeID[2], "Cluster2"][,1],
  Differentiation@meta.data[ Differentiation@meta.data$EdgeID == EdgeID.df$EdgeID[3], "Cluster2"][,1],
  Differentiation@meta.data[ Differentiation@meta.data$EdgeID == EdgeID.df$EdgeID[4], "Cluster2"][,1],
  Differentiation@meta.data[ Differentiation@meta.data$EdgeID == EdgeID.df$EdgeID[5], "Cluster2"][,1],
  Differentiation@meta.data[ Differentiation@meta.data$EdgeID == EdgeID.df$EdgeID[6], "Cluster2"][,1],
  Differentiation@meta.data[ Differentiation@meta.data$EdgeID == EdgeID.df$EdgeID[7], "Cluster2"][,1],
  Differentiation@meta.data[ Differentiation@meta.data$EdgeID == EdgeID.df$EdgeID[8], "Cluster2"][,1],
  Differentiation@meta.data[ Differentiation@meta.data$EdgeID == EdgeID.df$EdgeID[9], "Cluster2"][,1],
  Differentiation@meta.data[ Differentiation@meta.data$EdgeID == EdgeID.df$EdgeID[10], "Cluster2"][,1],
  Differentiation@meta.data[ Differentiation@meta.data$EdgeID == EdgeID.df$EdgeID[11], "Cluster2"][,1],
  Differentiation@meta.data[ Differentiation@meta.data$EdgeID == EdgeID.df$EdgeID[12], "Cluster2"][,1],
  Differentiation@meta.data[ Differentiation@meta.data$EdgeID == EdgeID.df$EdgeID[13], "Cluster2"][,1],
  Differentiation@meta.data[ Differentiation@meta.data$EdgeID == EdgeID.df$EdgeID[14], "Cluster2"][,1],
  Differentiation@meta.data[ Differentiation@meta.data$EdgeID == EdgeID.df$EdgeID[15], "Cluster2"][,1],
  Differentiation@meta.data[ Differentiation@meta.data$EdgeID == EdgeID.df$EdgeID[16], "Cluster2"][,1],
  outline=F, las=1, main = "Adipocyte signature", names = as.character(EdgeID.df$EdgeID), ylab="Scaled score")

rm(EdgeID, EdgeID.branchs, EdgeID.branchs, # PartStruct,ProjStruct,
   PHATE, PlotPG, sample.meta, sen.meta, tree.data,
   Tree_e2e, Tree_Graph, TreeEPG)



# Trajectory validation using additional tools -----------------------------------------------------
## Trajectory
# Set labels (adipocytes set to a single label to focus on differentiation)
DimPlot(Differentiation, reduction = "phate", group.by = "Label", label = TRUE, raster=FALSE)
table(Differentiation@meta.data[["Annotation"]] , Differentiation@meta.data[["seurat_clusters"]] )
table(Differentiation@meta.data[["Annotation"]] , Differentiation@meta.data[["Label"]] )
Differentiation$Adipogenesis.Label <- "NA"
Differentiation@meta.data[ Differentiation@meta.data$Annotation == "Adipocyte" & Differentiation@meta.data$seurat_clusters == "0", "Adipogenesis.Label"] <- 1
Differentiation@meta.data[ Differentiation@meta.data$Annotation == "Adipocyte" & Differentiation@meta.data$seurat_clusters == "3", "Adipogenesis.Label"] <- 1
Differentiation@meta.data[ Differentiation@meta.data$Annotation == "FAP" & Differentiation@meta.data$seurat_clusters == "3", "Adipogenesis.Label"] <- 2
Differentiation@meta.data[ Differentiation@meta.data$Annotation == "FAP" & Differentiation@meta.data$seurat_clusters == "1", "Adipogenesis.Label"] <- 3

# Extract pseudo-time from the original object
ElPiTime <- as.data.frame(Differentiation@meta.data[,c("Label",paste0(path, "_Pt"))])
ElPiTime$Cell <- rownames(ElPiTime)
colnames(ElPiTime)[2] <- "ElPiGraph"

## Test with different tools
# Slingshot
Test <- slingshot(Embeddings(Differentiation, "phate"), clusterLabels=Differentiation$Adipogenesis.Label)
SlingTime <- as.data.frame(slingPseudotime(Test))
SlingTime$Cell <- rownames(SlingTime)
colnames(SlingTime)[1] <- "Slingshot"

# TSCAN
Test <- exprmclust(as.data.frame(t(Embeddings(Differentiation, "phate"))), reduce=F,clusternum=2:50)
TSCANtime <- TSCANorder(Test, orderonly=F)
colnames(TSCANtime)[c(1,3)] <- c("Cell","TSCAN")

# DPT
DM <- DiffusionMap(Embeddings(Differentiation, "harmony"))
dpt <- DPT(DM)
DESTINY <- as.data.frame(DPT$DPT1942) # First cell in same trajectory as ElPiGraph.R
DESTINY$Cell <- rownames(Differentiation@meta.data)
colnames(DESTINY)[1] <- "Destiny"

# SCORPIUS
space <- reduce_dimensionality(Embeddings(Differentiation, "harmony"), dist = "pearson", ndim = 2)
traj <- infer_trajectory(space)
SCORPIUS <- as.data.frame(traj$time)
SCORPIUS$Cell <- rownames(SCORPIUS)
colnames(SCORPIUS)[1] <- "SCORPIUS"

# Monocle3
MD <- data.frame(gene_short_name = rownames(Differentiation))
rownames(MD) <- rownames(Differentiation)
cds <- new_cell_data_set(Differentiation@assays$originalexp@counts, cell_metadata = Differentiation@meta.data, gene_metadata = MD)
cds <- preprocess_cds(cds, num_dim = 50, verbose=T, use_genes = VariableFeatures(Differentiation))
reducedDims(cds)[["Aligned"]] <- as.matrix(Embeddings(Differentiation, "harmony"))
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds)
cds@clusters[["UMAP"]]$clusters <- as.factor(Differentiation$Adipogenesis.Label)
cds <- learn_graph(cds)
cds <- order_cells(cds, root_cells = names(which.min(Differentiation$Pt)))
Mcl3 <- as.data.frame(cds@principal_graph_aux[["UMAP"]]$pseudotime)
Mcl3$Cell <- rownames(Mcl3)
colnames(Mcl3)[1] <- "Monocle3"

## Combine results
Pt <- merge(ElPiTime, SlingTime, by="Cell")
Pt <- merge(Pt, TSCANtime[,c(1,3)], by="Cell")
Pt <- merge(Pt, SCORPIUS, by="Cell")
Pt <- merge(Pt, DESTINY, by="Cell")
Pt <- merge(Pt, Mcl3, by="Cell")
Pt <- Pt[ !is.na(Pt$ElPiGraph),]

# Calculate and plot the correlation
barplot(abs(cor(Pt[,3:8], method="spearman")[c(2,3,5,4,6),1]), ylim=c(0,1), las=1)

## Meta-stable states
# Scorpius
H_R2 <- hist(Pt[ Pt$Cell %in% rownames(Differentiation@meta.data[ grep("R1",Differentiation@meta.data$Dataset),]),"SCORPIUS"], breaks=seq(0,1.0,by=0.05))
H_R1 <- hist(Pt[ Pt$Cell %in% rownames(Differentiation@meta.data[ grep("R2",Differentiation@meta.data$Dataset),]),"SCORPIUS"], breaks=seq(0,1.0,by=0.05))
Blt <- barplot(rowMeans(cbind(H_R1$counts/sum(H_R1$counts),H_R2$counts/sum(H_R2$counts))), ylim=c(0, 0.2))
arrows(Blt, rowMeans(cbind(H_R1$counts/sum(H_R1$counts),H_R2$counts/sum(H_R2$counts)))+(apply(cbind(H_R1$counts/sum(H_R1$counts),H_R2$counts/sum(H_R2$counts)),1,FUN="sd")/sqrt(2)),Blt, rowMeans(cbind(H_R1$counts/sum(H_R1$counts),H_R2$counts/sum(H_R2$counts)))-(apply(cbind(H_R1$counts/sum(H_R1$counts),H_R2$counts/sum(H_R2$counts)),1,FUN="sd")/sqrt(2)),code = 3, angle = 90, len=0.05)

# Destiny
H_R2 <- hist(Pt[ Pt$Cell %in% rownames(Differentiation@meta.data[ grep("R1",Differentiation@meta.data$Dataset),]),"Destiny"], breaks=seq(0,0.5,by=0.025))
H_R1 <- hist(Pt[ Pt$Cell %in% rownames(Differentiation@meta.data[ grep("R2",Differentiation@meta.data$Dataset),]),"Destiny"], breaks=seq(0,0.5,by=0.025))
Blt <- barplot(rowMeans(cbind(H_R1$counts/sum(H_R1$counts),H_R2$counts/sum(H_R2$counts))), ylim=c(0, 0.5))
arrows(Blt, rowMeans(cbind(H_R1$counts/sum(H_R1$counts),H_R2$counts/sum(H_R2$counts)))+(apply(cbind(H_R1$counts/sum(H_R1$counts),H_R2$counts/sum(H_R2$counts)),1,FUN="sd")/sqrt(2)),Blt, rowMeans(cbind(H_R1$counts/sum(H_R1$counts),H_R2$counts/sum(H_R2$counts)))-(apply(cbind(H_R1$counts/sum(H_R1$counts),H_R2$counts/sum(H_R2$counts)),1,FUN="sd")/sqrt(2)),code = 3, angle = 90, len=0.05)

# Monocle3
H_R2 <- hist(Pt[ Pt$Cell %in% rownames(Differentiation@meta.data[ grep("R1",Differentiation@meta.data$Dataset),]),"Monocle3"], breaks=seq(0,30,by=0.75))
H_R1 <- hist(Pt[ Pt$Cell %in% rownames(Differentiation@meta.data[ grep("R2",Differentiation@meta.data$Dataset),]),"Monocle3"], breaks=seq(0,30,by=0.75))
Blt <- barplot(rowMeans(cbind(H_R1$counts/sum(H_R1$counts),H_R2$counts/sum(H_R2$counts))), ylim=c(0, 0.2))
arrows(Blt, rowMeans(cbind(H_R1$counts/sum(H_R1$counts),H_R2$counts/sum(H_R2$counts)))+(apply(cbind(H_R1$counts/sum(H_R1$counts),H_R2$counts/sum(H_R2$counts)),1,FUN="sd")/sqrt(2)),Blt, rowMeans(cbind(H_R1$counts/sum(H_R1$counts),H_R2$counts/sum(H_R2$counts)))-(apply(cbind(H_R1$counts/sum(H_R1$counts),H_R2$counts/sum(H_R2$counts)),1,FUN="sd")/sqrt(2)),code = 3, angle = 90, len=0.05)




### Transcription factor waves ---------------------------------------------------------------------------------
# Import list of TFs in mice
TF <- read.delim("Mus_musculus_TF.txt", header=T)

# Subset smoothend expression of genes associated with pseudo-time in LFD to only TFs
SmoothTF <- Smooth[ rownames(Smooth) %in% TF[,2],]

# Seriate the result
SmoothTF <- SmoothTF[ get_order(seriate(SmoothTF, method="PCA_angle")),]

# Create heatmaps (manual reordering to make it look nice)
HTM1 <- Heatmap(SmoothTF[c(1,3,4,9,6,7,2,11,5,10,8,12,13),], cluster_columns=F, cluster_rows=F)
HTM2 <- Heatmap(SmoothTF[c(15,14,18,17,23,16,22,19,21,20,24),], cluster_columns=F, cluster_rows=F)
HTM3 <- Heatmap(SmoothTF[c(25,29,30,31,32,33,27,34,28,26,35,36,37,38,39,40),], cluster_columns=F, cluster_rows=F)
HTM4 <- Heatmap(SmoothTF[c(41,42,48,72,52,43,44,46,45,47,50,61,49,69,51,55,53,54,60,59,58,65,62,68,56,71,64,57,74,75,63,70,66,67,73,76,77),], cluster_columns=F, cluster_rows=F)

# Plot the heatmap
HTM3 %v% HTM4 %v% HTM1 %v% HTM2

## Differential association with pseudo-time in LFD and HFD ------------------------------------------------------------
## Pattern test
# Setup: R1
Differentiation_R1 <- subset(Differentiation, subset = Replicate == "R1")
Pseudotime <- data.frame(curve1 = Differentiation_R1$Pt, curve2 = Differentiation_R1$Pt)
Weights <- data.frame(LFD = rep(0, ncol(Differentiation_R1)),HFD = rep(0, ncol(Differentiation_R1)))
rownames(Weights) = colnames(Differentiation_R1)
Weights[ rownames(Weights) %in% rownames(Differentiation_R1@meta.data[Differentiation_R1@meta.data$Diet == "LFD",]),"LFD"] <- 1
Weights[ rownames(Weights) %in% rownames(Differentiation_R1@meta.data[Differentiation_R1@meta.data$Diet == "HFD",]),"HFD"] <- 1
Counts <- Differentiation_R1@assays$originalexp@counts
Counts <- Counts[ rownames(Counts) %in% expressed_genes,]
R1_sce <- fitGAM(counts = Counts, pseudotime = Pseudotime, cellWeights = Weights, nknots = 4, verbose = TRUE)

# Setup: R2
Differentiation_R1 <- subset(Differentiation, subset = Replicate == "R2")
Pseudotime <- data.frame(curve1 = Differentiation_R1$Pt, curve2 = Differentiation_R1$Pt)
Weights <- data.frame(LFD = rep(0, ncol(Differentiation_R1)),HFD = rep(0, ncol(Differentiation_R1)))
rownames(Weights) = colnames(Differentiation_R1)
Weights[ rownames(Weights) %in% rownames(Differentiation_R1@meta.data[Differentiation_R1@meta.data$Diet == "LFD",]),"LFD"] <- 1
Weights[ rownames(Weights) %in% rownames(Differentiation_R1@meta.data[Differentiation_R1@meta.data$Diet == "HFD",]),"HFD"] <- 1
Counts <- Differentiation_R1@assays$originalexp@counts
Counts <- Counts[ rownames(Counts) %in% expressed_genes,]
R2_sce <- fitGAM(counts = Counts, pseudotime = Pseudotime, cellWeights = Weights, nknots = 4, verbose = TRUE)

## Tests and filtering
# Perform tests
Pat_R1 <- as.data.frame(patternTest(R1_sce))
Pat_R2 <- as.data.frame(patternTest(R2_sce))

# Subset to genes associated with pseudo-time
Pat_R2 <- Pat_R2[ rownames(Pat_R2) %in% Ass,]
Pat_R1 <- Pat_R1[ rownames(Pat_R1) %in% Ass,]

# FDR-adjustment
Pat_R1$FDR <- p.adjust(Pat_R1$pvalue, method="fdr")
Pat_R2$FDR <- p.adjust(Pat_R2$pvalue, method="fdr")

# Subset to significant ones
Pat_R1 <- Pat_R1[ Pat_R1$FDR <= 0.05,]
Pat_R2 <- Pat_R2[ Pat_R2$FDR <= 0.05,]

# Remove NAs
Pat_R1 <- Pat_R1[ grep("NA", rownames(Pat_R1), invert=T),]
Pat_R2 <- Pat_R2[ grep("NA", rownames(Pat_R2), invert=T),]

# Combine the results
Pat <- rownames(Pat_R1[ Pat_R1$FDR <= 0.05 & rownames(Pat_R1) %in% rownames(Pat_R2[ Pat_R2$FDR <= 0.05,]),])

### Heatmap of differentially associated TFs -----------------------------------------
# Import list of mouse TFs
TF <- read.delim("Mus_musculus_TF.txt", header=T)

# Subset results to only TFs
Pat_TFs <- Pat[ Pat %in% TF[,2]]

# Predict smoothend expression
Smooth_R1 <- predictSmooth(R1_sce, gene = Pat_TFs, tidy=F, n=100)
Smooth_R2 <- predictSmooth(R2_sce, gene = Pat_TFs, tidy=F, n=100)

# Average smoothend expression across replicates and scale
Smooth <- Smooth_R1
for (i in 1:nrow(Smooth)) { Smooth[i,] <- colMeans(rbind(Smooth_R1[i,], Smooth_R2[i,])) }
Smooth <- t(scale(t(Smooth)))

# Reorder the data (manually to make it look good visually)
Smooth <- Smooth[ c(5,17,25,8,3,9,13,6,28,10,24,19,2,4,27,7,18,23,20,1,14,11,22,12,16,15,29,26,21),]

# Plot the heatmap
Heatmap( Smooth, cluster_columns=F, cluster_rows=F, right_annotation = rowAnnotation(bar = c("U","I","U","U","I","U","I","U","I","U","U","I","A","I","U","I","U","A","U","A","A","A","I","U","U","U","I","U","U")))

# Cell-cell signaling ------------------------------------------------------------------------
## Process data from CellPhoneDB
# Import gene lists and annotations (downloaded from here https://www.cellphonedb.org/downloads)
Genes_Annotations <- read.delim("gene_input.txt", header=T, sep=",")
Complex_Annotation <- read.delim("complex_input.txt", header=T, sep=",")
Interactions <- read.delim("interaction_curated.txt", sep = ",", header=T)
Complexes <- read.delim("complex_curated.txt", sep = ",")
Proteins <- read.delim("protein_curated.txt", sep = ",")

# Setup to loop
Pairs <- as.data.frame(matrix(ncol=7, nrow=1))
colnames(Pairs) <- c("ID","Receptor","Ligand","Receptor_Genes","Ligand_Genes", "Type_Receptor", "Type_Ligand")

# Loop through all interactions and save information about the partners and their interactions
Counter <- 1
for (i in 1:nrow(Interactions)) {
  # Check if partner a is a secreted factor
  if (nrow(Complexes[ Complexes$secreted == "TRUE" & Complexes[,1] == as.character(Interactions[i, "partner_a"]),]) == 1 | nrow(Proteins[ Proteins$uniprot == as.character(Interactions[i, "partner_a"]) & Proteins$secreted == "True",]) == 1) {
    # Check if partner b is a receptor
    if (nrow(Complexes[ Complexes$receptor == "TRUE" & Complexes[,1] == as.character(Interactions[i, "partner_b"]),]) == 1 | nrow(Proteins[ Proteins$uniprot == as.character(Interactions[i, "partner_b"]) & Proteins$receptor == "True",]) == 1) {
      # Save the interaction information
      Pairs[Counter,1] <- as.character(Interactions[i,1])
      if (nrow(Complexes[ Complexes$secreted == "TRUE" & Complexes[,1] == as.character(Interactions[i, "partner_a"]),]) == 1) {
        Pairs[Counter,"Ligand"] <- as.character(Complexes[ Complexes$secreted == "TRUE" & Complexes[,1] == as.character(Interactions[i, "partner_a"]),1])
        Pairs[Counter,"Ligand_Genes"] <- paste0(Genes_Annotations[ Genes_Annotations$uniprot %in% as.character(unlist(Complexes[ Complexes$secreted == "TRUE" & Complexes[,1] == as.character(Interactions[i, "partner_a"]),c(2,3,4,5)])),3], collapse = ",")
        Pairs[Counter,"Type_Ligand"] <- "Complex"
      }
      if (nrow(Proteins[ Proteins$uniprot == as.character(Interactions[i, "partner_a"]) & Proteins$secreted == "True",]) == 1) {
        Pairs[Counter,"Ligand"] <- paste0(as.character(Genes_Annotations[ Genes_Annotations$uniprot %in% Proteins[ Proteins$uniprot == as.character(Interactions[i, "partner_a"]) & Proteins$secreted == "True",1],3]), collapse = ",")
        Pairs[Counter,"Ligand_Genes"] <- paste0(as.character(Genes_Annotations[ Genes_Annotations$uniprot %in% Proteins[ Proteins$uniprot == as.character(Interactions[i, "partner_a"]) & Proteins$secreted == "True",1],3]), collapse = ",")
        Pairs[Counter,"Type_Ligand"] <- "Interaction"
      }
      
      if (nrow(Complexes[ Complexes$receptor == "TRUE" & Complexes[,1] == as.character(Interactions[i, "partner_b"]),]) == 1) {
        Pairs[Counter,"Receptor"] <- as.character(Complexes[ Complexes$receptor == "TRUE" & Complexes[,1] == as.character(Interactions[i, "partner_b"]),1])
        Pairs[Counter,"Receptor_Genes"] <- paste0(Genes_Annotations[ Genes_Annotations$uniprot %in% as.character(unlist(Complexes[ Complexes$receptor == "TRUE" & Complexes[,1] == as.character(Interactions[i, "partner_b"]),c(2,3,4,5)])),3], collapse = ",")
        Pairs[Counter,"Type_Receptor"] <- "Complex"
      }
      if (nrow(Proteins[ Proteins$uniprot == as.character(Interactions[i, "partner_b"]) & Proteins$receptor == "True",]) == 1) {
        Pairs[Counter,"Receptor"] <- paste0(as.character(Genes_Annotations[ Genes_Annotations$uniprot %in% Proteins[ Proteins$uniprot == as.character(Interactions[i, "partner_b"]) & Proteins$receptor == "True",1],3]), collapse = ",")
        Pairs[Counter,"Receptor_Genes"] <- paste0(as.character(Genes_Annotations[ Genes_Annotations$uniprot %in% Proteins[ Proteins$uniprot == as.character(Interactions[i, "partner_b"]) & Proteins$receptor == "True",1],3]), collapse = ",")
        Pairs[Counter,"Type_Receptor"] <- "Interaction"
      }
      Counter <- Counter + 1
    }
  }
  # Check if partner b is a secreted factor
  if (nrow(Complexes[ Complexes$secreted == "TRUE" & Complexes[,1] == as.character(Interactions[i, "partner_b"]),]) == 1 | nrow(Proteins[ Proteins$uniprot == as.character(Interactions[i, "partner_b"]) & Proteins$secreted == "True",]) == 1) {
    # Check if partner b is a receptor
    if (nrow(Complexes[ Complexes$receptor == "TRUE" & Complexes[,1] == as.character(Interactions[i, "partner_a"]),]) == 1 | nrow(Proteins[ Proteins$uniprot == as.character(Interactions[i, "partner_a"]) & Proteins$receptor == "True",]) == 1) {
      # Save the interaction information
      Pairs[Counter,1] <- as.character(Interactions[i,1])
      if (nrow(Complexes[ Complexes$secreted == "TRUE" & Complexes[,1] == as.character(Interactions[i, "partner_b"]),]) == 1) {
        Pairs[Counter,"Ligand"] <- as.character(Complexes[ Complexes$secreted == "TRUE" & Complexes[,1] == as.character(Interactions[i, "partner_b"]),1])
        Pairs[Counter,"Ligand_Genes"] <- paste0(Genes_Annotations[ Genes_Annotations$uniprot %in% as.character(unlist(Complexes[ Complexes$secreted == "TRUE" & Complexes[,1] == as.character(Interactions[i, "partner_b"]),c(2,3,4,5)])),3], collapse = ",")
        Pairs[Counter,"Type_Ligand"] <- "Complex"
      }
      if (nrow(Proteins[ Proteins$uniprot == as.character(Interactions[i, "partner_b"]) & Proteins$secreted == "True",]) == 1) {
        Pairs[Counter,"Ligand"] <- paste0(as.character(Genes_Annotations[ Genes_Annotations$uniprot %in% Proteins[ Proteins$uniprot == as.character(Interactions[i, "partner_b"]) & Proteins$secreted == "True",1],3]), collapse = ",")
        Pairs[Counter,"Ligand_Genes"] <- paste0(as.character(Genes_Annotations[ Genes_Annotations$uniprot %in% Proteins[ Proteins$uniprot == as.character(Interactions[i, "partner_b"]) & Proteins$secreted == "True",1],3]), collapse = ",")
        Pairs[Counter,"Type_Ligand"] <- "Interaction"
      }
      
      if (nrow(Complexes[ Complexes$receptor == "TRUE" & Complexes[,1] == as.character(Interactions[i, "partner_a"]),]) == 1) {
        Pairs[Counter,"Receptor"] <- as.character(Complexes[ Complexes$receptor == "TRUE" & Complexes[,1] == as.character(Interactions[i, "partner_a"]),1])
        Pairs[Counter,"Receptor_Genes"] <- paste0(Genes_Annotations[ Genes_Annotations$uniprot %in% as.character(unlist(Complexes[ Complexes$receptor == "TRUE" & Complexes[,1] == as.character(Interactions[i, "partner_a"]),c(2,3,4,5)])),3], collapse = ",")
        Pairs[Counter,"Type_Receptor"] <- "Complex"
      }
      if (nrow(Proteins[ Proteins$uniprot == as.character(Interactions[i, "partner_a"]) & Proteins$receptor == "True",]) == 1) {
        Pairs[Counter,"Receptor"] <- paste0(as.character(Genes_Annotations[ Genes_Annotations$uniprot %in% Proteins[ Proteins$uniprot == as.character(Interactions[i, "partner_a"]) & Proteins$receptor == "True",1],3]), collapse = ",")
        Pairs[Counter,"Receptor_Genes"] <- paste0(as.character(Genes_Annotations[ Genes_Annotations$uniprot %in% Proteins[ Proteins$uniprot == as.character(Interactions[i, "partner_a"]) & Proteins$receptor == "True",1],3]), collapse = ",")
        Pairs[Counter,"Type_Receptor"] <- "Interaction"
      }
      Counter <- Counter + 1
    }
  }
}

## Import gene mappings between human and mouse (from Ensembl biomart)
Convert <- read.delim("Human_Mouse.txt", header=T)

## First analyze expression of the receptors in early pre-adipocytes
# Subset
Receptors <- Pairs[ duplicated(Pairs$V2) == F,]

# Split into wether its a complex or single receptor.
Complexes <- Receptors[ Receptors$V6 == "Complex",]
Interactions <- Receptors[ Receptors$V6 != "Complex",]

# Process complexes. For each complex, find the lowest expressed member in each diet.
for (i in 1:nrow(Complexes)) {
  Tmp <- Convert[ Convert[,1] %in% unlist(strsplit(as.character(Complexes[i,4]),",")),]
  for (q in 1:nrow(Tmp)) {
    if (sum(rownames(Differentiation@assays$originalexp@data) %in% as.character(Tmp[q,2])) == 1) {
      Tmp[q,3] <- mean(Differentiation@assays$originalexp@data[ rownames(Differentiation@assays$originalexp@data) %in% as.character(Tmp[q,2]), colnames(Differentiation@assays$originalexp@data) %in% rownames(Differentiation@meta.data[ Differentiation@meta.data$Type == "Early_Pre" & Differentiation$Diet == "LFD",])])
      Tmp[q,4] <- mean(Differentiation@assays$originalexp@data[ rownames(Differentiation@assays$originalexp@data) %in% as.character(Tmp[q,2]), colnames(Differentiation@assays$originalexp@data) %in% rownames(Differentiation@meta.data[ Differentiation@meta.data$Type == "Early_Pre" & Differentiation$Diet == "HFD",])])
    } else {
      Tmp[q,3] <- 0
      Tmp[q,4] <- 0
    }
  }
  # LFD
  Tmp <- Tmp[ order(Tmp[,1], -Tmp[,3]),]
  Tmp2 <- Tmp[ duplicated(Tmp$Gene.name) == F,]
  Complexes[i,8] <- min(Tmp2[,3])
  # HFD
  Tmp <- Tmp[ order(Tmp[,1], -Tmp[,4]),]
  Tmp2 <- Tmp[ duplicated(Tmp$Gene.name) == F,]
  Complexes[i,9] <- min(Tmp2[,4])
}

# Process interactions.
for (i in 1:nrow(Interactions)) {
  if (sum(rownames(Differentiation@assays$originalexp@data) %in% as.character(Convert[ Convert[,1] %in% unlist(strsplit(as.character(Interactions[i,4]),",")),2])) == 0) {
    Interactions[i,8] <- 0
    Interactions[i,9] <- 0
  }
  if (sum(rownames(Differentiation@assays$originalexp@data) %in% as.character(Convert[ Convert[,1] %in% unlist(strsplit(as.character(Interactions[i,4]),",")),2])) == 1) {
    Interactions[i,8] <- mean(Differentiation@assays$originalexp@data[ rownames(Differentiation@assays$originalexp@data) %in% as.character(Convert[ Convert[,1] %in% unlist(strsplit(as.character(Interactions[i,4]),",")),2]), colnames(Differentiation@assays$originalexp@data) %in% rownames(Differentiation@meta.data[ Differentiation@meta.data$Type == "Early_Pre" & Differentiation$Diet == "LFD",])])
    Interactions[i,9] <- mean(Differentiation@assays$originalexp@data[ rownames(Differentiation@assays$originalexp@data) %in% as.character(Convert[ Convert[,1] %in% unlist(strsplit(as.character(Interactions[i,4]),",")),2]), colnames(Differentiation@assays$originalexp@data) %in% rownames(Differentiation@meta.data[ Differentiation@meta.data$Type == "Early_Pre" & Differentiation$Diet == "HFD",])])
  }
  if (sum(rownames(Differentiation@assays$originalexp@data) %in% as.character(Convert[ Convert[,1] %in% unlist(strsplit(as.character(Interactions[i,4]),",")),2])) > 1) {
    Interactions[i,8] <- Interactions[i,8] <- max(Matrix::rowMeans(Differentiation@assays$originalexp@data[ rownames(Differentiation@assays$originalexp@data) %in% as.character(Convert[ Convert[,1] %in% unlist(strsplit(as.character(Interactions[i,4]),",")),2]), colnames(Differentiation@assays$originalexp@data) %in% rownames(Differentiation@meta.data[ Differentiation@meta.data$Type == "Early_Pre" & Differentiation$Diet == "LFD",])]))
    Interactions[i,9] <- Interactions[i,8] <- max(Matrix::rowMeans(Differentiation@assays$originalexp@data[ rownames(Differentiation@assays$originalexp@data) %in% as.character(Convert[ Convert[,1] %in% unlist(strsplit(as.character(Interactions[i,4]),",")),2]), colnames(Differentiation@assays$originalexp@data) %in% rownames(Differentiation@meta.data[ Differentiation@meta.data$Type == "Early_Pre" & Differentiation$Diet == "HFD",])]))
  }
}

# Combine the information
Receptors <- rbind(Complexes, Interactions)

# Subset to only expressed receptors (higher than 0.1 in each LFD or HFD)
Receptors <- Receptors[ Receptors[,9] >= 0.1 | Receptors[,8] >= 0.1,] # 32 receptors

## Now analyze expression of the ligands in cell types endogenous to the adipose tissue
# Subset to ligands of expressed receptors
Ligands <- Pairs[ Pairs[,2] %in% Receptors[,2],]
Ligands <- Ligands[ duplicated(Ligands[,3]) == F,]
Ligands[,5] <- as.character(Ligands[,5])

# Loop through all ligands
AllLigands <- data.frame()
for (i in 2:nrow(Ligands)) {
  Tmp <- Convert[ Convert[,1] %in% unlist(strsplit(Ligands[i,5],",")),]
  Tmp$Adipocytes_LFD <- 0
  Tmp$Adipocytes_HFD <- 0
  Tmp$Mesothelial_LFD <- 0
  Tmp$Mesothelial_HFD <- 0
  Tmp$Endothelial_LFD <- 0
  Tmp$Endothelial_HFD <- 0
  Tmp$FAP_LFD <- 0
  Tmp$FAP_HFD <- 0
  Tmp$Immune_LFD <- 0
  Tmp$Immune_HFD <- 0
  # Loop over all genes homologous to the ligands and get their mean expression
  for (q in 1:nrow(Tmp)) {
    if (sum(rownames(eWAT@assays$originalexp@counts) %in% as.character(Tmp[q,2])) > 0) {
      Tmp[q,"Adipocytes_LFD"] <- mean(eWAT@assays$originalexp@counts[ rownames(eWAT@assays$originalexp@counts) %in% as.character(Tmp[q,2]), colnames(eWAT@assays$originalexp@counts) %in% rownames(eWAT@meta.data[ eWAT@meta.data$Annotation == "Adipocyte" & Data$Diet == "LFD",])])
      Tmp[q,"Adipocytes_HFD"] <- mean(eWAT@assays$originalexp@counts[ rownames(eWAT@assays$originalexp@counts) %in% as.character(Tmp[q,2]), colnames(eWAT@assays$originalexp@counts) %in% rownames(eWAT@meta.data[ eWAT@meta.data$Annotation == "Adipocyte" & Data$Diet == "HFD",])])
      Tmp[q,"Mesothelial_LFD"] <- mean(eWAT@assays$originalexp@counts[ rownames(eWAT@assays$originalexp@counts) %in% as.character(Tmp[q,2]), colnames(eWAT@assays$originalexp@counts) %in% rownames(eWAT@meta.data[ eWAT@meta.data$Annotation == "Mesothelial" & Data$Diet == "LFD",])])
      Tmp[q,"Mesothelial_HFD"] <- mean(eWAT@assays$originalexp@counts[ rownames(eWAT@assays$originalexp@counts) %in% as.character(Tmp[q,2]), colnames(eWAT@assays$originalexp@counts) %in% rownames(eWAT@meta.data[ eWAT@meta.data$Annotation == "Mesothelial" & Data$Diet == "HFD",])])
      Tmp[q,"Endothelial_LFD"] <- mean(eWAT@assays$originalexp@counts[ rownames(eWAT@assays$originalexp@counts) %in% as.character(Tmp[q,2]), colnames(eWAT@assays$originalexp@counts) %in% rownames(eWAT@meta.data[ eWAT@meta.data$Annotation == "Endothelial" & Data$Diet == "LFD",])])
      Tmp[q,"Endothelial_HFD"] <- mean(eWAT@assays$originalexp@counts[ rownames(eWAT@assays$originalexp@counts) %in% as.character(Tmp[q,2]), colnames(eWAT@assays$originalexp@counts) %in% rownames(eWAT@meta.data[ eWAT@meta.data$Annotation == "Endothelial" & Data$Diet == "HFD",])])
      Tmp[q,"FAP_LFD"] <- mean(eWAT@assays$originalexp@counts[ rownames(eWAT@assays$originalexp@counts) %in% as.character(Tmp[q,2]), colnames(eWAT@assays$originalexp@counts) %in% rownames(eWAT@meta.data[ eWAT@meta.data$Annotation == "FAP" & Data$Diet == "LFD",])])
      Tmp[q,"FAP_HFD"] <- mean(eWAT@assays$originalexp@counts[ rownames(eWAT@assays$originalexp@counts) %in% as.character(Tmp[q,2]), colnames(eWAT@assays$originalexp@counts) %in% rownames(eWAT@meta.data[ eWAT@meta.data$Annotation == "FAP" & Data$Diet == "HFD",])])
      Tmp[q,"Immune_LFD"] <- mean(eWAT@assays$originalexp@counts[ rownames(eWAT@assays$originalexp@counts) %in% as.character(Tmp[q,2]), colnames(eWAT@assays$originalexp@counts) %in% rownames(eWAT@meta.data[ eWAT@meta.data$Annotation == "Immune" & Data$Diet == "LFD",])])
      Tmp[q,"Immune_HFD"] <- mean(eWAT@assays$originalexp@counts[ rownames(eWAT@assays$originalexp@counts) %in% as.character(Tmp[q,2]), colnames(eWAT@assays$originalexp@counts) %in% rownames(eWAT@meta.data[ eWAT@meta.data$Annotation == "Immune" & Data$Diet == "HFD",])])
    }
  }
  AllLigands <- rbind(AllLigands, Tmp)
}

# Remove ligands without names
AllLigands <- AllLigands[ AllLigands[,2] != "",]

# Subset to expressed ligands (>= 0.1) across cell types and diets and merge results
AllLigands <- AllLigands[which(apply(AllLigands[,3:ncol(AllLigands)],1,FUN="max") >= 0.1),]
Ligands <- merge(Ligands, AllLigands[,c(1,3:ncol(AllLigands))], by.x=3, by.y=1) # 17 ligands

## Combine the result of the ligand and receptor analysis
Interactions <- Pairs[ Pairs[,2] %in% Receptors[,2] & Pairs[,3] %in% Ligands[,1],1:3]
colnames(Receptors)[8] <- "Receptor_LFD"
colnames(Receptors)[9] <- "Receptor_HFD"
Interactions <- merge(Interactions, Receptors[,c(2,8,9)], by.x=2, by.y=1)
Interactions <- merge(Interactions, Ligands[,c(1,8:17)], by.x=3, by.y=1) # 18 receptors, 17 ligands, 29 interactions

## Heatmap the expression of receptors in early pre-adipocytes
Receptors <- Interactions[ duplicated(Interactions[,2])==F,]
rownames(Receptors) <- Receptors[,2]
ComplexHeatmap::Heatmap(Receptors[c(5,6,3,4,1,7,8,9,10,11,15,12,13,14,2,17,16,18),4:5], cluster_rows=F, cluster_columns=F, col = colorRamp2(breaks = seq(0,1.65, length.out=11), rev(brewer.pal(11,"Spectral"))))

## Heatmap of the expression of ligands in cell types endogenous to the adipose tissue
Ligands <- Interactions[ duplicated(Interactions[,1])==F,]
rownames(Ligands) <- Ligands[,1]
ComplexHeatmap::Heatmap(Ligands[c(4,3,1,5,10,6,7,8,9,14,17,11,13,12,2,15,16),6:15], cluster_rows=F, cluster_columns=F, col = colorRamp2(breaks = seq(0,1.65, length.out=11), rev(brewer.pal(11,"Spectral"))))





