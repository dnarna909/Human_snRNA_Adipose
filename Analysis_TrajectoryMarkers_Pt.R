### DISCLAIMER
#   Some of the algorithms are non-deterministic making the results slightly different from run to run.
# Many of the algorithms are frequently updated with new releases making the results slightly different from version to version.
# Depending on your system, this code below may not produce an exact 1:1 copy of the results. 
# This adversely affects clustering, filtering etc. if the code is naively copy-pasted.
# For better reproducibility, please download quality-filtered Rds objects and work on downstream analyses. 
# <br>
#   change "D:/" to "/media/jianie/Extreme SSD1/" # if in LUNIX computer

### Load libraries -------------------------------------------------------------------------------------------------------------
#rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.
memory.limit(size = 1e+13)

# https://kstreet13.github.io/bioc2020trajectories/articles/workshopTrajectories.html
# source(paste0("/media/jianie/Extreme SSD1/2022-09-01 STARR Combine/Project Parameters.R"), local = knitr::knit_global())
source(paste0(Disk, Project.folder, "/", "Project Parameters_GCCRI.R"), local = knitr::knit_global())
library(tradeSeq)
suppressPackageStartupMessages({
  library(slingshot); library(SingleCellExperiment)
  library(RColorBrewer); library(scales)
  library(viridis); library(UpSetR)
  library(pheatmap); library(msigdbr)
  library(fgsea); library(knitr)
  library(ggplot2); library(gridExtra)
  library(tradeSeq)
})

# prepare folder
export.folder <- paste0("TrajectoryMarkers")
dir.create(file.path(paste0(Disk, Project.folder, "/", Rds.folder, "/"), export.folder), showWarnings = FALSE)
dir <- paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder, "/");dir

## Association tests ----------------------------------------------------
### NOTE: This step is non-deterministic! Results vary from run-to-run
# Setup
# # Differential expression along trajectory (with tradeSeq) --------------------------------------------------
# ## Find expressed genes (in at least 10% of the nuclei)
# sca <-  FromMatrix(as.matrix(Differentiation@assays$originalexp@data), Differentiation@meta.data)
# expressed_genes <- names(which(freq(sca) >= 0.1))
# if(exists("expressed_genes")){rm(sca)};gc()

seurat.file
path
Select.Group 
compare.group 
groups

if (exists("Differentiation") == F){
  Differentiation <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/",seurat.file) )
  Differentiation <- Seurat::SetIdent(Differentiation, value = Select.Group)
  Differentiation <- subset(Differentiation, idents = TRUE);gc()
} 
Counts <- Differentiation@assays$originalexp@counts; dim(Counts)
# Counts <- Counts[ rownames(Counts) %in% expressed_genes,]
BPPARAM <- BiocParallel::bpparam()
BPPARAM$workers <- 4 # use 4 cores
Pseudotime <- data.frame(matrix(data = Differentiation@meta.data[[paste0(path, "_Pt")]], 
                                ncol = 1, nrow = length(colnames(Differentiation))))
colnames(Pseudotime) <- paste0("Pt")
Weights <- data.frame(matrix(data = 1, ncol = 1, nrow = length(colnames(Differentiation))))
colnames(Weights) <- paste0("Weights")
rownames(Weights) = colnames(Differentiation)

set.seed(3)
icMat <- evaluateK(counts = Counts,
                   pseudotime = Pseudotime,
                   cellWeights = Weights,
                   conditions = factor(Differentiation@meta.data[[compare.group]]),
                   nGenes = 300,
                   k = 3:7)
sce <- tradeSeq::fitGAM(counts = as.matrix(Counts), pseudotime = as.matrix(Pseudotime), cellWeights = as.matrix(Weights), nknots = 5, verbose = TRUE,
                        conditions = factor(Differentiation@meta.data[[compare.group]]), 
                        parallel=TRUE, BPPARAM = BPPARAM );gc()# average 1 hr per subject

rowData(sce)$assocRes <- associationTest(sce, lineages = TRUE, l2fc = log2(2))
AssocRes <- as.data.frame(rowData(sce)$assocRes)
AssocRes$FDR <- p.adjust(AssocRes$pvalue, "fdr")

condRes <- as.data.frame(conditionTest(sce, l2fc = log2(2)))
condRes$FDR <- p.adjust(condRes$pvalue, "fdr")
mean(condRes$FDR <= 0.05, na.rm = TRUE)

save(AssocRes, condRes, 
     file = paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder,  "/", unlist(stringr::str_split(seurat.file, ".Rds"))[1], "_Pt.DEGs.RData" ))
saveRDS(sce, file = paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder,  "/", unlist(stringr::str_split(seurat.file, ".Rds"))[1], ".sce.Rds" ) )

rm(sce, AssocRes, condRes );gc()

rm(Counts, Differentiation, Differentiation, Pseudotime, Weights, 
   ls.names, ff, BPPARAM)
gc()

rm(dir)
gc()
sessionInfo()