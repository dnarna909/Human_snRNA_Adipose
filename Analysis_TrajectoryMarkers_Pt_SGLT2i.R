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

# extract each sample file name
if (file.exists(paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder,  "/", "Differentiation.ls.names.Rds" ) ) == F )  {
  Differentiation <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/",seurat.file) )
  Differentiation <- Seurat::SetIdent(Differentiation, value = Select.Group)
  Differentiation <- subset(Differentiation, idents = TRUE);gc()
  ls.names <- unique(Differentiation@meta.data[["patient_id"]])
  saveRDS(unique(Differentiation@meta.data[["patient_id"]]), paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder,  "/", "Differentiation.ls.names.Rds" ))
} else {
  ls.names = readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder,  "/", "Differentiation.ls.names.Rds") )
}

# calculate the unfinshed sample file name
finished.names <- list.files( paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder, "/"))[grep(paste0(".sce.Rds"), 
                                                                                          list.files( paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder, "/")))];
finished.names;
need.names <- setdiff(ls.names , stringr::str_split_fixed(finished.names, "[.]", n=2)[,1] );need.names;

if (length(need.names) > 0 ) {
  if (exists("Differentiation") == F){
    Differentiation <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/",seurat.file) )
    Differentiation <- Seurat::SetIdent(Differentiation, value = Select.Group)
    Differentiation <- subset(Differentiation, idents = TRUE);gc()
  } 
  for (dd in need.names) {
    # dd = need.names[1]; dd
    Differentiation.ls <- subset(Differentiation, subset = patient_id == dd);gc()
    # Pseudotime <- data.frame(matrix(data = Differentiation.ls@meta.data[[paste0(path, "_Pt")]], 
    #                                 ncol = length(groups), nrow = length(colnames(Differentiation.ls))))
    # colnames(Pseudotime) <- groups
    # Weights <- data.frame(matrix(data = 0, ncol = length(groups), nrow = length(colnames(Differentiation.ls))))
    # colnames(Weights) <- groups
    # rownames(Weights) = colnames(Differentiation.ls)
    # for (ii in groups) {
    #   Weights[ rownames(Weights) %in% rownames(Differentiation.ls@meta.data[Differentiation.ls@meta.data[[compare.group]] == ii,]), ii] <- 1
    # }
    Counts <- Differentiation.ls@assays$originalexp@counts; dim(Counts)
    # Counts <- Counts[ rownames(Counts) %in% expressed_genes,]
    # set.seed(3)
    # icMat <- evaluateK(counts = Counts,
    #                    pseudotime = Pseudotime,
    #                    cellWeights = Weights,
    #                    conditions = factor(Differentiation.ls@meta.data[[compare.group]]),
    #                    nGenes = 300,
    #                    k = 3:7)
    BPPARAM <- BiocParallel::bpparam()
    BPPARAM$workers <- 4 # use 4 cores
    # sce <- tradeSeq::fitGAM(counts = Counts, pseudotime = Pseudotime, cellWeights = Weights, nknots = 5, verbose = TRUE,
    #                         parallel=TRUE, BPPARAM = BPPARAM
    #                         );gc()# average 1 hr per subject
    Pseudotime <- data.frame(matrix(data = Differentiation.ls@meta.data[[paste0(path, "_Pt")]], 
                                    ncol = 1, nrow = length(colnames(Differentiation.ls))))
    colnames(Pseudotime) <- paste0("patient_id", "_", dd)
    Weights <- data.frame(matrix(data = 1, ncol = 1, nrow = length(colnames(Differentiation.ls))))
    colnames(Weights) <- paste0("patient_id", "_", dd)
    rownames(Weights) = colnames(Differentiation.ls)
    sce <- tradeSeq::fitGAM(counts = as.matrix(Counts), pseudotime = as.matrix(Pseudotime), cellWeights = as.matrix(Weights), nknots = 5, verbose = TRUE,
                               conditions = factor(Differentiation.ls@meta.data[[compare.group]]), 
                            parallel=TRUE, BPPARAM = BPPARAM );gc()# average 1 hr per subject
    
    rowData(sce)$assocRes <- associationTest(sce, lineages = TRUE, l2fc = log2(2))
    AssocRes <- as.data.frame(rowData(sce)$assocRes)
    AssocRes$FDR <- p.adjust(AssocRes$pvalue, "fdr")
    
    condRes <- as.data.frame(conditionTest(sce, l2fc = log2(2)))
    condRes$FDR <- p.adjust(condRes$pvalue, "fdr")
    mean(condRes$FDR <= 0.05, na.rm = TRUE)
    
    # PatRes <- as.data.frame(patternTest(sce, l2fc = log2(2))) #  fcMedian is the median fold change between all those points by points comparison
    # PatRes$FDR <- p.adjust(PatRes$pvalue, method="fdr")
    
    saveRDS(sce, file = paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder,  "/", dd, ".sce.Rds" ) )
    saveRDS(AssocRes, file = paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder,  "/", dd, ".AssocRes.Rds" ) )
    saveRDS(condRes, file = paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder,  "/", dd, ".condRes.Rds" ) )
    # saveRDS( PatRes, file = paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder,  "/", dd, ".PatRes.Rds" ) )
    rm(sce, AssocRes, condRes , PatRes);gc()
  }
} 

finished.names <- list.files( paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder, "/"))[grep(paste0(".sce.Rds"), 
                                                                                                           list.files( paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder, "/")))];
finished.names;
need.names <- setdiff(ls.names , stringr::str_split_fixed(finished.names, "[.]", n=2)[,1] );need.names;
rm(Counts, Differentiation, Differentiation.ls, Pseudotime, Weights, dd, 
   ls.names, ff, BPPARAM)
gc()

if (length(need.names) == 0 ) {
# Pat_ls <- list();
AssocRes <- list(); condRes<- list(); sce_ls <- list()
for (ff in stringr::str_split_fixed(finished.names, "[.]", n=2)[,1]  ) {
  sce_ls[[ff]] <- readRDS(file = paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder,  "/", ff, ".sce.Rds" ) )
  AssocRes[[ff]] <- readRDS(file = paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder,  "/", ff, ".AssocRes.Rds" ) )
  condRes[[ff]] <- readRDS(file = paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder,  "/", ff, ".condRes.Rds" ) )
  # PatRes[[ff]] <- readRDS(file = paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder,  "/", ff, ".PatRes.Rds" ) )
}
gc()
save(AssocRes, condRes, # PatRes, 
     file = paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder,  "/", unlist(stringr::str_split(seurat.file, ".Rds"))[1], "_Pt.DEGs.RData" ))
saveRDS(sce_ls, file = paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder,  "/", unlist(stringr::str_split(seurat.file, ".Rds"))[1], ".sce.Rds" ) )
}

if (file.exists(paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder,  "/", unlist(stringr::str_split(seurat.file, ".Rds"))[1], ".sce.Rds" ) ) == T &
    file.exists(paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder,  "/", unlist(stringr::str_split(seurat.file, ".Rds"))[1], "_Pt.DEGs.RData" ) ) == T )  {
  for (ff in stringr::str_split_fixed(finished.names, "[.]", n=2)[,1] ) {
    file.remove( paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder,  "/", ff, ".sce.Rds") )
    file.remove( paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder,  "/", ff, ".AssocRes.Rds") )
    file.remove( paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder,  "/", ff, ".condRes.Rds") )
    # file.remove( paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder,  "/", ff, ".PatRes.Rds") )
  }
}

rm( sce_ls,AssocRes, condRes )# PatRes,
rm(dir, finished.names,
   need.names, ff)
gc()
sessionInfo()