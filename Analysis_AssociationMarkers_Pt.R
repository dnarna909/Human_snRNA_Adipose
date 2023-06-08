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

# source(paste0("/media/jianie/Extreme SSD1/2022-09-01 STARR Combine/Project Parameters.R"), local = knitr::knit_global())
source(paste0(Disk, Project.folder, "/", "Project Parameters_GCCRI.R"), local = knitr::knit_global())
library(tradeSeq)


# prepare folder
export.folder <- paste0("AssocationMarkers")
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
if (file.exists(paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder,  "/", "Differentiation.ls.Rds" ) ) == F )  {
  Differentiation <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/",seurat.file) )
  
  Differentiation.ls <- list()
  for (dd in unique(Differentiation@meta.data$Dataset)) {
    # dd= unique(Differentiation@meta.data$Dataset)[1];dd
    Differentiation.ls[[dd]] <- subset(Differentiation, subset = Dataset == dd)
  }
  saveRDS(Differentiation.ls, file = paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder,  "/", "Differentiation.ls.Rds" ) )
  if(exists("Differentiation.ls")){ rm(Differentiation)};gc()
  
} else {
  Differentiation.ls <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder,  "/", "Differentiation.ls.Rds") )
}

if (file.exists(paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder,  "/", "sce_ls.Rds" ) ) == F )  {
  sce_ls <- list()
  for (dd in names(Differentiation.ls)) {
    Pseudotime <- data.frame(curve1 = Differentiation.ls[[dd]]@meta.data[[paste0(path, "_Pt")]])
    Weights <- data.frame( DD = rep(1, ncol(Differentiation.ls[[dd]])))
    colnames(Weights) <- dd
    rownames(Weights) = colnames(Differentiation.ls[[dd]])
    Counts <- Differentiation.ls[[dd]]@assays$originalexp@counts
    # Counts <- Counts[ rownames(Counts) %in% expressed_genes,]
    sce_ls[[dd]] <- tradeSeq::fitGAM(counts = Counts, pseudotime = Pseudotime, cellWeights = Weights, nknots = 4, verbose = TRUE)
  }
  saveRDS(sce_ls, file = paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder,  "/", "sce_ls.Rds" ) )
  gc()
  if(exists("sce_ls")){ rm(Differentiation.ls)};gc()
} else {
  sce_ls <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder,  "/", "sce_ls.Rds") )
  if(length(names(sce_ls)) < length(names(Differentiation.ls))){
    new.dd <- setdiff(names(Differentiation.ls), names(sce_ls))
    for (dd in new.dd ) {
      Pseudotime <- data.frame(curve1 = Differentiation.ls[[dd]]@meta.data[[paste0(path, "_Pt")]])
      Weights <- data.frame( DD = rep(1, ncol(Differentiation.ls[[dd]])))
      colnames(Weights) <- dd
      rownames(Weights) = colnames(Differentiation.ls[[dd]])
      Counts <- Differentiation.ls[[dd]]@assays$originalexp@counts
      # Counts <- Counts[ rownames(Counts) %in% expressed_genes,]
      sce_ls[[dd]] <- tradeSeq::fitGAM(counts = Counts, pseudotime = Pseudotime, cellWeights = Weights, nknots = 4, verbose = TRUE)
    }
    saveRDS(sce_ls, file = paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder,  "/", "sce_ls.Rds" ) )
    gc()
  }
  
}

Asso_ls <- list()
for (dd in names(sce_ls)) {
  Asso_ls[[dd]] <- as.data.frame(associationTest(sce_ls[[dd]]))
  Asso_ls[[dd]]$FDR <- p.adjust(Asso_ls[[dd]]$pvalue, method="fdr")
}
gc()
saveRDS(Asso_ls, file = paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder,  "/", unlist(stringr::str_split(seurat.file, ".Rds"))[1], "_Pt.DEGs.Rds" ))

if (file.exists(paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder,  "/", unlist(stringr::str_split(seurat.file, ".Rds"))[1], "_Pt.DEGs.Rds" ) ) == T )  {
  file.remove( paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder,  "/", "sce_ls.Rds") )
}

rm(Asso_ls, sce_ls)
gc()
