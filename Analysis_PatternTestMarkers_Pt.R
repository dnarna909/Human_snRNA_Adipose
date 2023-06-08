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
export.folder <- paste0("PatternMarkers")
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

if (file.exists(paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder,  "/", compare.group, ".sce_ls.Rds" ) ) == F )  {
  Differentiation <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/",seurat.file) )
  Differentiation <- Seurat::SetIdent(Differentiation, value = Select.Group)
  Differentiation <- subset(Differentiation, idents = TRUE);gc()
  Pseudotime <- data.frame(matrix(data = Differentiation@meta.data[[paste0(path, "_Pt")]], 
                                  ncol = length(groups), nrow = length(colnames(Differentiation))))
  Weights <- data.frame(matrix(data = 0, ncol = length(groups), nrow = length(colnames(Differentiation))))
  colnames(Weights) <- groups
  rownames(Weights) = colnames(Differentiation)
  for (ii in groups) {
    Weights[ rownames(Weights) %in% rownames(Differentiation@meta.data[Differentiation@meta.data[[compare.group]] == ii,]), ii] <- 1
  }
  Counts <- Differentiation@assays$originalexp@counts
  # Counts <- Counts[ rownames(Counts) %in% expressed_genes,]
  gc()
  sce_ls <- tradeSeq::fitGAM(counts = Counts, pseudotime = Pseudotime, cellWeights = Weights, nknots = 4, verbose = TRUE); gc()
  saveRDS(sce_ls, file = paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder,  "/", compare.group, ".sce_ls.Rds" ) )
  gc()
  if(exists("sce_ls")){ rm(Differentiation)};gc()
} else {
  sce_ls <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder,  "/", compare.group, ".sce_ls.Rds") )
}
gc()

Pat_ls <- as.data.frame(patternTest(sce_ls))
Pat_ls$FDR <- p.adjust(Pat_ls$pvalue, method="fdr")
gc()
saveRDS(Pat_ls, file = paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder,  "/", 
                              unlist(stringr::str_split(seurat.file, ".Rds"))[1], "_", compare.group, "_Pt.DEGs.Rds" ))

if (file.exists(paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder,  "/", 
                       unlist(stringr::str_split(seurat.file, ".Rds"))[1], "_", compare.group, "_Pt.DEGs.Rds" ) ) == T )  {
  file.remove( paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder,  "/", compare.group, ".sce_ls.Rds") )
}

rm(Pat_ls, sce_ls)
gc()
