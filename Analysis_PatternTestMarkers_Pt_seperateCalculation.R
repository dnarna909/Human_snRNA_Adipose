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
compare.group 
groups

# extract each sample file name
if (file.exists(paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder,  "/", "Differentiation.ls.names.Rds" ) ) == F )  {
  Differentiation <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/",seurat.file) )
  saveRDS(unique(Differentiation@meta.data$Dataset), paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder,  "/", "Differentiation.ls.names.Rds" ))
} else {
  ls.names = readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder,  "/", "Differentiation.ls.names.Rds") )
}

# calculate the unfinshed sample file name
finished.names <- list.files( paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder, "/"))[grep(paste0(".sce_ls.Rds"), 
                                                                                          list.files( paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder, "/")))];
finished.names;
need.names <- setdiff(ls.names , finished.names);need.names;

if (length(need.names) > 0 ) {
  if (exists("Differentiation") == F){
    Differentiation <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/",seurat.file) )
  } 
  for (dd in need.names) {
    Differentiation.ls <- subset(Differentiation, subset = Dataset == dd)
    Pseudotime <- data.frame(curve1 = Differentiation.ls@meta.data[[paste0(path, "_Pt")]])
    Weights <- data.frame(matrix(data = 0, ncol = length(groups), nrow = length(colnames(Differentiation.ls))))
    colnames(Weights) <- groups
    rownames(Weights) = colnames(Differentiation.ls)
    for (ii in groups) {
      Weights[ rownames(Weights) %in% rownames(Differentiation.ls@meta.data[Differentiation.ls@meta.data[[Compare.G]] == groups[ii],]), groups[ii]] <- 1
    }
    Counts <- Differentiation.ls@assays$originalexp@counts
    # Counts <- Counts[ rownames(Counts) %in% expressed_genes,]
    sce_ls <- tradeSeq::fitGAM(counts = Counts, pseudotime = Pseudotime, cellWeights = Weights, nknots = 4, verbose = TRUE)
    saveRDS(sce_ls, file = paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder,  "/", dd, "_sce_ls.Rds" ) )
    rm(sce_ls)
    
  }
  }
}

  for (dd in unique(Differentiation@meta.data$Dataset)) {
if (file.exists(paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder,  "/", dd, "Differentiation.ls.Rds" ) ) == F )  {
    # dd= unique(Differentiation@meta.data$Dataset)[1];dd
    Differentiation.ls <- subset(Differentiation, subset = Dataset == dd)
    
  }
  saveRDS(Differentiation.ls, file = paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder,  "/", "Differentiation.ls.Rds" ) )
  if(exists("Differentiation.ls")){ rm(Differentiation)};gc()
  
} else {
  Differentiation.ls <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder,  "/", "Differentiation.ls.Rds") )
}


for (dd in names(Differentiation.ls)) {
  if (file.exists(paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder,  "/", dd, "_sce_ls.Rds" ) ) == F )  {
    Pseudotime <- data.frame(curve1 = Differentiation.ls[[dd]]@meta.data[[paste0(path, "_Pt")]])
    Weights <- data.frame(matrix(data = 0, ncol = length(groups), nrow = length(colnames(Differentiation.ls[[dd]]))))
    colnames(Weights) <- groups
    rownames(Weights) = colnames(Differentiation.ls[[dd]])
    for (ii in groups) {
      Weights[ rownames(Weights) %in% rownames(Differentiation.ls[[dd]]@meta.data[Differentiation.ls[[dd]]@meta.data[[Compare.G]] == groups[ii],]), groups[ii]] <- 1
    }
    Counts <- Differentiation.ls[[dd]]@assays$originalexp@counts
    # Counts <- Counts[ rownames(Counts) %in% expressed_genes,]
    sce_ls <- tradeSeq::fitGAM(counts = Counts, pseudotime = Pseudotime, cellWeights = Weights, nknots = 4, verbose = TRUE)
    saveRDS(sce_ls, file = paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder,  "/", dd, "_sce_ls.Rds" ) )
    rm(sce_ls)
  }
  gc()
}

if (file.exists(paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder,  "/", "sce_ls.Rds" ) ) == F )  {
  
  
  sce_ls <- list()
  sce_ls[[dd]] <- 
    
    if(exists("sce_ls")){ rm(Differentiation.ls)};gc()
} else {
  sce_ls <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder,  "/", "sce_ls.Rds") )
  if(length(names(sce_ls)) < length(names(Differentiation.ls))){
    new.dd <- setdiff(names(Differentiation.ls), names(sce_ls))
    for (dd in new.dd ) {
      Pseudotime <- data.frame(curve1 = Differentiation.ls[[dd]]@meta.data[[paste0(path, "_Pt")]])
      Weights <- data.frame(matrix(data = 0, ncol = length(groups), nrow = length(colnames(Differentiation.ls[[dd]]))))
      colnames(Weights) <- groups
      rownames(Weights) = colnames(Differentiation.ls[[dd]])
      for (ii in groups) {
        Weights[ rownames(Weights) %in% rownames(Differentiation.ls[[dd]]@meta.data[Differentiation.ls[[dd]]@meta.data[[Compare.G]] == groups[ii],]), groups[ii]] <- 1
      }
      Counts <- Differentiation.ls[[dd]]@assays$originalexp@counts
      # Counts <- Counts[ rownames(Counts) %in% expressed_genes,]
      sce_ls[[dd]] <- tradeSeq::fitGAM(counts = Counts, pseudotime = Pseudotime, cellWeights = Weights, nknots = 4, verbose = TRUE)
    }
    saveRDS(sce_ls, file = paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder,  "/", "sce_ls.Rds" ) )
    gc()
  }
}

Pat_ls <- list()
for (dd in names(sce_ls)) {
  Pat_ls[[dd]] <- as.data.frame(patternTest(sce_ls[[dd]]))
  Pat_ls[[dd]]$FDR <- p.adjust(Pat_ls[[dd]]$pvalue, method="fdr")
}
gc()
saveRDS(Pat_ls, file = paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder,  "/", unlist(stringr::str_split(seurat.file, ".Rds"))[1], "_Pt.DEGs.Rds" ))

if (file.exists(paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder,  "/", unlist(stringr::str_split(seurat.file, ".Rds"))[1], "_Pt.DEGs.Rds" ) ) == T )  {
  file.remove( paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder,  "/", "sce_ls.Rds") )
}

rm(Pat_ls, sce_ls)
gc()
