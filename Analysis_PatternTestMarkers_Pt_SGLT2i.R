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
export.folder <- paste0("PatternMarkers_SGLT2i")
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
finished.names <- list.files( paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder, "/"))[grep(paste0(".sce_ls.Rds"), 
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
    Differentiation.ls <- subset(Differentiation, subset = patient_id == dd);gc()
    Pseudotime <- data.frame(matrix(data = Differentiation.ls@meta.data[[paste0(path, "_Pt")]], 
                                    ncol = length(groups), nrow = length(colnames(Differentiation.ls))))
    colnames(Pseudotime) <- groups
    Weights <- data.frame(matrix(data = 0, ncol = length(groups), nrow = length(colnames(Differentiation.ls))))
    colnames(Weights) <- groups
    rownames(Weights) = colnames(Differentiation.ls)
    for (ii in groups) {
      Weights[ rownames(Weights) %in% rownames(Differentiation.ls@meta.data[Differentiation.ls@meta.data[[compare.group]] == ii,]), ii] <- 1
    }
    Counts <- Differentiation.ls@assays$originalexp@counts; dim(Counts)
    # Counts <- Counts[ rownames(Counts) %in% expressed_genes,]
    sce_ls <- tradeSeq::fitGAM(counts = Counts, pseudotime = Pseudotime, cellWeights = Weights, nknots = 4, verbose = TRUE);gc()# average 1 hr per subject
    saveRDS(sce_ls, file = paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder,  "/", dd, ".sce_ls.Rds" ) )
    rm(sce_ls);gc()
  }
} 

finished.names <- list.files( paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder, "/"))[grep(paste0(".sce_ls.Rds"), 
                                                                                                           list.files( paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder, "/")))];
finished.names;
need.names <- setdiff(ls.names , stringr::str_split_fixed(finished.names, "[.]", n=2)[,1] );need.names;

if (length(need.names) == 0 ) {
Pat_ls <- list()
for (ff in stringr::str_split_fixed(finished.names, "[.]", n=2)[,1]  ) {
  sce_ls <- readRDS(file = paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder,  "/", ff ) )
  Pat_ls[[ff]] <- as.data.frame(patternTest(sce_ls)) #  fcMedian is the median fold change between all those points by points comparison
  Pat_ls[[ff]]$FDR <- p.adjust(Pat_ls[[ff]]$pvalue, method="fdr")
}
gc()
saveRDS(Pat_ls, file = paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder,  "/", unlist(stringr::str_split(seurat.file, ".Rds"))[1], "_Pt.DEGs.Rds" )) 
}

if (file.exists(paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder,  "/", unlist(stringr::str_split(seurat.file, ".Rds"))[1], "_Pt.DEGs.Rds" ) ) == T )  {
  for (ff in finished.names ) {
    file.remove( paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder,  "/", ff) )
  }
}

rm(Pat_ls, sce_ls)
rm(Counts, Differentiation, Differentiation.ls, Pseudotime, Weights, dd, dir, finished.names,
   ii, ls.names, need.names, ff, sce_ls, Pat_ls)
gc()
sessionInfo()