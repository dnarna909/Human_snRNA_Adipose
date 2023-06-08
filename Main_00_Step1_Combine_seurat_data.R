#' DISCLAIMER
#' Some of the algorithms are non-deterministic making the results slightly different from run to run.
#' Many of the algorithms are frequently updated with new releases making the results slightly different from version to version.
#' Depending on your system, this code below may not produce an exact 1:1 copy of the results. 
#' This adversely affects clustering, filtering etc. if the code is naively copy-pasted.
#' For better reproducibility, please download quality-filtered Rds objects and work on downstream analyses. 
#' change "D:/" to "/media/jianie/Extreme SSD1/" # if in LUNIX computer

# Load libraries -------------------------------------------------------------------------------------------------------------
# rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.
memory.limit(size = 1e+13) # 'memory.limit()' is Windows-specific

# set parameters
source(paste0(Disk, Project.folder, "/","Project Parameters.R"), local = knitr::knit_global())


## 2. Merge all the datasets -------------------------------------------------------------------------------------------
##' open files
sample.files = list.files(paste0(Disk, Project.folder, "/",  Rds.folder, "/seurat.list/"))

## merge
eWAT <- merge(x=readRDS(paste0(Disk, Project.folder, "/",  Rds.folder, "/seurat.list/", sample.files[1] )), 
              y= readRDS(paste0(Disk, Project.folder, "/",  Rds.folder, "/seurat.list/", sample.files[2] )))
saveRDS(eWAT, file = paste0(Disk, Project.folder, "/",  Rds.folder, "/SAT1.Rds"))

gc() #free up memrory and report the memory usage.
for (ii in sample.files[3:28])  {
  eWAT <- merge(x= readRDS(paste0(Disk, Project.folder, "/",  Rds.folder, "/SAT1.Rds")) , 
                y= readRDS(paste0(Disk, Project.folder, "/",  Rds.folder, "/seurat.list/", ii) ))
  saveRDS(eWAT, file = paste0(Disk, Project.folder, "/",  Rds.folder, "/SAT1.Rds"))
  rm(eWAT)
  print(ii)
  gc() #free up memrory and report the memory usage.
}
# crashed because of too big file, can only process 28 samples in Dell pc

eWAT <- merge(x=readRDS(paste0(Disk, Project.folder, "/",  Rds.folder, "/seurat.list/", sample.files[29] )), 
              y= readRDS(paste0(Disk, Project.folder, "/",  Rds.folder, "/seurat.list/", sample.files[30] )))
saveRDS(eWAT, file = paste0(Disk, Project.folder, "/",  Rds.folder, "/SAT2.Rds"))
for (ii in sample.files[31:length(sample.files)])  {
  eWAT <- merge(x= readRDS(paste0(Disk, Project.folder, "/",  Rds.folder, "/SAT2.Rds")) , 
                y= readRDS(paste0(Disk, Project.folder, "/",  Rds.folder, "/seurat.list/", ii) ))
  saveRDS(eWAT, file = paste0(Disk, Project.folder, "/",  Rds.folder, "/SAT2.Rds"))
  rm(eWAT)
  print(ii)
  gc() #free up memrory and report the memory usage.
}

gc() #free up memrory and report the memory usage.
eWAT <- merge(x=readRDS(paste0(Disk, Project.folder, "/",  Rds.folder, "/", "SAT1.Rds" )), 
              y= readRDS(paste0(Disk, Project.folder, "/",  Rds.folder, "/", "SAT2.Rds" )))
# eWAT <- merge(x=seurat.list[[1]] , y= seurat.list[2:length(seurat.list)])# long time step-1.5 hr
saveRDS(eWAT, file = paste0(Disk, Project.folder, "/",  Rds.folder, "/SAT.Rds"))


# 3. Merge and seperate STARR and SGLT2 --------------------------------------------------------------------------
## merge STARR
eWAT <- merge(x=readRDS(paste0(Disk, Project.folder, "/",  Rds.folder, "/seurat.list/", sample.files[1] )), 
              y= readRDS(paste0(Disk, Project.folder, "/",  Rds.folder, "/seurat.list/", sample.files[2] )))
saveRDS(eWAT, file = paste0(Disk, Project.folder, "/",  Rds.folder, "/SAT1.Rds"))

gc() #free up memrory and report the memory usage.
for (ii in sample.files[3:18])  {
  eWAT <- merge(x= readRDS(paste0(Disk, Project.folder, "/",  Rds.folder, "/SAT1.Rds")) , 
                y= readRDS(paste0(Disk, Project.folder, "/",  Rds.folder, "/seurat.list/", ii) ))
  saveRDS(eWAT, file = paste0(Disk, Project.folder, "/",  Rds.folder, "/SAT1.Rds"))
  rm(eWAT)
  print(ii)
  gc() #free up memrory and report the memory usage.
}

## merge SGLT2
gc() #free up memrory and report the memory usage.
eWAT <- merge(x=readRDS(paste0(Disk, Project.folder, "/",  Rds.folder, "/seurat.list/", sample.files[19] )), 
              y= readRDS(paste0(Disk, Project.folder, "/",  Rds.folder, "/seurat.list/", sample.files[20] )))
saveRDS(eWAT, file = paste0(Disk, Project.folder, "/",  Rds.folder, "/SAT2.Rds"))
for (ii in sample.files[21:length(sample.files)])  {
  eWAT <- merge(x= readRDS(paste0(Disk, Project.folder, "/",  Rds.folder, "/SAT2.Rds")) , 
                y= readRDS(paste0(Disk, Project.folder, "/",  Rds.folder, "/seurat.list/", ii) ))
  saveRDS(eWAT, file = paste0(Disk, Project.folder, "/",  Rds.folder, "/SAT2.Rds"))
  rm(eWAT)
  print(ii)
  gc() #free up memrory and report the memory usage.
}

## combine all
gc() #free up memrory and report the memory usage.
eWAT <- merge(x=readRDS(paste0(Disk, Project.folder, "/",  Rds.folder, "/", "SAT1.Rds" )), 
              y= readRDS(paste0(Disk, Project.folder, "/",  Rds.folder, "/", "SAT2.Rds" )))
# eWAT <- merge(x=seurat.list[[1]] , y= seurat.list[2:length(seurat.list)])# long time step-1.5 hr
saveRDS(eWAT, file = paste0(Disk, Project.folder, "/",  Rds.folder, "/SAT.Rds"))


