###' DISCLAIMER
#' Some of the algorithms are non-deterministic making the results slightly different from run to run.
#' Many of the algorithms are frequently updated with new releases making the results slightly different from version to version.
#' Depending on your system, this code below may not produce an exact 1:1 copy of the results. 
#' This adversely affects clustering, filtering etc. if the code is naively copy-pasted.
#' For better reproducibility, please download quality-filtered Rds objects and work on downstream analyses. 
#' change "D:/" to "/media/jianie/Extreme SSD1/" # if in LUNIX computer

##' Load libraries-------------------------------------------------------------------------------------------------------------
# rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memory and report the memory usage.
memory.limit(size = 1e+13) # 'memory.limit()' is Windows-specific
if (serve == "GCCRI"){
  source(paste0(Disk, Project.folder, "/", "Project Parameters_GCCRI.R"), local = knitr::knit_global())
} else {
  source(paste0(Disk, Project.folder, "/", "Project Parameters.R"), local = knitr::knit_global())
}
# https://satijalab.org/seurat/articles/visualization_vignette.html
library(Seurat)
# library(SeuratData)
library(ggplot2)
library(patchwork)

# import packages
# library("velocyto.R") # linux -specific

# Import file name -----------------------------------------------------------------------------------------------
# metadata.new <- readRDS(paste0(Disk, Project.folder, "/",  Rds.folder, "/", "SAT_Sample.meta.data.Rds")) ## This object is downloadable from Open Science Framework, with further annotations as subsequent scripts add to the file

# Import data and subset ---------------------------------------------------------------------------------------------------
export.folder <- "Features.plotting"
dir.create(file.path(paste0(Disk, Project.folder, "/", figures.folder, "/"), export.folder), showWarnings = FALSE)

# Import data and subset ---------------------------------------------------------------------------------------------------------------------------------
for (seurat.file in seurat.files) {
  # seurat.file = seurat.files[1]
  file.name <- stringr::str_split(seurat.file, ".Rds")[[1]][1]
  ident1 = ifelse(length(stringr::str_split( file.name, "_")[[1]]) ==4, "Subtype",
                  ifelse(length(stringr::str_split( file.name, "_")[[1]]) ==3, "Annotation"))
  seurat.data <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/", seurat.file))
  seurat.data <- SetIdent(seurat.data, value = ident1)
  
  if(exists("split.groups")) {
  seurat.data@meta.data$orig.ident <- seurat.data@meta.data[[split.groups]]
  seurat.data <- subset(seurat.data, subset = orig.ident %in% select.idents)
  }
  gc() #free up memory and report the memory usage.
  
  features = features[features %in% rownames(seurat.data) ] 
  
  if (length(features) >0){
  export.filename <- ifelse(length(features) >2, paste0(paste(features[1:2], collapse = "."), "_", length(features) ),
                          paste(features, collapse = ".") )

  gc() #free up memory and report the memory usage.
  dir.create(file.path(paste0(Disk, Project.folder, "/", figures.folder, "/"), export.folder, "/", export.filename), showWarnings = FALSE)
  dir <- paste0(Disk, Project.folder, "/", figures.folder, "/", export.folder, "/", export.filename, "/")
  
  # Import metadata
  sample.meta <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/", stringr::str_split(file.name, "_")[[1]][1], "_Sample.meta.data.Rds")) %>% tibble::rownames_to_column(var = "rowname")
  Subtype.meta <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/", file.name, "_", file.type)) %>% tibble::rownames_to_column(var = "rowname")
  metadata.all <- left_join(Subtype.meta, sample.meta, by = "rowname") %>% tibble::column_to_rownames(var = "rowname")
  seurat.data <- AddMetaData( object = seurat.data, metadata = metadata.all)
  seurat.data <- SetIdent(seurat.data, value = ident1)
  
  if (exists("Select.Group")){
    seurat.data <- SetIdent(seurat.data, value = Select.Group)
    seurat.data <- subset(seurat.data, idents = TRUE)
  }
  gc() #free up memory and report the memory usage.
  
  if(exists("split.groups")) {
  table(seurat.data@meta.data[[split.groups]])
  levels(seurat.data@meta.data[[split.groups]])
  seurat.data@meta.data[[split.groups]] <- factor(seurat.data@meta.data[[split.groups]]  , 
                                      levels = select.idents)
  levels(seurat.data@meta.data[[split.groups]])
  table(seurat.data@meta.data[[split.groups]])
  }
  
  # put into seurat.data
  print(seurat.file)
  g <-DimPlot(seurat.data, group.by = ident1, raster = F)
  print(g)
  
  if(exists("split.groups")) {
  g <-DimPlot(seurat.data, group.by = split.groups, raster = F)
  print(g)
  }
  
  # Five visualizations of marker feature expression ------------
  ## Ridge plots -------------
  #' from ggridges. Visualize single cell expression distributions in each cluster 
  g <-RidgePlot(seurat.data, features = features, ncol = 2, group.by = ident1)
  print(g)
  png(file=paste0(dir, file.name, "_", ident1, "_", export.filename,  ".RidgePlot.png"), 
      width= ifelse(length(features) >=2, 10, 5), 
      height=ifelse(length(features) >2, (ceiling(length(features)/2))*((length(unique(seurat.data@meta.data[[ident1]])))*0.1 + 4.4), 
                    ((length(unique(seurat.data@meta.data[[ident1]])))*0.1 + 4.4)), res = 300, units = "in")
  print(g)
  dev.off()
  
  if(exists("split.groups")) {
    g <-RidgePlot(seurat.data, features = features, ncol = 2, group.by = split.groups)
  print(g)
  png(file=paste0(dir, file.name, "_", split.groups, "_", export.filename,  ".RidgePlot.png"), 
      width= ifelse(length(features) >=2, 10, 5), 
      height=ifelse(length(features) >2, (ceiling(length(features)/2))*((length(unique(seurat.data@meta.data[[split.groups]])))*0.1 + 4.4), 
                    ((length(unique(seurat.data@meta.data[[split.groups]])))*0.1 + 4.4)), res = 300, units = "in")
  print(g)
  dev.off()
  }
  
  ## Violin plot -------------------------- 
  #' Visualize single cell expression distributions in each cluster 
  g <-VlnPlot(seurat.data, features = features, group.by = ident1, raster=FALSE,
              pt.size = 0, y.max = max(seurat.data@assays$originalexp@data[ features ,], na.rm = TRUE) *1.05,
              )+ 
    # stat_compare_means(aes(group= seurat.data@meta.data[[ident1]]), label = "p.signif", size =10)+ 
    theme(plot.title = element_blank(), 
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12)
    ) 
  print(g)
  png(file=paste0(dir, file.name, "_", ident1, "_", export.filename, ".VlnPlot.png"), 
      width=  ifelse(length(features) >=2, 2*((length(unique(seurat.data@meta.data[[ident1]])))*0.1 + 4.4), ((length(unique(seurat.data@meta.data[[ident1]])))*0.1 + 4.4) ), 
      height= ifelse(length(features) > 2, (ceiling(length(features)/2))*((length(unique(seurat.data@meta.data[[ident1]])))*0.1 + 4.4), 
                     ((length(unique(seurat.data@meta.data[[ident1]])))*0.1 + 4.4) ) , res = 300, units = "in")
  print(g)
  dev.off()
  
  if(exists("split.groups")) {
    g <-VlnPlot(seurat.data, features = features, group.by = split.groups)
  print(g)
  png(file=paste0(dir, file.name, "_", split.groups, "_", export.filename, ".VlnPlot.png"), 
      width=  ifelse(length(features) >=2, 2*((length(unique(seurat.data@meta.data[[split.groups]])))*0.1 + 4.4), ((length(unique(seurat.data@meta.data[[split.groups]])))*0.1 + 4.4) ), 
      height= ifelse(length(features) > 2, (ceiling(length(features)/2))*((length(unique(seurat.data@meta.data[[split.groups]])))*0.1 + 4.4), 
                     ((length(unique(seurat.data@meta.data[[split.groups]])))*0.1 + 4.4) ) , res = 300, units = "in")
  print(g)
  dev.off()
  }

  ## Feature plot ------------- 
  #' visualize feature expression in low-dimensional space
  g <-FeaturePlot(seurat.data, features = features, raster = F, ncol = 2)
  print(g)
  png(file=paste0(dir, file.name, "_",  export.filename, ".FeaturePlot.png"), 
      width= ifelse(length(features) >=2, 10, 5), 
      height=ifelse(length(features) >=2, (ceiling(length(features)/2))*5, 5), res = 300, units = "in")
  print(g)
  dev.off()
  
  if(exists("split.groups")) {
  g <-FeaturePlot(seurat.data, features = features, raster = F, ncol = 2, split.by = split.groups )
  print(g)
  png(file=paste0(dir, file.name, "_", split.groups, "_" ,export.filename, ".FeaturePlot.png"), 
      width= ((length(unique(seurat.data@meta.data[[split.groups]])))*5 + 1), 
      height=length(features)*5, res = 300, units = "in")
  print(g)
  dev.off()
  }
  
  FeaturePlot(seurat.data, features = features, min.cutoff = "q10", max.cutoff = "q90", raster = F)
  
  ## Dot plots -----------------------
  #' the size of the dot corresponds to the percentage of cells expressing the
  #' feature in each cluster. The color represents the average expression level
  g <-DotPlot(seurat.data, features = features, group.by = ident1) + RotatedAxis()
  print(g)
  png(file=paste0(dir, file.name, "_", export.filename, ".DotPlot.png"), 
      width=  ((length(features))*0.2 + 4.4) , 
      height= ((length(unique(seurat.data@meta.data[[ident1]])))*0.1 + 4.4) , res = 300, units = "in")
  print(g)
  dev.off()
  
  ## Single cell heatmap of feature expression ---------------
  # DoHeatmap(subset(seurat.data, downsample = 100), features = features, size = 3)
  # DoHeatmap(seurat.data, features = VariableFeatures(seurat.data)[1:100], cells = 1:500, size = 4,angle = 90) + NoLegend()
  
  if(exists("split.groups")) {
    g <-DimPlot(seurat.data, group.by = ident1, split.by = split.groups, raster = F)
    print(g)
    
    ## split by groups ----------------
    # VlnPlot(seurat.data, features = "DoubletScore", split.by = split.groups)
    # VlnPlot(seurat.data, features = "subsets_Mito_percent", split.by = split.groups)
    #' Split visualization to view expression by groups (replaces FeatureHeatmap)
    seurat.data@meta.data$new.ident <- paste0(seurat.data@meta.data[[ident1]], "_", seurat.data@meta.data[[split.groups]])
    seurat.data <- SetIdent(seurat.data, value = "new.ident")
    
    g <-DotPlot(seurat.data, features = features, group.by = "new.ident") + RotatedAxis() + theme_bw()
    print(g)
    png(file=paste0(dir, file.name, "_", split.groups,"_",  ident1, "_",  export.filename,  ".DotPlot.png"), 
        width= ((length(features))*0.5 + 4.4) , 
        height= ((length(unique(seurat.data@meta.data[[ident1]])))*(length(unique(seurat.data@meta.data[[ident1]])))*0.03 + 2.4), res = 300, units = "in")
    print(g)
    dev.off()
    
    # cols = grDevices::colorRampPalette(brewer.pal(12, "Paired"))(length(unique(seurat.data@meta.data[["Pre_Post"]])))
    # seurat.data <- SetIdent(seurat.data, value = "Group")
    # DotPlot(seurat.data, features = features, split.by = "Pre_Post") + RotatedAxis() # , cols = cols
  }
  
  
  for (ff in features) {
    g <-  FeaturePlot(seurat.data, features = ff, min.cutoff = 1, max.cutoff = 3, raster=FALSE)
    print(g)
    png(file=paste0(dir, file.name, "_", ff, ".FeaturePlot.png"), 
        width= 5, 
        height=5, res = 300, units = "in")
    print(g)
    dev.off()
    
    ## Scatter plot of two features simultaneously 
    g <- FeatureScatter(seurat.data, feature1 = "Age", feature2 = ff, raster=FALSE, group.by = ident1)
    print(g)
    png(file=paste0(dir, file.name, "_", ident1, "_",  ff, ".FeatureScatter.Age.png"), 
        width= 15, 
        height=15, res = 300, units = "in")
    print(g)
    dev.off()
    
    if(exists("split.groups")) {
    g <- FeatureScatter(seurat.data, feature1 = "Age", feature2 = ff, raster=FALSE, group.by = split.groups)
    print(g)
    png(file=paste0(dir, file.name, "_", split.groups, "_",  ff, ".FeatureScatter.Age.png"), 
        width= 15, 
        height=15, res = 300, units = "in")
    print(g)
    dev.off()
    }
  }
  
  
  if(length(features) == 2){
    ## Visualize co-expression of two features simultaneously ---------------
    g <- FeaturePlot(seurat.data, features = c(features[1], features[2]), blend = TRUE)
    print(g)
    png(file=paste0(dir, file.name, "_", export.filename, ".FeaturePlot2.png"), 
        width= 20, 
        height=5, res = 300, units = "in")
    print(g)
    dev.off()
    if(exists("split.groups")) {
    ## Scatter plot of two features simultaneously 
    g <- FeatureScatter(seurat.data, feature1 = features[1], feature2 = features[2], raster=FALSE, group.by = split.groups)
    print(g)
    png(file=paste0(dir, file.name, "_",  split.groups, "_", export.filename, ".FeatureScatter.png"), 
        width= 15, 
        height=15, res = 300, units = "in")
    print(g)
    dev.off()
    }
    
    g <- FeatureScatter(seurat.data, feature1 = features[1], feature2 = features[2], raster=FALSE, group.by = ident1)
    print(g)
    png(file=paste0(dir, file.name, "_",  ident1, "_", export.filename, ".FeatureScatter.png"), 
        width= 15, 
        height=15, res = 300, units = "in")
    print(g)
    dev.off()
    
  }
  
  
  # DimPlot(seurat.data, group.by = "Label", split.by = "sample_id", raster = F)
  rm(seurat.data)
  gc() #free up memory and report the memory usage.
}
rm(g, metadata.all, sample.meta,
   Subtype.meta
   )

}
gc()

sessionInfo()
