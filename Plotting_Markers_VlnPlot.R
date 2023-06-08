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
source(paste0(Disk, Project.folder, "/", "Project Parameters.R"), local = knitr::knit_global())

# vlnplot for all markers
# https://github.com/ycl6/StackedVlnPlot
library(Seurat)
library(ggplot2)
library(cowplot)

# prepare folder
export.folder <- paste0("Markers")
dir.create(file.path(paste0(Disk, Project.folder, "/", figures.folder, "/"), export.folder), showWarnings = FALSE)
dir <- paste0(Disk, Project.folder, "/", figures.folder, "/", export.folder, "/");dir


for (seurat.file in seurat.files) {
  # markers.file = markers.files[1]; markers.file
  # seurat.file = seurat.files[1]
  
  # Load data.frame obj
  seurat.data <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/", seurat.file) )
  file.name <- stringr::str_split(seurat.file, ".Rds")[[1]][1]
  Annotation.file.name <- paste(stringr::str_split(file.name, "_")[[1]][c(1,2,3)], collapse = "_")
  sample.meta <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/",sample.file.type)) %>% tibble::rownames_to_column(var = "rowname")
  Subtype.meta <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/", Annotation.file.name, "_", file.type)) %>% tibble::rownames_to_column(var = "rowname")
  metadata <- left_join(Subtype.meta, sample.meta, by = "rowname") %>% tibble::column_to_rownames(var = "rowname")
  seurat.data <- AddMetaData(seurat.data, metadata = metadata)
  dim(seurat.data);
  if (exists("Select.Group")){
    seurat.data <- SetIdent(seurat.data, value = Select.Group)
    seurat.data <- subset(seurat.data, idents = TRUE)
  }
  dim(seurat.data)
  all.pbmc = seurat.data@assays[["originalexp"]]@data 
  df <- data.frame(Cell = rownames(seurat.data@meta.data),
                   Idents = seurat.data@meta.data[[group.x]] )
  rm(seurat.data); gc();
  
  # prepare folder
  # export.folder0 <- stringr::str_split(markers.file, "/")[[1]] ;export.folder0
  # dir1 <- paste0(Disk, Project.folder, "/", Rds.folder, "/")
  # for (ee in 1:(length(export.folder0) -1) ) {
  #   dir.create(file.path(dir1, export.folder0[ee]), showWarnings = FALSE)
  #   dir2 <- paste0(dir1, export.folder0[ee], "/")
  # }
  # dir2
  
  features.df <- rio::import(paste0(Disk, Project.folder, "/", Markers_file), sheet = sheet)
  features <-features.df$gene;features
  
  library(dplyr)
  pbmc = as.data.frame(t(as.matrix(all.pbmc[c(features, setdiff(rownames(all.pbmc) , features)[1]), ])))  %>% 
    dplyr::select(-one_of(setdiff(rownames(all.pbmc) , features)[1])) %>% 
    tibble::rownames_to_column(var = "Cell") %>%
    left_join( 
      df, 
      by ="Cell"
    )
  head(pbmc)[,]               
  
  
  if (exists("Compare.G")){
    pbmc.split <-split(pbmc, pbmc[, Compare.G])
    pbmc.split <- pbmc.split[sapply(pbmc.split, nrow)>0]
  } else {
    pbmc.split <- list(All = pbmc)
  }
  
  # Use melt to change data.frame format
  pbmc.melt <- lapply(pbmc.split, function(x) {
    reshape2::melt(x, id.vars = c("Cell","Idents"), measure.vars = features,
                   variable.name = "Feat", value.name = "Expr") %>% arrange(-Expr)
  })
  
  # pbmc <- reshape2::melt(pbmc, id.vars = c("Cell","Idents"), measure.vars = features,
  #                        variable.name = "Feat", value.name = "Expr")
  # head(pbmc, 10)
  
  # Identity on x-axis
  for (mm  in names(pbmc.melt)) {
    # mm  = names(pbmc.melt)[1];mm
    a <- ggplot(pbmc.melt[[mm]], aes(factor(Idents), Expr, fill = Idents)) +
      geom_violin(scale = "width", adjust = 1, trim = TRUE) +
      scale_y_continuous(expand = c(0, 0), position="left", labels = function(x)
        c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
      facet_grid(rows = vars(Feat), scales = "free") + # , switch = "y"
      theme_cowplot(font_size = 12) +
      theme(
        legend.position = "right", 
        legend.title=element_blank(),
        # legend.position = "none", 
        panel.spacing = unit(0, "lines"),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        axis.title.x=element_blank(),
        axis.text.x=element_text(angle = 45, hjust = 1),
        # strip.text.y.left = element_text(angle = 0),
        strip.text.y.right = element_text(angle = 0, hjust = 0)
      ) +
      ggtitle(paste0(sheet, " Markers plotted in ", mm, " Cells")) + 
      # xlab("Identity") + 
      ylab("Expression Level")
    print(a)
    
    png(filename = paste0(dir, 
                          stringr::str_split(seurat.file, ".Rds")[[1]][1], "_",
                          sheet, ".", mm, "_VlnPLot.png"), 
        width = (length(unique(pbmc.melt[[mm]]$Idents)) )+2.59, 
        height= (length(features)/4.48)+1.3429 , res = 300, units = "in")
    print(a)
    dev.off()
  }
  
}




rm(a, all.pbmc, df, metadata, pbmc, pbmc.melt,Annotation.file.name,
   pbmc.split, Result, sample.meta, Subtype.meta, dir, dir1,
   dir2, ee, export.folder, export.folder0, features, file.name,
    mm, n, nn,seurat.file, features.df,
   ff
   )
gc()
