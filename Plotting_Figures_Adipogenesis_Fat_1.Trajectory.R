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

dir.create(file.path(paste0(Disk, Project.folder, "/", figures.folder, "/",  "Trajectory.phate", "/")), showWarnings = FALSE)
dir0 <- paste0(Disk, Project.folder, "/", figures.folder, "/", "Trajectory.phate",  "/")
dir.create(file.path(paste0(Disk, Project.folder, "/", Rds.folder, "/",  "Trajectory.phate", "/")), showWarnings = FALSE)
dir1 <- paste0(Disk, Project.folder, "/", Rds.folder, "/",  "Trajectory.phate",  "/")

# Import data -----------------------------------------------------------------------------------------------------------------------------------------------------
for (nn in groups) {
  Adipogenesis <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/","Adipogenesis_PHATE_", nn, ".Rds"))
  meta_data_merge <- readRDS(file=paste0(Disk, Project.folder, "/", Rds.folder, "/","Adipogenesis.meta.data.PHATE_", nn, ".Pt.Rds"))
  load( file = paste0(Disk, Project.folder, "/", Rds.folder, "/", "Adipogenesis.PHATE_", nn, ".Branch.RData"))
  
  Adipogenesis <- AddMetaData( object = Adipogenesis, metadata = meta_data_merge)
  
  # plotting  ----------------------------------------------------------------------------------------------
  p <- DimPlot(Adipogenesis, reduction = "phate", group.by = "EdgeID", label = TRUE, raster=FALSE)
  # DimPlot(Adipogenesis, reduction = "phate", group.by = "EdgeID", split.by =  "Age_Group", label = T, raster=FALSE)
  # DimPlot(Adipogenesis, reduction = "phate", group.by = "EdgeID", split.by =  "Status", label = T, raster=FALSE)
  print(p)
  png(file=paste0(dir0, "Adipogenesis.PHATE_", nn, "_",  "EdgeID",  ".DimPlot.png"), 
      width= 7, height=6, res = 300, units = "in")
  print(p)
  dev.off()
  
  p <- DimPlot(Adipogenesis, reduction = "phate", group.by = "EdgeID", label = FALSE, raster=FALSE) +
          annotate("point", x = TreeEPG[[1]]$NodePositions[, "phate_1"], y=TreeEPG[[1]]$NodePositions[, "phate_2"], colour = "blue")
  print(p)
  png(file=paste0(dir0, "Adipogenesis.PHATE_", nn, "_",  "EdgeID.dot",  ".DimPlot.png"), 
      width= 7, height=6, res = 300, units = "in")
  print(p)
  dev.off()
  
  
  p <- DimPlot(Adipogenesis, reduction = "phate", group.by = col.groups, label = FALSE, raster=FALSE) +
          annotate("point", x = TreeEPG[[1]]$NodePositions[, "phate_1"], y=TreeEPG[[1]]$NodePositions[, "phate_2"], colour = "black",  fill="white", size = 3, shape =21) + 
          annotate("text", x = TreeEPG[[1]]$NodePositions[, "phate_1"]*0.9, y=TreeEPG[[1]]$NodePositions[, "phate_2"]*0.9, label = rownames(as.data.frame(TreeEPG[[1]]$NodePositions)) )
  print(p)
  png(file=paste0(dir0, "Adipogenesis.PHATE_", nn, "_",  col.groups, ".", "NodePositions",  ".DimPlot.png"), 
      width= 7, height=6, res = 300, units = "in")
  print(p)
  dev.off()
  
  
  Node.df <- as.data.frame(TreeEPG[[1]]$NodePositions) 
  Node.list <- list()
  for (aa in 1:nrow(allDistancesDf)){
    Branch.name <- names(sort(lengths(Tree_e2e), decreasing = TRUE)[aa])
    Node.list[[Branch.name]] <- na.omit(Node.df[as.numeric(allDistancesDf[Branch.name, ]), ]) %>%
      tibble::rownames_to_column(var = "Node") %>%
      mutate(Branch = Branch.name, Node.Numers = n(), Node.seq = 1:n())
    # Node.list[[Branch.name]][[Branch.name]] <- Branch.name
  }
  Node.list.df <- bind_rows(Node.list)

  library(RColorBrewer)
  library(circlize)
  coul <-colorRampPalette(brewer.pal(9, "Set1"))(length(unique(Node.list.df$Branch)))
  colors = coul[1: length(unique(Node.list.df$Branch))]
  pie(rep(1, length(coul)), col = coul , main="")
  names(colors) <- sort(unique(Node.list.df$Branch))
  Node.list.df <- Node.list.df %>% left_join(as.data.frame(colors) %>% tibble::rownames_to_column(var = "Branch") , by = "Branch")

  Node.list.df2 <- Node.list.df %>%
    group_by(Node) %>%
    dplyr::filter(n() > 1) %>%
    mutate( amount = 1/n(), radius= sqrt(amount)/400)

  my_cols2 <- scales::hue_pal()(length(unique(Adipogenesis@meta.data[[col.groups]])))
  scales::show_col(scales::hue_pal()(length(unique(Adipogenesis@meta.data[[col.groups]]))))
  # my_cols2 <- colorRampPalette(brewer.pal(12, "Paired"))(length(unique(Adipogenesis@meta.data[[col.groups]])))
  # pie(rep(1, length(my_cols2)), col = my_cols2, main="")
  names(my_cols2) <- sort(unique(Adipogenesis@meta.data[[col.groups]]))

  library(ggforce)
  p <-DimPlot(Adipogenesis, reduction = "phate", group.by = col.groups, label = TRUE, raster=FALSE, cols = my_cols2,  repel=TRUE) +
    # geom_line(data = Node.list.df, aes(x = phate_1, y=phate_2, group = Branch), linewidth = 1.5, linetype='dashed')  + # colour = "black", 
    geom_point(data = Node.list.df, aes(x = phate_1, y=phate_2, fill = Branch),colour = "black",  size = 4, shape =21) +# fill="white", 
    geom_arc_bar(data = Node.list.df2,  aes(x0=phate_1,y0=phate_2,r0=0,r=radius,amount=amount,fill=Branch),stat="pie", linewidth = 0.5)+
    coord_equal()+ # so one gets actually circles
    scale_fill_manual(values=c(colors))   +
    labs(fill="Branch", color="Group", group = "Branch")
  print(p)
  png(file=paste0(dir0, "Adipogenesis.PHATE_", nn, "_",  col.groups,  ".DimPlot.png"), 
      width= 12, height=6, res = 300, units = "in")
  print(p)
  dev.off()
  
  
  Adipogenesis <- SetIdent(Adipogenesis, value = col.groups)
  cell.types = unique(Adipogenesis@meta.data[[col.groups]])
  for (cell.type in cell.types){
    # cell.type <-"CA3" 
    hight.cells <- WhichCells(Adipogenesis,  idents =  cell.type )
    p <- DimPlot(Adipogenesis, reduction = "phate", group.by = col.groups, label = FALSE, raster=FALSE, 
                 cols= "grey", cells.highlight = list(hight.cells), cols.highlight = c("red"),  repel=TRUE) +
      # geom_line(data = Node.list.df, aes(x = phate_1, y=phate_2, group = Branch), linewidth = 1.5, linetype='dashed')  + # colour = "black", 
      geom_point(data = Node.list.df, aes(x = phate_1, y=phate_2, fill = Branch),colour = "black",  size = 4, shape =21) +# fill="white", 
      geom_arc_bar(data = Node.list.df2,  aes(x0=phate_1,y0=phate_2,r0=0,r=radius,amount=amount,fill=Branch),stat="pie", linewidth = 0.5)+
      coord_equal()+ # so one gets actually circles
      scale_fill_manual(values=c(colors))   +
      labs(fill="Branch",  group = "Branch", title = cell.type ) +
      guides(fill = guide_legend(order = 1), color = "none")
    print(p)
    png(file=paste0(dir0, "Adipogenesis.PHATE_", nn, "_",  col.groups, "_", cell.type,  ".DimPlot.png"), 
        width= 7, height=6, res = 300, units = "in")
    print(p)
    dev.off()
  }
  
  
  
  
  # DimPlot(Adipogenesis, reduction = "phate", group.by = col.groups, label = FALSE, raster=FALSE) + 
  #   geom_line(data = Node.df.MaxBranch, aes(x = phate_1, y=phate_2),colour = "black", linewidth = 1.5, linetype='dashed') +
  #   geom_point(data = Node.df.MaxBranch, aes(x = phate_1, y=phate_2),colour = "black",  fill="white", size = 6, shape =21) +
  #   geom_text(data = Node.df.MaxBranch, aes(x = .data[["phate_1"]]*0.9, y=.data[["phate_2"]]*0.9, label = rownames(Node.df.MaxBranch)))
  # 
  # DimPlot(Adipogenesis, reduction = "phate", group.by = col.groups, label = FALSE, raster=FALSE) + 
  #   geom_line(data = Node.df.MaxBranch, aes(x = phate_1, y=phate_2),colour = "red", linewidth = 1.5, linetype='dashed') +
  #   geom_point(data = Node.df.MaxBranch, aes(x = phate_1, y=phate_2),colour = "black",  fill="white", size = 6, shape =21) +
  #   geom_text(data = Node.df.MaxBranch, aes(x = .data[["phate_1"]]*0.9, y=.data[["phate_2"]]*0.9, label = rownames(Node.df.MaxBranch)))+ 
  #   
  #   geom_line(data = Node.df.MaxBranch2, aes(x = phate_1, y=phate_2),colour = "green", linewidth = 1.5, linetype='dashed') +
  #   geom_point(data = Node.df.MaxBranch2, aes(x = phate_1, y=phate_2),colour = "black",  fill="white", size = 6, shape =21) +
  #   geom_text(data = Node.df.MaxBranch2, aes(x = .data[["phate_1"]]*0.9, y=.data[["phate_2"]]*0.9, label = rownames(Node.df.MaxBranch2)))+ 
  #   
  #   geom_line(data = Node.df.MaxBranch3, aes(x = phate_1, y=phate_2),colour = "black", linewidth = 1.5, linetype='dashed') +
  #   geom_point(data = Node.df.MaxBranch3, aes(x = phate_1, y=phate_2),colour = "black",  fill="white", size = 6, shape =21) +
  #   geom_text(data = Node.df.MaxBranch3, aes(x = .data[["phate_1"]]*0.9, y=.data[["phate_2"]]*0.9, label = rownames(Node.df.MaxBranch3)))
  # 
  # DimPlot(Adipogenesis, reduction = "phate", group.by = col.groups, label = TRUE, raster=FALSE) + 
  #   geom_line(data = Node.df.MaxBranch, aes(x = phate_1, y=phate_2),colour = "red", linewidth = 1.5, linetype='dashed') +
  #   geom_point(data = Node.df.MaxBranch, aes(x = phate_1, y=phate_2),colour = "black",  fill="white", size = 6, shape =21) +
  #   
  #   geom_line(data = Node.df.MaxBranch2, aes(x = phate_1, y=phate_2),colour = "green", linewidth = 1.5, linetype='dashed') +
  #   geom_point(data = Node.df.MaxBranch2, aes(x = phate_1, y=phate_2),colour = "black",  fill="white", size = 6, shape =21) +
  #   
  #   geom_line(data = Node.df.MaxBranch3, aes(x = phate_1, y=phate_2),colour = "black", linewidth = 1.5, linetype='dashed') +
  #   geom_point(data = Node.df.MaxBranch3, aes(x = phate_1, y=phate_2),colour = "black",  fill="white", size = 6, shape =21) +
  #   
  #   geom_line(data = Node.df.MaxBranch4, aes(x = phate_1, y=phate_2),colour = "blue", linewidth = 1.5, linetype='dashed') +
  #   geom_point(data = Node.df.MaxBranch4, aes(x = phate_1, y=phate_2),colour = "black",  fill="white", size = 6, shape =21) +
  #   
  #   geom_line(data = Node.df.MaxBranch5, aes(x = phate_1, y=phate_2),colour = "purple", linewidth = 1.5, linetype='dashed') +
  #   geom_point(data = Node.df.MaxBranch5, aes(x = phate_1, y=phate_2),colour = "black",  fill="white", size = 6, shape =21) +
  #   
  #   geom_line(data = Node.df.MaxBranch6, aes(x = phate_1, y=phate_2),colour = "yellow", linewidth = 1.5, linetype='dashed') +
  #   geom_point(data = Node.df.MaxBranch6, aes(x = phate_1, y=phate_2),colour = "black",  fill="white", size = 6, shape =21) +
  #   
  #   geom_line(data = Node.df.MaxBranch7, aes(x = phate_1, y=phate_2),colour = "orange", linewidth = 1.5, linetype='dashed') +
  #   geom_point(data = Node.df.MaxBranch7, aes(x = phate_1, y=phate_2),colour = "black",  fill="white", size = 6, shape =21) +
  #   
  #   geom_line(data = Node.df.MaxBranch8, aes(x = phate_1, y=phate_2),colour = "pink", linewidth = 1.5, linetype='dashed') +
  #   geom_point(data = Node.df.MaxBranch8, aes(x = phate_1, y=phate_2),colour = "black",  fill="white", size = 6, shape =21) +
  #   
  #   geom_line(data = Node.df.MaxBranch9, aes(x = phate_1, y=phate_2),colour = "darkgray", linewidth = 1.5, linetype='dashed') +
  #   geom_point(data = Node.df.MaxBranch9, aes(x = phate_1, y=phate_2),colour = "black",  fill="white", size = 6, shape =21) + 
  #   theme(legend.position="right")
  
  
  
  # Get only pseudotime and put pseudotime(pt) to metadata (the differentiation branch connecting FAPs and adipocytes)
  ## extract pseudotime 
  gc() #free up memrory and report the memory usage.
  
  rm(meta_data_merge, 
     tree_data, TreeEPG, 
     Tree_Graph, Tree_e2e,
     PartStruct, ProjStruct, Node.df, Node.list.df, Node.list.df2,p,PlotPG,
     allDistancesDf, sort.BranchNames,
     Pt.list, EdgeID.branchs, Node.list, my_cols2, nn, aa, dir1, dir0,Branch.name, colors, coul, cell.type, cell.types, hight.cells, MaxBranch,MaxBranch2, MaxBranch3,maxLength, 
     Adipogenesis,  meta_data_merge )
  gc() #free up memrory and report the memory usage.
}
gc() #free up memrory and report the memory usage.



