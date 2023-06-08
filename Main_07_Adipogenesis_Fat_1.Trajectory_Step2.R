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
  
  # Extract phate coordinates and train a tree ----------------------------------------------------------------------------------------------
  tree_data <- Embeddings(Adipogenesis, "phate")
  TreeEPG <- computeElasticPrincipalTree(X = tree_data, NumNodes = 30, Lambda = .005, Mu = .001, drawAccuracyComplexity = FALSE, drawPCAView = FALSE, drawEnergy = FALSE)
  
  # plotting
  PlotPG(X = tree_data, TargetPG = TreeEPG[[1]], DimToPlot = 1:2, Do_PCA = FALSE) # tree map
  # save figure
  png(filename = paste0(dir0, "Adipogenesis.PHATE_", nn, ".tree..png"),
      width = 8, height = 5, units = "in", res = 300)
  PlotPG(X = tree_data, TargetPG = TreeEPG[[1]], DimToPlot = 1:2, Do_PCA = FALSE)
  dev.off()
  
  #' ### Infer pseudotime 
  Tree_Graph <- ConstructGraph(TreeEPG[[1]])
  Tree_e2e <- GetSubGraph(Net = Tree_Graph, Structure = 'branches')
  PartStruct <- PartitionData(X = tree_data, NodePositions = TreeEPG[[1]]$NodePositions)
  ProjStruct <- project_point_onto_graph(X = tree_data,
                                         NodePositions = TreeEPG[[1]]$NodePositions,
                                         Edges = TreeEPG[[1]]$Edges$Edges,
                                         Partition = PartStruct$Partition)
  plot(tree_data, pch = 21, bg = "lightblue", col = "blue",  cex = 1 ) 
  plot(TreeEPG[[1]]$NodePositions)
  text((TreeEPG[[1]]$NodePositions)*0.9, rownames(as.data.frame(TreeEPG[[1]]$NodePositions)), family = "sans")
  gc() #free up memrory and report the memory usage.
  
  # save figure
  png(filename = paste0(dir0, "Adipogenesis.tree_data.", nn, ".png"), 
      width = 6, height = 5, units = "in", res = 300)
  plot(tree_data, pch = 20, bg = "lightblue", col = "gray",  cex = 0.3 ) 
  dev.off()
  
  png(filename = paste0(dir0, "Adipogenesis.TreeEPG.NodePositions.", nn, ".png"), 
      width = 6, height = 5, units = "in", res = 300)
  plot(TreeEPG[[1]]$NodePositions)
  text((TreeEPG[[1]]$NodePositions)*0.9, rownames(as.data.frame(TreeEPG[[1]]$NodePositions)), family = "sans")
  dev.off()
  
  # put Node to seurat data
  Adipogenesis@meta.data$EdgeID <- ProjStruct$EdgeID
  print(DimPlot(Adipogenesis, reduction = "phate", group.by = "EdgeID", label = TRUE, raster=FALSE))
  # DimPlot(Adipogenesis, reduction = "phate", group.by = "EdgeID", split.by =  "Age_Group", label = T, raster=FALSE)
  # DimPlot(Adipogenesis, reduction = "phate", group.by = "EdgeID", split.by =  "Status", label = T, raster=FALSE)
  
  print(DimPlot(Adipogenesis, reduction = "phate", group.by = "EdgeID", label = FALSE, raster=FALSE) +
          annotate("point", x = TreeEPG[[1]]$NodePositions[, "phate_1"], y=TreeEPG[[1]]$NodePositions[, "phate_2"], colour = "blue"))
  
  print(DimPlot(Adipogenesis, reduction = "phate", group.by = "Annotation", label = FALSE, raster=FALSE) +
          annotate("point", x = TreeEPG[[1]]$NodePositions[, "phate_1"], y=TreeEPG[[1]]$NodePositions[, "phate_2"], colour = "blue"))
  
  print(DimPlot(Adipogenesis, reduction = "phate", group.by = "Subtype", label = FALSE, raster=FALSE) +
          annotate("point", x = TreeEPG[[1]]$NodePositions[, "phate_1"], y=TreeEPG[[1]]$NodePositions[, "phate_2"], colour = "black",  fill="white", size = 3, shape =21) + 
          annotate("text", x = TreeEPG[[1]]$NodePositions[, "phate_1"]*0.9, y=TreeEPG[[1]]$NodePositions[, "phate_2"]*0.9, label = rownames(as.data.frame(TreeEPG[[1]]$NodePositions)) ) )
  
  # Get the max length of those lists #
  maxLength <- max(lengths(Tree_e2e))
  
  # generating a dataframe from the nested list, making all lengths equal
  allDistancesDf <- as.data.frame(do.call(rbind, lapply(Tree_e2e, `length<-`, maxLength)))
  
  sort.BranchNames <- names(sort(lengths(Tree_e2e), decreasing = TRUE))
  MaxBranch <- names(Tree_e2e[lengths(Tree_e2e) == maxLength])
  MaxBranch2 <- names(sort(lengths(Tree_e2e), decreasing = TRUE)[2])
  MaxBranch3 <- names(sort(lengths(Tree_e2e), decreasing = TRUE)[3])
  
  # plot NodePositions in DimPlot
  Node.df <- as.data.frame(TreeEPG[[1]]$NodePositions) 
  # Node.df.MaxBranch <- na.omit(Node.df[as.numeric(allDistancesDf[MaxBranch, ]), ]) %>% tibble::rownames_to_column(var = "Node") 
  # Node.df.MaxBranch2 <- na.omit(Node.df[as.numeric(allDistancesDf[MaxBranch2, ]), ])
  # Node.df.MaxBranch3 <- na.omit(Node.df[as.numeric(allDistancesDf[MaxBranch3, ]), ])
  # Node.df.MaxBranch4 <- na.omit(Node.df[as.numeric(allDistancesDf[names(sort(lengths(Tree_e2e), decreasing = TRUE)[4]), ]), ])
  # Node.df.MaxBranch5 <- na.omit(Node.df[as.numeric(allDistancesDf[names(sort(lengths(Tree_e2e), decreasing = TRUE)[5]), ]), ])
  # Node.df.MaxBranch6 <- na.omit(Node.df[as.numeric(allDistancesDf[names(sort(lengths(Tree_e2e), decreasing = TRUE)[6]), ]), ])
  # Node.df.MaxBranch7 <- na.omit(Node.df[as.numeric(allDistancesDf[names(sort(lengths(Tree_e2e), decreasing = TRUE)[7]), ]), ])
  # Node.df.MaxBranch8 <- na.omit(Node.df[as.numeric(allDistancesDf[names(sort(lengths(Tree_e2e), decreasing = TRUE)[8]), ]), ])
  # Node.df.MaxBranch9 <- na.omit(Node.df[as.numeric(allDistancesDf[names(sort(lengths(Tree_e2e), decreasing = TRUE)[9]), ]), ])
  
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
    mutate( amount = 1/n(), radius= sqrt(amount)/500)
  
  my_cols2 <- scales::hue_pal()(length(unique(Adipogenesis@meta.data[["Subtype"]])))
  scales::show_col(scales::hue_pal()(length(unique(Adipogenesis@meta.data[["Subtype"]]))))
  # my_cols2 <- colorRampPalette(brewer.pal(12, "Paired"))(length(unique(Adipogenesis@meta.data[["Subtype"]])))
  # pie(rep(1, length(my_cols2)), col = my_cols2, main="") 
  names(my_cols2) <- sort(unique(Adipogenesis@meta.data[["Subtype"]]))
  
  library(ggforce)
  print(DimPlot(Adipogenesis, reduction = "phate", group.by = "Subtype", label = TRUE, raster=FALSE, cols = my_cols2,  repel=TRUE) +
          geom_line(data = Node.list.df, aes(x = phate_1, y=phate_2, group = Branch), linewidth = 1.5, linetype='dashed')  + # colour = "black", 
          geom_point(data = Node.list.df, aes(x = phate_1, y=phate_2, fill = Branch),colour = "black",  size = 6, shape =21) +# fill="white", 
          geom_arc_bar(data = Node.list.df2,  aes(x0=phate_1,y0=phate_2,r0=0,r=radius,amount=amount,fill=Branch),stat="pie", linewidth = 0.5)+
          coord_equal()+ # so one gets actually circles
          scale_fill_manual(values=c(colors))   +
          labs(fill="Branch", color="Group", group = "Branch"))
  
  print(DimPlot(Adipogenesis, reduction = "phate", group.by = "EdgeID", label = TRUE, raster=FALSE,  repel=TRUE) +
          geom_line(data = Node.list.df, aes(x = phate_1, y=phate_2, group = Branch), linewidth = 1.5, linetype='dashed')  + # colour = "black", 
          geom_point(data = Node.list.df, aes(x = phate_1, y=phate_2, fill = Branch),colour = "black",  size = 6, shape =21) +# fill="white", 
          geom_arc_bar(data = Node.list.df2,  aes(x0=phate_1,y0=phate_2,r0=0,r=radius,amount=amount,fill=Branch),stat="pie", linewidth = 0.5)+
          coord_equal()+ # so one gets actually circles
          scale_fill_manual(values=c(colors))   +
          labs(fill="Branch", color="Group", group = "Branch") )
  
  # Get only pseudotime and put pseudotime(pt) to metadata (the differentiation branch connecting FAPs and adipocytes)
  ## extract pseudotime 
  Pt.list <- list()
  for (i in names(Tree_e2e) ) {
    Pt.list[[i]] <- getPseudotime(ProjStruct = ProjStruct, NodeSeq = rev(Tree_e2e[[i]])) 
  }
  
  ## put pseudotime to metadata
  meta_data_pt <- Adipogenesis@meta.data[, c("Label", "Annotation", "Subtype")] %>% tibble::rownames_to_column(var = "rowname")
  # EdgeID.branchs <- data.frame(matrix(nrow = length(unique(Adipogenesis$EdgeID)), ncol = 0))
  EdgeID.branchs <- data.frame()
  for (ii in names(Pt.list) ) {
    meta_data_pt[[paste0(ii,"_Pt")]] <-  Pt.list[[ii]]$Pt
    meta_data_pt[[paste0(ii)]] <-  ifelse(meta_data_pt[[paste0(ii,"_Pt")]] != "NA", ii, NA)
    EdgeID.branchs <- rbind(EdgeID.branchs, as.data.frame(table(Adipogenesis$EdgeID, meta_data_pt[[ii]]))) %>% dplyr::filter(Freq != 0)
  }
  colnames(EdgeID.branchs) <- c("EdgeID", "Branch", "CellNumbers")
  gc() #free up memrory and report the memory usage.
  
  
  ## add metadata to seurat data
  metadata <- Adipogenesis@meta.data %>%  tibble::rownames_to_column(var = "rowname")
  intersect(colnames(meta_data_pt), colnames(metadata))
  meta_data_merge <- metadata %>% full_join(meta_data_pt, by = c(intersect(colnames(meta_data_pt), colnames(metadata)))) %>%  tibble::column_to_rownames("rowname")
  # Adipogenesis <- AddMetaData( object = Adipogenesis, metadata = meta_data_merge)
  
  ## save
  saveRDS(meta_data_merge, file=paste0(Disk, Project.folder, "/", Rds.folder, "/","Adipogenesis.meta.data.PHATE_", nn, ".Pt.Rds"))
  save(tree_data, TreeEPG, Tree_Graph, Tree_e2e, PartStruct, ProjStruct,  allDistancesDf, sort.BranchNames,
       maxLength, MaxBranch, MaxBranch2, MaxBranch3, Pt.list, EdgeID.branchs, 
       file = paste0(Disk, Project.folder, "/", Rds.folder, "/", "Adipogenesis.PHATE_", nn, ".Branch.RData"))
  # saveRDS(Adipogenesis, paste0(dir1,"Adipogenesis.", nn, ".Branch.Rds"))
  gc() #free up memrory and report the memory usage.
  
  rm(meta_data_merge, 
     tree_data, TreeEPG, 
     Tree_Graph, Tree_e2e,
     PartStruct, ProjStruct, Node.df, Node.list.df, Node.list.df2,
     allDistancesDf, sort.BranchNames,
     Pt.list, EdgeID.branchs, Node.list, PlotPG, Branch.name, colors, coul, dir1, dir0, i, ii, MaxBranch,MaxBranch2, MaxBranch3, maxLength, my_cols2, nn, aa, 
     Adipogenesis, meta_data_pt, metadata, meta_data_merge )
  gc() #free up memrory and report the memory usage.
}
gc() #free up memrory and report the memory usage.



