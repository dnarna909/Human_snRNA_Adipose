  ### DISCLAIMER
#   Some of the algorithms are non-deterministic making the results slightly different from run to run.
# Many of the algorithms are frequently updated with new releases making the results slightly different from version to version.
# Depending on your system, this code below may not produce an exact 1:1 copy of the results. 
# This adversely affects clustering, filtering etc. if the code is naively copy-pasted.
# For better reproducibility, please download quality-filtered Rds objects and work on downstream analyses. 
# <br>
#   change "D:/" to "/media/jianie/Extreme SSD1/" # if in LUNIX computer

### Load libraries -------------------------------------------------------------------------------------------------------------
rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.
memory.limit(size = 1e+13)


# source(paste0("/media/jianie/Extreme SSD1/2022-09-01 STARR_SGLT2 Combine/Project Parameters.R"), local = knitr::knit_global())
source(paste0("D:/2022-09-01 STARR_SGLT2 Combine/Project Parameters.R"), local = knitr::knit_global())

# source(paste0(Disk, "00_Functions_Refs/Functions_Human Fat snRNA.R"), local = knitr::knit_global())
# or sys.source("your-script.R", envir = knitr::knit_global())
source(paste0(Disk, "00_Functions_Refs/Functions_marker.genes.R"), local = knitr::knit_global())

# import data -------------------------------------
new.meta.data <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/","CellCycle.Senescence.stage.metadata.Rds"))  %>%  
  tibble::rownames_to_column(var = "rowname")

seurat.data <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/", "SAT_Annotated_Subtype.Rds"))
metadata <- seurat.data@meta.data  %>% select(-one_of(c(intersect(colnames(seurat.data@meta.data), colnames(new.meta.data)) )))%>%  
  tibble::rownames_to_column(var = "rowname") 
# intersect(colnames(metadata), colnames(new.meta.data))
meta_data_merge <- metadata %>% left_join(new.meta.data, by = c("rowname"))  %>%  tibble::column_to_rownames("rowname")
seurat.data <- AddMetaData( object = seurat.data, metadata = meta_data_merge)



# FindMarkers and plotting ----------------------------------------------------------------------------
library(dplyr)
organism = "Homo sapiens"
# Prepare for WikiPathways analysis
#gmt <- rWikiPathways::downloadPathwayArchive(organism = "Homo sapiens", format = "gmt")
gmt <- rWikiPathways::downloadPathwayArchive(date="20210510", organism = "Homo sapiens", format = "gmt")
# gmt <- rWikiPathways::downloadPathwayArchive(organism = organism, format = "gmt")
wp2gene <- clusterProfiler::read.gmt(gmt)
#wp2gene <- read.gmt(paste0("/media/jianie/Extreme SSD1/2021-01-11 Figures for Grant resubmission/wikipathways-20210210-gmt-Homo_sapiens.gmt"))
wp2gene <- wp2gene %>% tidyr::separate(term, c("name", "version", "wpid", "org"), "%")
wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) # TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid, name) # TERM2NAME
OrgDb = org.Hs.eg.db


Results <- list()

# Analysis 1 --------------
Annotation.type=NULL; ident.group = "Senescence.stage"; ident1 = "sens_4"; ident2 = NULL;
Idents(object = seurat.data) <- ident.group;
Result1 <- FindMarkers(seurat.data, 
                       ident.1 = ident1, ident.2 = ident2, logfc.threshold = 0.25,
                       test.use = "wilcox",
                       min.pct = 0.2) ;
Result1$pct.compare <- Result1$pct.1/Result1$pct.2
file.name <- ifelse(is.null(Annotation.type), paste0(ident1, ifelse(is.null(ident2), "", paste0("_", ident2))), 
                    paste0(Annotation.type, "_",ident1,ifelse(is.null(ident2), "", paste0("_", ident2)))
)
saveRDS(Result1, file = paste0(Disk, Project.folder, "/", Rds.folder, "/",
                               file.name, "_MarkersPathwaysResult.Rds"))
Results[["sens_4"]] <- Find.Markers.Plot(Result = Result1, ident.group.= "Senescence.stage",  ident1. = "sens_4")
DimPlot(seurat.data, group.by = ident.group, raster=FALSE )
FeaturePlot(seurat.data, features = Result1 %>% arrange(avg_log2FC, p_val_adj, pct.compare) %>% top_n(3) %>% rownames(), raster=FALSE )
FeaturePlot(seurat.data, features = c("COL1A1", "GPC6", "TEAD1", "COL3A1"), raster=FALSE )
VlnPlot(seurat.data, features = c("COL1A1"), raster=FALSE, group.by = "Senescence.stage" )



# Analysis 2 ------------------------
Annotation.type = "FAP"; ident.group = "Senescence.stage"; ident1 = "sens_4"; ident2 = NULL;
Idents(object = seurat.data) <- "Annotation";
Result1 <- FindMarkers(seurat.data, 
                       ident.1 = ident1, ident.2 = ident2,
                       group.by = ident.group, subset.ident = Annotation.type,
                       test.use = "wilcox", logfc.threshold = 0.25,
                       min.pct = 0.2) ;
file.name <- ifelse(is.null(Annotation.type), paste0(ident1, ifelse(is.null(ident2), "", paste0("_", ident2))), 
                    paste0(Annotation.type, "_",ident1,ifelse(is.null(ident2), "", paste0("_", ident2)))
)
Result1$pct.compare <- Result1$pct.1/Result1$pct.2
saveRDS(Result1, file = paste0(Disk, Project.folder, "/", Rds.folder, "/",
                               file.name, "_MarkersPathwaysResult.Rds"))
Results[["FAP.sens_4"]] <- Find.Markers.Plot(Result = Result1, Annotation.type. = "FAP", ident.group.= "Senescence.stage",  ident1. = "sens_4")

Annotation.type = "FAP"; ident.group = "Senescence.stage"; ident1 = "sens_4"; ident2 = "sens_1";
Idents(object = seurat.data) <- "Annotation";
Result1 <- FindMarkers(seurat.data, 
                       ident.1 = ident1, ident.2 = ident2, 
                       group.by = ident.group, subset.ident = Annotation.type,
                       logfc.threshold = 0.25,
                       test.use = "wilcox",
                       min.pct = 0.2) ;
file.name <- ifelse(is.null(Annotation.type), paste0(ident1, ifelse(is.null(ident2), "", paste0("_", ident2))), 
                    paste0(Annotation.type, "_",ident1,ifelse(is.null(ident2), "", paste0("_", ident2)))
)
Result1$pct.compare <- Result1$pct.1/Result1$pct.2
saveRDS(Result1, file = paste0(Disk, Project.folder, "/", Rds.folder, "/",
                               file.name, "_MarkersPathwaysResult.Rds"))
Results[["FAP.sens_4.sens_1"]] <- Find.Markers.Plot(Result = Result1, Annotation.type = "FAP", ident.group= "Senescence.stage",  
                                                           ident1. = "sens_1", ident2. = "sens_4" )

saveRDS(Results, file = paste0(Disk, Project.folder, "/", Rds.folder, "/",
                               "All_MarkersPathwaysResult.Rds"))

# Analysis 3 ----------------------
Annotation.type = "Endothelial"; ident.group = "Senescence.stage"; ident1 = "sens_4"; ident2 = NULL;
Idents(object = seurat.data) <-  "Annotation";
Result1 <- FindMarkers(seurat.data, 
                       ident.1 = ident1, ident.2 = ident2, 
                       group.by = ident.group, subset.ident = Annotation.type,
                       logfc.threshold = 0.25,
                       test.use = "wilcox",
                       min.pct = 0.2) ;
file.name <- ifelse(is.null(Annotation.type), paste0(ident1, ifelse(is.null(ident2), "", paste0("_", ident2))), 
                    paste0(Annotation.type, "_",ident1,ifelse(is.null(ident2), "", paste0("_", ident2)))
)
Result1$pct.compare <- Result1$pct.1/Result1$pct.2
saveRDS(Result1, file = paste0(Disk, Project.folder, "/", Rds.folder, "/",
                               file.name, "_MarkersPathwaysResult.Rds"))
Results[["Endothelial.sens_4"]] <- Find.Markers.Plot(Result = Result1, Annotation.type = "Endothelial", ident.group= "Senescence.stage",  ident1 = "sens_4")


# Analysis 4 --------------------
Annotation.type = "Adipocytes"; ident.group = "Senescence.stage"; ident1 = "sens_4"; ident.2 = NULL;
Idents(object = seurat.data) <- "Annotation";
Result1 <- FindMarkers(seurat.data, 
                       ident.1 = ident1, ident.2 = ident2, 
                       group.by = ident.group, subset.ident = Annotation.type,
                       logfc.threshold = 0.25,
                       test.use = "wilcox",
                       min.pct = 0.2) ;
file.name <- ifelse(is.null(Annotation.type), paste0(ident1, ifelse(is.null(ident2), "", paste0("_", ident2))), 
                    paste0(Annotation.type, "_",ident1,ifelse(is.null(ident2), "", paste0("_", ident2)))
)
Result1$pct.compare <- Result1$pct.1/Result1$pct.2
saveRDS(Result1, file = paste0(Disk, Project.folder, "/", Rds.folder, "/",
                               file.name, "_MarkersPathwaysResult.Rds"))
Results[["Adipocytes.sens_4"]] <- Find.Markers.Plot(Result = Result1, Annotation.type = "Adipocytes", ident.group= "Senescence.stage",  
                                                           ident1 = "sens_4" )

saveRDS(Results, file = paste0(Disk, Project.folder, "/", Rds.folder, "/",
                               "All_MarkersPathwaysResult.Rds"))

# Analysis 5 ---------------------------
Annotation.type = "Adipocytes"; ident.group = "Senescence.stage"; ident1 = "sens_1"; ident2 = "sens_3";
Idents(object = seurat.data) <- "Annotation";
Result1 <- FindMarkers(seurat.data, 
                       ident.1 = ident1, ident.2 = ident2, 
                       group.by = ident.group, subset.ident = Annotation.type,
                       logfc.threshold = 0.25,
                       test.use = "wilcox",
                       min.pct = 0.2) ;
file.name <- ifelse(is.null(Annotation.type), paste0(ident1, ifelse(is.null(ident2), "", paste0("_", ident2))), 
                    paste0(Annotation.type, "_",ident1,ifelse(is.null(ident2), "", paste0("_", ident2)))
)
Result1$pct.compare <- Result1$pct.1/Result1$pct.2
saveRDS(Result1, file = paste0(Disk, Project.folder, "/", Rds.folder, "/",
                               file.name, "_MarkersPathwaysResult.Rds"))
Results[["Adipocytes.sens_1.sens_3"]] <- Find.Markers.Plot(Result = Result1, Annotation.type = "Adipocytes", ident.group= "Senescence.stage",  
                                                           ident1 = "sens_1", ident2 = "sens_3" )


# Analysis 6 ---------------------------------
Annotation.type = "Adipocytes"; ident.group = "Senescence.stage"; ident1 = "sens_1"; ident2 = "sens_4";
Idents(object = seurat.data) <- "Annotation";
Result1 <- FindMarkers(seurat.data, 
                       ident.1 = ident1, ident.2 = ident2, 
                       group.by = ident.group, subset.ident = Annotation.type,
                       logfc.threshold = 0.25,
                       test.use = "wilcox",
                       min.pct = 0.2) ;
file.name <- ifelse(is.null(Annotation.type), paste0(ident1, ifelse(is.null(ident2), "", paste0("_", ident2))), 
                    paste0(Annotation.type, "_",ident1,ifelse(is.null(ident2), "", paste0("_", ident2)))
)
Result1$pct.compare <- Result1$pct.1/Result1$pct.2
saveRDS(Result1, file = paste0(Disk, Project.folder, "/", Rds.folder, "/",
                               file.name, "_MarkersPathwaysResult.Rds"))
Results[["Adipocytes.sens_1.sens_4"]] <- Find.Markers.Plot(Result = Result1, Annotation.type = "Adipocytes", ident.group= "Senescence.stage",  
                                                           ident1. = "sens_1", ident2. = "sens_4" )

saveRDS(Results, file = paste0(Disk, Project.folder, "/", Rds.folder, "/",
                               "All_MarkersPathwaysResult.Rds"))


# Analysis 7 -----------------------------------
Annotation.type = "Adipocytes"; ident.group = "Senescence.stage"; ident1 = "sens_2"; ident2 = "sens_3";
Idents(object = seurat.data) <- "Annotation";
Result1 <- FindMarkers(seurat.data, 
                       ident.1 = ident1, ident.2 = ident2, 
                       group.by = ident.group, subset.ident = Annotation.type,
                       logfc.threshold = 0.25,
                       test.use = "wilcox",
                       min.pct = 0.2) ;
file.name <- ifelse(is.null(Annotation.type), paste0(ident1, ifelse(is.null(ident2), "", paste0("_", ident2))), 
                    paste0(Annotation.type, "_",ident1,ifelse(is.null(ident2), "", paste0("_", ident2)))
)
Result1$pct.compare <- Result1$pct.1/Result1$pct.2
saveRDS(Result1, file = paste0(Disk, Project.folder, "/", Rds.folder, "/",
                               file.name, "_MarkersPathwaysResult.Rds"))
Results[["Adipocytes.sens_2.sens_3"]] <- Find.Markers.Plot(Result = Result1, Annotation.type = "Adipocytes", ident.group= "Senescence.stage",  
                                                           ident1 = "sens_2", ident2 = "sens_3" )


# Analysis 8 ------------------------------------
Annotation.type = "Adipocytes"; ident.group = "Senescence.stage"; ident1 = "sens_3"; ident2 = "sens_4";
Idents(object = seurat.data) <- "Annotation";
Result1 <- FindMarkers(seurat.data, 
                       ident.1 = ident1, ident.2 = ident2, 
                       group.by = ident.group, subset.ident = Annotation.type,
                       logfc.threshold = 0.25,
                       test.use = "wilcox",
                       min.pct = 0.2) ;
file.name <- ifelse(is.null(Annotation.type), paste0(ident1, ifelse(is.null(ident2), "", paste0("_", ident2))), 
                    paste0(Annotation.type, "_",ident1,ifelse(is.null(ident2), "", paste0("_", ident2)))
)
Result1$pct.compare <- Result1$pct.1/Result1$pct.2
saveRDS(Result1, file = paste0(Disk, Project.folder, "/", Rds.folder, "/",
                               file.name, "_MarkersPathwaysResult.Rds"))
Results[["Adipocytes.sens_3.sens_4"]] <- Find.Markers.Plot(Result = Result1, Annotation.type = "Adipocytes", ident.group= "Senescence.stage",  
                                                           ident1 = "sens_3", ident2 = "sens_4" )


# Analysis 9 --------------------


