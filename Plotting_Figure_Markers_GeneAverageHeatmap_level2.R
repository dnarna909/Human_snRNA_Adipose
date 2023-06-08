#'  change "D:/" to "/media/jianie/Extreme SSD/" # if in LUNIX computer
# rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.
memory.limit(size = 1e+13)

# import parameters and packages --------------------------------------------------------------------------------------------
# source(paste0("/media/jianie/Extreme SSD1/2022-09-01 STARR_SGLT2 Combine/Project Parameters.R"), local = knitr::knit_global())
if (serve == "GCCRI"){
  source(paste0(Disk, Project.folder, "/", "Project Parameters_GCCRI.R"), local = knitr::knit_global())
} else {
  source(paste0(Disk, Project.folder, "/", "Project Parameters.R"), local = knitr::knit_global())
}
library(Seurat)
library(ggpubr)
source(paste0(Disk, "00_Functions_Refs/Functions_Human Fat snRNA.R"), local = knitr::knit_global())

# prepare folder
export.folder0 <- paste0("Heatmap")
dir.create(file.path(paste0(Disk, Project.folder, "/", figures.folder, "/"), export.folder0), showWarnings = FALSE)
dir.create(file.path(paste0(Disk, Project.folder, "/", Rds.folder, "/"), export.folder0), showWarnings = FALSE)
dir <- paste0(Disk, Project.folder, "/", figures.folder, "/", export.folder0, "/")
dir1 <- paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder0, "/" )

  # Import data -----------------------------------------------------------------------------------------------------------------------------------------------------
  file.name <- stringr::str_split(seurat.file, ".Rds")[[1]][1]
  Annotation.file.name <- paste(stringr::str_split(file.name, "_")[[1]][c(1,2,3)], collapse = "_")
  Subtype.file.name <- stringr::str_split(file.name, paste0(Annotation.file.name, "_"))[[1]][2]
  if (is.na(Subtype.file.name)) {Subtype.file.name = ""}
  
  file.name = paste0(Annotation.file.name, "_", if (exists("Select.Group")){paste0(Select.Group, "_")}, 
                     meta.group1, "_",  meta.group2, "_",  
                     "GeneAverages")
  
if(!file.exists(paste0(dir1, file.name, ".rds"))) {
  # Import metadata
  sample.meta <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/", sample.file)) %>% tibble::rownames_to_column(var = "rowname")
  Subtype.meta <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/", Annotation.file.name, "_", file.type)) %>% tibble::rownames_to_column(var = "rowname")
  metadata.all <- left_join(Subtype.meta, sample.meta, by = "rowname") %>% tibble::column_to_rownames(var = "rowname")
  # metadata.all <- Subtype.meta %>% tibble::column_to_rownames(var = "rowname")
  
  # Import seurat data
  if (exists("sub.folder")){
    seurat.data <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/", sub.folder, "/", seurat.file)) 
  }
  if (!exists("sub.folder")){
    seurat.data <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/", seurat.file)) 
  }
  
  seurat.data <- AddMetaData(seurat.data, metadata = metadata.all)
  if (exists("Select.Group")){
    seurat.data <- SetIdent(seurat.data, value = Select.Group)
    seurat.data <- subset(seurat.data, idents = TRUE)
  }
  gc() #free up memory and report the memory usage.
  
  
  # calculate average of genes
  Genes <- rio::import(paste0(Disk, Project.folder, "/", Markers_file), sheet = sheet) %>% 
    mutate(List = 1:n()) %>% dplyr::filter(!is.na(gene))
  
  Genes <- Genes[Genes$gene %in% rownames(seurat.data), ] # match to seurat data
  Averages <- Maker.gene.heatmap.Label.c1.c2(seurat.data, Genes, 
                                             condition1 = meta.group1, 
                                             condition2 = meta.group2, groups.in.condition2 = unique(seurat.data@meta.data[[meta.group2]]) ) 
  Averages <- Averages %>% arrange(type)
  rownames(Averages) <- paste0(1:nrow(Averages), paste0(".", Averages$gene))
  saveRDS(Averages, paste0(dir1,  file.name, ".rds"))
} else {
  Averages <- rio::import(paste0(dir1,  file.name, ".rds"))
  
  rm(
     sample.meta, Subtype.meta, metadata.all, Genes)
}


# Plot the heatmap
library(RColorBrewer)
library(circlize)
coul <-colorRampPalette(brewer.pal(12, "Paired"))(length(unique(Averages$type)))
colors = coul[1: length(unique(Averages$type))]
pie(rep(1, length(coul)), col = coul , main="") 
names(colors) <- sort(unique(Averages$type))

row_ha = rowAnnotation(
  df = data.frame(type = Averages$type),
  col = list(type = colors), 
  show_legend = F , show_annotation_name = F,
  annotation_legend_param = list(
    type = list( 
      title_gp = gpar(fontsize = 10, 
                      fontface = "bold")
      , ncol = 2,  
      labels_gp = gpar(fontsize = 8))
    ,
  annotation_name_side = "left" )
)
ComplexHeatmap::draw(row_ha)
legend.col = colorRamp2(seq(-2, 2, length = 3), c("blue", "#EEEEEE", "red"))
row_labels = structure(Averages[["gene"]], names = rownames(Averages))

ht_opt( message = FALSE)
ht_list = 
  Heatmap(
    as.matrix(Averages %>% dplyr::select(!one_of(c("List"))) %>% dplyr::select_if(is.numeric)), name = "mat",
    col = legend.col , 
    border = TRUE,
    border_gp = gpar(col = "black", lty = 1),
    
    cluster_columns = F, cluster_rows = F, 
    cluster_column_slices = FALSE,
    
    # column_title = paste0(export.folder, "_GeneAverages"), 
    column_title_side = c("bottom"),
    column_title_gp = gpar(fontsize = 10, fontface = "plain"),
    #column_title_gp = gpar(fontsize = 17, fontface = "bold"),
    column_title_rot = 0,
    
    column_split= stringr::str_split_fixed(colnames(Averages %>% dplyr::select(!one_of(c("List"))) %>% dplyr::select_if(is.numeric)), ":", 2)[, 2],
    column_gap = unit(0, "mm"),
    
    show_column_names = FALSE,
    column_names_gp = gpar(fontsize = 12), 
    column_names_rot = 45, 
    column_names_max_height = unit(6, "cm"),
    
    row_names_side = "left", 
    show_row_names = TRUE,
    row_names_max_width = unit(6, "cm"),
    row_names_gp = gpar(fontsize = 7, fontface = "italic"),
    row_names_rot = 0,
    
    left_annotation = row_ha, 
    
    heatmap_legend_param = list(direction = "vertical", legend_height = unit(3, "cm"), 
                                direction = "vertical",
                                title = "Scaled average expression", 
                                title_gp = gpar(fontsize = 10, fontfface = "plain"),
                                title_position = "leftcenter-rot"),
    # show_heatmap_legend = FALSE,
    
    row_dend_side = NULL,
    na_col = "black",
    row_labels = row_labels[rownames(Averages)],
    use_raster = T
  )

labels <- colnames(Averages %>% dplyr::select(!one_of(c("List"))) %>% dplyr::select_if(is.numeric))
png(file = paste0(dir, file.name, "_Heatmap.png")
    ,  units ="in", res = 300
    , width = (ncol(Averages)/25.7)+1.8
    ,  height = (nrow(Averages)/6)+0.7
    
)
ComplexHeatmap::draw(ht_list , 
                     heatmap_legend_side = "right"
)
dev.off()


rm(file.name, Annotation.file.name, Subtype.file.name, seurat.file,
   coul, colors, row_ha,legend.col,row_labels, Averages,
   export.folder0, dir)
