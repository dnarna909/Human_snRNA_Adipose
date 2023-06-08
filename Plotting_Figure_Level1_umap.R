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

# Import data -----------------------------------------------------------------------------------------------------------------------------------------------------
file.name <- stringr::str_split(seurat.file, ".Rds")[[1]][1]
Annotation.file.name <- paste(stringr::str_split(file.name, "_")[[1]][c(1,2,3)], collapse = "_")
Subtype.file.name <- stringr::str_split(file.name, paste0(Annotation.file.name, "_"))[[1]][2]
if (is.na(Subtype.file.name)) {Subtype.file.name = ""}

# Import metadata
sample.meta <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/", "SAT_Sample.meta.data.Rds")) %>% tibble::rownames_to_column(var = "rowname")
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

# prepare folder
export.folder0 <- paste0(Annotation.file.name)
dir.create(file.path(paste0(Disk, Project.folder, "/", figures.folder, "/"), export.folder0), showWarnings = FALSE)
dir <- paste0(Disk, Project.folder, "/", figures.folder, "/", export.folder0, "/")

# plot 1 ------------------------------------------------
set.seed(100)
p <- DimPlot(seurat.data, label = T, group.by = meta.group, pt.size = 0.01, raster=FALSE) + 
  # NoLegend() +
  labs(title = file.name, 
       tag = paste0("n = ", format(dim(seurat.data)[2], big.mark=",", scientific=FALSE) )) +
  theme(plot.title = element_text(color="black", size=14, face="bold", hjust=0.5, vjust=0),
        plot.tag = element_text(color="black", size=12, face="plain"),
        plot.tag.position = c(0.15, 0.16),
        plot.margin = unit(c(0,0,0,0), "cm"),
        panel.border = element_rect(linetype = "solid", 
                                    fill = NA, colour = "black"))
p
png(file=paste0(dir, file.name, "_", meta.group, "_umap_labeled_legend.png"), width=6, height=4, 
    res = 300, units = "in")
print(p)
dev.off()

p <- DimPlot(seurat.data, label = T, group.by = meta.group, pt.size = 0.01, raster=FALSE) + 
  NoLegend() +
  labs(title = file.name, 
       tag = paste0("n = ", format(dim(seurat.data)[2], big.mark=",", scientific=FALSE) )) +
  theme(plot.title = element_text(color="black", size=14, face="bold", hjust=0.5, vjust=0),
        plot.tag = element_text(color="black", size=12, face="plain"),
        plot.tag.position = c(0.15, 0.16),
        plot.margin = unit(c(0,0,0,0), "cm"),
        panel.border = element_rect(linetype = "solid", 
                                    fill = NA, colour = "black"))
p
png(file=paste0(dir, file.name, "_", meta.group, "_umap_labeled.png"), width=6, height=4, 
    res = 300, units = "in")
print(p)
dev.off()

p <- DimPlot(seurat.data, label = F, group.by = meta.group, pt.size = 0.01, raster=FALSE) + NoLegend() +
  labs(title = file.name, 
       tag = paste0("n = ", format(dim(seurat.data)[2], big.mark=",", scientific=FALSE) )) +
  theme(plot.title = element_text(color="black", size=14, face="bold", hjust=0.5, vjust=0),
        plot.tag = element_text(color="black", size=12, face="plain"),
        plot.tag.position = c(0.15, 0.16),
        plot.margin = unit(c(0,0,0,0), "cm"),
        panel.border = element_rect(linetype = "solid", 
                                    fill = NA, colour = "black"))
p
png(file=paste0(dir, file.name, "_", meta.group, "_umap_noLabel.png"), 
    width=4, height=4, 
    res = 300, units = "in")
print(p)
dev.off()


# plot 2 ------------------------------------------------
# seurat.data$Treatment_Group <- factor(seurat.data$Treatment_Group, 
#                                       levels = c("C_Visit 5", "C_Visit 11", "D_Visit 5", "D_Visit 11"))
for (ll in names(groups.list)) {
  split.group = groups.list[[ll]][["split.group"]]
  seurat.data <- SetIdent(seurat.data, value = groups.list[[ll]][["subset.group"]])
  seurat.data.1 <- subset(seurat.data, idents = TRUE)
  
  cells = as.data.frame(table(seurat.data.1@meta.data[[split.group]])) %>% dplyr::filter(Freq > 0)
  cells$numbers = paste0("n = ", format(cells$Freq, big.mark=",", scientific=FALSE))
  subtitle = paste(paste(cells$numbers[1:nrow(cells)],  collapse = "                                "))
  set.seed(100)
  p2 <- DimPlot(seurat.data.1, label = FALSE, group.by = meta.group, pt.size = 0.01, 
                split.by = split.group, ncol = nrow(cells), raster=FALSE) + NoLegend() + NoAxes() +
    labs(title = split.group, tag = subtitle) +
    theme(plot.title = element_text(color="black", size=14, face="bold", hjust=0.5, vjust=5),
          plot.tag = element_text(color="black", size=10, face="bold"),
          plot.tag.position = c(0.082, 0.79),# plot.tag.position = c(0.08, 0.815) for above the Subject line
          plot.margin = unit(c(0.5,0,0,0), "cm")
    )
  p2
  png(file=paste0(dir, file.name, "_", meta.group, "_umap_", split.group, ".png"), 
      width= nrow(cells)*2, height=3, res = 300, units = "in")
  print(p2)
  dev.off()
}

# plot 3 ------------------------------------------------
Sample.group = "Dataset"
cells = as.data.frame(table(seurat.data@meta.data[[Sample.group]]))
cells$numbers = paste0("n = ", format(cells$Freq, big.mark=",", scientific=FALSE))
nrow2 = round(nrow(cells)/6) 
nrow1 = floor(nrow(cells)/6)

if (nrow2 == 1 ){
  subtitle = paste(paste(cells$numbers[1:nrow(cells)],  
                         collapse = "                                                             "))
}

if (nrow2 == 2 ){
  subtitle = paste(paste(cells$numbers[1:6],  
                         collapse = "                                                             "), 
                   paste(cells$numbers[7:nrow(cells)], 
                         collapse = "                                                             "),  
                   sep="\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n")
}

if (nrow2 == 3 ){
  subtitle = paste(paste(cells$numbers[1:6],  
                         collapse = "                                                             "), 
                   paste(cells$numbers[7:12], 
                         collapse = "                                                             "), 
                   paste(cells$numbers[13:nrow(cells)], 
                         collapse = "                                                             "), 
                   sep="\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n")
}

if (nrow2 == 4 ){
  subtitle = paste(paste(cells$numbers[1:6],  
                         collapse = "                                 "), 
                   paste(cells$numbers[7:12], 
                         collapse = "                                 "), 
                   paste(cells$numbers[13:18], 
                         collapse = "                                 "), 
                   paste(cells$numbers[19:nrow(cells)], 
                         collapse = "                                 "), 
                   sep="\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n")
}
set.seed(100)
p3 <- DimPlot(seurat.data, label = FALSE, group.by = meta.group, 
              split.by = Sample.group, ncol = ifelse(nrow(cells) < 6, 1, 6), raster=FALSE) + NoLegend() + NoAxes() +
  labs(title = Sample.group, tag = subtitle) +
  theme(plot.title = element_text(color="black", size=14, face="bold", hjust=0.5, vjust=5),
        plot.tag = element_text(color="black", size=10, face="bold"),
        plot.tag.position = c(0.065, 0.73),# plot.tag.position = c(0.08, 0.815) for above the Subject line
        plot.margin = unit(c(0.5,0,0,0), "cm")
  )
p3
png(file=paste0(dir, file.name, "_", meta.group, "_umap_", Sample.group, ".png"), 
    width= ifelse(nrow(cells) < 6, nrow(cells)*2, 12), height=3*nrow2, res = 300, units = "in")
p3
dev.off()



rm(p,p2,p3)
rm(file.name, Annotation.file.name, Subtype.file.name,seurat.file,
   sample.meta, Subtype.meta, metadata.all, seurat.data,ll, seurat.data.1,split.group,
   export.folder0, dir, cells, subtitle, Sample.group, nrow2, nrow1)
