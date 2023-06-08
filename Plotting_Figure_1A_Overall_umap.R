#'  change "D:/" to "/media/jianie/Extreme SSD/" # if in LUNIX computer
# rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.
memory.limit(size = 1e+13)

# import parameters and packages --------------------------------------------------------------------------------------------
# source(paste0("/media/jianie/Extreme SSD1/2022-09-01 STARR_SGLT2 Combine/Project Parameters.R"), local = knitr::knit_global())
source(paste0(paste0(Disk, Project.folder, "/", "Project Parameters.R")), local = knitr::knit_global())
library(Seurat)
library(ggpubr)

# Import data -----------------------------------------------------------------------------------------------------------------------------------------------------
file.name <- stringr::str_split(seurat.file, ".Rds")[[1]][1]
Annotation.file.name <- paste(stringr::str_split(file.name, "_")[[1]][c(1,2,3)], collapse = "_")
Subtype.file.name <- stringr::str_split(file.name, paste0(Annotation.file.name, "_"))[[1]][2]
if (is.na(Subtype.file.name)) {Subtype.file.name = ""}

# Import metadata
sample.meta <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/", "SAT_Sample.meta.data.Rds")) %>% tibble::rownames_to_column(var = "rowname")
Subtype.meta <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/", file.name, "_", file.type)) %>% tibble::rownames_to_column(var = "rowname")
metadata.all <- left_join(Subtype.meta, sample.meta, by = "rowname") %>% tibble::column_to_rownames(var = "rowname")
# metadata.all <- Subtype.meta %>% tibble::column_to_rownames(var = "rowname")

# Import seurat data
seurat.data <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/", seurat.file))
seurat.data <- AddMetaData(seurat.data, metadata = metadata.all)
seurat.data <- SetIdent(seurat.data, value = "Select_Group")
seurat.data <- subset(seurat.data, idents = TRUE)
gc() #free up memory and report the memory usage.

# prepare folder
export.folder0 <- paste0(file.name)
dir.create(file.path(paste0(Disk, Project.folder, "/", figures.folder, "/"), export.folder0), showWarnings = FALSE)
dir <- paste0(Disk, Project.folder, "/", figures.folder, "/", export.folder0, "/")

# plot ------------------------------------------------
set.seed(100)
p <- DimPlot(seurat.data, label = T, group.by = meta.group, pt.size = 0.01, raster=FALSE) + NoLegend() +
  labs(title = file.name, 
       tag = paste0("n = ", format(dim(seurat.data)[2], big.mark=",", scientific=FALSE) )) +
  theme(plot.title = element_text(color="black", size=14, face="bold", hjust=0.5, vjust=0),
        plot.tag = element_text(color="black", size=12, face="plain"),
        plot.tag.position = c(0.15, 0.16),
        plot.margin = unit(c(0,0,0,0), "cm"),
        panel.border = element_rect(linetype = "solid", 
                                    fill = NA, colour = "black"))
p
png(file=paste0(dir, file.name, "_", meta.group, "_umap.png"), width=4, height=4, 
    res = 300, units = "in")
print(p)
dev.off()

# seurat.data$Treatment_Group <- factor(seurat.data$Treatment_Group, 
#                                       levels = c("C_Visit 5", "C_Visit 11", "D_Visit 5", "D_Visit 11"))
for (ll in names(groups.list)) {
  split.group = groups.list[[ll]][["split.group"]]
  seurat.data <- SetIdent(seurat.data, value = groups.list[[ll]][["subset.group"]])
  seurat.data.1 <- subset(seurat.data, idents = TRUE)
  
  cells = as.data.frame(table(seurat.data.1@meta.data[[split.group]]))
  cells$numbers = paste0("n = ", format(cells$Freq, big.mark=",", scientific=FALSE))
  subtitle = paste(paste(cells$numbers[1:nrow(cells)],  collapse = "                                "))
  set.seed(100)
  p2 <- DimPlot(seurat.data.1, label = FALSE, group.by = meta.group, pt.size = 0.01, 
                split.by = split.group, ncol = nrow(cells), raster=FALSE) + NoLegend() + NoAxes() +
    labs(title = split.group, tag = subtitle) +
    theme(plot.title = element_text(color="black", size=14, face="bold", hjust=0.5, vjust=5),
          plot.tag = element_text(color="black", size=10, face="bold"),
          plot.tag.position = c(0.082, 0.80),# plot.tag.position = c(0.08, 0.815) for above the Subject line
          plot.margin = unit(c(0.5,0,0,0), "cm")
    )
  p2
  png(file=paste0(dir, file.name, "_", meta.group, "_umap_", split.group, ".png"), 
      width= nrow(cells)*2, height=3, res = 300, units = "in")
  print(p2)
  dev.off()
}


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
                         collapse = "                                                             "), 
                   paste(cells$numbers[7:12], 
                         collapse = "                                                             "), 
                   paste(cells$numbers[13:18], 
                         collapse = "                                                             "), 
                   paste(cells$numbers[19:nrow(cells)], 
                         collapse = "                                                             "), 
                   sep="\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n")
}

set.seed(100)
p3 <- DimPlot(seurat.data, label = FALSE, group.by = meta.group, 
              split.by = Sample.group, ncol = ifelse(nrow(cells) < 6, 1, 6), raster=FALSE) + NoLegend() + NoAxes() +
  labs(title = Sample.group, tag = subtitle) +
  theme(plot.title = element_text(color="black", size=14, face="bold", hjust=0.5, vjust=5),
        plot.tag = element_text(color="black", size=10, face="bold"),
        plot.tag.position = c(0.065, 0.790),# plot.tag.position = c(0.08, 0.815) for above the Subject line
        plot.margin = unit(c(0.5,0,0,0), "cm")
  )
p3
png(file=paste0(dir, file.name, "_", meta.group, "_umap_", Sample.group, ".png"), 
    width= ifelse(nrow(cells) < 6, nrow(cells)*2, 12), height=3*nrow2, res = 300, units = "in")
p3
dev.off()

rm(p,p2,p3)

# new.meta.data <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/","CellCycle.Senescence.stage.metadata.Rds"))  %>%  tibble::rownames_to_column(var = "rowname")
# metadata <- seurat.data@meta.data  %>% select(-one_of(c(intersect(colnames(seurat.data@meta.data), colnames(new.meta.data)) )))%>%  tibble::rownames_to_column(var = "rowname") 
# # intersect(colnames(metadata), colnames(new.meta.data))
# meta_data_merge <- metadata %>% left_join(new.meta.data, by = c("rowname"))  %>%  tibble::column_to_rownames("rowname")
# seurat.data <- AddMetaData( object = seurat.data, metadata = meta_data_merge)
# p <- DimPlot(seurat.data, label = FALSE, group.by = "Annotation", split.by = "Phase", pt.size = 0.01, raster=FALSE) + # NoLegend() +
#   labs(title = "Cell types in SAT", 
#        # tag = paste0("n = ", format(dim(seurat.data)[2], big.mark=",", scientific=FALSE) )
#        ) +
#   theme(plot.title = element_text(color="black", size=14, face="bold", hjust=0.5, vjust=0),
#         plot.tag = element_text(color="black", size=12, face="plain"),
#         plot.tag.position = c(0.2, 0.16),
#         plot.margin = unit(c(0,0,0,0), "cm"),
#         panel.border = element_rect(linetype = "solid", 
#                                     fill = NA, colour = "black"))
# p
# png(file=paste0(Disk, Project.folder, "/", figures.folder, "/", "Figure 1A_Overall_umap_CellCyclePhase.png"), width=13, height=4,
#     res = 300, units = "in")
# p
# dev.off()


# meta_data_merge <- meta_data_merge %>% mutate(Branch.in = ifelse(Branch %in% c("Branch_1", "Branch_2", "Branch_5", "Branch_7", 
#                                                                                "Branch_9", "Branch_10", "Branch_11", "Branch_12"), "Others", paste0(Branch)) )
# table(meta_data_merge$Branch.in)
# seurat.data <- AddMetaData( object = seurat.data, metadata = meta_data_merge)
# p <- DimPlot(seurat.data, label = FALSE, group.by = "Annotation", split.by = "Branch.in", pt.size = 0.01, raster=FALSE) + # NoLegend() +
#   labs(title = "Cell types in SAT", 
#        # tag = paste0("n = ", format(dim(seurat.data)[2], big.mark=",", scientific=FALSE) )
#        ) +
#   theme(plot.title = element_text(color="black", size=14, face="bold", hjust=0.5, vjust=0),
#         plot.tag = element_text(color="black", size=12, face="plain"),
#         plot.tag.position = c(0.2, 0.16),
#         plot.margin = unit(c(0,0,0,0), "cm"),
#         panel.border = element_rect(linetype = "solid", 
#                                     fill = NA, colour = "black"))
# p
# png(file=paste0(Disk, Project.folder, "/", figures.folder, "/", "Figure 1A_Overall_umap_Branch.png"), width=17, height=4,
#     res = 300, units = "in")
# p
# dev.off()
# 
# p <- DimPlot(seurat.data, label = FALSE, group.by = "Annotation", split.by = "Senescence.stage", pt.size = 0.01, raster=FALSE) + # NoLegend() +
#   labs(title = "Cell types in SAT", 
#        # tag = paste0("n = ", format(dim(seurat.data)[2], big.mark=",", scientific=FALSE) )
#   ) +
#   theme(plot.title = element_text(color="black", size=14, face="bold", hjust=0.5, vjust=0),
#         plot.tag = element_text(color="black", size=12, face="plain"),
#         plot.tag.position = c(0.2, 0.16),
#         plot.margin = unit(c(0,0,0,0), "cm"),
#         panel.border = element_rect(linetype = "solid", 
#                                     fill = NA, colour = "black"))
# p
# png(file=paste0(Disk, Project.folder, "/", figures.folder, "/", "Figure 1A_Overall_umap_Senescence.stage.png"), width=25, height=4,
#     res = 300, units = "in")
# p
# dev.off()
# 
# 
# set.seed(100)
# p <- DimPlot(seurat.data, label = FALSE, group.by = "Dataset", pt.size = 0.01, raster=FALSE)  +
#   labs(title = "Cell types in SAT", 
#        tag = paste0("n = ", format(dim(seurat.data)[2], big.mark=",", scientific=FALSE) )) +
#   theme(plot.title = element_text(color="black", size=14, face="bold", hjust=0.5, vjust=0),
#         plot.tag = element_text(color="black", size=12, face="plain"),
#         plot.tag.position = c(0.2, 0.16),
#         plot.margin = unit(c(0,0,0,0), "cm"),
#         panel.border = element_rect(linetype = "solid", 
#                                     fill = NA, colour = "black"))
# p
# # pdf(file=paste0("Figure 1X_umap.pdf"), width=5.2, height=4)
# p
# # dev.off()


# seurat.data@meta.data$Subjects <- factor(seurat.data@meta.data$Subjects, levels = paste("Subject", 1:length(unique(seurat.data@meta.data$Subjects))))




# seurat.data$Status <- factor(seurat.data$Status, levels = c("Lean", "Overweight", "Obese"))
# cells = as.data.frame(table(seurat.data$Status))
# cells$numbers = paste0("n = ", format(cells$Freq, big.mark=",", scientific=FALSE))
# subtitle = paste(paste(cells$numbers[1:3],  collapse = "                                  "),  
#                  sep="\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n")
# set.seed(100)
# p2 <- DimPlot(seurat.data, label = FALSE, group.by = "Annotation", pt.size = 0.01, 
#               split.by = "Status", ncol = 3, raster=FALSE) + NoLegend() + NoAxes() +
#   labs(title = "Groups in BMI", tag = subtitle) +
#   theme(plot.title = element_text(color="black", size=14, face="bold", hjust=0.5, vjust=5),
#         plot.tag = element_text(color="black", size=10, face="bold"),
#         plot.tag.position = c(0.11, 0.80),# plot.tag.position = c(0.08, 0.815) for above the Subject line
#         plot.margin = unit(c(0.5,0,0,0), "cm")
#   )
# p2
# png(file=paste0(Disk, Project.folder, "/", figures.folder, "/", "Figure 1A_Overall_umap_split_Status.png"), width=6, height=3, res = 300, units = "in")
# p2
# dev.off()


# cells = as.data.frame(table(seurat.data$Age_Group))
# cells$numbers = paste0("n = ", format(cells$Freq, big.mark=",", scientific=FALSE))
# subtitle = paste(paste(cells$numbers[1:2],  collapse = "                                  "),  
#                  sep="\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n")
# set.seed(100)
# p4 <- DimPlot(seurat.data, label = FALSE, group.by = "Annotation", pt.size = 0.01, 
#               split.by = "Age_Group", ncol = 2, raster=FALSE) + NoLegend() + NoAxes() +
#   labs(title = "Groups in Age", tag = subtitle) +
#   theme(plot.title = element_text(color="black", size=14, face="bold", hjust=0.5, vjust=5),
#         plot.tag = element_text(color="black", size=10, face="bold"),
#         plot.tag.position = c(0.165, 0.80),# plot.tag.position = c(0.08, 0.815) for above the Subject line
#         plot.margin = unit(c(0.5,0,0,0), "cm")
#   )
# p4
# png(file=paste0(Disk, Project.folder, "/", figures.folder, "/", "Figure 1A_Overall_umap_split_Age_Group.png"), width=4, height=3, res = 300, units = "in")
# p4
# dev.off()
