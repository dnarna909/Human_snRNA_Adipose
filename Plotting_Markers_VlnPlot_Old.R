# vlnplot for all markers
# https://github.com/ycl6/StackedVlnPlot
library(Seurat)
library(ggplot2)
library(cowplot)

# Load data.frame obj
features <- c("Adgre1", "Fcgr1", "Cd68", "Col1a1", "Gpm6b", "Svep1", "Cdh5", "Pecam1",
              "Flt1", "Cd3d", "Trbc2", 
              "Cd79a", "Ms4a1", "Cd79b", "Msln", "Lrrn4",
              "Upk3b", "Pdgfrb", "Cspg4", "Hdc", "Cd274", "Sept3", "Xcr1",  "Ptprc", "Ncr1", "Gzma")

features.df <- rio::import("/media/jianie/Extreme SSD/2022-09-01 STARR_SGLT2 Combine/Fat_Subtype.xlsx", sheet = "Annotation_Markers")
features <-features.df$gene
seurat.file <- readRDS("/media/jianie/Extreme SSD/2022-09-01 STARR_SGLT2 Combine/Rds files_Aggr_all/SAT_Label_2000.Rds")
metadata <- readRDS("/media/jianie/Extreme SSD/2022-09-01 STARR_SGLT2 Combine/Rds files_Aggr_all/SAT_Label_2000_Subtype.metadata.Rds")
dim(seurat.file)
seurat.file <- AddMetaData(seurat.file, metadata = metadata)
group.x = "Annotation"

library(dplyr)
df <-  data.frame(Cell = rownames(seurat.file@meta.data),
                  Idents = seurat.file@meta.data[[group.x]] )
pbmc = as.data.frame(t(as.matrix(seurat.file@assays[["originalexp"]]@data[features, ]))) %>% 
  tibble::rownames_to_column(var = "Cell")
pbmc <- pbmc %>%
  left_join( 
    df, 
    by ="Cell"
  )
head(pbmc)[,1:6]                
head(pbmc)               

# Use melt to change data.frame format
pbmc <- reshape2::melt(pbmc, id.vars = c("Cell","Idents"), measure.vars = features,
                       variable.name = "Feat", value.name = "Expr")
head(pbmc, 10)

# Identity on x-axis
a <- ggplot(pbmc, aes(factor(Idents), Expr, fill = Idents)) +
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
  # ggtitle("Identity on x-axis") + 
  # xlab("Identity") + 
  ylab("Expression Level")
print(a)

png(filename = paste0(path1, figures.folder, "/", "All.markers_", group.x, "_VlnPlot.All.png"),
    width = 10, height = 8, units = "in", res = 300)
print(a)
dev.off()