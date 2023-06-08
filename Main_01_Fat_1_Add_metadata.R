# DISCLAIMER
# Some of the algorithms are non-deterministic making the results slightly different from run to run.
# Many of the algorithms are frequently updated with new releases making the results slightly different from version to version.
# Depending on your system, this code below may not produce an exact 1:1 copy of the results. 
# This adversely affects clustering, filtering etc. if the code is naively copy-pasted.
# For better reproducibility, please download quality-filtered Rds objects and work on downstream analyses. 
# change "D:/" to "/media/jianie/Extreme SSD1/" # if in LUNIX computer

rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.
memory.limit(size = 1e+13)

source(paste0("/media/jianie/Extreme SSD1/2022-09-01 STARR_SGLT2 Combine/Project Parameters.R"), local = knitr::knit_global())
# source(paste0("D:/2022-09-01 STARR_SGLT2 Combine/Project Parameters.R"), local = knitr::knit_global())

# Load libraries ------------------------------------------------------------------------
# library("velocyto.R") # linux -specific

# Import data and subset ------------------------------------------------------------------------------------------------
eWAT <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/","SAT_Label.Rds"))
metadata <- eWAT@meta.data  %>%  tibble::rownames_to_column(var = "rowname")

# add samples information -------------------------------------------------------------------------------------------------------------------------------------------
## import data
samples.df <- readRDS(file = paste0(sample.info.disk, sample.info.path, "all.patients.rds")) 
colnames(samples.df)
select.cols <- c("Dataset", "sample_id", "SingleNuclei", "record_id", "recruited_from", 
                 "Gender", 
                 "Age" ,"Age_Group",
                 "Height (inches)" , "Weight (pounds)", "BMI", "Status",
                 "Glucose" ,  "A1c Value", "fasting", # "dbc4_diabetes",
                 "Estimated Number of Cells",  "Median Genes per Cell", "Sequencing Saturation", "Median UMI Counts per Cell")
samples.df <- samples.df %>% select(all_of(select.cols)) # %>% mutate(Replicate = "R1")
samples.df.metadata <- metadata  %>% dplyr::select(one_of(c("rowname", "Dataset"))) %>% full_join(samples.df, by = "Dataset") # %>%  tibble::column_to_rownames("rowname")

# add annotation -----------------------------------------------------------------------------------------------------------------------------------------
Annoation.metadata <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/","Annotation.metadata.Rds"))

# add subtype
Subtype.data <- as.data.frame()
for (nn in unique(Annoation.metadata$Annoation) ) {
  if(file.exists(paste0(Disk, Project.folder, "/", Rds.folder, "/", nn, ".Subtype.metadata.Rds"))) {
  metadata.df <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/", nn, ".Subtype.metadata.Rds")) 
  Subtype.data <- rbind(Subtype.data, metadata.df)
  print(nn)
  } 
}
Subtype.metadata <- metadata %>% dplyr::select(one_of(c("rowname"))) %>% full_join(Subtype.data, by = "rowname")

# combine all data ---------------------------------------------------------------------------------------------------------------------------------------
all.metadata <- metadata %>% 
  full_join(samples.df.metadata , by = "rowname")  %>% 
  full_join(Annoation.metadata, by = "rowname") %>%
  full_join(Subtype.metadata, by = "rowname") %>%  
  tibble::column_to_rownames("rowname")
saveRDS(all.metadata, paste0(Disk, Project.folder, "/", Rds.folder, "/","Annotation.samples.Subtype.metadata.Rds"))

# add to meta data
eWAT <- AddMetaData( object = eWAT, metadata = all.metadata)
table(eWAT$Dataset)
table(eWAT$Status)
table(eWAT$Age_Group)
table(eWAT$Dataset)


## Import data and subset ------------------------------------------------------------------------------------------------
if(file.exists(paste0(Disk, Project.folder, "/", Rds.folder, "/","SAT_Annotated_Subtype.Rds"))) {eWAT <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/","SAT_Annotated_Subtype.Rds")) } else {
  eWAT <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/","SAT_Annotated.Rds"))
} 
metadata <- eWAT@meta.data

# import senescence meta data
metadata <- eWAT@meta.data  %>%  tibble::rownames_to_column(var = "rowname")
sens.meta.data <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/","sens.meta.data.Rds"))  %>%  tibble::rownames_to_column(var = "rowname")
new.meta.data <- sens.meta.data 
intersect(colnames(metadata ), colnames(new.meta.data))
meta_data_merge <- metadata %>% full_join(new.meta.data, by = c(intersect(colnames(metadata), colnames(new.meta.data)))) 

# import sens.meta.data.PHATE.Pt
sens.meta.data.PHATE.Pt <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/","sens.meta.data.PHATE.Pt.Rds"))  %>%  tibble::rownames_to_column(var = "rowname")
new.meta.data <- sens.meta.data.PHATE.Pt
intersect(colnames(meta_data_merge), colnames(new.meta.data ) )
meta_data_merge <- meta_data_merge %>% full_join(new.meta.data, by = c(intersect(colnames(meta_data_merge), colnames(new.meta.data) ))) 

# final
meta_data_merge <- meta_data_merge %>%  tibble::column_to_rownames("rowname")
eWAT <- AddMetaData( object = eWAT, metadata = meta_data_merge)

#
FeaturePlot(eWAT, label = F, raster=FALSE, features = "SenMayo")
FeaturePlot(eWAT, label = F, raster=FALSE, features = c("SenMayo", "Senescence.GO.score"), blend = T, order = T, min.cutoff = 0)
DimPlot(eWAT, label = T, raster=FALSE, group.by = "Subtype")
DimPlot(eWAT, label = F, raster=FALSE, group.by = "Senescence.stage")
DimPlot(eWAT, label = F, raster=FALSE, group.by = "Annotation", split.by = "Senescence.stage")
FeatureScatter(eWAT,raster=FALSE, feature1 = c("SenMayo"), feature2 = c("Senescence.GO.score")) 

# save picture 
png(filename = paste0(Disk, Project.folder, "/", figures.folder, "/", "SenMayo_Senescence.GO.score.FeatureScatter.png"), 
    width = 8, height = 6, res = 300, unit="in")
FeatureScatter(eWAT,raster=FALSE, feature1 = c("SenMayo"), feature2 = c("Senescence.GO.score")) 
dev.off()

png(filename = paste0(Disk, Project.folder, "/", figures.folder, "/", "Senescence.stage.DimPlot.png"), 
    width = 8, height = 6, res = 300, unit="in")
DimPlot(eWAT, label = F, raster=FALSE, group.by = "Senescence.stage")
dev.off()

png(filename = paste0(Disk, Project.folder, "/", figures.folder, "/", "Senescence.stage.Annotation.DimPlot.png"), 
    width = 24, height = 4, res = 300, unit="in")
DimPlot(eWAT, label = F, raster=FALSE, group.by = "Annotation", split.by = "Senescence.stage")
dev.off()

# Save results
saveRDS(eWAT, paste0(Disk, Project.folder, "/", Rds.folder, "/","SAT_Annotated_Subtype.Rds")) ## This object is downloadable from Open Science Framework with additional annotations added by subclustering each major cell type

