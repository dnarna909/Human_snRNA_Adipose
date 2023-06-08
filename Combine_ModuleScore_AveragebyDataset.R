
# rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.
memory.limit(size = 1e+13)

# import parameters and packages --------------------------------------------------------------------------------------------
# source(paste0("/media/jianie/Extreme SSD1/2022-09-01 STARR_SGLT2 Combine/Project Parameters.R"), local = knitr::knit_global())
source(paste0(paste0(Disk, Project.folder, "/", "Project Parameters.R")), local = knitr::knit_global())
library("RColorBrewer")
library("ggpubr")

# Import data -----------------------------------------------------------------------------------------------------------------------------------------------------
# export.folder <- "Compositions"
# dir.create(file.path(paste0(Disk, Project.folder, "/", Rds.folder, "/"), export.folder), showWarnings = FALSE)


metadata.ave <- list()
for (seurat.file in seurat.files) {
  # seurat.file = seurat.files[1]
  print(seurat.file)
  file.name <- stringr::str_split(seurat.file, ".Rds")[[1]][1];
  print(file.name )
  Annotation.file.name <- paste(stringr::str_split(file.name, "_")[[1]][c(1,2,3)], collapse = "_");
  print(Annotation.file.name)
  Subtype.file.name <- stringr::str_split(file.name, paste0(Annotation.file.name, "_"))[[1]][2]
  if (is.na(Subtype.file.name)) {Subtype.file.name = ""}
  print(Subtype.file.name)
  if( length(stringr::str_split(Subtype.file.name, paste0("_"))[[1]]) > 1){
    celltype.name <- stringr::str_split(Subtype.file.name, paste0("_"))[[1]][1]
    #label.name <- stringr::str_split(Subtype.file.name, paste0("_"))[[1]][2]
  } else {
    celltype.name <- Annotation.file.name
    # label.name <- stringr::str_split(Subtype.file.name, paste0("_"))[[1]][1]
  }
  print(celltype.name)
  
  # Import sample data
  sample.df <- readRDS(paste0(Disk, Project.folder, "/",  Rds.folder, "/", "Samples.df.Rds")) 
  
  # Import metadata
  sample.meta <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/", 
                                stringr::str_split(file.name, "_")[[1]][1], "_Sample.meta.data.Rds")) %>% 
    tibble::rownames_to_column(var = "rowname")
  
  Subtype.meta <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/", 
                                 file.name, "_", file.type)) %>% 
    tibble::rownames_to_column(var = "rowname")
  
  Subtype.meta1 <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/",  
                                  file.type1)) %>% 
    tibble::rownames_to_column(var = "rowname")
  
  metadata.all <- left_join(Subtype.meta, sample.meta, by = "rowname") %>% 
    left_join(Subtype.meta1, by = "rowname")
  
  row.groups.new <- paste(row.groups, col.groups, sep = "_")
  for (row.group.new in row.groups.new) {
    metadata.all[[row.group.new]] <-
      paste(metadata.all[[stringr::str_split(row.group.new, "_")[[1]][1]]], 
            metadata.all[[stringr::str_split(row.group.new, "_")[[1]][2]]], sep = "_")
  }
  
  for (col.group in col.groups) {
    # col.group = col.groups[1]
    for (col.group1 in col.groups1) { 
      # col.group1 = col.groups1[2]
      
      all.ave <- metadata.all %>% 
        arrange(vars(one_of(c( col.group, "Dataset")))) %>% 
        group_by_at(vars(one_of(c(col.group, "Dataset")))) %>%
        dplyr::filter(get(col.group1) > 0) %>% # , !is.na(get(col.group1))
        summarise(Avg = mean(get(col.group1), na.rm = TRUE), .groups = 'drop') %>% 
        spread(as.name(col.group), Avg) %>% 
        tibble::column_to_rownames(var = "Dataset")
      colnames(all.ave) <- paste0(col.group1, "_", colnames(all.ave))
      
      metadata.all.ave <- all.ave %>% 
        tibble::rownames_to_column(var = "Dataset") %>% 
        left_join(sample.df, by = "Dataset") %>%
        dplyr::filter(Dataset %in% unique(metadata.all$Dataset))
      
      all.ave <- all.ave %>% tibble::rownames_to_column(var = "Dataset")
      metadata.ave[[paste0(file.name, "_", col.group, "_" , col.group1)]] <- all.ave
      
      print(paste(col.group, col.group1, sep = ":"))
    }
  }
}
score.average <- Reduce(full_join, metadata.ave)
score.average <- score.average[,colSums(is.na(score.average))<nrow(score.average)] # Remove columns from dataframe where ALL values are NA
  
saveRDS(score.average, file = paste0(Disk, Project.folder, "/", Rds.folder, "/", Annotation.file.name, "_" ,stringr::str_split(file.type1, "_")[[1]][1], "_", "ModuleScoreAve.df.Rds") )


rm(dir1,  seurat.file , df, 
   file.name, Annotation.file.name, Subtype.file.name,
   sample.df,  sample.meta, Subtype.meta, Subtype.meta1, 
   metadata.all, row.groups.new, row.group.new, 
   col.group, col.group1,
   all.ave, metadata.all.ave, metadata.ave, score.average,
   celltype.name 
   
)
gc() #free up memrory and report the memory usage.








