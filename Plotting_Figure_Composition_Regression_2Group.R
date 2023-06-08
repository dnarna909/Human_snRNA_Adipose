
# rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.
memory.limit(size = 1e+13)

# import parameters and packages --------------------------------------------------------------------------------------------
# source(paste0("/media/jianie/Extreme SSD1/2022-09-01 STARR_SGLT2 Combine/Project Parameters.R"), local = knitr::knit_global())
source(paste0(paste0(Disk, Project.folder, "/", "Project Parameters.R")), local = knitr::knit_global())
source(paste0(paste0(Disk, "00_Functions_Refs/", "Functions_Human Fat snRNA.R")), local = knitr::knit_global() )
library("RColorBrewer")
library("ggpubr")
library(rstatix)
# Import data -----------------------------------------------------------------------------------------------------------------------------------------------------
# export.folder <- "Compositions"
# dir.create(file.path(paste0(Disk, Project.folder, "/", Rds.folder, "/"), export.folder), showWarnings = FALSE)
dir0 <- paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder, "/")
files <- list.files(dir0)

dir.create(file.path(paste0(Disk, Project.folder, "/", figures.folder, "/", export.folder, ".Regre.plotting")), showWarnings = FALSE)
dir.create(file.path(paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder, ".VariableSelection", "/")), showWarnings = FALSE)
dir1 <- paste0(Disk, Project.folder, "/", figures.folder, "/", export.folder, ".Regre.plotting", "/")


for (file in files) {
  # file = files[1]
  load(paste0(dir0, file))
  
  # prepare file name
  file.name <- stringr::str_split(file, ".metadata_Compositions.RData")[[1]][1]
  Annotation.file.name <- paste(stringr::str_split(file.name, "_")[[1]][c(1,2,3)], collapse = "_")
  Subtype.file.name <- stringr::str_split(file.name, paste0(Annotation.file.name, "_"))[[1]][2]
  if( length(stringr::str_split(Subtype.file.name, paste0("_"))[[1]]) > 1){
    celltype.name <- stringr::str_split(Subtype.file.name, paste0("_"))[[1]][1]
    #label.name <- stringr::str_split(Subtype.file.name, paste0("_"))[[1]][2]
  } else {
    celltype.name <- Annotation.file.name
    # label.name <- stringr::str_split(Subtype.file.name, paste0("_"))[[1]][1]
  }
  
  # import data
  Multi.Regre.Result.list <- list()
  for (list.name in names(Composition.lists)) {
    # list.name = names(Composition.lists)[1]
    label.name <- stringr::str_split(list.name, paste0("_"))[[1]][length(stringr::str_split(list.name, paste0("_"))[[1]])] # the last one in list name
    x.name <- stringr::str_split(list.name, paste0("_"))[[1]][2] # the second one in list name
    Composition_3 <- Composition.lists[[list.name]][["Composition_3All"]] %>% 
      dplyr::filter(
        !!sym(Select_Group) %in% Select_Group.names #"Middle_Lean", "Middle_Overweight" 
      )  
    Composition_3[[Select_Group]] <- factor(Composition_3[[Select_Group]], levels = Select_Group.names) # "Middle_Lean", "Middle_Overweight"
    Composition_2 <- Composition.lists[[list.name]][["Composition_2All"]] %>% dplyr::filter( !!sym("Dataset") %in% Composition_3[["Dataset"]] ) %>% 
      select_if(~sum(!is.na(.)) > 0) %>%
      select(which(!plyr::numcolwise(sum)(., na.rm=TRUE) %in% 0))
    Composition <- Composition.lists[[list.name]][["Composition_All"]] %>% tibble::rownames_to_column(var = "rowname") %>%
      dplyr::filter( !!sym("rowname") %in% Composition_3[["Dataset"]] ) %>% tibble::column_to_rownames(var = "rowname") %>% 
      select(which(!colSums(., na.rm=TRUE) %in% 0))
    clusters = colnames(Composition)
    # display.brewer.pal(length(clusters), "Set1")
    # Composition_3$x.label <- factor(Composition_3[[label.name]], levels=clusters)


# plot regression figure
library(ggplot2)
    df <- Composition_2
p <- ggplot(df, aes(x= .data[[column1]], y=.data[[column2]], group=.data[[Select_Group]])) +
  geom_point(aes(shape=.data[[Select_Group]], color=.data[[Select_Group]]))
print(p)
png(file=paste0(dir1, file.name, "_", x.name, "_", "Composition_", analysis.group , ".ggpairs.png"), 
    width=ncol(Composition), height=ncol(Composition), res = 300, units = "in")
print(p)
dev.off()


