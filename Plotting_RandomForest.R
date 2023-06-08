

# rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.
memory.limit(size = 1e+13)

# import parameters and packages --------------------------------------------------------------------------------------------
# source(paste0("/media/jianie/Extreme SSD1/2022-09-01 STARR_SGLT2 Combine/Project Parameters.R"), local = knitr::knit_global())
source(paste0(paste0(Disk, Project.folder, "/", "Project Parameters.R")), local = knitr::knit_global())
library("RColorBrewer")
library("ggpubr")


# Import data -----------------------------------------------------------------------------------------------------------------------------------------------------
dir0 <- paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder, "/")
files <- list.files(dir0)

dir.create(file.path(paste0(Disk, Project.folder, "/", figures.folder, "/", export.folder)), showWarnings = FALSE)
dir1 <- paste0(Disk, Project.folder, "/", figures.folder, "/", export.folder, "/")

for (file in files) {
  # file = files[1]
  load(paste0(dir0, file))
  
  file.name <- stringr::str_split(file, ".Outcomes.ML.RData")[[1]][1]
  
  for (cc in names(Test.list)){
    for (tt in names(Test.list[[cc]]) ) {
      if ( tt %in% c("randomForest", "randomForest2", "bagging")){
        print(paste(cc, tt))
      varImpPlot(Test.list[[cc]][[tt]])
      png(filename = paste0(dir1,  file.name , "_", cc, "_", tt, ".png"), 
          width = 12, 
          height = 8, res = 300, units = "in")
      varImpPlot(Test.list[[cc]][[tt]])
      dev.off()
      print(paste(cc, tt))
      }
    }
  }
}

rm(data, Import.list, mse.df, Test.list, 
   cc,tt, dir0, dir1, file, files,
   file.name)
