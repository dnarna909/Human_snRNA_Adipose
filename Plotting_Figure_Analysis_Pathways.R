### DISCLAIMER
#   Some of the algorithms are non-deterministic making the results slightly different from run to run.
# Many of the algorithms are frequently updated with new releases making the results slightly different from version to version.
# Depending on your system, this code below may not produce an exact 1:1 copy of the results. 
# This adversely affects clustering, filtering etc. if the code is naively copy-pasted.
# For better reproducibility, please download quality-filtered Rds objects and work on downstream analyses. 
# <br>
#   change "D:/" to "/media/jianie/Extreme SSD1/" # if in LUNIX computer

### Load libraries -------------------------------------------------------------------------------------------------------------
#rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.
memory.limit(size = 1e+13)

# source(paste0("/media/jianie/Extreme SSD1/2022-09-01 STARR Combine/Project Parameters.R"), local = knitr::knit_global())
source(paste0(Disk, Project.folder, "/", "Project Parameters.R"), local = knitr::knit_global())

for (pathways.file in pathways.files) {
  # pathways.file = pathways.files[2]
  rm(Result)
  Result <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/", pathways.file))  
  
  # prepare folder
  export.folder0 <- paste0("Pathways")
  dir.create(file.path(paste0(Disk, Project.folder, "/", figures.folder, "/"), export.folder0), showWarnings = FALSE)
  dir <- paste0(Disk, Project.folder, "/", figures.folder, "/", export.folder0, "/")
  
  for (nn in names(Result) ) {
    
    if (nrow(Result[[nn]]) >0 ) {
      # p.adjust or pvalue
      df <- Result[[nn]] %>% dplyr::filter( pvalue < 0.05
                                            # ,Description != "Cytoplasmic ribosomal proteins"
      ) %>% 
        mutate(y = ifelse(Direction == "Up", -log10(pvalue), 
                          ifelse(Direction == "Down", log10(pvalue), NA)),
               Direction = factor(Direction, levels = c("Up", "Down"))
        )
      
      
      
      ## plot ggboxplot for all data
      p <- ggbarplot(df, x = "Description", y = "y", 
                     fill = "Direction",               # change fill color by cyl
                     color = "white",            # Set bar border colors to white
                     palette = c("red3", "skyblue3", "lightgray"),            # jco journal color palett. see ?ggpar
                     width = 1,
                     #main = paste0("Top ", "significant pathways"), 
                     #font.main = c(14, "bold", "black"), 
                     xlab = FALSE,
                     font.x = c(14, "bold", "black"),
                     ylab = "-log10(Pvalue) for Induced / log10(Pvalue) for Repressed", # "-log10(AdjPvalue) for Induced / log10(AdjPvalue) for Repressed", 
                     # font.y = c(0, "bold", "white"), # error:Error in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  : polygon edge not found
                     # facet.by = "clusters", 
                     scales = "free",  ncol = 1,
                     #panel.labs = list(clusters = c("","","","","")),
                     
                     #sort.val = "asc",          # Sort the value in ascending order
                     #sort.by.groups = TRUE,     # Don't sort inside each group
                     
                     legend.title = "Group",
                     font.legend = c(20, "plain", "black"),
                     legend ="top",
                     # x.text.angle = 90,           # Rotate vertically x axis texts
                     rotate = TRUE,
                     #orientation = "horizontal",
                     ggtheme = theme_minimal() )  +
        geom_hline(yintercept = -log10(0.05), color = "black", linetype="dashed", size=0.5) + 
        geom_hline(yintercept = log10(0.05), color = "black", linetype="dashed", size=0.5) +
        border() +
        theme(axis.text.y = element_text(size = 15, colour = "black", face = "bold"),
              strip.background = element_blank(),
              strip.text.x = element_blank(),
              panel.spacing.y = unit(0, "lines")) 
      p
      png(filename = paste0(dir,  sub("_[^_]+$", "", stringr::str_split(pathways.file, "/")[[1]][2]), "_",
                            stringr::str_split(nn, ":")[[1]][1], ".", "bar.plot.png"), 
          width = 18, 
          height= nrow(df) * 0.22 + 1 , res = 300, units = "in")
      print(p)
      dev.off()
    }
  }
}
rm(Result, export.folder0, dir, nn, df,
   p)

