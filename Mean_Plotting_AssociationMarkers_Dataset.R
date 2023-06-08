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

# prepare folder
# export.folder <- paste0("AssociationMarkers")
dir.create(file.path(paste0(Disk, Project.folder, "/", figures.folder, "/"), export.folder), showWarnings = FALSE)
dir <- paste0(Disk, Project.folder, "/", figures.folder, "/", export.folder, "/");dir

for (markers.file in markers.files) {
  # markers.file = markers.files[1]; markers.file
  load(paste0(Disk, Project.folder, "/", Rds.folder, "/", markers.file))  
  names(comb.df); names(comb.ls);
  
  # prepare folder
  export.folder0 <- stringr::str_split(markers.file, "/")[[1]] ;export.folder0
  dir1 <- paste0(Disk, Project.folder, "/", Rds.folder, "/")
  for (ee in 1:(length(export.folder0) -1) ) {
    dir.create(file.path(dir1, export.folder0[ee]), showWarnings = FALSE)
    dir2 <- paste0(dir1, export.folder0[ee], "/")
  }
  dir2
  
  # keep only genes in more than 2 people
  Result <- list()
  Result.comb <- list()

  if ( length(grep("meanLogFC", colnames(comb.ls[[1]]), ignore.case =TRUE)) > 0 ){
    FC.name ="meanLogFC"} else {
      FC.name = "fcMedian"}
  FC.name;
  
  Result.comb <- lapply(comb.ls, function(x) {
    if (length(grep(FC.name, colnames(x), ignore.case =TRUE)) > 1 ) {
      x <- x 
    } 
  })
  Result.comb <- Filter(Negate(is.null), Result.comb)
  
  # select significant genes 
  Result <- lapply(Result.comb, function(x) {
    bind_rows( x %>%  as.data.frame() %>% 
                 tibble::rownames_to_column(var = "ID") %>%  
                 dplyr::select(-ID) %>%  
                 tibble::column_to_rownames(var = "rowname") %>%  
                 filter_at(vars(contains("FDR") ), all_vars(. <= 0.05)) %>% # "pvalue"
                 filter_at(vars(contains(FC.name)), all_vars(. > 0 )) %>% 
                 tibble::rownames_to_column(var = "rowname")  %>%  
                 mutate( Direction = "Up"),
               x %>%  as.data.frame() %>% 
                 tibble::rownames_to_column(var = "ID") %>%  
                 dplyr::select(-ID) %>%  
                 tibble::column_to_rownames(var = "rowname") %>%  
                 filter_at(vars(contains("FDR")  ), all_vars(. <= 0.05)) %>% # "pvalue"
                 filter_at(vars(contains(FC.name)), all_vars(. < 0 )) %>% 
                 tibble::rownames_to_column(var = "rowname")  %>%  
                 mutate( Direction = "Down") 
    ) 
  })
  Result <- lapply(Result, function(x) {
    if (nrow(x) != 0 ) {
      x <- x 
    } 
  })
  Result <- Filter(Negate(is.null), Result)
  
  Result <- lapply(Result, function(x) {
    x  %>% 
      mutate(  
        Mean_waldSTAT = rowMeans(dplyr::select(., contains("waldStat")), na.rm = TRUE),
        SD_waldSTAT = rowSds(as.matrix(dplyr::select(., contains("waldStat")), na.rm = TRUE)),
        SEM_waldSTAT = rowSds(as.matrix(dplyr::select(., contains("waldStat")), na.rm = TRUE))/sqrt(length(grep("waldStat", colnames(x), ignore.case =TRUE))),
        
        Mean_AVG_logFC = rowMeans(dplyr::select(., contains(FC.name)), na.rm = TRUE),
        SD_AVG_logFC = rowSds(as.matrix(dplyr::select(., contains(FC.name)), na.rm = TRUE)),
        SEM_AVG_logFC = rowSds(as.matrix(dplyr::select(., contains(FC.name)), na.rm = TRUE))/sqrt(length(grep(FC.name, colnames(x), ignore.case =TRUE))),
        
        Mean_P_Val = rowMeans(dplyr::select(., contains("pvalue") ), na.rm = TRUE),# 
        SD_P_Val = rowSds(as.matrix(dplyr::select(., contains("pvalue") ), na.rm = TRUE)),
        SEM_P_Val = rowSds(as.matrix(dplyr::select(., contains("pvalue") ), na.rm = TRUE))/sqrt(length(grep("pvalue", colnames(x), ignore.case =TRUE))),
        
        Mean_FDR = rowMeans(dplyr::select(., contains("FDR") ), na.rm = TRUE),# 
        SD_FDR = rowSds(as.matrix(dplyr::select(., contains("FDR") ), na.rm = TRUE)),
        SEM_FDR = rowSds(as.matrix(dplyr::select(., contains("FDR") ), na.rm = TRUE))/sqrt(length(grep("FDR", colnames(x), ignore.case =TRUE))),
        
      ) %>%
      arrange(-Mean_AVG_logFC) %>%
      mutate( Rank.ID.logFC = 1:n()) %>%
      arrange(-Mean_waldSTAT) %>%
      mutate( Rank.ID.waldSTAT = 1:n()) %>%
      arrange(Mean_FDR) %>%
      mutate( Rank.ID.FDR = 1:n()) %>%
      mutate( Rank.ID = Rank.ID.FDR + Rank.ID.logFC + Rank.ID.waldSTAT) %>%
      arrange(Rank.ID)
  })
  
  for (nn in names(Result)) {
    if (nrow(Result[[nn]]) > 0 ) {
      # nn = names(Result)[1];nn
      df <- Result[[nn]]
      n = 10
      labels = unique(c(df %>% dplyr::filter(Direction == "Up") %>% dplyr::top_n(-n, Rank.ID) %>% pull(rowname),
                        df %>% dplyr::filter(Direction == "Down") %>% dplyr::top_n(n, Rank.ID) %>% pull(rowname)));labels
      # labels = unique(c(df %>% dplyr::filter(Direction == "Up", mLog10FDR == Inf) %>% pull(rowname),
      #                   df %>% dplyr::filter(Direction == "Down", mLog10FDR == Inf) %>% pull(rowname)));labels; length(labels)
      
      df <- df %>% mutate(
        Plotname = ifelse(rowname %in% labels, rowname, ""),
        Direction = factor(Direction, levels = c("Up", "Down", "non")),
        mLog10FDR = -log10(Mean_FDR)
      )
      options(ggrepel.max.overlaps = Inf) 
      p <- ggscatter(df, x = "Mean_AVG_logFC", y = "mLog10FDR", 
                     color = "Direction",palette = c("red3", "skyblue3", "lightgray"),
                     # color = "Sig.overall", palette = c("red3",  "orange", "skyblue3","lightgray", "snow"),
                     #label = "Row.names", repel = TRUE, label.select = labels,
                     shape = 20, 
                     size = ifelse(df[["Direction"]] == "non", 1, 3), # Points color, shape and size
                     xlab = "avg_log2FoldChange", 
                     ylab = "-log10(FDR)",
                     title = paste0(nn)
      ) +
        # ylim(0, max(df$mLog10FDR)) +
        # geom_abline(intercept = 0, slope = 1, color = "darkgray", linetype="solid", linewidth=1) +
        geom_vline(xintercept = 0, color = "black", linetype="dashed", linewidth=0.2) +
        # geom_hline(yintercept = 0, color = "gray", linetype="dashed", linewidth=0.2) 
        theme_bw() +
        theme(legend.position = "none" ) +
        # geom_text(aes(label = .data[["Plotname"]]), size = 4, color ="black") +
        # ggrepel::geom_text_repel(aes(label = .data[["Plotname"]]),
        #                          min.segment.length = -0, seed = 42, box.padding = 0.8
        #                          , max.overlaps = Inf
        # )  +
        geom_errorbar(data = df,
                      aes(xmin = Mean_AVG_logFC - SEM_AVG_logFC, xmax = Mean_AVG_logFC + SEM_AVG_logFC),
                      width = 0.2, color= "black",
                      # position = position_dodge(0.05)
        )+
        ggrepel::geom_text_repel(aes(label = .data[["Plotname"]]), size= 3,
                                 nudge_x = .15,
                                 box.padding = 0.5,
                                 nudge_y = 1,
                                 segment.curvature = -0.1,
                                 segment.ncp = 1,
                                 segment.angle = 1
        ) 
      print(p)
      png(filename = paste0(dir, 
                            sub("_[^_]+$", "", stringr::str_split(markers.file, "/")[[1]][2]), 
                            "_", paste0(nn), "_Pt_logFC_ggscatter.png"), 
          width = 6, 
          height= 6 , res = 300, units = "in")
      print(p)
      dev.off()
      
      
      p <- ggscatter(df, x = "Mean_waldSTAT", y = "mLog10FDR", 
                     color = "Direction",palette = c("red3", "skyblue3", "lightgray"),
                     # color = "Sig.overall", palette = c("red3",  "orange", "skyblue3","lightgray", "snow"),
                     #label = "Row.names", repel = TRUE, label.select = labels,
                     shape = 20, 
                     size = ifelse(df[["Direction"]] == "non", 1, 3), # Points color, shape and size
                     xlab = "Mean of wald statistic", 
                     ylab = "-log10(FDR)",
                     title = paste0(nn)
      ) +
        # geom_abline(intercept = 0, slope = 1, color = "darkgray", linetype="solid", linewidth=1) +
        geom_vline(xintercept = 0, color = "black", linetype="dashed", linewidth=0.2) +
        # geom_hline(yintercept = 0, color = "gray", linetype="dashed", linewidth=0.2) 
        theme_bw() +
        theme(legend.position = "none" ) +
        # geom_text(aes(label = .data[["Plotname"]]), size = 4, color ="black") +
        # ggrepel::geom_text_repel(aes(label = .data[["Plotname"]]),
        #                          min.segment.length = -0, seed = 42, box.padding = 0.8
        #                          , max.overlaps = Inf
        # )  
        geom_errorbar(data = df,
                      aes(xmin = Mean_waldSTAT - SEM_waldSTAT, xmax = Mean_waldSTAT + SEM_waldSTAT),
                      width = 0.2, color= "black",
                      # position = position_dodge(0.05)
        )+
        ggrepel::geom_text_repel(aes(label = .data[["Plotname"]]),size= 3,
                                 nudge_x = .15,
                                 box.padding = 0.5,
                                 nudge_y = 1,
                                 segment.curvature = -0.1,
                                 segment.ncp = 1,
                                 segment.angle = 1
        )
      print(p)
      png(filename = paste0(dir, 
                            sub("_[^_]+$", "", stringr::str_split(markers.file, "/")[[1]][2]), "_",
                            paste0(nn), "_Pt_waldSTAT_ggscatter.png"), 
          width = 6, 
          height= 6 , res = 300, units = "in")
      print(p)
      dev.off()
    }
    
    
  }
  
  saveRDS(Result, 
          paste0(Disk, Project.folder, "/", Rds.folder, "/", stringr::str_split(markers.file, ".RData")[[1]][1],
                 ".Mean.Rds") 
  )
  print(markers.file)
}
rm(comb.df, comb.ls, Result, export.folder0, dir, nn, df, 
   markers.file, p, labels, n , Result.comb,
   ee, ll, export.folder
)
gc()
