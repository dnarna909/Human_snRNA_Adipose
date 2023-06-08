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
export.folder <- paste0("Markers")
dir.create(file.path(paste0(Disk, Project.folder, "/", figures.folder, "/"), export.folder), showWarnings = FALSE)
dir <- paste0(Disk, Project.folder, "/", figures.folder, "/", export.folder, "/");dir

for (markers.file in markers.files) {
  # markers.file = markers.files[3]; markers.file
  load(paste0(Disk, Project.folder, "/", Rds.folder, "/", markers.file))  
  names(Result.comb.df); names(Result.comb.list);
  
  # prepare folder
  export.folder0 <- stringr::str_split(markers.file, "/")[[1]] ;export.folder0
  dir1 <- paste0(Disk, Project.folder, "/", Rds.folder, "/")
  for (ee in 1:(length(export.folder0) -1) ) {
    dir.create(file.path(dir1, export.folder0[ee]), showWarnings = FALSE)
    dir2 <- paste0(dir1, export.folder0[ee], "/")
  }
  dir2
  
  # keep only genes in more than 2 people
  Result.sig <- list()
  Result <- list()
  Result.comb <- list()
  for (ll in names(Result.comb.list)) {
    # ll = names(Result.comb.list)[1];ll
    Result.comb[[ll]] <- lapply(Result.comb.list[[ll]], function(x) {
      if (length(grep("avg_log2FC", colnames(x), ignore.case =TRUE)) > 1 ) {
        x <- x 
      } 
    })
    Result.comb[[ll]] <- Filter(Negate(is.null), Result.comb[[ll]])
    
    # select significant genes 
    Result.sig[[ll]] <- lapply(Result.comb[[ll]], function(x) {
      bind_rows( x %>%  as.data.frame() %>% 
                   tibble::rownames_to_column(var = "ID") %>%  
                   dplyr::select(-ID) %>%  
                   tibble::column_to_rownames(var = "rowname") %>%  
                   filter_at(vars(contains("p_val") & !contains("p_val_adj") ), all_vars(. < 0.05)) %>% # 
                   filter_at(vars(contains("avg_log2FC")), all_vars(. > 0 )) %>% 
                   tibble::rownames_to_column(var = "rowname")  %>%  
                   mutate( Direction = "Up"),
                 x %>%  as.data.frame() %>% 
                   tibble::rownames_to_column(var = "ID") %>%  
                   dplyr::select(-ID) %>%  
                   tibble::column_to_rownames(var = "rowname") %>%  
                   filter_at(vars(contains("p_val") & !contains("p_val_adj") ), all_vars(. < 0.05)) %>% # 
                   filter_at(vars(contains("avg_log2FC")), all_vars(. < 0 )) %>% 
                   tibble::rownames_to_column(var = "rowname")  %>%  
                   mutate( Direction = "Down") 
      ) 
    })
    Result.sig[[ll]] <- Filter(Negate(is.null), Result.sig[[ll]])
    
    
    Result[[ll]] <- lapply(Result.comb[[ll]], function(x) {
      x  %>% 
        mutate(  
          Mean_AVG_log2FC = rowMeans(dplyr::select(., contains("avg_log2FC")), na.rm = TRUE),
          SD_AVG_log2FC = rowSds(as.matrix(dplyr::select(., contains("avg_log2FC")), na.rm = TRUE)),
          SEM_AVG_log2FC = rowSds(as.matrix(dplyr::select(., contains("avg_log2FC")), na.rm = TRUE))/sqrt(length(grep("avg_log2FC", colnames(x), ignore.case =TRUE))),
          
          Mean_PCT.dif = rowMeans(dplyr::select(., contains("pct.dif")), na.rm = TRUE),
          SD_PCT.dif = rowSds(as.matrix(dplyr::select(., contains("pct.dif")), na.rm = TRUE)),
          SEM_PCT.dif = rowSds(as.matrix(dplyr::select(., contains("pct.dif")), na.rm = TRUE))/sqrt(length(grep("pct.dif", colnames(x), ignore.case =TRUE))),
          
          Mean_PCT.1 = rowMeans(dplyr::select(., contains("pct.1")), na.rm = TRUE),
          SD_PCT.1 = rowSds(as.matrix(dplyr::select(., contains("pct.1")), na.rm = TRUE)),
          SEM_PCT.1 = rowSds(as.matrix(dplyr::select(., contains("pct.1")), na.rm = TRUE))/sqrt(length(grep("pct.1", colnames(x), ignore.case =TRUE))),
          
          Mean_PCT.2 = rowMeans(dplyr::select(., contains("pct.2")), na.rm = TRUE),
          SD_PCT.2 = rowSds(as.matrix(dplyr::select(., contains("pct.2")), na.rm = TRUE)),
          SEM_PCT.2 = rowSds(as.matrix(dplyr::select(., contains("pct.2")), na.rm = TRUE))/sqrt(length(grep("pct.2", colnames(x), ignore.case =TRUE))),
          
          Mean_P_Val = rowMeans(dplyr::select(., contains("p_val")&!contains("p_val_adj") ), na.rm = TRUE),# 
          SD_P_Val = rowSds(as.matrix(dplyr::select(., contains("p_val")&!contains("p_val_adj")  ), na.rm = TRUE)),
          SEM_P_Val = rowSds(as.matrix(dplyr::select(., contains("p_val")&!contains("p_val_adj") ), na.rm = TRUE))/sqrt(length(grep("avg_log2FC", colnames(x), ignore.case =TRUE))),
          
        ) %>%
        mutate(Mean_mlog10_P = -log10(Mean_P_Val))  %>%
        arrange(-Mean_mlog10_P) %>% mutate( Rank.ID.mlog10_P = 1:n()) %>%
        arrange(-Mean_PCT.dif) %>% mutate( Rank.pos.ID.pct = 1:n()) %>%
        arrange(Mean_PCT.dif) %>% mutate( Rank.neg.ID.pct = 1:n()) %>%
        arrange(-Mean_AVG_log2FC) %>%  mutate( Rank.ID.pos.log2FC = 1:n()) %>%
        arrange(Mean_AVG_log2FC) %>% mutate( Rank.ID.neg.log2FC = 1:n()) %>%
        mutate( Rank.ID = Rank.pos.ID.pct + Rank.ID.pos.log2FC + Rank.ID.mlog10_P,
                Rank.ID.neg = Rank.neg.ID.pct + Rank.ID.neg.log2FC + Rank.ID.mlog10_P
        ) %>% # + Rank.ID.mlog10_P
        arrange(Rank.ID) %>% mutate(Rank.ID.New = 1:n()) %>% 
        arrange(Rank.ID.neg) %>% mutate(Rank.ID.neg.New = 1:n())  %>% 
        arrange(Rank.ID)
    })
    
    
    for (nn in names(Result[[ll]]) ) {
      # nn = names(Result[[ll]])[1];nn
          Result[[ll]][[nn]] <- Result[[ll]][[nn]] %>% left_join(Result.sig[[ll]][[nn]] %>%
                                                                   dplyr::select(one_of(c("Direction", "rowname"))), by = "rowname" )
          colnames( Result[[ll]][[nn]])
      
      if (nrow(Result[[ll]][[nn]]%>% dplyr::filter(Direction %in% c("Up", "Down")))  > 0 ) {
        if (!is.null(Result.sig[[ll]][[nn]])){
          df <- Result[[ll]][[nn]] %>% dplyr::filter(Direction %in% c("Up", "Down")) 
          
          n = 10
          labels = unique(c(
            df %>% dplyr::filter(Direction == "Up") %>% dplyr::top_n(-n, Rank.ID) %>% pull(rowname),
            df %>% dplyr::filter(Direction == "Down") %>% dplyr::top_n(-n, Rank.ID.neg) %>% pull(rowname)));labels
          df <- df %>% mutate(
            Plotname = ifelse(rowname %in% labels, rowname, ""),
            Direction = factor(Direction, levels = c("Up", "Down", "non")),
            mLog10PVal = -log10(Mean_P_Val)
          )
          options(ggrepel.max.overlaps = Inf) 
          p <- ggscatter(df, x = "Mean_AVG_log2FC", y = "mLog10PVal", 
                         color = "Direction",palette = c("red3", "skyblue3", "lightgray"),
                         # color = "Sig.overall", palette = c("red3",  "orange", "skyblue3","lightgray", "snow"),
                         #label = "Row.names", repel = TRUE, label.select = labels,
                         shape = 20, 
                         size = ifelse(df[["Direction"]] == "non", 1, 3), # Points color, shape and size
                         xlab = "avg_log2FoldChange", 
                         ylab = "-log10(pvalue)",
                         title = paste0(ll, ", ",nn)
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
            # )  +
            geom_errorbar(data = df,
                          aes(xmin = Mean_AVG_log2FC - SEM_AVG_log2FC, xmax = Mean_AVG_log2FC + SEM_AVG_log2FC),
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
                                sub("_[^_]+$", "", stringr::str_split(markers.file, "/")[[1]][3]), 
                                "_", paste0(ll, ".",nn), "_l2_log2FC_ggscatter.png"), 
              width = 6, 
              height= 6 , res = 300, units = "in")
          print(p)
          dev.off()
          
          
          p <- ggscatter(df, x = "Mean_PCT.dif", y = "mLog10PVal", 
                         color = "Direction",palette = c("red3", "skyblue3", "lightgray"),
                         # color = "Sig.overall", palette = c("red3",  "orange", "skyblue3","lightgray", "snow"),
                         #label = "Row.names", repel = TRUE, label.select = labels,
                         shape = 20, 
                         size = ifelse(df[["Direction"]] == "non", 1, 3), # Points color, shape and size
                         xlab = "Difference of Cell Percentage", 
                         ylab = "-log10(pvalue)",
                         title = paste0(ll, ", ",nn)
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
                          aes(xmin = Mean_PCT.dif - SEM_PCT.dif, xmax = Mean_PCT.dif + SEM_PCT.dif),
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
                                sub("_[^_]+$", "", stringr::str_split(markers.file, "/")[[1]][3]), "_",
                                paste0(ll, ".",nn), "_l2_pct.dif_ggscatter.png"), 
              width = 6, 
              height= 6 , res = 300, units = "in")
          print(p)
          dev.off()
          
        } 
      }
    }
  }
  saveRDS(Result, 
          paste0(Disk, Project.folder, "/", Rds.folder, "/", stringr::str_split(markers.file, ".RData")[[1]][1],
                 ".Mean.Rds") 
  )
  print(markers.file)
}
rm(Result.comb.df, Result.comb.list, Result, export.folder0, dir, nn, df, 
   markers.file, p, labels, n , Result.comb,
   ee, ll, export.folder, Result.sig
)
gc()
