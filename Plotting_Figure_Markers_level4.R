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

for (markers.file in markers.files) {
  # markers.file = markers.files[1]
  Result.combine.list <- list()
  
  Result <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/", markers.file))  
  Result <- lapply(Result, function(x) {
    mutate(x,
           Sig= factor(p_val_adj <= 0.001 , levels = c("TRUE", "FALSE")),
           type = factor(type, levels = c("Up", "Down")),
           pct.dif = pct.1 - pct.2) %>% # & avg_log2FC > 0.2
      tibble::rownames_to_column(var = "rowname")
  })
  
  # prepare folder
  export.folder0 <- paste0("Markers")
  dir.create(file.path(paste0(Disk, Project.folder, "/", figures.folder, "/"), export.folder0), showWarnings = FALSE)
  dir <- paste0(Disk, Project.folder, "/", figures.folder, "/", export.folder0, "/")
  
  for (nn in names(Result) ) {
    if (nrow(Result[[nn]]) >0 ) {
      # nn2 = names(Result)[grepl(stringr::str_split(nn, "_")[[1]][1], names(Result), fixed = TRUE)]
      nn2 = names(Result)[startsWith(names(Result), paste0(unique(Result[[nn]]$Group1)) )] # stringr::str_split(nn, "_")[[1]][1]
      nn2.1 = names(Result)[startsWith(names(Result), paste0(unique(Result[[nn]]$Group1), "_",unique(Result[[nn]]$Group2)) )]
      nn2.2 = nn2[!nn2 %in% nn2.1]
      
      if (!rlang::is_empty(nn2.2)){
        
        
        # unique(Result[[nn]]$Group3)
        compare.groups <- c(unique(Result[[nn]]$ident.1), unique(Result[[nn]]$ident.2))
        list.name1 = paste0(unique(Result[[nn]]$Group1), "_", 
                            unique(Result[[nn]]$Group2), "_",
                            paste(compare.groups, collapse = ":")) 
        list.name2 = paste0(unique(Result[[nn]]$Group1), "_", 
                            stringr::str_split(stringr::str_split(nn2.2, paste0(unique(Result[[nn]]$Group1), "_") )[[1]], "_")[[2]][1], "_",
                            paste(compare.groups, collapse = ":")) 
        t.test.group = c(Reduce(setdiff, strsplit(c(list.name1, list.name2), split = "_")), 
                         Reduce(setdiff, strsplit(c(list.name2, list.name1), split = "_")) )
        
        ids1 <- c()
        for (ii in nn2.1){
          id1 <- c(Reduce(setdiff, strsplit(c(ii, list.name1), split = "_")))
          ids1 <- c(id1, ids1)
          rm(ii, id1)
        }
        ids2 <- c()
        for (ii in nn2.2){
          id1 <- c(Reduce(setdiff, strsplit(c(ii, list.name2), split = "_")))
          ids2 <- c(id1, ids2)
          rm(ii, id1)
        }
        
        if(!(list.name1 %in% names(Result.combine.list)) ){
          df1 <- Result[nn2.1] %>% purrr::imap(function(x, y) x %>% 
                                                 rename_with(~paste(., y, sep = '.'), -c(rowname, Group1, Group2, Compare.Group, ident.1, ident.2))) %>% 
            purrr::reduce(full_join, by = c('rowname', 'Group1', 'Group2', 'Compare.Group', 'ident.1', 'ident.2'))
          colnames(df1)[2]
          df2 <- Result[nn2.2] %>% purrr::imap(function(x, y) x %>% rename_with(~paste(., y, sep = '.'), -c(rowname, Group1, Group2, Group3, Compare.Group, ident.1, ident.2))) %>% 
            purrr::reduce(full_join, by = c('rowname', 'Group1', 'Group2', 'Compare.Group', 'ident.1', 'ident.2'))
          colnames(df2)[2]
          
          df.list <- list()
          df.list[[paste0(list.name1, "_mean")]] <-data.frame(
            rowname = df1$rowname,
            p_val.mean = rowMeans(subset(df1, select = colnames(df1)[startsWith(colnames(df1), "p_val.")]), na.rm = TRUE),
            avg_log2FC.mean = rowMeans(subset(df1, select = colnames(df1)[startsWith(colnames(df1), "avg_log2FC.")]), na.rm = TRUE),
            pct.1.mean = rowMeans(subset(df1, select = colnames(df1)[startsWith(colnames(df1), "pct.1.")]), na.rm = TRUE),
            pct.2.mean = rowMeans(subset(df1, select = colnames(df1)[startsWith(colnames(df1), "pct.2.")]), na.rm = TRUE),
            p_val_adj.mean = rowMeans(subset(df1, select = colnames(df1)[startsWith(colnames(df1), "p_val_adj.")]), na.rm = TRUE),
            pct.dif.mean = rowMeans(log2(subset(df1, select = colnames(df1)[startsWith(colnames(df1), "pct.dif.")])), na.rm = TRUE), 
            
            p_val.sd = apply(subset(df1, select = colnames(df1)[startsWith(colnames(df1), "p_val.")]), 1, sd, na.rm = TRUE),
            avg_log2FC.sd = apply(subset(df1, select = colnames(df1)[startsWith(colnames(df1), "avg_log2FC.")]), 1, sd, na.rm = TRUE),
            pct.1.sd = apply(subset(df1, select = colnames(df1)[startsWith(colnames(df1), "pct.1.")]), 1, sd, na.rm = TRUE),
            pct.2.sd = apply(subset(df1, select = colnames(df1)[startsWith(colnames(df1), "pct.2.")]), 1, sd, na.rm = TRUE),
            p_val_adj.sd = apply(subset(df1, select = colnames(df1)[startsWith(colnames(df1), "p_val_adj.")]), 1, sd, na.rm = TRUE),
            pct.dif.sd = apply(log2(subset(df1, select = colnames(df1)[startsWith(colnames(df1), "pct.dif.")])),  1, sd, na.rm = TRUE),
            
            n = length(colnames(df1)[startsWith(colnames(df1), "p_val.")]),
            n.ids = paste(ids1, collapse = ", ")
            
          ) %>% 
            full_join(df1 %>% dplyr::select(c(rowname, Group1, Group2, Compare.Group, ident.1, ident.2)), by = 'rowname') %>%
            mutate(
              p_val.sem = p_val.sd/sqrt(n),
              avg_log2FC.sem = avg_log2FC.sd/sqrt(n),
              pct.1.sem = pct.1.sd/sqrt(n),
              pct.2.sem = pct.2.sd/sqrt(n),
              p_val_adj.sem = p_val_adj.sd/sqrt(n),
              pct.dif.sem = pct.dif.sd/sqrt(n)
            )
          
          df.list[[paste0(list.name2, "_mean")]] <-data.frame(
            rowname = df2$rowname,
            p_val.mean = rowMeans(subset(df2, select = colnames(df2)[startsWith(colnames(df2), "p_val.")]), na.rm = TRUE),
            avg_log2FC.mean = rowMeans(subset(df2, select = colnames(df2)[startsWith(colnames(df2), "avg_log2FC.")]), na.rm = TRUE),
            pct.1.mean = rowMeans(subset(df2, select = colnames(df2)[startsWith(colnames(df2), "pct.1.")]), na.rm = TRUE),
            pct.2.mean = rowMeans(subset(df2, select = colnames(df2)[startsWith(colnames(df2), "pct.2.")]), na.rm = TRUE),
            p_val_adj.mean = rowMeans(subset(df2, select = colnames(df2)[startsWith(colnames(df2), "p_val_adj.")]), na.rm = TRUE),
            pct.dif.mean = rowMeans(log2(subset(df2, select = colnames(df2)[startsWith(colnames(df2), "pct.dif.")])), na.rm = TRUE), 
            
            p_val.sd = apply(subset(df2, select = colnames(df2)[startsWith(colnames(df2), "p_val.")]), 1, sd, na.rm = TRUE),
            avg_log2FC.sd = apply(subset(df2, select = colnames(df2)[startsWith(colnames(df2), "avg_log2FC.")]), 1, sd, na.rm = TRUE),
            pct.1.sd = apply(subset(df2, select = colnames(df2)[startsWith(colnames(df2), "pct.1.")]), 1, sd, na.rm = TRUE),
            pct.2.sd = apply(subset(df2, select = colnames(df2)[startsWith(colnames(df2), "pct.2.")]), 1, sd, na.rm = TRUE),
            p_val_adj.sd = apply(subset(df2, select = colnames(df2)[startsWith(colnames(df2), "p_val_adj.")]), 1, sd, na.rm = TRUE),
            pct.dif.sd = apply(log2(subset(df2, select = colnames(df2)[startsWith(colnames(df2), "pct.dif.")])),  1, sd, na.rm = TRUE),
            
            n = length(colnames(df2)[startsWith(colnames(df2), "p_val.")]),
            n.ids = paste(ids2, collapse = ", ")
            
          ) %>% 
            full_join(df2 %>% dplyr::select(c(rowname, Group1, Group2, Compare.Group, ident.1, ident.2)), by = 'rowname') %>%
            mutate(
              p_val.sem = p_val.sd/sqrt(n),
              avg_log2FC.sem = avg_log2FC.sd/sqrt(n),
              pct.1.sem = pct.1.sd/sqrt(n),
              pct.2.sem = pct.2.sd/sqrt(n),
              p_val_adj.sem = p_val_adj.sd/sqrt(n),
              pct.dif.sem = pct.dif.sd/sqrt(n)
            )
          
          df.list[[paste0(list.name1, "_data")]] <-df1
          df.list[[paste0(list.name2, "_data")]] <-df2
          
          for (tt in c("avg_log2FC.", "pct.1.", "pct.2.", "pct.dif.")) {
            stats <- as.data.frame(do.call(rbind, lapply(1:nrow(df1), function(i){
              (matrixTests::row_t_welch(df1[i, c(colnames(df1)[startsWith(colnames(df1), tt)]) ], 
                                        df2[i, c(colnames(df2)[startsWith(colnames(df2), tt)]) ] ))
            })))
            stats <- stats %>% mutate(rowname = df1$rowname, 
                                      Direction = ifelse( mean.x > mean.y & pvalue < 0.05, "Up", 
                                                          ifelse( mean.x < mean.y & pvalue < 0.05, "Down", "non")),
                                      Sig.overall = ifelse( pvalue < 0.05, "Sig", "non.sig"))
            rownames(stats) = df1$rowname
            df.list[[paste0(unique(Result[[nn]]$Group1), "_", tt, paste(t.test.group, collapse = ":") ) ]] <- stats
            print(tt)
          }
          rm(stats, tt)
          Result.combine.list[[paste0(list.name1)]] <- df.list
          
          n = 5
          for (tt in c("avg_log2FC.", "pct.dif.")) {
            df <- df.list[[paste0(unique(Result[[nn]]$Group1), "_", tt, paste(t.test.group, collapse = ":") ) ]]
            labels = c(df %>% dplyr::filter(Sig.overall == "Sig") %>% dplyr::top_n(-n, pvalue) %>% pull(rowname))
            df <- df %>% mutate(
              Plotname = ifelse(rowname %in% labels, rowname, ""),
              Direction = factor(Direction, levels = c("Up", "Down", "non")),
              Sig.overall = factor(Sig.overall, levels = c("Sig", "non-sig")))
            
            p <- ggscatter(df, x = "mean.x", y = "mean.y", 
                           color = "Direction",palette = c("red3", "skyblue3", "lightgray"),
                           # color = "Sig.overall", palette = c("red3",  "orange", "skyblue3","lightgray", "snow"),
                           #label = "Row.names", repel = TRUE, label.select = labels,
                           shape = 20, 
                           size = ifelse(df[["Direction"]] == "non", 1, 3), # Points color, shape and size
                           # add = "reg.line",  # Add regression line
                           # add.params = list(color = "blue", fill = "lightgray", linetype="dashed", size=0.2), # Customize reg. line
                           # conf.int = TRUE, # Add confidence interval
                           # cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                           # cor.method = "pearson",
                           # cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n"),
                           xlab = paste0(tt, "(", paste(compare.groups, collapse = ":"), ")"," in Group_", t.test.group[1]), 
                           ylab = paste0(tt, "(", paste(compare.groups, collapse = ":"), ")"," in Group_", t.test.group[2]),
                           title = unique(Result[[nn]]$Group1)
                           
            ) +
              geom_abline(intercept = 0, slope = 1, color = "darkgray", linetype="solid", linewidth=1) +
              # geom_vline(xintercept = 0, color = "gray", linetype="dashed", size=0.2) +
              # geom_hline(yintercept = 0, color = "gray", linetype="dashed", size=0.2) 
              theme_bw() +
              theme(legend.position = "none" ) +
              #geom_text(aes(label = .data[["Plotname"]]), size = 4, color ="black") +
              ggrepel::geom_text_repel(aes(label = .data[["Plotname"]]),
                                       min.segment.length = -0, seed = 42, box.padding = 0.8
                                       , max.overlaps = Inf
              )  
            # # geom_errorbar(data = df, 
            #               aes(ymin = mean.y - sqrt(var.y)/sqrt(obs.y), ymax = mean.y + sqrt(var.y)/sqrt(obs.y)),
            #               width = 0.2,
            #               position = position_dodge(0.05)) 
            # +ggrepel::geom_text_repel(aes(label = .data[["Plotname"]]),
            #   nudge_x = .15,
            #   box.padding = 0.5,
            #   nudge_y = 1,
            #   segment.curvature = -0.1,
            #   segment.ncp = 1,
            #   segment.angle = 1
            # ) +
            # geom_smooth(method="lm", formula=y~x+0) # force to pass 0
            # stat_function(fun = function(x) predict(fit, newdata = data.frame(avg_log2FC.x = x)),
            #                 color = "blue",linetype="dashed", size=0.2) +
            ###   annotate(label = sprintf("y = %.3f x\nR = %.2f", coef(fit), summary(fit)$r.squared),
            #            geom = "text" 
            #            , x = (min(df[["avg_log2FC.x"]])) + (0.9*(min(df[["avg_log2FC.x"]]))), 
            #            y = max(df[["avg_log2FC.y"]])
            #            #, size = 12
            #            )  # force to pass 0
            
            print(p)
            
            # anova(lm(avg_log2FC.y ~ avg_log2FC.x*avg_log2FC.y, data = df))
            # ###' if ":" Pr(F) < 0.05, 
            # ###' Interaction is not significant, so the slope across groups
            # ###' is not different.
            # anova(lm(avg_log2FC.y ~ avg_log2FC.x + avg_log2FC.y, data = df))
            # ###' The category variable (Species) is significant,
            # ###' so the intercepts among groups are different
            # 
            # library(strucchange)
            # sctest(df$avg_log2FC.y ~ df$avg_log2FC.x, type = "Chow", point = 10)
            # ###' We can reject the null hypothesis of the test because the p-value is less than 0.05. 
            # ###' This means we have enough evidence to conclude that the data contains a structural 
            # ###' breakpoint.
            # ###' To put it another way, two regression lines can match the data pattern better than 
            # ###' a single regression line.
            # 
            
            png(filename = paste0(dir, 
                                  sub("_[^_]+$", "", stringr::str_split(markers.file, "/")[[1]][3]), "_",
                                  stringr::str_split(list.name1, ":")[[1]][1], "_", tt, "_level4_ggscatter.png"), 
                width = 6, 
                height= 6 , res = 300, units = "in")
            print(p)
            dev.off()
          }
        }
      }
    }
  }
  saveRDS(Result.combine.list, 
          paste0(Disk, Project.folder, "/", Rds.folder, "/", stringr::str_split(markers.file, ".Rds")[[1]][1],
                 ".combined.Rds") )
  print(markers.file)
}
rm(Result, export.folder0, dir, nn, df, tt,
   df1, df2, df.list, ids1, ids2, list.name1, list.name2, nn2, nn2.1, nn2.2, 
   markers.file,compare.groups,t.test.group, 
   p, Result.combine.list, labels, n)

