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
types = c("avg_log2FC", "pct.dif")

for (markers.file in markers.files) {
  # markers.file = markers.files[3]; markers.file
  load(paste0(Disk, Project.folder, "/", Rds.folder, "/", markers.file), verbose = TRUE)  
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
  Result.anova <- list()
  Result.t.test <- list()
  Result.t.test2 <- list()
  Result.summary <- list()
  Result <- list()
  Result.df_long <- list()
  Result.all <- list()
  for (ll in names(Result.comb.df)) {
    # ll = names(Result.comb.df)[1];ll
    Result.all <- Result.comb.df
    compare.groups = names(Result.comb.list[[ll]])
    for (type in types) {
      # type ="avg_log2FC";
      
      # Prepare data ------
      cols.FC = grep(type, colnames(Result.comb.df[[ll]]), ignore.case =TRUE);cols.FC
      cols.gene = grep("rowname", colnames(Result.comb.df[[ll]]), ignore.case =TRUE);cols.gene
      df <- Result.comb.df[[ll]][, c(cols.gene, cols.FC)] 
      all.genes <- Result.comb.df[[ll]][, c(cols.gene, cols.FC)][["rowname"]];head(all.genes)
      sample.name <- colnames(Result.comb.df[[ll]][, c(cols.FC)] );head(sample.name)
      df_long <- df %>%
        pivot_longer(cols = 2:ncol(.)) %>%
        pivot_wider(names_from = rowname, values_from = value) %>%
        separate(name, into = c("trt.sample"), sep = paste0("_", type)) %>%
        mutate(group = stringr::str_split_fixed(trt.sample, "_sample_", 2)[, 1] ,
               sample_id = stringr::str_split_fixed(trt.sample, "_sample_", 2)[, 2]) %>% 
        mutate(group = factor(group, levels = compare.groups)) %>% 
        relocate(group)
      table(df_long[["group"]]); unique(df_long[["group"]])
      head(df_long)[1:6, c("group", "sample_id")]
      Result.df_long[[ll]] <- df_long 
      df_long <- df_long %>%
        dplyr::select(-trt.sample, -sample_id) 
      head(df_long)[1:6, 1:6]
      head(df_long)[1:6, c(1:6,ncol(df_long))]
      
      # 1 way ANOVA: One-way ANOVA analysis by row -----
      fit_aov <- function(col) {
        aov(col ~ group, data = df_long)
      }
      # anovas <- map(df_long[, 1:(ncol(df_long)-1)], fit_aov) # group column in the last col
      anovas <- map(df_long[, 2:(ncol(df_long))], fit_aov) # group column in the first col
      summary(anovas[[1]]);
      names(anovas)[1];summary(anovas[[1]])[[1]][1,5];
      
      # save result2: add p value to Result.all[[ll]]
      for (nn in names(anovas) ) {
        Result.all[[ll]][which(Result.all[[ll]]$rowname == nn), paste0("anova.", type, ".pval")] <- summary(anovas[[nn]])[[1]][1,5] 
      }
      Result.all[[ll]] <- Result.all[[ll]] %>% arrange((!!sym(paste0("anova.", type, ".pval"))))
      
      # visualize the genes
      Result.all[[ll]][1,c("rowname", paste0("anova.", type, ".pval"))]
      ggplot(df_long, aes(x = group, y = .data[[Result.all[[ll]][1,"rowname"]]], colour = group)) + 
        geom_point()
      
      # 1 way ANOVA and save as data frame -----
      library(rstatix)
      anovas2_ls <- list()
      for (cc in colnames(df_long[, 2:(ncol(df_long))])) {
        # cc = colnames(df_long[, 2:(ncol(df_long))])[1];cc
        test.df <- df_long %>% 
          dplyr::select(one_of(c("group", cc)))
        colnames(test.df) <- c("group", "gene")
        if(sum(test.df[["gene"]]) > 0){
          anovas2_ls[[cc]] <-
            rstatix::anova_test(data = test.df, formula= formula(paste("gene", "~", "group"))) %>%
            adjust_pvalue(method = "bonferroni") %>%
            add_significance("p.adj") %>% 
            add_significance("p") %>%
            ungroup() %>% as_tibble() %>% 
            mutate(Y.value = cc,
                   group = Select.Group,
                   sig = ifelse(.data[["p"]] < 0.05 & !is.na(.data[["p"]]) , "Yes", "ns"),
                   Test.type = "1ANOVA"
            )
        }
      }
      anovas2 = do.call(rbind, anovas2_ls)
      anovas2 <- anovas2 %>% 
        rename(rowname = Y.value)%>% 
        arrange(p.adj)
      head(anovas2)
      Result.anova[[ll]][[type]] <- anovas2
      

      # t.test in 2 groups and save as data frame -----
      if (length(compare.groups)  == 2){
        df_long2 <- df_long # %>% dplyr::filter(group %in% c("Older_Lean" , "Older_Overweight"))
        fit_t.test2 <- function(col) {
          tidy(t.test(col ~ group, data = df_long2))
        }
        t.test2 <- map(df_long2, fit_t.test2) # group column in the first col
        t.test2[[1]];
        names(anovas)[1];
        t.test2[[1]]$p.value;
        
        # save result1
        res.t.test = do.call(rbind, t.test2)
        res.t.test <- res.t.test %>% 
          tibble::rownames_to_column(var = "rowname") %>% 
          arrange(p.value)
        head(res.t.test)
        Result.t.test2[[ll]][[type]] <- res.t.test
        
        # save result2: add p value to df
        for (nn in names(t.test2) ) {
          Result.all[[ll]][which(Result.all[[ll]]$rowname == nn), paste0("t.test2.", type, ".pval") ] <- t.test2[[nn]]$p.value 
        }
        Result.all[[ll]] <- Result.all[[ll]] %>% arrange((!!sym(paste0("t.test2.", type, ".pval"))))
        
        # visualize the genes
        Result.all[[ll]][1,c("rowname", paste0("t.test2.", type, ".pval"))]
        ggplot(df_long2, aes(x = group, y = .data[[Result.all[[ll]][1,"rowname"]]], colour = group)) + 
          geom_point()
        rm(df_long2, t.test2, res.t.test)
      }
      
      
      # t_test in multiple groups and save as data frame -----
      library(rstatix)
      t_test_ls <- list()
      for (cc in colnames(df_long[, 2:(ncol(df_long))])) {
        # cc = colnames(df_long[, 2:(ncol(df_long))])[1];cc
        test.df <- df_long %>% 
          dplyr::select(one_of(c("group", cc)))
        colnames(test.df) <- c("group", "gene")
        t_test_ls[[cc]] <-
          rstatix::t_test(data = test.df, formula= formula(paste("gene", "~", "group"))) %>%
          adjust_pvalue(method = "bonferroni") %>%
          add_significance("p.adj") %>% 
          add_significance("p") %>%
          ungroup() %>% as_tibble() %>% 
          # add_xy_position(x = group) %>% 
          mutate(Y.value = cc,
                 group = Select.Group,
                 sig = ifelse(.data[["p"]] < 0.05 & !is.na(.data[["p"]]) , "Yes", "ns"),
                 Test.type = "T-test"
          ) 
      }
      t_test2 = do.call(rbind, t_test_ls)
      t_test2 <- t_test2 %>% 
        rename(rowname = Y.value)%>% 
        arrange(p.adj)
      head(t_test2)
      Result.t.test[[ll]][[type]] <- t_test2
      
      # summary ------
      df.summary <- df_long %>% 
        as.data.frame() %>% 
        group_by(group) %>%
        summarise(
          across(where(is.numeric), 
                 .fns = 
                   list(
                     min =  ~min(., na.rm = TRUE)# 
                     , count = ~n()
                     , median =  ~median(., na.rm = TRUE)
                     , mean =  ~mean(., na.rm = TRUE)
                     , stdev = ~sd(., na.rm = TRUE)
                     , q25 = ~quantile(., 0.25)
                     , q75 = ~quantile(., 0.75)
                     , max =  ~max(., na.rm = TRUE)
                     ))) 
      df.summary <-df.summary %>%
        pivot_longer(where(is.numeric), names_sep='_', names_to=c('rowname', '.value'))
      df.summary <-df.summary %>%
        mutate(sem = stdev/sqrt(count) ) %>%
        arrange((!!sym("rowname")))
      head(df.summary)
      Result.summary[[ll]][[type]] <- df.summary
      
      # plot DEGs in heatmap ---------------------------
      Averages <- Result.comb.df[[ll]][, c(cols.gene, cols.FC)] %>%
        dplyr::filter( rowname %in% c(Result.all[[ll]] %>% dplyr::filter((!!sym(paste0("anova.", type, ".pval"))) < 0.05) %>% pull(rowname)) ) %>%
        tibble::rownames_to_column(var = "id") %>%
        dplyr::select(-id)%>%
        tibble::column_to_rownames(var = "rowname"); head(Averages[, 1:2])
      sample.df <- data.frame(sample.name = colnames(Result.comb.df[[ll]][, c(cols.FC)] ) )
      sample.df <- sample.df %>% 
        mutate(group = stringr::str_split_fixed(sample.name, "_sample_", 2)[, 1],
               sample_id2 = stringr::str_split_fixed(sample.name, "_sample_", 2)[, 2]) %>% 
        mutate(sample_id = stringr::str_split_fixed(sample_id2, paste0("_", type), 2)[, 1] ) %>% 
        mutate(group = factor(group, levels = compare.groups))
      head(sample.df)
      
      library(RColorBrewer)
      library(circlize)
      coul <-colorRampPalette(brewer.pal(12, "Paired"))(length(unique(sample.df$group)))
      colors = coul[1: length(unique(sample.df$group))]
      pie(rep(1, length(coul)), col = coul , main="") 
      names(colors) <- sort(unique(sample.df$group))
      
      col_ha = columnAnnotation(
        df = data.frame(type = sample.df$group),
        col = list(type = colors), 
        show_legend = F , show_annotation_name = F,
        annotation_legend_param = list(
          group = list( 
            title_gp = gpar(fontsize = 10, 
                            fontface = "bold")
            , ncol = 2,  
            labels_gp = gpar(fontsize = 8))
          ,
          annotation_name_side = "top" )
      )
      ComplexHeatmap::draw(col_ha)
      hist(as.matrix(Averages))
      # legend.col = colorRamp2(seq(-0.5, 0.5, length = 3), c("blue", "#EEEEEE", "red"))
      # legend.col = colorRamp2(seq(min(as.matrix(Averages)), max(as.matrix(Averages)), length = 3), c("blue", "#EEEEEE", "red"))
      legend.col = colorRamp2(seq(-0.4, 0.4, length = 3), c("blue", "#EEEEEE", "red"))
      row_labels = rownames(Averages)# structure(Averages[["gene"]], names = rownames(Averages))
      
      cl = kmeans(as.matrix(Averages), centers = 4, nstart = 100)$cluster
      dend = as.dendrogram(hclust(dist(as.matrix(Averages))))
      dend = dendextend::color_branches(dend, k = 2)
      
      ht_opt( message = FALSE)
      set.seed(123)
      ht_list = 
        Heatmap(
          as.matrix(Averages), name = "mat",
          col = legend.col , 
          border = TRUE,
          border_gp = gpar(col = "black", lty = 1),
          
          cluster_columns = F, 
          cluster_column_slices = FALSE,
          
          cluster_rows = FALSE,
          # cluster_rows = TRUE, row_split = 2, # Split by dendrogram
          # cluster_rows = dend, row_split = 2, # Split by dendrogram
          # row_km = 4, row_km_repeats = 100, # Split by k-means clustering 
          row_split = cl, # Split by k-means clustering 
          # row_split = rep(c("A", "B"), 9), # Split by categorical variables
          # cluster_row_slices = FALSE,
          
          
          # column_title = paste0(export.folder, "_GeneAverages"), 
          column_title_side = c("bottom"),
          column_title_gp = gpar(fontsize = 10, fontface = "plain"),
          #column_title_gp = gpar(fontsize = 17, fontface = "bold"),
          column_title_rot = 0,
          
          column_split= sample.df$group,
          column_gap = unit(1, "mm"),
          
          show_column_names = FALSE,
          column_names_gp = gpar(fontsize = 12), 
          column_names_rot = 45, 
          column_names_max_height = unit(6, "cm"),
          
          row_names_side = "left", 
          show_row_names = FALSE,
          row_names_max_width = unit(6, "cm"),
          row_names_gp = gpar(fontsize = 7, fontface = "italic"),
          row_names_rot = 0,
          
          bottom_annotation = col_ha, 
          
          heatmap_legend_param = list(direction = "vertical", legend_height = unit(3, "cm"), 
                                      direction = "vertical",
                                      title = type , # "Log2 Fold Change of Markers", 
                                      title_gp = gpar(fontsize = 10, fontfface = "plain"),
                                      title_position = "leftcenter-rot"),
          # show_heatmap_legend = FALSE,
          
          row_dend_side = NULL,
          na_col = "black",
          row_labels = row_labels[rownames(Averages)]
        )
      ComplexHeatmap::draw(ht_list , 
                           heatmap_legend_side = "right"
      )
      Averages.df <- as.data.frame(Averages) 
      Averages.df$cluster <- cl
      Averages.df$rowname <- rownames(Averages.df)
      table(Averages.df$cluster);head(Averages.df[, 1:2])
      df <- df %>% left_join(Averages.df %>%
                               dplyr::select(one_of(c("rowname", "cluster"))), by = "rowname")
      colnames(df)
      table(df$cluster); unique(df$cluster);
      Result[[ll]][[type]] <- df
      Result.all[[ll]] <- Result.all[[ll]] %>% left_join(Averages.df %>%
                                             dplyr::select(one_of(c("rowname", "cluster"))), by = "rowname")
      colnames(Result.all[[ll]])
      table(Result.all[[ll]]$cluster); unique(Result.all[[ll]]$cluster);
      
      dir;
      png(filename = paste0(dir, 
                            sub("_[^_]+$", "", stringr::str_split(markers.file, "/")[[1]][3]), 
                            "_", paste0(ll, ".", type), "_l2_ANOVA_Heatmap.png"), 
          width = (ncol(Averages)/6)+1.8, 
          height= (nrow(Averages)/120)+0.7 , res = 300, units = "in")
      ComplexHeatmap::draw(ht_list , 
                           heatmap_legend_side = "right")
      dev.off()
    }
    rm(compare.groups, type, cols.FC, cols.gene,
       df, all.genes, sample.name, df_long, anovas,anovas2_ls,
       cc, nn, test.df, anovas2,
       t_test_ls, t_test2,
       df.summary, Averages, sample.df, 
       coul, colors, legend.col, row_labels, 
       cl, dend, ht_list, Averages.df
       )
    gc()
  }
  
  save(Result.all,
       Result.anova ,
       Result.t.test ,
       Result.t.test2 ,
       Result.summary ,
       Result ,
       Result.df_long , 
       file = paste0(dir2, 
                     # sub("_[^_]+$", "", stringr::str_split(markers.file, "/")[[1]][3]), 
                     unlist(stringr::str_split(stringr::str_split(markers.file, "/")[[1]][3], "[.]"))[1], 
                     "_", paste0(Select.Group, "_", Compare.G),
                     ".Mean.ANOVA.RData") 
  )
  rm(
     Result.all, Result.anova, Result.t.test ,
     Result.t.test2 ,
     Result.summary ,
     Result ,
     Result.df_long)
  gc()
  print(markers.file)
}
rm(Result.comb.df, Result.comb.list, 
   export.folder, export.folder0,
   dir,
   dir1, dir2, types,col_ha, 
)
gc()
