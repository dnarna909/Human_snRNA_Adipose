
# rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.
memory.limit(size = 1e+13)

# import parameters and packages --------------------------------------------------------------------------------------------
# source(paste0("/media/jianie/Extreme SSD1/2022-09-01 STARR_SGLT2 Combine/Project Parameters.R"), local = knitr::knit_global())
source(paste0(paste0(Disk, Project.folder, "/", "Project Parameters.R")), local = knitr::knit_global())
library("RColorBrewer")
library("ggpubr")

# Import data -----------------------------------------------------------------------------------------------------------------------------------------------------
dir.create(file.path(paste0(Disk, Project.folder, "/", figures.folder, "/", export.folder, ".plotting")), showWarnings = FALSE)
dir1 <- paste0(Disk, Project.folder, "/", figures.folder, "/", export.folder, ".plotting", "/")
dir.create(file.path(paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder)), showWarnings = FALSE)
dir0 <- paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder,  "/")

# set up values -----------------
file.name <- stringr::str_split(file.variables, ".Rds")[[1]][1]
file.name
Annotation.file.name <- paste(stringr::str_split(file.name, "_")[[1]][c(1,2,3)], collapse = "_")
Annotation.file.name 
Type.file.name <- paste(stringr::str_split(file.name, "_")[[1]][c(4)], collapse = "_")
Type.file.name 

# Import metadata ------------------
sample.meta <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/", Sample.meta.file)) %>% 
  tibble::rownames_to_column(var = "rowname")
sample.meta <- sample.meta[ duplicated(sample.meta[["Dataset"]])==F,] 
sample.meta <- sample.meta[ duplicated(sample.meta[["sample_id"]])==F,] 
sample.meta <- sample.meta %>% 
  dplyr::filter((!!sym(Select.Group)) == T)
sample.meta <- sample.meta[,colSums(is.na(sample.meta))<nrow(sample.meta)] # Remove columns from dataframe where ALL values are NA


# Import Variable data ------------------
var.df <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/", file.variables)) # contain "Dataset"
var.df <- var.df[,colSums(is.na(var.df))<nrow(var.df)] 
colnames(var.df)
var.df[["Dataset"]]

# combine ---------
all.df <- sample.meta %>% left_join(var.df, by = "Dataset")
colnames(all.df)
all.df.n <- all.df %>% select(where(is.numeric))%>% 
  select(which(!colSums(., na.rm=TRUE) %in% 0))
numeric.cols = intersect(setdiff(colnames(var.df), "Dataset"), colnames(all.df.n))
# numeric.cols = colnames(all.df.n)[colnames(all.df) %in% setdiff(colnames(var.df), "Dataset")]
numeric.cols


ano.results <- data.frame()
t.test.results <- data.frame()
summary.results <- data.frame()

for (nn in numeric.cols){
  # nn = numeric.cols[1]
  # Comparisons by one way ANOVA
  library(rstatix)
  table(all.df[[compare.group]])
  
  df_test <- all.df[!is.na(all.df[[nn]]), ] %>% data.frame() %>% 
    group_by_at(vars(one_of(c(compare.group))))   %>%
    # nest() %>%
    dplyr::filter(n() > 1 ) %>%
    ungroup() %>%  
    dplyr::filter(n_distinct((!!sym(compare.group))) >1) %>% as.data.frame() 
  table(df_test[[compare.group]])
  df_test[[compare.group]] <- factor(df_test[[compare.group]]  , 
                                     levels = levels(df_test[[compare.group]])[ levels(df_test[[compare.group]]) %in% unique(df_test[[compare.group]])])
  table(df_test[[compare.group]])
  
  
  if(nrow(df_test) > 2){
      my_comparisons <- list()
    if( length(levels(df_test[[compare.group]])) > 2) {
      comparisons = combn(levels(df_test[[compare.group]]),2)
      for (cc in 1:ncol(comparisons)) {
        my_comparisons[[cc]] <- as.vector(c(comparisons[1, cc], comparisons[2, cc]))
      }
    } else {
            my_comparisons[[1]] <- as.vector(c(levels(df_test[[compare.group]])))
    }
    my_comparisons 
    
    
    # summary
    summary.df <-df_test %>% 
      group_by((!!sym(compare.group))) %>%
      get_summary_stats(as.name(nn), type = "common")
    summary.results <- rbind(summary.results, summary.df)
    
    # simple t.test
    test.result <- df_test%>%
      # nest() %>%
      # pluck("data", 46) %>% 
      rstatix::t_test(formula(paste(nn, "~", compare.group))) %>%
      adjust_pvalue(method = "bonferroni") %>%
      add_significance("p.adj") %>% 
      add_significance("p") %>%
      mutate(method = "T-test") %>% 
      add_xy_position(x = compare.group) 
    test.result <- test.result %>%
      mutate(y.position = (max(df_test[[nn]], na.rm = T))* seq(1,1.5, 0.5/(length(my_comparisons)-1)) ) 
    
    t.test.results <- rbind(t.test.results, test.result )
    # rm(test.result)
    
    # 1 ANOVA
    ano.result <- df_test %>%
      # as.data.frame() %>%
      rstatix::anova_test(formula(paste(nn, "~", compare.group))) %>%
      adjust_pvalue(method = "bonferroni") %>%
      add_significance("p.adj") %>% 
      add_significance("p") %>%
      ungroup() %>% as_tibble() %>% 
      mutate(Y.value = nn,
             sig = ifelse(.data[["p"]] < 0.05 & !is.na(.data[["p"]]) , "Yes", "ns"),
             Test.type = "1ANOVA",
             file.type.1 = file.name,
             file.type.2 = Annotation.file.name,
             file.type.3 = Type.file.name
      )
    ano.results <- rbind(ano.results, ano.result)
    
    if (TRUE %in% c(ano.result$p <= 0.06) | TRUE %in% c(test.result$p <= 0.06) ) {
      p1 <- ggbarplot(df_test, x = compare.group, y = nn,
                      add = c("mean_se", "jitter"),
                      position = position_dodge(),
                      color = compare.group, 
                      # id = "patient_id", line.color = "gray", line.size = 0.4,
                      palette = "npg",
                      # facet.by = c( "Gender"), short.panel.labs = T,
                      # ncol = length(unique(df_test[["Gender"]])),
                      xlab = "", ylab = y.lab,
                      ylim = c(0, (max(df_test[[nn]], na.rm = T)*1.5))
      )+labs(color=NULL, title = nn)+
        # stat_compare_means(method = "anova", label.y = 90 ) +
        # stat_pvalue_manual(stat.anova, label = "p.adj" , size = 3)+ 
        stat_compare_means(aes_string(group= compare.group),
                           method = "anova",
                           # paired = T,
                           # label = "p.format", size = 4,
                           hide.ns = TRUE,  show.legend = F,
                           label.y = (max(df_test[[nn]], na.rm = T)*1.5))+
        # stat_compare_means(comparisons = my_comparisons, 
        #                    mapping=aes(label = format.pval(..p.adj.., digits = 1)),
        #                    method = "t.test", 
        #                    # paired = T,
        #                    # label = "p.signif", size = 4, 
        #                    label = "p.format", size = 2, 
        #                    hide.ns = TRUE,  show.legend = F, 
        #                    step.increase = ((max(df_test[[vv]], na.rm = T))*0.5)/500,
        #                    # label.y = (max(df_test[[vv]], na.rm = T)*1) 
        #                    ) +
        stat_pvalue_manual( test.result  #%>% dplyr::filter(p.adj <= 0.05)
                            # , label = "p = {p.adj}"
                            , label = "{p.adj}{p.adj.signif}"
                            # , label = "p.adj"
                            , hide.ns = TRUE
                            , size = 3
        )+ 
        theme(
          legend.position="right", # legend.justification="middle",
          axis.text.x=element_blank(),
          #axis.text.x = element_text(angle = 0, hjust = 0.5),
          axis.ticks.x=element_blank()
        ) + guides(color=guide_legend(nrow=length(levels(df_test[[compare.group]])),byrow=TRUE))
      print(p1)
      
      png(filename = paste0(dir1, stringr::str_split(file.variables, "[.]")[[1]][1],  "_", Select.Group, ".", nn, ".barplot.png"), 
          width = (length(levels(df_test[[compare.group]]))/1.5)+ (max(nchar(levels(df_test[[compare.group]])))/8), 
          height = 5, res = 300, units = "in")
      print(p1)
      dev.off()
    }
  }
}

save(ano.results,  t.test.results, summary.results,
     file = paste0(dir0, 
                   stringr::str_split(file.variables, ".Rds")[[1]][1], "_", Select.Group,  "_Results.RData"))

rm(dir1, dir0,
   file.name, Annotation.file.name, Type.file.name,
   nn, numeric.cols,
   all.df, ano.result, ano.results, df_test, p1, sample.meta, summary.df, summary.results,
   t.test.results, test.result, var.df
)








