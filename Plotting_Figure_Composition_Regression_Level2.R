
# rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.
memory.limit(size = 1e+13)

# import parameters and packages --------------------------------------------------------------------------------------------
# source(paste0("/media/jianie/Extreme SSD1/2022-09-01 STARR_SGLT2 Combine/Project Parameters.R"), local = knitr::knit_global())
source(paste0(paste0(Disk, Project.folder, "/", codes.folder, "/", "Project Parameters.R")), local = knitr::knit_global())
source(paste0(paste0(Disk, "00_Functions_Refs/", "Functions_Human Fat snRNA.R")), local = knitr::knit_global() )
library("RColorBrewer")
library("ggpubr")
library(rstatix)

export.folder; Select_Group; Select_Group.names; x.factor; Compare.group ; analysis.group
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
      select_if(~sum(!is.na(.)) > 0)  
    Composition_2 <- Composition_2 %>%
      select( where(is.character) ) %>%
      full_join(Composition_2 %>%
                  tibble::column_to_rownames(var = "Dataset")%>%
                  select( where(is.numeric) ) %>%
                  select(which(!(plyr::numcolwise(sum)(., na.rm=TRUE) %in% 0)))%>%
                  tibble::rownames_to_column(var = "Dataset"), by = "Dataset")
    Composition <- Composition.lists[[list.name]][["Composition_All"]] %>% tibble::rownames_to_column(var = "rowname") %>%
      dplyr::filter( !!sym("rowname") %in% Composition_3[["Dataset"]] ) %>% tibble::column_to_rownames(var = "rowname") %>% 
      select(which(!colSums(., na.rm=TRUE) %in% 0))
    clusters = colnames(Composition)
    # display.brewer.pal(length(clusters), "Set1")
    # Composition_3$x.label <- factor(Composition_3[[label.name]], levels=clusters)
    
    # correlation in Composition_2  ----------
    ## analyze all pairs of correlation by pearson
    Corr.Result.pearson.pairs <- as.data.frame(cor(Composition_2%>% dplyr::select(where(is.numeric)))) # Building a correlation matrix, "pearson" (default)
    Corr.Result.pearson.pairs <- Corr.Result.pearson.pairs[ , colSums(is.na(Corr.Result.pearson.pairs))<(nrow(Corr.Result.pearson.pairs)-1 )]
    Corr.Result.pearson.pairs <- Corr.Result.pearson.pairs[rowSums(is.na(Corr.Result.pearson.pairs))<(ncol(Corr.Result.pearson.pairs))  , ]
    numeric.variables = setdiff(colnames(Corr.Result.pearson.pairs), clusters)
    Multi.Regre.Result.list[[list.name]][[x.name]][["cor"]] <- Corr.Result.pearson.pairs
    
    ## plot all pairs of correlation by pearson
    g <- Composition %>% 
      GGally::ggpairs(., 
                      upper = list(continuous = my_custom_cor),
                      diag = list(continuous = "densityDiag"),
                      lower = list(continuous = my_custom_smooth_color),
                      axisLabels = "none")   # all possible  pairwise plots
    print(g)
    png(file=paste0(dir1, file.name, "_", x.name, "_", "Composition_", analysis.group , ".ggpairs.png"), 
        width=ncol(Composition), height=ncol(Composition), res = 300, units = "in")
    print(g)
    dev.off()
    
    ## plot all pairs of correlation by pearson with clinical variables
    # g <- Composition_2  %>% dplyr::select(colnames(Corr.Result.pearson.pairs)) %>% # clusters, Compare.group
    #   GGally::ggpairs(., 
    #                   upper = list(continuous = my_custom_cor),
    #                   diag = list(continuous = "densityDiag"),
    #                   lower = list(continuous = my_custom_smooth_color),
    #                   axisLabels = "none")   # all possible  pairwise plots
    # print(g)
    # png(file=paste0(dir1, file.name, "_", x.name, "_", "Composition_Clinicals_", analysis.group , ".ggpairs.png"), 
    #     width=length(colnames(Corr.Result.pearson.pairs)), height=length(colnames(Corr.Result.pearson.pairs)), res = 300, units = "in")# c(clusters, Compare.group)
    # print(g)
    # dev.off()
    g <- Composition_2 %>% dplyr::select(clusters, Compare.group) %>% # clusters, Compare.group
      GGally::ggpairs(.,
                      upper = list(continuous = my_custom_cor),
                      diag = list(continuous = "densityDiag"),
                      lower = list(continuous = my_custom_smooth_color),
                      axisLabels = "none")   # all possible  pairwise plots
    print(g)
    png(file=paste0(dir1, file.name, "_", x.name, "_", "Composition_Compare.group_", analysis.group , ".ggpairs.png"),
        width=length(c(clusters, Compare.group)), height=length(c(clusters, Compare.group)), res = 300, units = "in")# c(clusters, Compare.group)
    print(g)
    dev.off()
    
    # plot multiple regression in Composition_2 ---------
    data <- Composition_2
    library(jtools)
    
    if (x.name == label.name){
      for (cell in clusters) {
        overall_p <- function(my_model) {
          f <- summary(my_model)$fstatistic
          p <- pf(f[1],f[2],f[3],lower.tail=F)
          attributes(p) <- NULL
          return(p)
        }
        
        formula2 <- formula(paste0(cell, "~", paste(Compare.group, collapse = " + ") ))
        model2 <- lm(formula2, data = data)
        m2 <- summary(model2)
        #extract overall p-value of model
        overall_p(model2)
        if(overall_p(model2) <= 0.05){
          RP = substitute("R"^"2"*" = "*a *", p"*" = "*p, list(a = round(m2$r.squared, 2), p = round(overall_p(model2), 2)))
          # R2 = substitute("R"^"2"*" = "*a, list(a = round(m2$r.squared, 2)))
          # p2 = substitute("p"*" = "*p, list(p = round(overall_p(model2), 2) ))
          # tag = substitute(" y = "*b*"+("*b1*")", 
          #                  list(b0=round(m2$coefficients[1], 2), 
          #                       b1= paste(paste(round(m2$coefficients[2:(length(Compare.group)+1),1], 2), rownames(m2$coefficients)[2:(length(Compare.group)+1)], sep = "*"), collapse = ")+(") ))
          
          # effect_plot
          gg_color_hue <- function(n) {
            hues = seq(15, 375, length = n +1 )
            hcl(h = hues, l = 65, c = 100)[1:n]
          }
          n = length(model2$residuals)
          cols = gg_color_hue(n)
          
          graphics.off()
          effect_plot(model2, pred = !!x.factor, 
                      interval = TRUE, plot.points = TRUE, colors = "darkgray", robust = "HC2", 
                      point.color = cols, point.alpha =1, point.size = 3, y.label = "Average nuclei percentage") +
            labs(title=cell, tag = RP ) +
            theme(panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
                  plot.title = element_text(color = "black", size = 14, face = "bold", hjust = 0.5),
                  axis.title.x = element_text(color = "black", size = 13, face = "plain"),
                  axis.title.y = element_text(color = "black", size = 13, face = "plain"),
                  plot.tag = element_text(color = "black", size = 12, face = "plain"),
                  plot.tag.position = c(0.8, 0.2) )
          ggsave(filename = paste0(dir1, file.name, "_", x.name, "_", cell, "_Composition_", x.factor, ".", length(Compare.group), "_",analysis.group ,"_effect_plot.png"), width = 4, height = 4)
          
          
          
          # plot_summs and plot_coefs
          ## Plot coefficient uncertainty as normal distributions
          library("broom.mixed")
          print("model details")
          jtools::plot_summs(model2, scale = TRUE, plot.distributions = TRUE, inner_ci_level = .9) #
          summ(model2, scale = TRUE) # Standardized/scaled coefficients
          summ(model2, center = TRUE) # Mean-centered variables
          summ(model2, confint = TRUE, digits = 3) # Confidence intervals
          summ(model2, vifs = TRUE, scale = TRUE, center = TRUE)
          Publish::publish(model2)
          
          ## Comparing model coefficients visually
          graphics.off()
          plot_summs(model2, scale = TRUE, robust = "HC2", model.names = c(paste(Compare.group, collapse = " + ") ))+
            theme(panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
                  legend.position = "bottom")
          # plot_coefs(model1, model2, scale = TRUE, model.names = c("model1", "model2"))
          ggsave(filename = paste0(dir1, file.name, "_", x.name, "_", cell, "_Composition_", length(Compare.group), "_",analysis.group , ".Multi.Corr.Result.png"), width = 5, height = length(Compare.group)/2)
          
          # Table output for Word and RMarkdown documents
          library(huxtable)
          library(officer)
          library(flextable)
          set_summ_defaults(digits = 3, pvals = TRUE, robust = "HC2", confint = TRUE)
          export_summs(model2, scale = TRUE, robust = "HC2", to.file = "docx", 
                       model.names = c( paste(Compare.group, collapse = " + ")),
                       statistics = c("N" = "nobs", 
                                      "R2" = "r.squared",
                                      "Adjusted R2" = "adj.r.squared", 
                                      "Residual standard error (RSE)" = "sigma",
                                      "F-statistic" = "statistic", 
                                      "P-value" = "p.value"),
                       file.name = paste0(dir1, file.name, "_", x.name, "_", cell, "_Composition_", length(Compare.group), "_",analysis.group , ".Multi.Corr.Result.docx")) # , error_format = "[{conf.low}, {conf.high}]"
          
          export_summs(model2, scale = TRUE, robust = "HC2", to.file = "HTML", 
                       model.names = c( paste(Compare.group, collapse = " + ")),
                       statistics = c("N" = "nobs", 
                                      "R2" = "r.squared",
                                      "Adjusted R2" = "adj.r.squared",
                                      "Residual standard error (RSE)" = "sigma",
                                      "F-statistic" = "statistic", 
                                      "P-value" = "p.value"),
                       file.name = paste0(dir1, file.name, "_", x.name, "_", cell, "_Composition_",  length(Compare.group), "_",analysis.group , ".Multi.Corr.Result.HTML")) # , error_format = "[{conf.low}, {conf.high}]"
          # export_summs(model1, model2, scale = TRUE, error_format = "[{conf.low}, {conf.high}]", to.file = "pdf", file.name = "Myeloid.Multi.Correlation.pdf")
          
        }
      }
    }
    
   
    }
}

# rm(dir1, dir0,  files, file,
#    file.name, Annotation.file.name, Subtype.file.name,
#    celltype.name, list.name, label.name, x.name,
#    Composition_3, Composition,
#    clusters.lists, Composition.lists, clusters
# )
# rm(bartfit, Composition_2, Corr.Result.pearson.pairs, data,
#    g, Import.list, m2, model2, mse.df,
#                   Multi.Regre.Result.list, per.var.df, rf.nhanes,
#                   Test.list, x, xtest, xtrain,
#    cell, cols, cols.names, formula2, n, ntest, numeric.variables,
#    ord, RP, testid, trainid, uu, y, y.pred, yhat.bart, ytest, ytrain
#    )
gc()
sessionInfo()






