library("RColorBrewer")
library("ggpubr")

# rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.
memory.limit(size = 1e+13)

# import parameters and packages --------------------------------------------------------------------------------------------
# source(paste0("/media/jianie/Extreme SSD1/2022-09-01 STARR_SGLT2 Combine/Project Parameters.R"), local = knitr::knit_global())
source(paste0(paste0(Disk, Project.folder, "/", "Project Parameters.R")), local = knitr::knit_global())

# Import data -----------------------------------------------------------------------------------------------------------------------------------------------------
# export.folder <- "Compositions"
# dir.create(file.path(paste0(Disk, Project.folder, "/", Rds.folder, "/"), export.folder), showWarnings = FALSE)
dir0 <- paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder, "/")
files <- list.files(dir0)

dir.create(file.path(paste0(Disk, Project.folder, "/", figures.folder, "/", export.folder, ".plotting")), showWarnings = FALSE)
dir1 <- paste0(Disk, Project.folder, "/", figures.folder, "/", export.folder, ".plotting", "/")


# Figure 1I
load(paste0(Disk, Project.folder, "/", Rds.folder, "/","Figure 1D_Composition_list.RData"))
# Composition_3$Dataset <- factor(Composition_3$Dataset, levels = paste("Subject", 1:nrow(Composition_3)))
Composition <- Composition.lists$SAT_Annotation$Composition_All
Composition_3 <- Composition.lists$SAT_Annotation$Composition_3All

# plot table numbers
library(readxl)
cell_types_file <- read_excel(paste0(Disk, Project.folder, "/", "Fat_Subtype.xlsx"), sheet = "SAT")

colnames(Composition)
clusters = cell_types_file$Annotation
display.brewer.pal(length(clusters), "Set1")
Composition_3$Annotation <- factor(Composition_3$Annotation, levels=clusters)

# plotting
p1 <- ggplot(Composition_3, aes(x = Annotation, y = composition, group = Annotation)) + 
  geom_boxplot( size=0.5) +
  geom_point(aes(color=Dataset), size=2) +
  labs(title="Cell type distribution", y = "Average percentage of all nuclei") +
  ylim(0, max(Composition_3$composition)) +
  theme_classic()+
  theme(panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "black", size = 13, face = "plain"),
        axis.text.x = element_text(color = "black", size = 11.5, face = "plain", angle = 45, hjust = 1),
        plot.title = element_text(color = "black", size = 14, face = "bold", hjust = 0.5),
        #legend.text = element_text(color = "black", size = 10.5, face = "plain"),
        #legend.position = c(0.9, 0.8),
        #legend.position = c(.99, .99),
        legend.position = "none",
        #legend.justification = c("right", "top"),
        #legend.box.just = "right",
        legend.spacing.y = unit(0.15, 'cm') ,
        legend.key.size = unit(0.45, "cm"),
        legend.background = element_rect( fill = "grey98", color = "grey98", size = 1)
  ) 
p1
png(file=paste0(Disk, Project.folder, "/", figures.folder, "/", "Statistic_Plotting_Overall/Figure 1D_Overall_Statistic_Composition_Box.png"), 
    width=3.7, height=4, res = 300, units = "in")
p1
dev.off()

 

# plotting Multiple Linear Regression 
load(paste0(Disk, Project.folder, "/", Rds.folder, "/","Composition_All.RData"))

library(ggrepel)
all <- data.frame()
for (i in 1:length(Multi.Regre.Result.all[["Multi.Result.clean.list"]])) { all <- rbind(all, Multi.Regre.Result.all[["Multi.Result.clean.list"]][[i]]) }

# plotting Multiple Linear Regression for all
library(dplyr)
List <- data.frame(`F-Statistic` = sort(unique(all$`F-Statistic`)), "List" = length(unique(all$`F-Statistic`)):1 )
all2 <- all
for (r in 1:nrow(all)) {all2$decile_rank[r] <- List[which(List$F.Statistic == all$`F-Statistic`[r]), "List"]} 
highlight <- all2[all2$'overall p.value' < 0.05 & all2$`R2` > 0.8, ] %>% group_by(Outcome_variable_y) %>% dplyr::filter(`F-Statistic` %in% c(max(`F-Statistic`), max( `F-Statistic`[`F-Statistic`!=max(`F-Statistic`)] )))
# highlight <- all2[all2$'overall p.value' < 0.05 & all2$`RSE error rate` < 0.5, ] %>% group_by(Outcome_variable_y) %>% dplyr::filter(`F-Statistic` %in% c(max(`F-Statistic`), max( `F-Statistic`[`F-Statistic`!=max(`F-Statistic`)] )))
library(ggrepel)
p2 <- ggscatter(all, x= "RSE error rate", y= "F-Statistic",
                color = "Outcome_variable_y", shape = 20, 
                size= c(ifelse((all$'overall p.value' < 0.05), 4, 1)), # Points color, shape and size
                #title = paste(names(Multi.Result.clean.list), collapse = "_"),
                title = "Average fraction of all nuclei",
                xlab = "Error Rate",
                ylab = paste0("F-Statistic")
) + xlim(0, max(all$`RSE error rate`)) + ylim(0, max(all$`F-Statistic`) )+ labs(color = "Cell Types")+
  geom_point(data=highlight, aes(x=`RSE error rate`, y=`F-Statistic`, colour= `Outcome_variable_y`), 
             size= c(ifelse((highlight$p.value < 0.05), 1.8, 0.1))) +
  geom_text_repel(data=highlight, aes(x=`RSE error rate`, y=`F-Statistic`, 
                                      label = .data[["Predictor_Compare.group"]],
                                      color = `Outcome_variable_y`,
                                      #colour= c(ifelse((highlight$p.value < 0.05), "Significant", `Outcome_variable_y`)) 
                                      ), 
                  size = 4.2, fontface="bold", box.padding = 1, max.overlaps = Inf ) +
  # geom_hline(yintercept = 5.5, linetype = "dashed", size=0.2) +
  geom_vline(xintercept = all[all$`F-Statistic` == max(all$`F-Statistic`),]$`RSE error rate`, linetype = "dashed", size=0.2) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
        legend.text = element_text(color = "black", size = 12, face = "plain"),
        legend.position = c(.79, .62),
        legend.background = element_rect( fill = "grey98", color = "grey98", size = 1),
        legend.spacing.y = unit(0.19, 'cm') ,
        legend.key.size = unit(0.45, "cm"),
        #legend.position = "right"
        ) 
p2
#ggsave(filename = paste0(Disk, Project.folder, "/", figures.folder, "/", "Statistic_Plotting_Overall/", "Overall.Multi.Correlation_scatter.png"), plot = p2, width = 4.6224, height = 4  )

# highlight <- all2[all2$'overall p.value' < 0.05 & all2$`R2` > 0.75, ] %>% group_by(Outcome_variable_y) %>% dplyr::filter(`mlogpvalue` %in% c(max(`mlogpvalue`)) | `R2` %in% c(max(`R2`)) )
library(ggrepel)
p3 <- ggscatter(all, x= "RSE error rate", y= "adj.R2",
                color = "Outcome_variable_y", shape = 20, 
                size= c(ifelse((all$'overall p.value' < 0.05), 4, 1)), # Points color, shape and size
                #title = paste(names(Multi.Result.clean.list), collapse = "_"),
                title = "Average fraction of all nuclei",
                xlab = "Error Rate",
                ylab = paste0("adj.R2") # expression(-log[10] ~ pValue)
) + xlim(0, max(all$`RSE error rate`)) + ylim(0, max(all$`adj.R2`) )+ labs(color = "Cell Types")+
  geom_point(data=highlight, aes(x=`RSE error rate`, y=`adj.R2`, colour= `Outcome_variable_y`), 
             size= c(ifelse((highlight$p.value < 0.05), 1.8, 0.1))) +
  geom_text_repel(data=highlight, aes(x=`RSE error rate`, y=`adj.R2`, 
                                      label = .data[["Predictor_Compare.group"]],
                                      color = `Outcome_variable_y`,
                                      #colour= c(ifelse((highlight$p.value < 0.05), "Significant", `Outcome_variable_y`)) 
  ), 
  size = 4.2, fontface="bold", box.padding = 1, max.overlaps = Inf ) +
  # geom_hline(yintercept = -log10(0.05), linetype = "dashed", size=0.2) +
  geom_vline(xintercept = all[all$`adj.R2` == max(all$`adj.R2`),]$`RSE error rate`, linetype = "dashed", size=0.2) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
        legend.text = element_text(color = "black", size = 12, face = "plain"),
        legend.position = c(.29, .62),
        legend.background = element_rect( fill = "grey98", color = "grey98", size = 1),
        legend.spacing.y = unit(0.19, 'cm') ,
        legend.key.size = unit(0.45, "cm"),
        #legend.position = "right"
  ) 
p3
ggsave(filename = paste0(Disk, Project.folder, "/", figures.folder, "/", "Statistic_Plotting_Overall/", "Overall.Multi.Correlation_scatter.png"), plot = p3, width = 4.6224, height = 4  )
write.csv(as.data.frame(all2), file = paste0(Disk, Project.folder, "/", figures.folder, "/", "Statistic_Plotting_Overall/All2.csv"), quote = FALSE)


# Export multiple regression table
data <- Composition_2
library(jtools)
cells = unique(highlight$Outcome_variable_y)

for (cell in cells) {
  highlight1 = highlight[as.character(highlight$Outcome_variable_y) == as.character(cell), ]
  F.values = as.character(sort(unique(highlight1$`F-Statistic`)))
  
  f = F.values[1]
  df.highlight = highlight1[as.character(highlight1$`F-Statistic`) == as.character(f), ]
  factors2 = df.highlight[order(df.highlight$p.value), ]$Predictor_variables_x
  formula2 <- formula(paste0(cell, "~", paste(factors2, collapse = " + ") ))
  
  model2 <- lm(formula2, data = data)
  summ(model2, scale = TRUE) # Standardized/scaled coefficients
  summ(model2, center = TRUE) # Mean-centered variables
  summ(model2, confint = TRUE, digits = 3) # Confidence intervals
  m2 <- summary(model2)
  # summ(model2, vifs = TRUE, scale = TRUE, center = TRUE)
  R2 = substitute("R"^"2"*" = "*a, list(a = round(m2$r.squared, 2)))
  tag = substitute(" y = "*b*"+("*b1*")", 
                   list(b0=round(m2$coefficients[1], 2), 
                        b1= paste(paste(round(m2$coefficients[2:(nrow(df.highlight)+1),1], 2), rownames(m2$coefficients)[2:(nrow(df.highlight)+1)], sep = "*"), collapse = ")+(") ))
  
  # effect_plot
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n +1 )
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  n = length(model2$residuals)
  cols = gg_color_hue(n)
  for (b in 1:length(factors2)) {
    effect_plot(model2, pred = !!factors2[b], interval = TRUE, plot.points = TRUE, colors = "darkgray", robust = "HC2", 
                point.color = cols, point.alpha =1, point.size = 3, y.label = "Average percentage of all nuclei") +
      labs(title=cell, tag = R2 ) +
      theme(panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
            plot.title = element_text(color = "black", size = 14, face = "bold", hjust = 0.5),
            axis.title.x = element_text(color = "black", size = 13, face = "plain"),
            axis.title.y = element_text(color = "black", size = 13, face = "plain"),
            plot.tag = element_text(color = "black", size = 12, face = "plain"),
            plot.tag.position = c(0.8, 0.3) )
    ggsave(filename = paste0(Disk, Project.folder, "/", figures.folder, "/", "Statistic_Plotting_Overall/", cell, "_2.", factors2[b], ".Multi.Correlation.effect_plot.png"), width = 4, height = 4)
  }
  
  # plot_summs and plot_coefs
  ## Plot coefficient uncertainty as normal distributions
  jtools::plot_summs(model2, scale = TRUE, plot.distributions = TRUE, inner_ci_level = .9) #
  plot_summs(model2, scale = TRUE, plot.distributions = TRUE, inner_ci_level = .9) 
  plot_summs(model2, scale = TRUE, plot.distributions = TRUE)
  plot_summs(model2, scale = TRUE, plot.distributions = TRUE)
  ## Comparing model coefficients visually
  plot_summs(model2, scale = TRUE, robust = "HC2", model.names = c(paste(factors2, collapse = " + ") ))+
    theme(panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
          legend.position = "bottom")
  # plot_coefs(model1, model2, scale = TRUE, model.names = c("model1", "model2"))
  ggsave(filename = paste0(Disk, Project.folder, "/", figures.folder, "/", "Statistic_Plotting_Overall/", cell, ".", length(factors2), ".Multi.Correlation.coefficient.plotting.png"), width = 5, height = 2.5)
  
  # Table output for Word and RMarkdown documents
  set_summ_defaults(digits = 3, pvals = TRUE, robust = "HC2", confint = TRUE)
  export_summs(model2, scale = TRUE, robust = "HC2", to.file = "docx", 
               model.names = c( paste(factors2, collapse = " + ")),
               statistics = c("N" = "nobs", 
                              "R2" = "r.squared",
                              "Adjusted R2" = "adj.r.squared", 
                              "Residual standard error (RSE)" = "sigma",
                              "F-statistic" = "statistic", 
                              "P-value" = "p.value"),
               file.name = paste0(Disk, Project.folder, "/", figures.folder, "/", "Statistic_Plotting_Overall/", cell, ".", length(factors2), ".Multi.Correlation.docx")) # , error_format = "[{conf.low}, {conf.high}]"
  
  export_summs(model2, scale = TRUE, robust = "HC2", to.file = "HTML", 
               model.names = c( paste(factors2, collapse = " + ")),
               statistics = c("N" = "nobs", 
                              "R2" = "r.squared",
                              "Adjusted R2" = "adj.r.squared",
                              "Residual standard error (RSE)" = "sigma",
                              "F-statistic" = "statistic", 
                              "P-value" = "p.value"),
               file.name = paste0(Disk, Project.folder, "/", figures.folder, "/", "Statistic_Plotting_Overall/", cell, ".", length(factors2), ".Multi.Correlation.HTML")) # , error_format = "[{conf.low}, {conf.high}]"
  # export_summs(model1, model2, scale = TRUE, error_format = "[{conf.low}, {conf.high}]", to.file = "pdf", file.name = "Myeloid.Multi.Correlation.png")
  
  if (length(F.values) > 1 ){
    f = F.values[2]
    df.highlight = highlight1[as.character(highlight1$`F-Statistic`) == as.character(f), ]
    factors1 = df.highlight[order(df.highlight$p.value), ]$Predictor_variables_x 
    formula1 <- formula(paste0(cell, "~", paste(factors1, collapse = " + ") ))
    model1 <- lm(formula1, data = data)
    summ(model1, scale = TRUE)
    summ(model1, scale = TRUE, center = TRUE) 
    m1 <- summary(model1)
    R2 = substitute("R"^"2"*" = "*a, list(a = round(m1$r.squared, 2)))
    tag = substitute(" y = "*b*"+("*b1*")", 
                     list(b0=round(m1$coefficients[1], 2), 
                          b1= paste(paste(round(m1$coefficients[2:(nrow(df.highlight)+1),1], 2), rownames(m1$coefficients)[2:(nrow(df.highlight)+1)], sep = "*"), collapse = ")+(") ))
    
    # effect_plot
    n = length(model1$residuals)
    cols = gg_color_hue(n)
    for (b in 1:length(factors1)) {
      effect_plot(model1, pred = !!factors1[b], interval = TRUE, plot.points = TRUE, colors = "darkgray", robust = "HC2", 
                  point.color = cols, point.alpha =1, point.size = 3, y.label = "Average percentage of all nuclei") +
        labs(title=cell, tag = R2 ) +
        theme(panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
              plot.title = element_text(color = "black", size = 14, face = "bold", hjust = 0.5),
              axis.title.x = element_text(color = "black", size = 13, face = "plain"),
              axis.title.y = element_text(color = "black", size = 13, face = "plain"),
              plot.tag = element_text(color = "black", size = 12, face = "plain"),
              plot.tag.position = c(0.8, 0.3) )
      ggplot2::ggsave(filename = paste0(Disk, Project.folder, "/", figures.folder, "/", "Statistic_Plotting_Overall/", cell, "_1.", factors1[b], ".Multi.Correlation.effect_plot.png"), width = 4, height = 4)
    }
    
    # plot_summs and plot_coefs
    ## Plot coefficient uncertainty as normal distributions
    jtools::plot_summs(model1, scale = TRUE, plot.distributions = TRUE, inner_ci_level = .9) #
    plot_summs(model2, scale = TRUE, plot.distributions = TRUE, inner_ci_level = .9) 
    plot_summs(model2, scale = TRUE, plot.distributions = TRUE)
    plot_summs(model1, model2, scale = TRUE, plot.distributions = TRUE)
    ## Comparing model coefficients visually
    plot_summs(model1, model2, scale = TRUE, robust = "HC2", model.names = c(paste(factors1, collapse = " + "), paste(factors2, collapse = " + ") ))+
      theme(panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
            legend.position = "bottom") + guides(color=guide_legend(nrow=2, byrow=TRUE))
    # plot_coefs(model1, model2, scale = TRUE, model.names = c("model1", "model2"))
    ggsave(filename = paste0(Disk, Project.folder, "/", figures.folder, "/", "Statistic_Plotting_Overall/", cell, ".", nrow(highlight), ".Multi.Correlation.coefficient.plotting.png"), width = 5, height = 2.5)
    
    # Table output for Word and RMarkdown documents
    set_summ_defaults(digits = 3, pvals = TRUE, robust = "HC2", confint = TRUE)
    export_summs(model1, model2, scale = TRUE, robust = "HC2", to.file = "docx", 
                 model.names = c(paste(factors1, collapse = " + "), paste(factors2, collapse = " + ")),
                 statistics = c("N" = "nobs", 
                                "R2" = "r.squared",
                                "Adjusted R2" = "adj.r.squared", 
                                "Residual standard error (RSE)" = "sigma",
                                "F-statistic" = "statistic", 
                                "P-value" = "p.value"),
                 file.name = paste0(Disk, Project.folder, "/", figures.folder, "/", "Statistic_Plotting_Overall/", cell, ".", nrow(highlight), ".Multi.Correlation.docx")) # , error_format = "[{conf.low}, {conf.high}]"
    
    export_summs(model1, model2, scale = TRUE, robust = "HC2", to.file = "HTML", 
                 model.names = c(paste(factors1, collapse = " + "), paste(factors2, collapse = " + ")),
                 statistics = c("N" = "nobs", 
                                "R2" = "r.squared",
                                "Adjusted R2" = "adj.r.squared",
                                "Residual standard error (RSE)" = "sigma",
                                "F-statistic" = "statistic", 
                                "P-value" = "p.value"),
                 file.name = paste0(Disk, Project.folder, "/", figures.folder, "/", "Statistic_Plotting_Overall/", cell, ".", nrow(highlight), ".Multi.Correlation.HTML")) # , error_format = "[{conf.low}, {conf.high}]"
    # export_summs(model1, model2, scale = TRUE, error_format = "[{conf.low}, {conf.high}]", to.file = "pdf", file.name = "Myeloid.Multi.Correlation.png")
  }
}

  

# plotting Simple Linear Regression for all
library(ggrepel)
Corr.Result.clean <- Corr.Result[!Corr.Result$Predictor_variables_x == "(Intercept)", ]
all <- Corr.Result.clean

library(dplyr)
List <- data.frame(`mlogpvalue` = sort(unique(all$`mlogpvalue`)), "List" = length(unique(all$`mlogpvalue`)):1 )
all2 <- all
for (r in 1:nrow(all)) {all2$decile_rank[r] <- List[which(List$mlogpvalue == all$`mlogpvalue`[r]), "List"]} 
highlight <- all2[all2$'overall p.value' < 0.05, ] # %>% group_by(Outcome_variable_y) %>% dplyr::filter(`F-Statistic` %in% c(max(`F-Statistic`), max( `F-Statistic`[`F-Statistic`!=max(`F-Statistic`)] )))
library(ggrepel)
p4 <- ggscatter(all2, x= "RSE error rate", y= "F-Statistic",
                color = "Outcome_variable_y", shape = 20, 
                size= c(ifelse((all2$'overall p.value' < 0.05), 4, 1)), # Points color, shape and size
                #title = paste(names(Multi.Result.clean.list), collapse = "_"),
                title = "Average fraction of all nuclei",
                xlab = "Error rate",
                ylab = paste0("F-Statistic")
) + xlim(0, max(all2$`RSE error rate`)) + ylim(0, max(all2$`F-Statistic`) )+ labs(color = "Cell Types")+
  geom_point(data=highlight, aes(x=`RSE error rate`, y=`F-Statistic`, colour= `Outcome_variable_y`), 
             size= c(ifelse((highlight$p.value < 0.05), 1.8, 0.1))) +
  geom_text_repel(data=highlight, aes(x=`RSE error rate`, y=`F-Statistic`, 
                                      label = .data[["Predictor_variables_x"]],
                                      color = `Outcome_variable_y`,
                                      #colour= c(ifelse((highlight$p.value < 0.05), "Significant", `Outcome_variable_y`)) 
  ), 
  size = 4.2, fontface="bold", box.padding = 1, max.overlaps = Inf ) +
  # geom_hline(yintercept = 5.5, linetype = "dashed", size=0.2) +
  geom_vline(xintercept = all2[all2$`F-Statistic` == max(all2$`F-Statistic`),]$`RSE error rate`, linetype = "dashed", size=0.2) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
        legend.text = element_text(color = "black", size = 12, face = "plain"),
        legend.position = c(.29, .62),
        legend.background = element_rect( fill = "grey98", color = "grey98", size = 1),
        legend.spacing.y = unit(0.19, 'cm') ,
        legend.key.size = unit(0.45, "cm"),
        #legend.position = "right"
  ) 
p4

highlight <- all2[all2$'overall p.value' < 0.05 & all2$`R2` > 0.6, ] %>% group_by(Outcome_variable_y) %>% dplyr::filter(`mlogpvalue` %in% c(max(`mlogpvalue`)) | `R2` %in% c(max(`R2`)) )
p5 <- ggscatter(all2, x= "R2", y= "mlogpvalue",
                color = "Outcome_variable_y", shape = 20, 
                size= c(ifelse((all2$'overall p.value' < 0.05), 4, 1)), # Points color, shape and size
                #title = paste(names(Multi.Result.clean.list), collapse = "_"),
                title = "Average fraction of all nuclei",
                xlab = "R2",
                ylab = expression(-log[10] ~ pValue)
) + xlim(0, max(all2$`R2`)) + ylim(0, max(all2$`mlogpvalue`) )+ labs(color = "Cell Types")+
  geom_point(data=highlight, aes(x=`R2`, y=`mlogpvalue`, colour= `Outcome_variable_y`), 
             size= c(ifelse((highlight$p.value < 0.05), 1.8, 0.1))) +
  geom_text_repel(data=highlight, aes(x=`R2`, y=`mlogpvalue`, 
                                      label = .data[["Predictor_variables_x"]],
                                      color = `Outcome_variable_y`,
                                      #colour= c(ifelse((highlight$p.value < 0.05), "Significant", `Outcome_variable_y`)) 
  ), 
  size = 4.2, fontface="bold", box.padding = 1, max.overlaps = Inf ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", size=0.2) +
  geom_vline(xintercept = all2[all2$`mlogpvalue` == max(all2$`mlogpvalue`),]$`R2`, linetype = "dashed", size=0.2) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
        legend.text = element_text(color = "black", size = 12, face = "plain"),
        legend.position = c(.29, .62),
        legend.background = element_rect( fill = "grey98", color = "grey98", size = 1),
        legend.spacing.y = unit(0.19, 'cm') ,
        legend.key.size = unit(0.45, "cm"),
        #legend.position = "right"
  ) 
p5
ggsave(filename = paste0(Disk, Project.folder, "/", figures.folder, "/", "Statistic_Plotting_Overall/", "Overall.Simple.Correlation_scatter.png"), plot = p5, width = 4.6224, height = 4  )
write.csv(as.data.frame(all2), file = paste0(Disk, Project.folder, "/", figures.folder, "/", "Statistic_Plotting_Overall/All2.Simple..csv"), quote = FALSE)



# plot box plotting for gender effect
p <- ggplot(Composition_3, aes(x = Annotation, y = composition, fill=Gender, 
                               color = Gender, group = interaction(Annotation, Gender))) + 
  geom_boxplot(alpha = 0.1, width=0.75, size=0.5) +
  geom_point(aes(color=Gender), size=2, position = position_dodge(width=0.75)) +
  labs(title="SAT cell type distribution", y = "Average fractioin of all nuclei") +
  ylim(0, max(Composition_3$composition)*1.1) +
  theme_classic()+
  theme(panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "black", size = 13, face = "plain"),
        axis.text.x = element_text(color = "black", size = 11.5, face = "plain", angle = 45, hjust = 1),
        plot.title = element_text(color = "black", size = 14, face = "bold", hjust = 0.5),
        legend.text = element_text(color = "black", size = 10.5, face = "plain"),
        #legend.position = c(0.9, 0.8),
        legend.position = c(.99, .99),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.spacing.y = unit(0.15, 'cm') ,
        legend.key.size = unit(0.45, "cm"),
        legend.background = element_rect( fill = "grey98", color = "grey98", size = 1)
  ) 
lable = Composition_3 %>% group_by(Annotation) %>% select(Annotation,composition) %>% summarise_if(is.numeric, max) %>% as.data.frame() %>% pull(composition)
p
p + stat_compare_means(show.legend=FALSE, aes(label = ..p.signif..), method = "t.test",
                       label.y = lable, hide.ns = TRUE, size = 10)
png(file=paste0(Disk, Project.folder, "/", figures.folder, "/", "Statistic_Plotting_Overall/Figure 1D_Overall_Gender_Composition.png"), 
    width=4, height=4)
p 
dev.off()


# plot table numbers
library(readxl)
cell_types_file <- read_excel("D:/2021-01-11 Figures for Grant resubmission/Fat_Subtype.xlsx", sheet = "Overall")


# plot table numbers
#Plot your table with table Grob in the library(gridExtra)
library(gridExtra)
df <- Cell_numbers
df$Sum <- rowSums(df)
n = nrow(df)+1
df[n, ] <- colSums(df)
rownames(df)[which(rownames(df) == n)] <- "Total"
df <- format(df, big.mark=",", scientific=FALSE)
df$Dataset <- rownames(df)
colnames(df)
df <- df[, c("Dataset", cell_types_file$Annotation, "Sum")]
ss <- tableGrob(df, rows = NULL)
n1 <- nchar(paste(colnames(df), collapse = ""))
n2 <- max(nchar(rownames(df)))
n3 <- nrow(df) + 1
plot(ss)
plot <- cowplot::plot_grid(ss)
plot
png(file=paste0(Disk, Project.folder, "/", figures.folder, "/", "Statistic_Plotting_Overall/Figure 1D_Overall_Cell Numbers.png"), 
    width= (0.09*n1)+1.25+0.35, height= n3/3.5)
plot
dev.off()

# df <- round(Composition*100, 1)
df <- round(Composition, 1)
# df$Sum <- rowSums(df)
df$Sum <- 100
n = nrow(df)+1
df[n, ] <- colSums(df)
rownames(df)[which(rownames(df) == n)] <- "Total"
df <- round(df, 1)
df$Dataset <- rownames(df)
colnames(df)
df <- df[, c("Dataset", cell_types_file$Annotation, "Sum")]
ss <- tableGrob(df, rows = NULL)
n1 <- nchar(paste(colnames(df), collapse = ""))
n2 <- max(nchar(rownames(df)))
n3 <- nrow(df) + 1
plot(ss)
plot <- cowplot::plot_grid(ss)
plot
png(file=paste0(Disk, Project.folder, "/", figures.folder, "/", "Statistic_Plotting_Overall/Figure 1D_Overall_Composition table.png"), 
    width= (0.09*n1)+1.25+0.12, height= n3/3.5)
plot
dev.off()

df <- round(p.Corr.Result_p,3)
df$'Cell Types' <- rownames(df)
colnames(df)
df <- df[, c("Cell Types", "Age_R", "Age_p.value", "BMI_R", "BMI_p.value", 
             "FastingGlucose_R", "FastingGlucose_p.value", "HbA1c_R", "HbA1c_p.value",
             "Insulin_R", "Insulin_p.value")]
ss <- tableGrob(df, rows = NULL)
n1 <- nchar(paste(colnames(df), collapse = ""))
n2 <- max(nchar(rownames(df)))
n3 <- nrow(df) + 1
plot(ss)
plot <- cowplot::plot_grid(ss)
plot
png(file=paste0(Disk, Project.folder, "/", figures.folder, "/", "Statistic_Plotting_Overall/Figure 1D_Overall_Composition.Correlation.png"), 
    width= (0.09*n1)+1.25 +0.8, height= n3/3.5)
plot
dev.off()




# Senescence.stage ------------------------------------------------------
Composition_3 <- Composition.lists$SAT_Senescence.stage$Composition_3All

new.meta.data <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/","CellCycle.Senescence.stage.metadata.Rds"))  %>%  
  tibble::rownames_to_column(var = "rowname")
new.meta.data1 <- new.meta.data %>% select(-one_of(c(intersect(colnames(new.meta.data), colnames(Composition_3)) ))) %>% mutate(subject_id = new.meta.data$subject_id)
new.meta.data.sample <- new.meta.data1[!duplicated(new.meta.data1$subject_id),]
meta_data_merge <- Composition_3 %>% left_join(new.meta.data.sample, by = c("subject_id")) 

p1 <- ggplot(meta_data_merge, aes(x = Senescence.stage, y = composition, group = Senescence.stage)) + 
  geom_boxplot(aes(x = Senescence.stage), size=0.5) +
  geom_point(aes(color=Dataset), size=2) +
  labs(title="", y = "Average percentage of all nuclei") +
  ylim(0, max(meta_data_merge$composition)) +
  theme_classic()+
  theme(panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "black", size = 13, face = "plain"),
        axis.text.x = element_text(color = "black", size = 11.5, face = "plain", angle = 45, hjust = 1),
        plot.title = element_text(color = "black", size = 14, face = "bold", hjust = 0.5),
        #legend.text = element_text(color = "black", size = 10.5, face = "plain"),
        #legend.position = c(0.9, 0.8),
        #legend.position = c(.99, .99),
        legend.position = "none",
        #legend.justification = c("right", "top"),
        #legend.box.just = "right",
        legend.spacing.y = unit(0.15, 'cm') ,
        legend.key.size = unit(0.45, "cm"),
        legend.background = element_rect( fill = "grey98", color = "grey98", size = 1)
  ) 
p1
png(file=paste0(Disk, Project.folder, "/", figures.folder, "/", "Figure 1D_Overall_Senescence.stage_Composition_Box.png"), 
    width=3.7, height=4, res = 300, units = "in")
p1
dev.off()









