
# rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.
memory.limit(size = 1e+13)

# import parameters and packages --------------------------------------------------------------------------------------------
# source(paste0("/media/jianie/Extreme SSD1/2022-09-01 STARR_SGLT2 Combine/Project Parameters.R"), local = knitr::knit_global())
source(paste0(paste0(Disk, Project.folder, "/", "Project Parameters.R")), local = knitr::knit_global())
library("RColorBrewer")
library("ggpubr")

# Import data -----------------------------------------------------------------------------------------------------------------------------------------------------
dir.create(file.path(paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder)), showWarnings = FALSE)
dir0 <- paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder,  "/")


# Import metadata ------------------
sample.meta <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/", Sample.meta.file)) %>% 
  tibble::rownames_to_column(var = "rowname")
sample.meta <- sample.meta[ duplicated(sample.meta[["Dataset"]])==F,] 
sample.meta <- sample.meta[ duplicated(sample.meta[["sample_id"]])==F,] 
sample.meta <- sample.meta %>% 
  dplyr::filter((!!sym(Select.Group)) == T)
sample.meta <- sample.meta[,colSums(is.na(sample.meta))<nrow(sample.meta)] # Remove columns from dataframe where ALL values are NA
colnames(sample.meta)
sample.meta[["Dataset"]]
numeric.cols.c <- sample.meta %>% dplyr::select(where(is.numeric)) %>% 
  dplyr::select(which(!colSums(., na.rm=TRUE) %in% 0)) %>% colnames()

# Import Variable data 1 ------------------
file.name <- stringr::str_split(file.variables, ".Rds")[[1]][1]
file.name
Annotation.file.name <- paste(stringr::str_split(file.name, "_")[[1]][c(1,2,3)], collapse = "_")
Annotation.file.name 
Type.file.name <- paste(stringr::str_split(file.name, "_")[[1]][c(4)], collapse = "_")
Type.file.name 
var.df <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/", file.variables))  # contain "Dataset"
var.df <- var.df[,colSums(is.na(var.df))<nrow(var.df)] 
colnames(var.df)
var.df[["Dataset"]]

# Import Variable data 2 ------------------
var.df2 <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/", file.variables2))  # contain "Dataset"
var.df2 <- var.df2[,colSums(is.na(var.df2))<nrow(var.df2)] 
colnames(var.df2)
var.df2[["Dataset"]]


# combine ---------
all.df <- sample.meta %>% left_join(var.df, by = "Dataset") %>% 
  left_join(var.df2, by = "Dataset") %>% 
  tibble::column_to_rownames(var = "Dataset")
colnames(all.df)
all.df.n <- all.df %>% dplyr::select(where(is.numeric)) %>% 
  dplyr::select(which(!colSums(., na.rm=TRUE) %in% 0))
numeric.cols1 = intersect(setdiff(c(colnames(var.df), colnames(var.df2)), "Dataset"), colnames(all.df.n))
# numeric.cols = colnames(all.df.n)[colnames(all.df) %in% setdiff(colnames(var.df), "Dataset")]
numeric.cols1

if (Select.Group %in% c("Compare_Group1", "Compare_Group3", "Select_Group")){
  numeric.cols2 =  c( 
    "Blood.Urea.Nitrogen..BUN.", "Creatinine", "Sodium", "Potassium", 
    "Calcium", "Protein..Total", "Albumin", "Alkaline.Phosphatase", "AST", "ALT",  
    "Red.Blood.Cell...RBC", "Hemoglobin", "Hematocrit", "MCV", "MCH", 
    "Platelets", "White.Blood.Cell...WBC", "HDL", "LDL", "Weight_pounds"
    ,
    "Age","A1c_Value", "Glucose", "Height_inches", "BMI", 
    "Blood.Pressure.Systolic", "Blood.Pressure.Diastolic", 
    "Insulin", "FFA",# "HOMA_IR","Weight_kg","Triglycerides", 
    
    "Cholesterol..Total",   "Bilirubin..Total",
    "TNF.alpha", "IL.6", "IL.8.CXCL8", "VEGF", "IL.7", "CCL3.MIP.1alpha", 
    "CCL4.MIP.1", "RAGE.AGER", "SOST.Sclerostin", "ADAMTS13", "Osteropontin.OPN", "ICAM.1.CD54", 
    "IL.15", "GDF.15", "TNF.RI.TNFRSF1A", "Fas.TNFRSF6.CD95", "GDNF", "ActivinA" 
  )
}

if (Select.Group %in% c("Compare_Group2")){
  numeric.cols2 =  c( 
    "Blood.Urea.Nitrogen..BUN.", "Creatinine", "Sodium", "Potassium", 
    "Calcium", "Protein..Total", "Albumin", "Alkaline.Phosphatase", "AST", "ALT", 
    "Chloride", "CO2", "White.Blood.Cell...WBC", "Red.Blood.Cell...RBC", "Hemoglobin", "Hematocrit",
    "MCV", "MCH", "Platelets", 
    "Monocyte", "Neutrophils", "Eosinophils", 
    "MCHC", "RDW","Triglycerides", "Bilirubin..Total",
    "Cholesterol..Total",  "HDL", "LDL", "VLDL", 
    "A1c_Value", "Glucose", 
    "Tissue...Fat.", "Region...Fat.", "Fat..g.","Lean..g.", "Fat_Lean_Ratio", 
    "Tissue..g.",  "BMC..g.", "Fat.Free..g.", "T.Mass..kg.", "Weight_kg",
    "BMI", "Blood.Pressure.Systolic", "Blood.Pressure.Diastolic",
    "Heart.Rate..bpm.", "Respirations.per.minute",
    "HOMA_IR",  "Insulin", "FFA", "TNF.alpha", "IL.6", "IL.8.CXCL8", "VEGF", "IL.7", "CCL3.MIP.1alpha", 
    "CCL4.MIP.1", "RAGE.AGER", "SOST.Sclerostin", "ADAMTS13", "Osteropontin.OPN", "ICAM.1.CD54", 
    "IL.15", "GDF.15", "TNF.RI.TNFRSF1A", "Fas.TNFRSF6.CD95", "GDNF", "ActivinA" )
}

numeric.cols = c(numeric.cols1, numeric.cols2)
numeric.variables =c(numeric.cols, charac.variables)

mse.df <- data.frame()
Import.list <- list()
Test.list <- list()
per.var.df <- data.frame()

for (uu in  numeric.cols1  ){
  # uu = numeric.cols1[1]
  data <- all.df[!is.na(all.df[[uu]]), ] %>% 
    dplyr::select(one_of(c(uu, numeric.cols2, charac.variables))) 
  data <-  data[,colSums(is.na(data)) == 0] 
  
  if (length(unique(data[[uu]])) > 5 ){
    n <- nrow(data)
    set.seed (13)
    ntest <- trunc(n / 2)
    trainid <- sample (1:n, ntest)
    testid <- (1:n)[-trainid]
    if (length(unique(data[trainid,][[uu]])) > 5 ){
      library("randomForest")
      library(caret)
      
      print(uu)
      # bagging -------------------------
      #' bagging is simply a special case of a random forest with m = p
      rf.nhanes <- randomForest(formula(paste0(uu, "~.")), data = data, subset = trainid, 
                                mtry = (ncol(data) -1)  ,
                                importance = TRUE)
      rf.nhanes
      per.var.df[uu, "bagging.per.var"] <- round(100 * rf.nhanes$rsq[length(rf.nhanes$rsq)], digits = 2)
      
      # varImpPlot(rf.nhanes) 
      y.pred <- predict(rf.nhanes, newdata = data[testid,])
      plot(x = y.pred, y= data[testid, uu]) 
      abline (0 , 1)
      mse.df[uu, "bagging.MSE"] <- mean (( y.pred - data[testid, uu] ) ^2)
      cols.names <- colnames(as.data.frame(importance ( rf.nhanes )))
      Import.list[[uu]] <-          as.data.frame(importance ( rf.nhanes )) %>% 
        `colnames<-`(paste0(cols.names, ".bagging")) %>%
        tibble::rownames_to_column(var = "rowname")
      Test.list[[uu]][["bagging"]] <- rf.nhanes
      
      # confusionMatrix( as.factor(predict(rf.nhanes, data[testid,])), as.factor(data[testid, "composition"]))
      
      # randomForest ------------------------
      # Summary: 
      # The random forest was built to predict HbA1c.
      # The default number of trees is 1000
      # 4 variables were considered at each internal node by default when building the trees
      
      # Step 3.2: Plot the training MSE by number of trees to select optional number of trees.**
      # 
      rf.nhanes <- randomForest(formula(paste0(uu, "~.")), data = data, subset = trainid,
                                mtry = (ncol(data) -1) ,
                                ntree = 1000)
      rf.nhanes
      y.pred <- predict(rf.nhanes, newdata = data[testid,])
      mean (( y.pred - data[testid, uu] ) ^2)
      plot(x = y.pred, y= data[testid, uu]) 
      abline (0 , 1)
      per.var.df[uu, "randomForest.per.var"] <- round(100 * rf.nhanes$rsq[length(rf.nhanes$rsq)], digits = 2)
      mse.df[uu, "randomForest.MSE"] <- mean (( y.pred - data[testid, uu] ) ^2)
      cols.names <- colnames(as.data.frame(importance ( rf.nhanes )))
      Import.list[[uu]] <-  full_join(Import.list[[uu]], 
                                      as.data.frame(importance ( rf.nhanes )) %>% 
                                        `colnames<-`(paste0(cols.names, ".randomForest"))%>%
                                        tibble::rownames_to_column(var = "rowname"), by = "rowname")
      # rf.mse <- data.frame(Trees = rep(1:1000), mse = c(rf.nhanes$mse))
      # ggplot(data = rf.mse, aes(x=Trees, y=mse)) + geom_point() + geom_line() + 
      #   ggtitle("Changes in training MSE by number of trees, p = 4") + theme(plot.title = element_text(hjust = 0.5))
      Test.list[[uu]][["randomForest"]] <- rf.nhanes
      
      rf.nhanes <- randomForest(formula(paste0(uu, "~.")), data = data, subset = trainid,
                                mtry = (ncol(data) -1) /2 ,
                                importance = TRUE)
      rf.nhanes
      y.pred <- predict(rf.nhanes, newdata = data[testid,])
      mean (( y.pred - data[testid, uu] ) ^2)
      plot(x = y.pred, y= data[testid, uu]) 
      abline (0 , 1)
      per.var.df[uu, "randomForest2.per.var"] <- round(100 * rf.nhanes$rsq[length(rf.nhanes$rsq)], digits = 2)
      mse.df[uu, "randomForest2.MSE"] <- mean (( y.pred - data[testid, uu] ) ^2)
      # varImpPlot ( rf.nhanes ) 
      cols.names <- colnames(as.data.frame(importance ( rf.nhanes )))
      Import.list[[uu]] <-   full_join(Import.list[[uu]], 
                                       as.data.frame(importance ( rf.nhanes )) %>% 
                                         `colnames<-`(paste0(cols.names, ".randomForest2"))%>%
                                         tibble::rownames_to_column(var = "rowname"), by = "rowname")
      #' The first is based upon
      #' the mean decrease of accuracy in predictions on the out of bag samples when
      #' a given variable is permuted.
      #' 
      #' The second is a measure of the total decrease
      #' in node impurity that results from splits over that variable, averaged over
      #' all trees
      Test.list[[uu]][["randomForest2"]] <- rf.nhanes
      
      # boosting --------------------------------------
      library (gbm)
      set.seed (1)
      # boost.data <- gbm ( composition ~ . , data = data[trainid,],
      #                     distribution = "gaussian" , n.trees = 5000 ,
      #                     interaction.depth = 4)
      # summary(boost.data)
      # # plot ( boost.data , i = "GDF.15" )
      # # plot ( boost.data , i = "Blood.Pressure.Diastolic" )
      # yhat.boost <- predict ( boost.data ,
      #                         newdata = data[testid,] , n.trees = 5000)
      # mean (( yhat.boost - data[testid, "composition"] ) ^2)
      # mse.df[uu, "boosting.MSE"] <- mean (( yhat.boost - data[testid, "composition"] ) ^2)
      # cols.names <- colnames(as.data.frame(summary ( boost.data )))
      # Import.list[[uu]] <- full_join(Import.list[[uu]],
      #                                as.data.frame(summary (boost.data )) %>%
      #                                  `colnames<-`(paste0(cols.names, ".boosting"))%>%
      #                                  tibble::rownames_to_column(var = "rowname"), by =  "rowname")
      # Test.list[[uu]][["boosting"]] <-boost.data
      # 
      # 
      # boost.data <- gbm ( composition ~ . , data = data[trainid,],
      #                     distribution = "gaussian" , n.trees = 5000 ,
      #                     interaction.depth = 4, shrinkage = 0.2 , verbose =F)
      # summary(boost.data)
      # # plot ( boost.data , i = "GDF.15" )
      # # plot ( boost.data , i = "Blood.Pressure.Diastolic" )
      # yhat.boost <- predict ( boost.data ,
      #                         newdata = data[testid,] , n.trees = 5000)
      # mean (( yhat.boost - data[testid, "composition"] ) ^2)
      # mse.df[uu, "boosting2.MSE"] <- mean (( yhat.boost - data[testid, "composition"] ) ^2)
      # cols.names <- colnames(as.data.frame(summary ( boost.data )))
      # Import.list[[uu]] <- full_join(Import.list[[uu]], 
      #                                as.data.frame(summary (boost.data )) %>% 
      #                                  `colnames<-`(paste0(cols.names, ".boosting2"))%>%
      #                                  tibble::rownames_to_column(var = "rowname"), by =  "rowname")
      # Test.list[[uu]][["boosting2"]] <-boost.data
      
      # Bayesian Additive Regression Trees --------------------------------------------
      library (BART)
      x <- data [ , -which(names(data) %in% c(uu, charac.variables))]
      y <- data [ , uu ]
      xtrain <- x [ trainid, ]
      ytrain <- y [ trainid ]
      xtest <- x [ testid , ]
      ytest <- y [ testid ]
      set.seed (1)
      bartfit <- gbart ( xtrain , ytrain , x.test = xtest )
      yhat.bart <- bartfit$yhat.test.mean
      mean (( ytest - yhat.bart ) ^2)
      ord <- order ( bartfit $varcount.mean , decreasing = T )
      bartfit$varcount.mean[ord]
      mse.df[uu, "BART.MSE"] <- mean (( ytest - yhat.bart ) ^2)
      cols.names <- colnames(as.data.frame(bartfit$varcount.mean[ord]))
      Import.list[[uu]] <- full_join(Import.list[[uu]], 
                                     as.data.frame(bartfit$varcount.mean[ord]) %>% 
                                       `colnames<-`(paste0(cols.names, ".BART"))%>%
                                       tibble::rownames_to_column(var = "rowname"), by =  "rowname")
      Test.list[[uu]][["BART"]] <-bartfit
      print(uu)
      
    }
  }
} 
Import.result <- bind_rows(Import.list, .id = "Basic.variables") 
save(mse.df,  Import.list , Test.list, all.df, trainid,   testid , per.var.df, Import.result, 
     file = paste0(dir0,  Select.Group , "_Basic.variables", "_", "Clinical.ML.RData")
)



