
# rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.
memory.limit(size = 1e+13)

# import parameters and packages --------------------------------------------------------------------------------------------
# source(paste0("/media/jianie/Extreme SSD1/2022-09-01 STARR_SGLT2 Combine/Project Parameters.R"), local = knitr::knit_global())
source(paste0(paste0(Disk, Project.folder, "/", codes.folder, "/","Project Parameters.R")), local = knitr::knit_global())
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
if (exists("Select_Group.names")) {
  sample.meta <- sample.meta %>% 
    dplyr::filter(
      !!sym(Select_Group) %in% Select_Group.names  
    )  
} else {
  sample.meta <- sample.meta %>% 
    dplyr::filter((!!sym(Select_Group)) == T)
}
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

if (Select_Group %in% c("Compare_Group1", "Compare_Group3", "Select_Group")){
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

if (Select_Group %in% c("Compare_Group2")){
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

Rsquare.df <- data.frame()
mse.df <- data.frame()
Import.list <- list()
Test.list <- list()
per.var.df <- data.frame()
data.list <-list()
Test.id.list <-list()  
Train.id.list <-list()  

for (uu in  numeric.cols1  ){
  # uu = numeric.cols1[1]
  data <- all.df[!is.na(all.df[[uu]]), ] %>% 
    dplyr::select(one_of(c(uu, numeric.cols2, charac.variables))) 
  data <-  data[,colSums(is.na(data)) == 0] 
  table(unlist(is.na(data)))
  data.list[[uu]] <- data
  
  n <- nrow(data)
  name.data <- data[, sapply(data, class) %in% c('factor', 'character')]
  for (cc in 1:ncol(name.data)) {
    data <- data %>% data.frame() %>% 
      tibble::rownames_to_column(var = "rowname")%>% 
      group_by_at(vars(one_of(c(colnames(name.data)[cc]))))   %>%
      # nest() %>%
      dplyr::filter(n() > 1 ) %>%
      ungroup() %>%  
      dplyr::filter(n_distinct((!!sym(colnames(name.data)[cc]))) >1) %>% 
      as.data.frame() %>% 
      tibble::column_to_rownames(var = "rowname")
    table(data[[colnames(name.data)[cc]]])
    # table(data[, colnames(name.data)[cc]])
    if (is.factor(data[[colnames(name.data)[cc]]])){
      data[[colnames(name.data)[cc]]] <- factor(data[[colnames(name.data)[cc]]]  , 
                                                levels = levels(data[[colnames(name.data)[cc]]])[ levels(data[[colnames(name.data)[cc]]]) %in% unique(data[[colnames(name.data)[cc]]])])
      table(data[[colnames(name.data)[cc]]])
    }
    if (is.character(data[[colnames(name.data)[cc]]])){
      data[[colnames(name.data)[cc]]] <- factor(data[[colnames(name.data)[cc]]]  , 
                                                levels = unique(data[[colnames(name.data)[cc]]]))
    }
  }
  random.sample <- function(x) {
    success <- FALSE
    while (!success) {
      # do something
      n <- nrow(data)
      name.data <- data[, sapply(data, class) %in% c('character', 'factor')]
      ntest <- trunc(n*2 / 3)
      trainid <- sample (1:n, ntest)
      check.data <- x[trainid,][, sapply(x[trainid,], class) %in% c('character', 'factor')]
      check.result <- c()
      for (cc in 1:ncol(check.data)) {
        check.result <- c(check.result,  length(unique(check.data[, cc])) == length(unique(name.data[, cc])) )
        n_occur <- data.frame(table(check.data[, cc]))
        check.result <- c(check.result,  length(unique(n_occur[n_occur$Freq > 1,"Var1"])) == length(unique(name.data[, cc])) )
      }; check.result
      # check for success
      success <- !(FALSE %in% check.result)
    }
    return(trainid)
  }
  trainid <- random.sample(data);trainid
  check.data <- data[trainid,][, sapply(data[trainid,], class) %in% c('character', 'factor')]
  check.result <- c()
  for (cc in 1:ncol(check.data)) {
    check.result <- c(check.result,  length(unique(check.data[, cc])) == length(unique(name.data[, cc])) )
    n_occur <- data.frame(table(check.data[, cc]))
    check.result <- c(check.result,  length(unique(n_occur[n_occur$Freq > 1,"Var1"])) == length(unique(name.data[, cc])) )
  }; check.result
  testid <- (1:n)[-trainid]
  Test.id.list[[uu]] <- testid 
  Train.id.list[[uu]] <-trainid 
  
  
  if (length(unique(data[[uu]])) > 5 ){
    if (length(unique(data[trainid,][[uu]])) > 5 ){
      print(uu)
      
      # 1. bagging -------------------------
      #' bagging is simply a special case of a random forest with m = p
      library("randomForest")
      library(caret)
      rf.nhanes <- NULL
      rf.nhanes <- randomForest(formula(paste0(uu, "~.")), data = data, subset = trainid, 
                                mtry = (ncol(data) -1)  ,
                                importance = TRUE)
      rf.nhanes
      if(!is.null(rf.nhanes)){
        per.var.df[uu, "bagging.per.var"] <- round(100 * rf.nhanes$rsq[length(rf.nhanes$rsq)], digits = 2)
        
        # varImpPlot(rf.nhanes) 
        y.pred <- predict(rf.nhanes, newdata = data[testid,])
        plot(x = y.pred, y= data[testid, uu]) 
        abline (0 , 1)
        mse.df[uu, "bagging.MSE"] <- mean(( y.pred - data[testid, uu] ) ^2)
        Rsquare.df[uu, "bagging.Rsquare"] <- caret::R2(y.pred, data[testid, uu])
        cols.names <- colnames(as.data.frame(importance ( rf.nhanes )))
        Import.list[[uu]] <-  as.data.frame(importance ( rf.nhanes )) %>% 
          `colnames<-`(paste0(cols.names, ".bagging")) %>%
          tibble::rownames_to_column(var = "rowname")
        Test.list[[uu]][["bagging"]] <- rf.nhanes
        
        # confusionMatrix( as.factor(predict(rf.nhanes, data[testid,])), as.factor(data[testid, "composition"]))
        
      }
      
      # 2. randomForest ------------------------
      # Summary: 
      # The random forest was built to predict HbA1c.
      # The default number of trees is 1000
      # 4 variables were considered at each internal node by default when building the trees
      
      # Step 3.2: Plot the training MSE by number of trees to select optional number of trees.**
      # 
      rf.nhanes <- NULL
      rf.nhanes <- randomForest(formula(paste0(uu, "~.")), data = data, subset = trainid,
                                mtry = (ncol(data) -1) ,
                                ntree = 1000)
      rf.nhanes
      if(!is.null(rf.nhanes)){
        y.pred <- predict(rf.nhanes, newdata = data[testid,])
        mean (( y.pred - data[testid, uu] ) ^2)
        plot(x = y.pred, y= data[testid, uu]) 
        abline (0 , 1)
        per.var.df[uu, "randomForest.per.var"] <- round(100 * rf.nhanes$rsq[length(rf.nhanes$rsq)], digits = 2)
        mse.df[uu, "randomForest.MSE"] <- mean (( y.pred - data[testid, uu] ) ^2)
        Rsquare.df[uu, "randomForest.Rsquare"] <- caret::R2(y.pred, data[testid, uu])
        cols.names <- colnames(as.data.frame(importance ( rf.nhanes )))
        Import.list[[uu]] <-  full_join(Import.list[[uu]], 
                                        as.data.frame(importance ( rf.nhanes )) %>% 
                                          `colnames<-`(paste0(cols.names, ".randomForest"))%>%
                                          tibble::rownames_to_column(var = "rowname"), by = "rowname")
        # rf.mse <- data.frame(Trees = rep(1:1000), mse = c(rf.nhanes$mse))
        # ggplot(data = rf.mse, aes(x=Trees, y=mse)) + geom_point() + geom_line() + 
        #   ggtitle("Changes in training MSE by number of trees, p = 4") + theme(plot.title = element_text(hjust = 0.5))
        Test.list[[uu]][["randomForest"]] <- rf.nhanes
      }
      
      rf.nhanes <- NULL
      rf.nhanes <- randomForest(formula(paste0(uu, "~.")), data = data, subset = trainid,
                                mtry = (ncol(data) -1) /2 ,
                                importance = TRUE)
      rf.nhanes
      if(!is.null(rf.nhanes)){
        y.pred <- predict(rf.nhanes, newdata = data[testid,])
        mean (( y.pred - data[testid, uu] ) ^2)
        plot(x = y.pred, y= data[testid, uu]) 
        abline (0 , 1)
        per.var.df[uu, "randomForest2.per.var"] <- round(100 * rf.nhanes$rsq[length(rf.nhanes$rsq)], digits = 2)
        Rsquare.df[uu, "randomForest2.Rsquare"] <- caret::R2(y.pred, data[testid, uu])
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
      }
      
      # 3. boosting --------------------------------------
      library (gbm)
      set.seed (1)
      # boost.data <- gbm ( formula(paste0(uu, "~.")) , data = data[trainid,],
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
      
      # 4. Bayesian Additive Regression Trees --------------------------------------------
      library (BART)
      x <- data [ , -which(names(data) %in% c(uu, charac.variables))]
      y <- data [ , uu ]
      xtrain <- x [ trainid, ]
      ytrain <- y [ trainid ]
      xtest <- x [ testid , ]
      ytest <- y [ testid ]
      set.seed (1)
      bartfit <- NULL
      bartfit <- gbart ( xtrain , ytrain , x.test = xtest )
      if(!is.null(bartfit)){
        yhat.bart <- bartfit$yhat.test.mean
        mean (( ytest - yhat.bart ) ^2)
        ord <- order( bartfit$varcount.mean , decreasing = T )
        bartfit$varcount.mean[ord]
        mse.df[uu, "BART.MSE"] <- mean (( ytest - yhat.bart ) ^2)
        Rsquare.df[uu, "BART.Rsquare"] <- caret::R2(yhat.bart, ytest)
        cols.names <- colnames(as.data.frame(bartfit$varcount.mean[ord]))
        Import.list[[uu]] <- full_join(Import.list[[uu]], 
                                       as.data.frame(bartfit$varcount.mean[ord]) %>% 
                                         `colnames<-`(paste0(cols.names, ".BART"))%>%
                                         tibble::rownames_to_column(var = "rowname"), by =  "rowname")
        Test.list[[uu]][["BART"]] <-bartfit
      }
      
      # 5. Compute ridge regression: -----------
      # Setup a grid range of lambda values:
      library(tidyverse)
      library(caret)
      library(glmnet)
      lambda <- 10^seq(-3, 3, length = 100)
      
      # Compute ridge regression:
      set.seed(123)
      ridge <- NULL
      ridge <- train(
        formula(paste0(uu, "~.")), data = data[trainid, ], method = "glmnet",
        trControl = trainControl("cv", number = 10),
        tuneGrid = expand.grid(alpha = 0, lambda = lambda)
      )
      if (!is.null(ridge)){
        # Model coefficients
        coef(ridge$finalModel, ridge$bestTune$lambda)
        # Make predictions
        predictions <- ridge %>% predict(data[testid, ])
        # Model prediction performance
        data.frame(
          RMSE = caret::RMSE(predictions, data[testid, uu]),
          Rsquare = caret::R2(predictions, data[testid, uu])
        )
        mse.df[uu, "ridge.MSE"] <- caret::RMSE(predictions, data[testid, uu])
        Rsquare.df[uu, "ridge.Rsquare"] <- caret::R2(predictions, data[testid, uu])
        cols.names <- "coef"
        # Import.list[[uu]] <-  as.data.frame( coef(ridge$finalModel, ridge$bestTune$lambda)) %>% 
        #   `colnames<-`(paste0(cols.names, ".ridge")) %>%
        #   tibble::rownames_to_column(var = "rowname")
        Import.list[[uu]] <-  full_join(Import.list[[uu]], 
                                        as.data.frame(coef(ridge$finalModel, ridge$bestTune$lambda)) %>% 
                                          `colnames<-`(paste0(cols.names, ".ridge"))%>%
                                          tibble::rownames_to_column(var = "rowname"), by = "rowname")
        Test.list[[uu]][["ridge"]] <- ridge
      }
      
      # 6. Compute lasso regression: --------
      # Build the model
      set.seed(123)
      lasso <- NULL
      lasso <- train(
        formula(paste0(uu, "~.")), data = data[trainid, ], method = "glmnet",
        trControl = trainControl("cv", number = 10),
        tuneGrid = expand.grid(alpha = 1, lambda = lambda)
      )
      if (!is.null(lasso)){
        # Model coefficients
        coef(lasso$finalModel, lasso$bestTune$lambda)
        # Make predictions
        predictions <- lasso %>% predict(data[testid, ])
        # Model prediction performance
        data.frame(
          RMSE = caret::RMSE(predictions, data[testid, uu]),
          Rsquare = caret::R2(predictions, data[testid, uu])
        )
        mse.df[uu, "lasso.MSE"] <- caret::RMSE(predictions, data[testid, uu])
        Rsquare.df[uu, "lasso.Rsquare"] <- caret::R2(predictions, data[testid, uu])
        cols.names <- "coef"
        # Import.list[[uu]] <-  as.data.frame( coef(lasso$finalModel, lasso$bestTune$lambda)) %>% 
        #   `colnames<-`(paste0(cols.names, ".lasso")) %>%
        #   tibble::rownames_to_column(var = "rowname")
        Import.list[[uu]] <-  full_join(Import.list[[uu]], 
                                        as.data.frame(coef(lasso$finalModel, lasso$bestTune$lambda)) %>% 
                                          `colnames<-`(paste0(cols.names, ".lasso"))%>%
                                          tibble::rownames_to_column(var = "rowname"), by = "rowname")
        Test.list[[uu]][["lasso"]] <- lasso
      }
      
      # 7. Elastic net regression: -----------
      # Build the model
      set.seed(123)
      elastic <- NULL
      elastic <- train(
        formula(paste0(uu, "~.")), data = data[trainid, ],  method = "glmnet",
        trControl = trainControl("cv", number = 10),
        tuneLength = 10
      )
      
      if (!is.null(elastic)){
        # Model coefficients
        coef(elastic$finalModel, elastic$bestTune$lambda)
        # Make predictions
        predictions <- elastic %>% predict(data[testid, ])
        # Model prediction performance
        data.frame(
          RMSE = caret::RMSE(predictions, data[testid, uu]),
          Rsquare = caret::R2(predictions, data[testid, uu])
        )
        mse.df[uu, "elastic.MSE"] <- caret::RMSE(predictions, data[testid, uu])
        Rsquare.df[uu, "elastic.Rsquare"] <- caret::R2(predictions, data[testid, uu])
        cols.names <- "coef"
        # Import.list[[uu]] <-  as.data.frame( coef(lasso$finalModel, lasso$bestTune$lambda)) %>% 
        #   `colnames<-`(paste0(cols.names, ".lasso")) %>%
        #   tibble::rownames_to_column(var = "rowname")
        Import.list[[uu]] <-  full_join(Import.list[[uu]], 
                                        as.data.frame(coef(elastic$finalModel, elastic$bestTune$lambda)) %>% 
                                          `colnames<-`(paste0(cols.names, ".elastic"))%>%
                                          tibble::rownames_to_column(var = "rowname"), by = "rowname")
        Test.list[[uu]][["elastic"]] <- elastic
        # resamples(Test.list[[uu]][c("ridge",   "lasso",   "elastic")]) %>% summary( metric = "RMSE")
      }
      rm(elastic)
      
      # 8. Computing principal component regression ------------
      library(tidyverse)
      library(caret)
      library(pls)
      # Build the model on training set
      set.seed(123)
      
      model <- NULL
      check.result <- c()
      for (cc in 1:ncol(check.data)) {
        check.result <- c(check.result,  length(unique(check.data[, cc])) == length(unique(name.data[, cc])) )
        n_occur <- data.frame(table(check.data[, cc]))
        check.result <- c(check.result,  length(unique(n_occur[n_occur$Freq > 2,"Var1"])) == length(unique(name.data[, cc])) )
      }; check.result
      if(!(FALSE %in% check.result)){
        model <- train(
          formula(paste0(uu, "~.")), data = data[trainid, ],  method = "pcr",
          scale = TRUE,preProc = c("center", "scale"),
          trControl = trainControl("cv", number = 10),
          # trControl = trainControl(method = "repeatedcv", repeats = 5),
          tuneLength = 5
        )
        if (!is.null(model)){
          coef(model$finalModel, model$bestTune$ncomp)
          # Plot model RMSE vs different values of components
          plot(model)
          # Print the best tuning parameter ncomp that
          # minimize the cross-validation error, RMSE
          model$bestTune
          # Make predictions
          predictions <- model %>% predict(data[testid, ])
          # Model performance metrics
          data.frame(
            RMSE = caret::RMSE(predictions, data[testid, uu]),
            Rsquare = caret::R2(predictions, data[testid, uu])
          )
          mse.df[uu, "pca.MSE"] <- caret::RMSE(predictions, data[testid, uu])
          Rsquare.df[uu, "pca.Rsquare"] <- caret::R2(predictions, data[testid, uu])
          cols.names <- "coef"
          # Import.list[[uu]] <-  as.data.frame( coef(lasso$finalModel, lasso$bestTune$lambda)) %>% 
          #   `colnames<-`(paste0(cols.names, ".lasso")) %>%
          #   tibble::rownames_to_column(var = "rowname")
          Import.list[[uu]] <-  full_join(Import.list[[uu]], 
                                          as.data.frame(coef(model$finalModel, model$bestTune$ncomp)) %>% 
                                            `colnames<-`(paste0(cols.names, ".pca"))%>%
                                            tibble::rownames_to_column(var = "rowname"), by = "rowname")
          Test.list[[uu]][["pca"]] <- model
        }
      }
      
      # 9. Computing partial least squares ---------------
      # Build the model on training set
      set.seed(123)
      model <- NULL
      if(!(FALSE %in% check.result)){
        model <- train(
          formula(paste0(uu, "~.")), data = data[trainid, ],  method = "pls",
          scale = TRUE,
          trControl = trainControl("cv", number = 10),
          tuneLength = 5
        )
        if (!is.null(model)){
          coef(model$finalModel, model$bestTune$ncomp)
          # Plot model RMSE vs different values of components
          plot(model)
          # Print the best tuning parameter ncomp that
          # minimize the cross-validation error, RMSE
          model$bestTune
          # Make predictions
          predictions <- model %>% predict(data[testid, ])
          # Model performance metrics
          data.frame(
            RMSE = caret::RMSE(predictions, data[testid, uu]),
            Rsquare = caret::R2(predictions, data[testid, uu])
          )
          mse.df[uu, "pls.MSE"] <- caret::RMSE(predictions, data[testid, uu])
          Rsquare.df[uu, "pls.Rsquare"] <- caret::R2(predictions, data[testid, uu])
          cols.names <- "coef"
          Import.list[[uu]] <-  full_join(Import.list[[uu]], 
                                          as.data.frame(coef(model$finalModel, model$bestTune$ncomp)) %>% 
                                            `colnames<-`(paste0(cols.names, ".pls"))%>%
                                            tibble::rownames_to_column(var = "rowname"), by = "rowname")
          Test.list[[uu]][["pls"]] <- model
        }
      }
      
      print(uu)
    }
  }
} 
Import.list <- lapply(Import.list, function(x){
  x <- x %>%
    as.data.frame() %>%
    mutate(
      tree.MSE.sum = rowSums(dplyr::select(., contains("IncMSE")&!contains("rowname")&!contains("coef") &!contains("sum") ), na.rm = TRUE),
      tree.Purity.sum = rowSums(dplyr::select(., contains("IncNodePurity")&!contains("rowname")&!contains("coef") &!contains("sum") ), na.rm = TRUE),
      tree.sum = rowSums(dplyr::select(., !contains("rowname")&!contains("coef") &!contains("sum") ), na.rm = TRUE),
      reg.sum = rowSums(dplyr::select(., !contains("rowname")&contains("coef")&!contains("sum") ), na.rm = TRUE)
      #, all.sum = rowSums(dplyr::select(., !contains("rowname") &!contains("sum") ), na.rm = TRUE)
    )
})

Import.result <- bind_rows(Import.list, .id = "Basic.variables") %>%
  dplyr::filter(!(rowname == "(Intercept)") )
if (exists("analysis.group")) {
  file.name = paste0(dir0,  analysis.group , "_Basic.variables", "_", "Clinical.ML.RData")
} else {
file.name = paste0(dir0,  Select_Group , "_Basic.variables", "_", "Clinical.ML.RData")
}
save(mse.df,  Import.list , Test.list, all.df, 
     per.var.df, 
     Rsquare.df,
     data.list,
     Test.id.list , 
     Train.id.list, Import.result ,
     file = file.name
)

rm(mse.df,  Import.list , Test.list, all.df, 
   per.var.df, 
   Rsquare.df,
   data.list,
   Test.id.list , 
   Train.id.list, Import.result)
rm(all.df.n, bartfit, data, lasso, model, n_occur, name.data,
   rf.nhanes, ridge, sample.meta, var.df, var.df2,
   x, xtest, xtrain, cc, check.result, cols.names, dir0, lambda,
   n, numeric.cols, numeric.cols.c, numeric.cols1, numeric.cols2,
   numeric.variables, ord, predictions, testid, trainid,
   uu, y, y.pred, yhat.bart, ytest, ytrain, random.sample,
   check.data, Annotation.file.name , file.name, Type.file.name )
gc()
sessionInfo()
