
# rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.
memory.limit(size = 1e+13)

# import parameters and packages --------------------------------------------------------------------------------------------
# source(paste0("/media/jianie/Extreme SSD1/2022-09-01 STARR_SGLT2 Combine/Project Parameters.R"), local = knitr::knit_global())
source(paste0(paste0(Disk, Project.folder, "/", "Project Parameters.R")), local = knitr::knit_global())
library("RColorBrewer")
library("ggpubr")

# Import data -----------------------------------------------------------------------------------------------------------------------------------------------------
# export.folder <- "Compositions"
# dir.create(file.path(paste0(Disk, Project.folder, "/", Rds.folder, "/"), export.folder), showWarnings = FALSE)
dir0 <- paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder, "/")
files <- list.files(dir0)

dir.create(file.path(paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder, ".VariableSelection")), showWarnings = FALSE)
dir1 <- paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder, ".VariableSelection", "/")
 # cat(paste(shQuote(colnames(Composition_3)[sapply(Composition_3, is.numeric)], type = "cmd"), collapse = ", "))
  numeric.variables =  c( "Blood.Urea.Nitrogen..BUN.", "Creatinine", "Sodium", "Potassium", 
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
  for (list.name in names(Composition.lists)) {
    # list.name = names(Composition.lists)[2]
    label.name <- stringr::str_split(list.name, paste0("_"))[[1]][length(stringr::str_split(list.name, paste0("_"))[[1]])]
    x.name <- stringr::str_split(list.name, paste0("_"))[[1]][2]
    Composition_3 <- Composition.lists[[list.name]][["Composition_3All"]] %>% 
      dplyr::filter(Treatment_Group %in% c("Notreatment", "NUT_Pre", "SGLT2i_Pre"),
                    BA_Group %in% c("Middle_Lean","Middle_Overweight", "Older_Lean", "Older_Overweight")#
      ) %>% 
      mutate(BA_Group = factor(BA_Group, levels = c("Middle_Lean", "Middle_Overweight", "Older_Lean", "Older_Overweight")))# 
    Composition <- Composition.lists[[list.name]][["Composition_All"]] 
    
    clusters = colnames(Composition)
    # display.brewer.pal(length(clusters), "Set1")
    # Composition_3$x.label <- factor(Composition_3[[label.name]], levels=clusters)
    mse.df <- data.frame()
    Import.list <- list()
    Test.list <- list()
    per.var.df <- data.frame()
    
    if (x.name == label.name){
      for (uu in  unique(Composition_3[[x.name]])){
        data <- Composition_3 %>%
          dplyr::filter(!!sym(x.name) == uu ) %>%
          dplyr::select(one_of(c("composition",numeric.variables))) %>%
          dplyr::filter(complete.cases(.)) 
        
        if (length(unique(data$composition)) > 5 ){
          n <- nrow(data)
          set.seed (13)
          ntest <- trunc(n / 2)
          trainid <- sample (1:n, ntest)
          testid <- (1:n)[-trainid]
          if (length(unique(data[trainid,]$composition)) > 5 ){
            library("randomForest")
            library(caret)
            
            print(uu)
            # bagging -------------------------
            #' bagging is simply a special case of a random forest with m = p
            rf.nhanes <- randomForest(composition~., data = data, subset = trainid, 
                                      mtry = length(numeric.variables) ,
                                      importance = TRUE)
            rf.nhanes
            per.var.df[uu, "bagging.per.var"] <- round(100 * rf.nhanes$rsq[length(rf.nhanes$rsq)], digits = 2)
            
            # varImpPlot(rf.nhanes) 
            y.pred <- predict(rf.nhanes, newdata = data[testid,])
            plot(x = y.pred, y= data[testid, "composition"]) 
            abline (0 , 1)
            mse.df[uu, "bagging.MSE"] <- mean (( y.pred - data[testid, "composition"] ) ^2)
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
            rf.nhanes <- randomForest(composition~., data = data, subset = trainid,
                                      mtry = length(numeric.variables) ,
                                      ntree = 1000)
            rf.nhanes
            y.pred <- predict(rf.nhanes, newdata = data[testid,])
            mean (( y.pred - data[testid, "composition"] ) ^2)
            plot(x = y.pred, y= data[testid, "composition"]) 
            abline (0 , 1)
            per.var.df[uu, "randomForest.per.var"] <- round(100 * rf.nhanes$rsq[length(rf.nhanes$rsq)], digits = 2)
            mse.df[uu, "randomForest.MSE"] <- mean (( y.pred - data[testid, "composition"] ) ^2)
            cols.names <- colnames(as.data.frame(importance ( rf.nhanes )))
            Import.list[[uu]] <-  full_join(Import.list[[uu]], 
                                            as.data.frame(importance ( rf.nhanes )) %>% 
                                              `colnames<-`(paste0(cols.names, ".randomForest"))%>%
                                              tibble::rownames_to_column(var = "rowname"), by = "rowname")
            # rf.mse <- data.frame(Trees = rep(1:1000), mse = c(rf.nhanes$mse))
            # ggplot(data = rf.mse, aes(x=Trees, y=mse)) + geom_point() + geom_line() + 
            #   ggtitle("Changes in training MSE by number of trees, p = 4") + theme(plot.title = element_text(hjust = 0.5))
            Test.list[[uu]][["randomForest"]] <- rf.nhanes
            
            rf.nhanes <- randomForest(composition~., data = data, subset = trainid,
                                      mtry = length(numeric.variables)/2 ,
                                      importance = TRUE)
            rf.nhanes
            y.pred <- predict(rf.nhanes, newdata = data[testid,])
            mean (( y.pred - data[testid, "composition"] ) ^2)
            plot(x = y.pred, y= data[testid, "composition"]) 
            abline (0 , 1)
            per.var.df[uu, "randomForest2.per.var"] <- round(100 * rf.nhanes$rsq[length(rf.nhanes$rsq)], digits = 2)
            mse.df[uu, "randomForest2.MSE"] <- mean (( y.pred - data[testid, "composition"] ) ^2)
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
            x <- data [ , -which(names(data) %in% c("composition"))]
            y <- data [ , "composition" ]
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
      save(mse.df,  Import.list , Test.list, data, trainid,   testid , per.var.df,
           file = paste0(dir1, file.name, "_", x.name, "_", "Variable.RData")
      )
    }
    
    
    if (x.name != label.name){
      for (ll in unique(Composition_3[[label.name]]) ) {
      }
    }
    
  }
}


  rm(dir1, dir0,  files, file,
   file.name, Annotation.file.name, Subtype.file.name,
   celltype.name, list.name, label.name, x.name,
   Composition_3, Composition,
   clusters.lists, Composition.lists, clusters,
   bartfit, data, Import.list, mse.df, rf.nhanes,
   Test.list, x, xtrain, xtest, cols.names,
   n, ntest, numeric.variables, ord, testid, trainid,
   uu, y, y.pred, yhat.bart, ytest, ytrain,
   per.var.df
)
gc() #free up memrory and report the memory usage.













