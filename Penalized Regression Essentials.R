# http://www.sthda.com/english/articles/37-model-selection-essentials-in-r/153-penalized-regression-essentials-ridge-lasso-elastic-net/#computing-ridge-regression
library(tidyverse)
library(caret)
library(glmnet)
sample.meta <- sample.meta %>% 
  select_if(~sum(!is.na(.)) > 0) %>% 
  select(where(~n_distinct(.) > 1)) %>% #  remove columns with same value in R
  select(-Dataset )# -sample_id_ori, -sample_id, -molecule_h5, -Metrics, -Samples, -I7_Index.ID, -I7_Index_Sequence, -FastqFolder 
sample.meta <- sample.meta[ , colSums(is.na(sample.meta)) == 0] # Remove columns from dataframe where some of values are NA
sample.meta$Adipocytes_avg_log2FC.cl.c
colnames(sample.meta)
set.seed(123)
training.samples <- sample.meta$Adipocytes_avg_log2FC.cl.c %>%
  createDataPartition(p = 0.8, list = FALSE)
train.data  <- sample.meta[training.samples, ]
test.data <- sample.meta[-training.samples, ]


# Setup a grid range of lambda values:
  lambda <- 10^seq(-3, 3, length = 100)

# Compute ridge regression:
set.seed(123)
ridge <- train(
  Adipocytes_avg_log2FC.cl.c ~., data = train.data, method = "glmnet",
  trControl = trainControl("cv", number = 10),
  tuneGrid = expand.grid(alpha = 0, lambda = lambda)
)
# Model coefficients
coef(ridge$finalModel, ridge$bestTune$lambda)
# Make predictions
predictions <- ridge %>% predict(test.data)
# Model prediction performance
data.frame(
  RMSE = RMSE(predictions, test.data$Adipocytes_avg_log2FC.cl.c),
  Rsquare = R2(predictions, test.data$Adipocytes_avg_log2FC.cl.c)
)

# Compute lasso regression:
  
  # Build the model
  set.seed(123)
lasso <- train(
  Adipocytes_avg_log2FC.cl.c ~., data = train.data, method = "glmnet",
  trControl = trainControl("cv", number = 10),
  tuneGrid = expand.grid(alpha = 1, lambda = lambda)
)
# Model coefficients
coef(lasso$finalModel, lasso$bestTune$lambda)
# Make predictions
predictions <- lasso %>% predict(test.data)
# Model prediction performance
data.frame(
  RMSE = RMSE(predictions, test.data$Adipocytes_avg_log2FC.cl.c),
  Rsquare = R2(predictions, test.data$Adipocytes_avg_log2FC.cl.c)
)


# Elastic net regression:
  
  # Build the model
  set.seed(123)
elastic <- train(
  Adipocytes_avg_log2FC.cl.c ~., data = train.data, method = "glmnet",
  trControl = trainControl("cv", number = 10),
  tuneLength = 10
)
# Model coefficients
coef(elastic$finalModel, elastic$bestTune$lambda)
# Make predictions
predictions <- elastic %>% predict(test.data)
# Model prediction performance
data.frame(
  RMSE = RMSE(predictions, test.data$Adipocytes_avg_log2FC.cl.c),
  Rsquare = R2(predictions, test.data$Adipocytes_avg_log2FC.cl.c)
)

# Comparing models performance:
models <- list(ridge = ridge, lasso = lasso, elastic = elastic)
resamples(models) %>% summary( metric = "RMSE")


res <- as.matrix(coef(elastic$finalModel, elastic$bestTune$lambda))











# Computing principal component regression
library(tidyverse)
library(caret)
library(pls)
# Build the model on training set
set.seed(123)
model <- train(
  Adipocytes_avg_log2FC.cl.c~., data = train.data, method = "pcr",
  scale = TRUE,
  trControl = trainControl("cv", number = 10),
  tuneLength = 10
)
# Plot model RMSE vs different values of components
plot(model)
# Print the best tuning parameter ncomp that
# minimize the cross-validation error, RMSE
model$bestTune
# Make predictions
predictions <- model %>% predict(test.data)
# Model performance metrics
data.frame(
  RMSE = caret::RMSE(predictions, test.data$Adipocytes_avg_log2FC.cl.c),
  Rsquare = caret::R2(predictions, test.data$Adipocytes_avg_log2FC.cl.c)
)

# Computing partial least squares
# Build the model on training set
set.seed(123)
model <- train(
  Adipocytes_avg_log2FC.cl.c~., data = train.data, method = "pls",
  scale = TRUE,
  trControl = trainControl("cv", number = 10),
  tuneLength = 10
)
# Plot model RMSE vs different values of components
plot(model)
# Print the best tuning parameter ncomp that
# minimize the cross-validation error, RMSE
model$bestTune


