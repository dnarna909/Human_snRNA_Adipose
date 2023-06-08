library(parallel)
detectCores();# physical and logical cores
detectCores(logical = FALSE) # return the number of physical cores.

s <- system.time({
  df.summary <- df_long[, 1:1000] %>% 
    group_by(group) %>%
    summarise(across(where(is.numeric), .fns = 
                       list(
                         min =  ~ min(., na.rm = TRUE),
                         count = ~ n(),
                         median =  ~ median(., na.rm = TRUE),
                         mean =  ~ mean(., na.rm = TRUE),
                         stdev = sd(., na.rm = TRUE),  
                         q25 = ~quantile(., 0.25),
                         q75 = ~quantile(., 0.75),
                         max =  ~ max(., na.rm = TRUE)))) 
});s

r <- mcmapply(  function(df_long) {
  df_long %>% 
    group_by(group) %>%
    summarise(across(where(is.numeric), .fns = 
                       list(
                         min =  ~ min(., na.rm = TRUE),
                         count = ~ n(),
                         median =  ~ median(., na.rm = TRUE),
                         mean =  ~ mean(., na.rm = TRUE),
                         stdev = sd(., na.rm = TRUE),  
                         q25 = ~quantile(., 0.25),
                         q75 = ~quantile(., 0.75),
                         max =  ~ max(., na.rm = TRUE)))) 
}, mc.cores = 4)







s <- system.time({
  mn <- lapply(Result, function(df) {
    lapply(df, function(dx) {
      quantile(dx[,2], 0.9, na.rm = TRUE)
    }) 
  })
}); s;
mn

s <- system.time({
  mn <- mclapply(Result, function(df) {
    lapply(df, function(dx) {
      quantile(dx[,2], 0.9, na.rm = TRUE)
    }) 
  }, 
  mc.cores = 4)
}); s; 
mn



f <- function(i) {
  map(df_long[, 2:(ncol(df_long))], fit_aov)
}
system.time(save1 <- lapply(1:2, f))
# user  system elapsed 
# 327.250   0.096 327.298 
system.time(save1 <- mclapply(1:2, f, mc.cores = 4))
# user  system elapsed 
# 244.239   1.139 257.616 
system.time(save1 <- mclapply(1:2, f, mc.cores = 12))
# user  system elapsed 
# 246.198   1.199 246.527 
system.time(save1 <- mclapply(1:2, f, mc.cores = 20))
# user  system elapsed 
# 249.626   1.167 251.088 

