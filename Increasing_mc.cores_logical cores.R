parallel::detectCores(logical = FALSE);parallel::detectCores()
message(paste("Number of cores available:", parallelly::availableCores()))

# https://stackoverflow.com/questions/61115433/increasing-mc-cores-beyond-the-number-of-logical-cores
sleepy <- function(i) {
  start <- Sys.time()
  Sys.sleep(i)
  as.numeric(Sys.time() - start)
}


mc.cores <- 100L
ntasks   <- 10000L

start <- Sys.time()
out <- parallel::mclapply(2/ntasks*runif(ntasks), sleepy, mc.cores = mc.cores)

real_duration <- as.numeric(Sys.time() - start)
cpu_duration <- sum(unlist(out))

data.frame(logical.cores = parallel::detectCores(),
           mc.cores      = mc.cores,
           speedup       = cpu_duration/real_duration)

Result <- list()
for (nn in 1: parallel::detectCores()){
  mc.cores <- nn
 
  start <- Sys.time()
  out <- parallel::mclapply(2/ntasks*runif(ntasks), sleepy, mc.cores = mc.cores)
  
  real_duration <- as.numeric(Sys.time() - start)
  cpu_duration <- sum(unlist(out))
  
  Result[[nn]] <- data.frame(logical.cores = parallel::detectCores(),
             mc.cores      = mc.cores,
             speedup       = cpu_duration/real_duration)
}
Result.df <- bind_rows(Result)
plot(Result.df$mc.cores, Result.df$speedup)
