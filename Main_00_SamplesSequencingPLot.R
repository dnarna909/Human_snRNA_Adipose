



##' 2.7 Histograms of quality measures for each dataset ----------------------------------------------------------------------------------
hist.Mito_percent <- lapply(sce.list, function(xx) { 
  ggplot(as.data.frame(xx@colData), aes(x=subsets_Mito_percent))+geom_histogram(binwidth=0.1)
} )
gridExtra::grid.arrange(grobs = hist.Mito_percent, cols = 4)


hist.Mito_percent <- lapply(seq_along(sce.list), function(i) { 
  paste(i)
  ggplot(as.data.frame(sce.list[[as.numeric(i)]]@colData), aes(x=subsets_Mito_percent))+geom_histogram(binwidth=1) + 
    geom_vline(xintercept = 15, linetype="solid", color = "red", size=1)  + 
    labs(title = names(sce.list)[as.numeric(i)])} )
gridExtra::grid.arrange(grobs = hist.Mito_percent, cols = 4)


hist.log.total <- lapply(seq_along(sce.list), function(i) { 
  ggplot(as.data.frame(sce.list[[as.numeric(i)]]@colData), aes(x=log10(total)))+geom_histogram(binwidth=0.5) + 
    xlim(1, 6) + 
    labs(title = names(sce.list)[as.numeric(i)])} )
gridExtra::grid.arrange(grobs = hist.log.total, cols = 4)

hist.total <- lapply(seq_along(sce.list), function(i) { 
  ggplot(as.data.frame(sce.list[[as.numeric(i)]]@colData), aes(x=total))+geom_histogram(binwidth=10000) + 
    xlim(0, 1e+05) + geom_vline(xintercept = 1000, linetype="solid", color = "red", size=1.5) + 
    labs(title = names(sce.list)[as.numeric(i)])  } )
gridExtra::grid.arrange(grobs = hist.total, cols = 4)


hist.log.detected <- lapply(seq_along(sce.list), function(i) { 
  ggplot(as.data.frame(sce.list[[as.numeric(i)]]@colData), aes(x=log10(detected)))+geom_histogram(binwidth=0.5) + 
    xlim(1, 6) + 
    labs(title = names(sce.list)[as.numeric(i)])} )
gridExtra::grid.arrange(grobs = hist.log.detected, cols = 4)


hist.detected <- lapply(seq_along(sce.list), function(i) { 
  ggplot(as.data.frame(sce.list[[as.numeric(i)]]@colData), aes(x=(detected)))+geom_histogram(binwidth=200) +
    xlim(0, 10000) + geom_vline(xintercept = 500, linetype="solid", color = "red", size=1.5) + 
    labs(title = names(sce.list)[as.numeric(i)])} )
gridExtra::grid.arrange(grobs = hist.detected, cols = 4)


hist.total.detected <- lapply(seq_along(sce.list), function(i) { 
  ggplot(as.data.frame(sce.list[[as.numeric(i)]]@colData), aes(x=(total/detected)))+geom_histogram(binwidth=0.5) + 
    xlim(0.5, 10) + 
    labs(title = names(sce.list)[as.numeric(i)])} )
gridExtra::grid.arrange(grobs = hist.total.detected, cols = 4)


hist.log.total.detected <- lapply(seq_along(sce.list), function(i) { 
  ggplot(as.data.frame(sce.list[[as.numeric(i)]]@colData), aes(x=log10(total)/log10(detected)))+geom_histogram(binwidth=0.02) + 
    xlim(0.9, 1.4) + geom_vline(xintercept = 1.25, linetype="solid", 
                                color = "red", size=1.5) + 
    labs(title = names(sce.list)[as.numeric(i)])} )
gridExtra::grid.arrange(grobs = hist.log.total.detected, cols = 4)


##' 2.8 save sce.list for plotting purpose--------------------------------------------------------------------------------------------------
save(hist.detected, hist.log.detected, hist.log.total, hist.log.total.detected,
     hist.Mito_percent, hist.total, hist.total.detected, file = paste0(Disk, Project.folder, "/", Rds.folder, "/QC hist_files.Rds"))

