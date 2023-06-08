#!/usr/bin/env Rscript
# https://github.com/czbiohub/tabula-muris-senis/blob/master/2_aging_signature/job.DGE_analysis/DGE_analysis.R
### Load libraries -------------------------------------------------------------------------------------------------------------
#rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.
memory.limit(size = 1e+13)

# source(paste0("/media/jianie/Extreme SSD1/2022-09-01 STARR Combine/Project Parameters.R"), local = knitr::knit_global())
source(paste0(Disk, Project.folder, "/", "Project Parameters_GCCRI.R"), local = knitr::knit_global())


start_time =proc.time()

library(MAST)

# Load parameters 
# args = commandArgs(trailingOnly=TRUE)
output_folder = args[1]
data_name = args[2]
str_n_gene = args[3]
analyte = args[4]
keep_sex = args[5]
flag_M_O = args[6]
tissue1 = args[7]; # unlist(strsplit(analyte, "[.]"))[1]
tissue2 = args[8]; # unlist(strsplit(analyte, "[.]"))[1]
celltype =  unlist(strsplit(unlist(strsplit(data_name, "[.Rds]"))[1], "_"))[4]
if (is.na(celltype) ){
  celltype = "All"
}
if(length(args) == 9){ others = args[9]}

print(paste0('output_folder: ', output_folder))
print(paste0('data_name: ', data_name))
print(paste0('str_n_gene: ', str_n_gene))
print(paste0('analyte: ', analyte))
print(paste0('keep_sex: ', keep_sex))
print(paste0('flag_M_O: ', flag_M_O))
print(paste0('tissue: ', tissue1, ", ",tissue2))
print(paste0('celltype: ', celltype))
print(paste0('others variables: ', others))

dir.create(file.path(paste0(Disk, Project.folder, "/", Rds.folder, "/"), output_folder), showWarnings = FALSE)

# Load data 
adata_temp = readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/", 
                            data_name))
sample.meta <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/", sample.file.type)) 
colnames(sample.meta)
adata_temp <- AddMetaData(adata_temp, metadata = sample.meta)

# Filter data
print('Before filtering')
print(dim(adata_temp))

# if (is.na(celltype)==FALSE){
#   print('tc analysis')
#   ind_select = (adata_temp$cell_ontology_class==celltype)
#   adata_temp = adata_temp[,ind_select]
# } else {
#   print('tissue analysis')
# }

if (is.na(tissue1)==FALSE){
  if (is.na(tissue2)==FALSE){
    print('tc analysis')
    table( (adata_temp@meta.data[, tissue2] == TRUE),  (adata_temp@meta.data[, tissue1] == TRUE))
    ind_select  = (adata_temp@meta.data[, tissue1] == TRUE) | (adata_temp@meta.data[, tissue2] == TRUE)
    table(ind_select)
    adata_temp = adata_temp[,ind_select]
  } else{
    ind_select  = (adata_temp@meta.data[, tissue1] == TRUE) 
    adata_temp = adata_temp[,ind_select]
  }
} else {
  print('tissue all analysis')
}

if (is.na(keep_sex)==FALSE){
  if (keep_sex=='M'){
    print('keep only male')
    ind_select = (adata_temp$Gender=='Male')
    adata_temp = adata_temp[,ind_select]
  }
  if (keep_sex=='F'){
    print('keep only female')
    ind_select = (adata_temp$Gender=='Female')
    adata_temp = adata_temp[,ind_select]
  } 
} else {
  print('keep both sexes')
}

if (is.na(flag_M_O)==FALSE){
  if (flag_M_O=='T'){
    print('keep only Middle and Older')
    ind_select = (adata_temp$Age_Group=='Middle')|(adata_temp$Age_Group=='Older')
    adata_temp = adata_temp[,ind_select]
  } else {
    print('no age filtering')
  }
} else {
  print('no age filtering')
}


print('After filtering')
print(dim(adata_temp))

if (dim(adata_temp)[2]>0){
  
  # Prepare sca object
  adata_temp <- Seurat::as.SingleCellExperiment(adata_temp) # converting the Seurat object to SingleCellExperiment
  sca <- SceToSingleCellAssay(adata_temp, class = "SingleCellAssay")
  if(exists("sca")){ rm(adata_temp)};gc()
  hist(colData(sca)$detected)
  colData(sca)$n_genes = scale(colData(sca)$detected) # colData(sca)$n_genes = scale(colData(sca)$n_genes) # n_gene (CDR)
  hist(colData(sca)$n_genes)
  sca_filt = sca[rowSums(assay(sca)) != 0, ]
  if(exists("sca_filt")){ rm(sca) }
  gc() #free up memrory and report the memory usage.
  
  # Set flags 
  if (length(unique(sca_filt$Gender))>1){
    flag_sex=TRUE
  } else {
    flag_sex=FALSE
  }
  print(paste0('flag_sex: ', flag_sex))
  if (str_n_gene=='T'){
    flag_n_gene=TRUE
  } else {
    flag_n_gene=FALSE
  }
  print(paste0('flag_n_gene: ', flag_n_gene))
  
  # DGE testing 
  covariate = ''
  if (flag_sex==TRUE){
    if (is.na(others) == FALSE){
      covariate = paste0(covariate, " + Gender", " + ", others)
    } else {
      covariate = paste0(covariate, " + Gender")
    }
  }
  if (flag_n_gene==TRUE){
    covariate = paste0(covariate, " + nFeature_originalexp")
  }
  
  print(paste0('covariate: ', covariate))
  
  gc() #free up memrory and report the memory usage.
  hist(colData(sca_filt)$Age)
  options(mc.cores = 4L)
  zlmCond <- zlm(formula = as.formula(paste0("~Age", covariate)), sca=sca_filt,  parallel = TRUE)
  rm(sca_filt)
  gc() #free up memrory and report the memory usage.
  
  summaryCond <- summary(zlmCond, doLRT="Age")
  # summaryCond_BMI <- summary(zlmCond, doLRT="BMI")
  if (flag_sex==TRUE){
    summaryCond_sex <- summary(zlmCond, doLRT='GenderMale')
  }
  summaryCond_others_ls <- list()
  if (is.na(others) == FALSE){
    for (vv in unlist(stringr::str_split(others, "[+]")) ) {
      summaryCond_others_ls[[vv]] <- summary(zlmCond, doLRT=vv)
    }
  }
  
  # Summarize results 
  summaryDt <- summaryCond$datatable
  dt1 = summaryDt[contrast=="Age" & component=="H", .(primerid, `Pr(>Chisq)`)]
  dt2 = summaryDt[contrast=="Age" & component=="logFC", .(primerid, coef, z)]
  de_res = merge(dt1, dt2, by="primerid")
  colnames(de_res) <- c("gene", "age.H_p", "age.logFC", 'age.logFC_z')
  de_res$age.H_fdr <- p.adjust(de_res$age.H_p, "fdr")
  
  
  if (flag_sex==TRUE){
    summaryDt <- summaryCond_sex$datatable
    dt_sex1 = summaryDt[contrast=="GenderMale" & component=="H", .(primerid, `Pr(>Chisq)`)]
    dt_sex2 = summaryDt[contrast=="GenderMale" & component=="logFC", .(primerid, coef, z)]
    de_res_sex = merge(dt_sex1, dt_sex2, by="primerid")
    colnames(de_res_sex) <- c("gene", "GenderMale.H_p", "GenderMale.logFC", 'GenderMale.logFC_z')
    de_res_sex$GenderMale.H_fdr <- p.adjust(de_res_sex$GenderMale.H_p, "fdr")
    de_res = merge(de_res, de_res_sex, by="gene")
  }    
  
  if (is.na(others) == FALSE){
    de_res_others_ls <- list()
    for (vv in unlist(stringr::str_split(others, "[+]")) ) {
      # Summarize results
      summaryDt <- summaryCond_others_ls[[vv]]$datatable
      dt1 = summaryDt[contrast==vv & component=="H", .(primerid, `Pr(>Chisq)`)]
      dt2 = summaryDt[contrast==vv & component=="logFC", .(primerid, coef, z)]
      de_res = merge(dt1, dt2, by="primerid")
      colnames(de_res) <- c("gene", paste0(vv, ".H_p"), paste0(vv, ".logFC"), paste0(vv, ".logFC_z"))
      de_res[[paste0(vv, ".H_fdr")]] <- p.adjust(de_res[[paste0(vv, ".H_p")]], "fdr")
      de_res_others_ls[[vv]] <- de_res
    }
  }
  de_res_others <- Reduce(de_res_others_ls, by='gene')
  
  # Write results
  saveRDS(de_res, paste0(Disk, Project.folder, "/", Rds.folder, "/", output_folder, "/", 
                         unlist(strsplit(data_name, "[.Rds]"))[1], "_", tissue, '_Age.RegressionDEGs.Rds'))
  write.csv(de_res,paste0(Disk, Project.folder, "/", Rds.folder, "/", output_folder, "/", 
                          unlist(strsplit(data_name, "[.Rds]"))[1], "_", tissue, '_Age.RegressionDEGs.csv'), row.names=FALSE)
}
gc() #free up memory and report the memory usage.

print('Finished')
print(proc.time() - start_time)
