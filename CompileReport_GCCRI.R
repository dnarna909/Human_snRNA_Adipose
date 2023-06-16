# 1. Create Parameters file --------------------------------
#' create `Project Parameters.R` files for each project
##' counts foler should be `ProjectName_aggr`
unlink(".RData")
rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memory and report the memory usage.
memory.limit(size = 1e+13)
# serve = "GCCRI"
serve = ""
if (serve == "GCCRI"){
  Disk <- c("/projects/f4fced47-44e8-4b7e-84ee-b14d0fce0693/")
  Project.folder <- c("")
} else {
  Disk <- c("/media/jianie/Extreme SSD/") 
  # Disk <- c("E:/")
  # Disk <- c("K:/")
  Project.folder <- c("2022-09-01 STARR_SGLT2 Combine")
  codes.folder <- c("Human_snRNA_Adipose")
  html.folder <- "html.result"
}


# 2. prerun ----------------------------------------------------------------------------------------------------------
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Step01_prerun_Creat_matrices.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), paste0('Step01_prerun_Creat_matrices.html'))) # 
#' generate: combine `project name_mat.rds`, feature.names.rds,  only need run once.

rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Step02_prerun_QC_Aggr_all.genes.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), paste0('Step02_prerun_QC_Aggr_all.genes.html'))) # 
#' generate: individual `sample number_mat.Rds`, combined project `dataset.group.Rds`, combined project `mat.Rds`(list), in the counts folder. only need run once.

# 3. QC ----------------------------------------------------------------------------------------------------------
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Step03_QC_Aggr_all.genes_step1.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), paste0('Step03_QC_Aggr_all.genes_step1.html'))) # 
#' QC, filter, and generate: `Ambient.Rds`, `QC files_Genes.RData`, `QC files_Cells.RData`, `Final.sce.list.Rds`, .only need run once.

rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Step03_QC_Aggr_all.genes_step2.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), paste0('Step03_QC_Aggr_all.genes_step2.html'))) # 
#' QC, filter, and generate: `sce.list folder`, `sce.list2 folder`, 
#' `Final.sce.list_part1.Rds`, `Final.sce.list_part2.Rds`, `Final.sce.list_part3.Rds`. only need run once.

# if samples <= 36
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Step04_QC_Aggr_all.genes_step3.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), paste0('Step04_QC_Aggr_all.genes_step3.html'))) # 
#' QC, filter, and generate: `rescaled.Rds`, `seurat.list folder`. only need run once.

# if samples > 36
# run  Step04_QC_Aggr_all.genes_step3.collab.R on colab:https://colab.research.google.com/#create=true&language=r
#' QC, filter, and generate:  `rescaled.Rds` or `rescaled_part1.Rds`, `rescaled_part2.Rds`, `rescaled_part3.Rds`,
#' `seurat.list.combine`: `seurat.list.rds`, or `seurat.list_part1.Rds`, `seurat.list_part2.Rds`,`seurat.list_part3.Rds`,`seurat.list_part4.Rds`. only need run once.
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Step04_QC_Aggr_all.genes_step4_after.collab.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), paste0('Step04_QC_Aggr_all.genes_step4_after.collab.html'))) # this step DONE in Lunix PC
#' seperate to `seurat.list folder`


rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Step05_Combine_seurat_data.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), paste0('Step05_Combine_seurat_data.html'))) # this step DONE in Lunix PC
#' `SAT folder`: combined `SAT1.Rds`, combined `SAT2.Rds`, or combined `SAT.Rds` (depending on sample number).only need run once.


rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Step06_Combine_seurat_data.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), paste0('Step06_Combine_seurat_data.html'))) # this step DONE in Lunix PC
#' `Seurat folder`: combined `SAT.Rds` (if sample number > 18 & <= 36).only need run once. 
#'run Step06_Combine_seurat_data.colab.R on colab:https://colab.research.google.com/#create=true&language=r
#' 25 min upload, 35 min running, 35 min download. 
#' save a copy of SAT.rds file (optional)

rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Step07_QC_Aggr_all.genes_step3.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), paste0('Step07_QC_Aggr_all.genes_step3.html'))) # 
#' import `SAT.Rds` for Seurat analysis. generate: `SAT.basic.meta.data.Rds`, 
#' `SAT_basic.meta.data.Rds`
#' `SAT_Label_1000.Rds`, `SAT_1000.meta.data.Rds`.
#' `SAT_Label_2000.Rds`, `SAT_2000.meta.data.Rds`.
#' `SAT_Label_3000.Rds`, `SAT_3000.meta.data.Rds`.

# Delete the agrregation balance samples
seurat.files = c("SAT_Label_1000_ori.Rds", "SAT_Label_2000_ori.Rds", "SAT_Label_3000_ori.Rds")
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Adjust_Balance_Samples_Aggregation.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Adjust_Balance_Samples_Aggregation', '.html')) )# 11∶24∶58 PM 
rm(seurat.files)

#' plot QC
# rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Step08_QC_Aggr_all.genes_step4.R"), output_format= "html_document",
#                   output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), paste0('Step08_QC_Aggr_all.genes_step4.html'))) # 

## Add sample metadata -------------------------------------
rmarkdown::render(input= paste0(Disk, Project.folder, "/","Metadata_Add_Samples.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Metadata_Add_Samples', '.html')) )# 11∶24∶58 PM 

## plot samples
export.folder = "Samples"
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Plotting_Samples.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Plotting_Samples', '.html')))
rm(export.folder)



seurat.files = c("SAT_Label_1000.Rds", "SAT_Label_2000.Rds", "SAT_Label_3000.Rds")
rmarkdown::render(input= paste0(Disk, Project.folder, "/","Metadata_Add_Samples_seurat_Plotting.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Metadata_Add_Samples_seurat_Plotting', '.html')) )# 11∶24∶58 PM 
rm(seurat.files)



# 4. Test Resolution --------------------------------------------------------------
## test all resolutions
seurat.file = "SAT_Label_3000.Rds"
res.range <- c(0.005, 0.01); res.by = 0.001
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_00_Step1_Test_Resolution.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_00_Step1_Test_Resolution', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 "_",  res.range[1], "_", res.range[2], '.html'))) # 

seurat.file = "SAT_Label_2000.Rds"
res.range <- c(0.005, 0.01); res.by = 0.001
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_00_Step1_Test_Resolution.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_00_Step1_Test_Resolution', stringr::str_split(seurat.file, ".Rds")[[1]][1],
                                                 "_",  res.range[1], "_", res.range[2], '.html'))) # 


# 5. cell type analysis -------------------------------------------------------------------------------------------------
# Overall: Annotation
#' check Resolution output files for best resolution and variable gene number
## compare clusters and select best resolution
seurat.file = "SAT_Label_2000.Rds"
res = "_0.005_0.01"  
select.res = "originalexp_snn_res.0.006"
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_01_Fat_Annotation_Compare.rmd"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_01_Human_Fat_Annotation_Compare_', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 stringr::str_split(select.res, "originalexp_snn_res")[[1]][2],  '.html')) ) # 11∶17∶45 PM
rm(seurat.file, res, select.res)

seurat.file = "SAT_Label_3000.Rds"
res = "_0.005_0.01"  
select.res = "originalexp_snn_res.0.005"
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_01_Fat_Annotation_Compare.rmd"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_01_Human_Fat_Annotation_Compare_', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 stringr::str_split(select.res, "originalexp_snn_res")[[1]][2],  '.html')) ) # 11∶17∶45 PM
rm(seurat.file, res, select.res)


## apply final resolution to seurat data
seurat.file = "SAT_Label_2000.Rds"
res = "_0.005_0.01"
select.res = "originalexp_snn_res.0.006"
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_00_Step2_Final_Resolution2.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_00_Step2_Final_Resolution_', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 stringr::str_split(select.res, "originalexp_snn_res")[[1]][2], '.html')) ) # 
rm(seurat.file, res, select.res)

## Add annotation metadata --------------------------------------------------------------------------
seurat.file = "SAT_Label_2000.Rds"
sample_annotation_file = "Fat_Subtype.xlsx"
sheets = c("SAT")
# need extra files: Fat_Subtype.xlsx, _Final_Res.metadata.Rds
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Metadata_Add_Annotation.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Metadata_Annotation_', stringr::str_split(seurat.file, ".Rds")[[1]][1], '.html')) )# 11∶24∶58 PM 
rm(seurat.file, sample_annotation_file, sheets)




# 6. cell subtype analysis -------------------------------------------------------------------------------------------------
## subset seurat ----
# seurat.files = c(
#   "SAT_Label_2000.Rds"
# )
# file.type = "Annotation.metadata.Rds"; # with Annotation.file.name # subset seurat.data
# subset.group = "Annotation" # subset seurat.data
# rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Subset_seurat.R"), output_format= "html_document",
#                   output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
#                                           paste0('Subset_seurat_', 
#                                                  stringr::str_split(seurat.files[1], "_")[[1]][1],"_", 
#                                                  stringr::str_split(file.type, ".")[[1]][1],'.html')))
# rm(seurat.files, file.type, subset.group)


## subset seurat by selection ----
seurat.files = c(
  "SAT_Label_2000.Rds"
)
file.type = "SAT_Sample.meta.data.Rds"; # with Annotation.file.name # subset seurat.data
subset.groups = c("Compare_Group1", "Compare_Group2", "Compare_Group3") # subset seurat.data
rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Subset_seurat_selection.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder),
                                          paste0('Subset_seurat_selection_',
                                                 stringr::str_split(seurat.files[1], "_")[[1]][1],"_",
                                                 stringr::str_split(file.type, ".")[[1]][1],'.html')))
rm(seurat.files, file.type, subset.group)

## FAP ------------------------------------------------------
## create seurat data
seurat.file = "SAT_Label_2000.Rds"
celltype.file =  "FAP.Rds"
rmarkdown::render(input= paste0(Disk, Project.folder, "/","Main_04_FAP_Fat_1_Step1.rmd"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_04_FAP_Fat_1_Step1_', stringr::str_split(seurat.file, ".Rds")[[1]][1], '.html'))) # 11∶32∶10 PM 

## test all resolutions
seurat.file = "SAT_Label_2000_FAP.Rds"
res.range <- c(0.02, 0.16); res.by = 0.02
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_00_Step1_Test_Resolution.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_00_Step1_Test_Resolution_', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 "_",  res.range[1], "_", res.range[2], '.html'))) # 
rm(seurat.file, res.range, res.by)

seurat.file = "SAT_Label_2000_FAP.Rds"
res.range <- c(0.1, 0.5); res.by = 0.05
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_00_Step1_Test_Resolution.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_00_Step1_Test_Resolution_', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 "_",  res.range[1], "_", res.range[2], '.html'))) # 
rm(seurat.file, res.range, res.by)

## compare clusters and select best resolution
seurat.file = "SAT_Label_2000_FAP.Rds"
res = "_0.02_0.16"  
select.res = "originalexp_snn_res.0.16"
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_04_Human_FAP_Subtype_Compare.rmd"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_04_Human_FAP_Subtype_Compare_', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 stringr::str_split(select.res, "originalexp_snn_res")[[1]][2], '.html'))) # 11∶36∶24 PM 
rm(seurat.file, res, select.res)

seurat.file = "SAT_Label_2000_FAP.Rds"
res = "_0.1_0.5"  
select.res = "originalexp_snn_res.0.45"
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_04_Human_FAP_Subtype_Compare.rmd"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_04_Human_FAP_Subtype_Compare_', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 stringr::str_split(select.res, "originalexp_snn_res")[[1]][2], '.html'))) # 11∶36∶24 PM 
rm(seurat.file, res, select.res)

seurat.file = "SAT_Label_2000_FAP.Rds"
res = "_0.1_0.5"  
select.res = "originalexp_snn_res.0.35"
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_04_Human_FAP_Subtype_Compare.rmd"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_04_Human_FAP_Subtype_Compare_', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 stringr::str_split(select.res, "originalexp_snn_res")[[1]][2], '.html'))) # 11∶36∶24 PM 
rm(seurat.file, res, select.res)


## combine average data
seurat.file = "SAT_Label_2000_FAP.Rds"
res = "_0.1_0.5"  
select.res = "originalexp_snn_res.0.45"
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "CombineAverageData.rmd"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('CombineAverageData_', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 stringr::str_split(select.res, "originalexp_snn_res")[[1]][2], '.html'))) # 11∶36∶24 PM 
rm(seurat.file, res, select.res)

seurat.file = "SAT_Label_2000_FAP.Rds"
res = "_0.1_0.5"  
select.res = "originalexp_snn_res.0.35"
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "CombineAverageData.rmd"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('CombineAverageData_', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 stringr::str_split(select.res, "originalexp_snn_res")[[1]][2], '.html'))) # 11∶36∶24 PM 
rm(seurat.file, res, select.res)

seurat.file = "SAT_Label_2000_FAP.Rds"
res = "_0.02_0.16"  
select.res = "originalexp_snn_res.0.16"
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "CombineAverageData.rmd"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('CombineAverageData_', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 stringr::str_split(select.res, "originalexp_snn_res")[[1]][2], '.html'))) # 11∶36∶24 PM 
rm(seurat.file, res, select.res)

## apply final resolution to seurat data
seurat.file = "SAT_Label_2000_FAP.Rds"
res = "_0.1_0.5"
select.res = "originalexp_snn_res.0.45"
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_00_Step2_Final_Resolution2.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_00_Step2_Final_Resolution_', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 stringr::str_split(select.res, "originalexp_snn_res")[[1]][2],  '.html')) ) # 
rm(seurat.file, res, select.res)

## Adipocytes ------------------------------------------------------------
## create seurat data
seurat.file = "SAT_Label_2000.Rds"
celltype.file =  "Adipocytes.Rds"
rmarkdown::render(input= paste0(Disk, Project.folder, "/","Main_05_Adipocytes_Fat_1_Step1.rmd"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_05_Adipocytes_Fat_1_Step1_', stringr::str_split(seurat.file, ".Rds")[[1]][1], '.html'))) # 11∶59∶03 PM

## test all resolutions
seurat.file = "SAT_Label_2000_Adipocytes.Rds"
res.range <- c(0.01, 0.05); res.by = 0.005
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_00_Step1_Test_Resolution.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_00_Step1_Test_Resolution', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 "_",  res.range[1], "_", res.range[2], '.html'))) # 
rm(seurat.file, res.range, res.by)

seurat.file = "SAT_Label_2000_Adipocytes.Rds"
res.range <- c(0.05, 0.5); res.by = 0.05
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_00_Step1_Test_Resolution.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_00_Step1_Test_Resolution', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 "_",  res.range[1], "_", res.range[2], '.html'))) # 
rm(seurat.file, res.range, res.by)

seurat.file = "SAT_Label_2000_Adipocytes.Rds"
res.range <- c(0.5, 1); res.by = 0.1
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_00_Step1_Test_Resolution.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_00_Step1_Test_Resolution', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 "_",  res.range[1], "_", res.range[2], '.html'))) # 
rm(seurat.file, res.range, res.by)

# compare clusters and select best resolution
seurat.file = "SAT_Label_2000_Adipocytes.Rds"
res = "_0.05_0.5"   
select.res = "originalexp_snn_res.0.05"
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_05_Human_Adipocytes_Subtype_Compare.rmd"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder),
                                          paste0('Main_05_Human_Adipocytes_Subtype_Compare_', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 stringr::str_split(select.res, "originalexp_snn_res")[[1]][2], '.html'))) # 12∶07∶55 AM
rm(seurat.file, res, select.res)

seurat.file = "SAT_Label_2000_Adipocytes.Rds"
res = "_0.05_0.5" 
select.res = "originalexp_snn_res.0.4"
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_05_Human_Adipocytes_Subtype_Compare.rmd"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder),
                                          paste0('Main_05_Human_Adipocytes_Subtype_Compare_', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 stringr::str_split(select.res, "originalexp_snn_res")[[1]][2], '.html'))) # 12∶07∶55 AM
rm(seurat.file, res, select.res)

seurat.file = "SAT_Label_2000_Adipocytes.Rds"
res = "_0.05_0.5" 
select.res = "originalexp_snn_res.0.25"
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_05_Human_Adipocytes_Subtype_Compare.rmd"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder),
                                          paste0('Main_05_Human_Adipocytes_Subtype_Compare_', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 stringr::str_split(select.res, "originalexp_snn_res")[[1]][2], '.html'))) # 12∶07∶55 AM
rm(seurat.file, res, select.res)

seurat.file = "SAT_Label_2000_Adipocytes.Rds"
res = "_0.05_0.5" 
select.res = "originalexp_snn_res.0.5"
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_05_Human_Adipocytes_Subtype_Compare.rmd"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder),
                                          paste0('Main_05_Human_Adipocytes_Subtype_Compare_', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 stringr::str_split(select.res, "originalexp_snn_res")[[1]][2], '.html'))) # 12∶07∶55 AM
rm(seurat.file, res, select.res)

## combine average data
seurat.file = "SAT_Label_2000_Adipocytes.Rds"
res = "_0.05_0.5"
select.res = "originalexp_snn_res.0.05"
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "CombineAverageData.rmd"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('CombineAverageData_', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 stringr::str_split(select.res, "originalexp_snn_res")[[1]][2], '.html'))) # 11∶36∶24 PM 
rm(seurat.file, res, select.res)

seurat.file = "SAT_Label_2000_Adipocytes.Rds"
res = "_0.05_0.5"
select.res = "originalexp_snn_res.0.4"
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "CombineAverageData.rmd"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('CombineAverageData_', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 stringr::str_split(select.res, "originalexp_snn_res")[[1]][2], '.html'))) # 11∶36∶24 PM 
rm(seurat.file, res, select.res)

seurat.file = "SAT_Label_2000_Adipocytes.Rds"
res = "_0.05_0.5"
select.res = "originalexp_snn_res.0.5"
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "CombineAverageData.rmd"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('CombineAverageData_', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 stringr::str_split(select.res, "originalexp_snn_res")[[1]][2], '.html'))) # 11∶36∶24 PM 
rm(seurat.file, res, select.res)

seurat.file = "SAT_Label_2000_Adipocytes.Rds"
res = "_0.05_0.5"
select.res = "originalexp_snn_res.0.25"
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "CombineAverageData.rmd"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('CombineAverageData_', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 stringr::str_split(select.res, "originalexp_snn_res")[[1]][2], '.html'))) # 11∶36∶24 PM 
rm(seurat.file, res, select.res)

## apply final resolution to seurat data
seurat.file = "SAT_Label_2000_Adipocytes.Rds"
res = "_0.05_0.5"
select.res = "originalexp_snn_res.0.25"
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_00_Step2_Final_Resolution2.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_00_Step2_Final_Resolution_', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 stringr::str_split(select.res, "originalexp_snn_res")[[1]][2],  '.html')) ) # 
rm(seurat.file, res, select.res)


## Endothelial --------------------------------------------------------------------------------------
## create seurat data
seurat.file = "SAT_Label_2000.Rds"
celltype.file =  "Endothelial.Rds"
rmarkdown::render(input= paste0(Disk, Project.folder, "/","Main_03_Endothelial_Fat_1_Step1.rmd"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_03_Endothelial_Fat_1_Step1_', stringr::str_split(seurat.file, ".Rds")[[1]][1], '.html'))) # 12∶36∶31 AM

## test all resolutions
seurat.file = "SAT_Label_2000_Endothelial.Rds"
res.range <- c(0.01, 0.1); res.by = 0.01
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_00_Step1_Test_Resolution.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_00_Step1_Test_Resolution', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 "_",  res.range[1], "_", res.range[2], '.html'))) # 
rm(seurat.file, res.range, res.by)

seurat.file = "SAT_Label_2000_Endothelial.Rds"
res.range <- c(0.1, 1);  res.by = 0.1
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_00_Step1_Test_Resolution.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_00_Step1_Test_Resolution', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 "_",  res.range[1], "_", res.range[2], '.html'))) # 
rm(seurat.file, res.range, res.by)

# compare clusters and select best resolution
seurat.file = "SAT_Label_2000_Endothelial.Rds"
res = "_0.01_0.1"
select.res = "originalexp_snn_res.0.01"
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_03_Human_Endothelial_Subtype_Compare.rmd"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder),
                                          paste0('Main_03_Human_Endothelial_Subtype_Compare_', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 stringr::str_split(select.res, "originalexp_snn_res")[[1]][2], '.html'))) #  12∶38∶06 AM
rm(seurat.file, res, select.res)

seurat.file = "SAT_Label_2000_Endothelial.Rds"
res = "_0.01_0.1"
select.res = "originalexp_snn_res.0.05"
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_03_Human_Endothelial_Subtype_Compare.rmd"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder),
                                          paste0('Main_03_Human_Endothelial_Subtype_Compare_', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 stringr::str_split(select.res, "originalexp_snn_res")[[1]][2], '.html'))) #  12∶38∶06 AM
rm(seurat.file, res, select.res)

seurat.file = "SAT_Label_2000_Endothelial.Rds"
res = "_0.1_1"
select.res = "originalexp_snn_res.0.4"
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_03_Human_Endothelial_Subtype_Compare.rmd"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder),
                                          paste0('Main_03_Human_Endothelial_Subtype_Compare_', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 stringr::str_split(select.res, "originalexp_snn_res")[[1]][2], '.html'))) #  12∶38∶06 AM
rm(seurat.file, res, select.res)

seurat.file = "SAT_Label_2000_Endothelial.Rds"
res = "_0.1_1"
select.res = "originalexp_snn_res.0.1"
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_03_Human_Endothelial_Subtype_Compare.rmd"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder),
                                          paste0('Main_03_Human_Endothelial_Subtype_Compare_', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 stringr::str_split(select.res, "originalexp_snn_res")[[1]][2], '.html'))) #  12∶38∶06 AM
rm(seurat.file, res, select.res)

## combine average data
seurat.file = "SAT_Label_2000_Endothelial.Rds"
res = "_0.01_0.1"
select.res = "originalexp_snn_res.0.01"
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "CombineAverageData.rmd"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('CombineAverageData_', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 stringr::str_split(select.res, "originalexp_snn_res")[[1]][2], '.html'))) # 11∶36∶24 PM 
rm(seurat.file, res, select.res)

seurat.file = "SAT_Label_2000_Endothelial.Rds"
res = "_0.01_0.1"
select.res = "originalexp_snn_res.0.05"
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "CombineAverageData.rmd"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('CombineAverageData_', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 stringr::str_split(select.res, "originalexp_snn_res")[[1]][2], '.html'))) # 11∶36∶24 PM 
rm(seurat.file, res, select.res)

seurat.file = "SAT_Label_2000_Endothelial.Rds"
res = "_0.01_0.1"
select.res = "originalexp_snn_res.0.1"
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "CombineAverageData.rmd"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('CombineAverageData_', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 stringr::str_split(select.res, "originalexp_snn_res")[[1]][2], '.html'))) # 11∶36∶24 PM 
rm(seurat.file, res, select.res)

seurat.file = "SAT_Label_2000_Endothelial.Rds"
res = "_0.1_1"
select.res = "originalexp_snn_res.0.4"
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "CombineAverageData.rmd"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('CombineAverageData_', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 stringr::str_split(select.res, "originalexp_snn_res")[[1]][2], '.html'))) # 11∶36∶24 PM 
rm(seurat.file, res, select.res)

## apply final resolution to seurat data
seurat.file = "SAT_Label_2000_Endothelial.Rds"
res = "_0.01_0.1"
select.res = "originalexp_snn_res.0.05"
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_00_Step2_Final_Resolution2.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_00_Step2_Final_Resolution_', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 stringr::str_split(select.res, "originalexp_snn_res")[[1]][2], '.html')) ) # 
rm(seurat.file, res, select.res)


## Immune ---------------------------------------------------------------------------------------------
## create seurat data
seurat.file = "SAT_Label_2000.Rds"
celltype.file =  "Immune.Rds"
rmarkdown::render(input= paste0(Disk, Project.folder, "/","Main_02_Immune_Fat_1_Step1.rmd"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_02_Immune_Fat_1_Step1_', stringr::str_split(seurat.file, ".Rds")[[1]][1], '.html'))) # 12∶22∶38 AM

## test all resolutions
seurat.file = "SAT_Label_2000_Immune.Rds"
res.range <- c(0.01, 0.1); res.by = 0.01
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_00_Step1_Test_Resolution.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_00_Step1_Test_Resolution', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 "_",  res.range[1], "_", res.range[2], '.html'))) # 
rm(seurat.file, res.range, res.by)

## compare clusters and select best resolution
seurat.file = "SAT_Label_2000_Immune.Rds"
res = "_0.01_0.1"
select.res = "originalexp_snn_res.0.01"
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_02_Human_Immune_Subtype_Compare.rmd"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_02_Human_Immune_Subtype_Compare_', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 stringr::str_split(select.res, "originalexp_snn_res")[[1]][2], '.html'))) # 12∶25∶47 AM
rm(seurat.file, res, select.res)

seurat.file = "SAT_Label_2000_Immune.Rds"
res = "_0.01_0.1"
select.res = "originalexp_snn_res.0.04"
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_02_Human_Immune_Subtype_Compare.rmd"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_02_Human_Immune_Subtype_Compare_', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 stringr::str_split(select.res, "originalexp_snn_res")[[1]][2], '.html'))) # 12∶25∶47 AM
rm(seurat.file, res, select.res)

seurat.file = "SAT_Label_2000_Immune.Rds"
res = "_0.01_0.1"
select.res = "originalexp_snn_res.0.1"
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_02_Human_Immune_Subtype_Compare.rmd"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_02_Human_Immune_Subtype_Compare_', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 stringr::str_split(select.res, "originalexp_snn_res")[[1]][2], '.html'))) # 12∶25∶47 AM
rm(seurat.file, res, select.res)

## combine average data
seurat.file = "SAT_Label_2000_Immune.Rds"
res = "_0.01_0.1"
select.res = "originalexp_snn_res.0.01"
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "CombineAverageData.rmd"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('CombineAverageData_', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 stringr::str_split(select.res, "originalexp_snn_res")[[1]][2], '.html'))) # 11∶36∶24 PM 
rm(seurat.file, res, select.res)

seurat.file = "SAT_Label_2000_Immune.Rds"
res = "_0.01_0.1"
select.res = "originalexp_snn_res.0.04"
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "CombineAverageData.rmd"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('CombineAverageData_', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 stringr::str_split(select.res, "originalexp_snn_res")[[1]][2], '.html'))) # 11∶36∶24 PM 
rm(seurat.file, res, select.res)

seurat.file = "SAT_Label_2000_Immune.Rds"
res = "_0.01_0.1"
select.res = "originalexp_snn_res.0.1"
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "CombineAverageData.rmd"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('CombineAverageData_', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 stringr::str_split(select.res, "originalexp_snn_res")[[1]][2], '.html'))) # 11∶36∶24 PM 
rm(seurat.file, res, select.res)

## apply final resolution to seurat data
seurat.file = "SAT_Label_2000_Immune.Rds"
res = "_0.01_0.1"
select.res = "originalexp_snn_res.0.01"
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_00_Step2_Final_Resolution2.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_00_Step2_Final_Resolution_', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 stringr::str_split(select.res, "originalexp_snn_res")[[1]][2], '.html')) ) # 
rm(seurat.file, res, select.res)



## FAP and Adipocytes combine --------------------------------------------------------


# 7. Add subtype metadata  ------------------------------------------------------------------------
seurat.files = c("SAT_Label_2000_FAP.Rds", "SAT_Label_2000_Endothelial.Rds",
                 "SAT_Label_2000_Immune.Rds", "SAT_Label_2000_Adipocytes.Rds")
# need extra files: Fat_Subtype.xlsx, Annotation.metadata.Rds, Final_Res.metadata.Rds,
for (seurat.file in seurat.files) {
  # seurat.file = "SAT_Label_2000_FAP.Rds"
  rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Metadata_Add_Subtype.R"), output_format= "html_document",
                    output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                            paste0('Metadata_Subtype_', stringr::str_split(seurat.file, ".Rds")[[1]][1], '.html')))
}


## combine
seurat.file = "SAT_Label_2000.Rds"
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Metadata_Combine_Subtype.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Metadata_Combine_Subtype_', stringr::str_split(seurat.file, ".Rds")[[1]][1], '.html')))
rm(seurat.file)



## subset seurat ----
seurat.files = c(
  "SAT_Label_2000_Immune.Rds"
  , "SAT_Label_2000_FAP.Rds"
  , "SAT_Label_2000_Adipocytes.Rds"
  , "SAT_Label_2000_Endothelial.Rds"
)
file.type = "Subtype.metadata.Rds"; # with Annotation.file.name # subset seurat.data
subset.group = "Subtype" # subset seurat.data
rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Subset_seurat.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Subset_seurat_', 
                                                 stringr::str_split(seurat.files[1], "_")[[1]][1],"_", 
                                                 stringr::str_split(file.type, ".")[[1]][1],'.html')))
rm(seurat.files, file.type, subset.group)


# 8. Senescence Score analysis -------------------------------------------------------------------------------------------------
seurat.file = "SAT_Label_2000.Rds"
file.type = "Subtype.metadata.Rds"
condition1 = "Annotation"; condition2 = "Subtype";
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_06_Senescence_GeneModule_Compare.rmd"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_06_Senescence_', stringr::str_split(seurat.file, ".Rds")[[1]][1], '.html')))
rm(seurat.file)

# combine values
seurat.files = c("SAT_Label_2000.Rds"
                 , "SAT_Label_2000_FAP.Rds"
                 , "SAT_Label_2000_Endothelial.Rds"
                 , "SAT_Label_2000_Immune.Rds"
                 , "SAT_Label_2000_Adipocytes.Rds"
);
col.groups1 = c("Senescence.GO.score", "Senescence.score", "SenMayo");
file.type1 = "Senescence_SAT_Label_2000/Senescence.Marker_ModuleScore_Label_metadata.rds";
col.groups = c("Annotation", "Subtype");
row.groups = "Dataset";
file.type = "Subtype.metadata.Rds";
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Combine_ModuleScore_AveragebyDataset.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Combine_ModuleScore_AveragebyDataset', 
                                                 '.html')))
rm(seurat.files, col.groups, row.groups, file.type, col.groups1, file.type1)


# plotting by grouups
seurat.files = c("SAT_Label_2000.Rds"
                 , "SAT_Label_2000_FAP.Rds"
                 , "SAT_Label_2000_Endothelial.Rds"
                 , "SAT_Label_2000_Immune.Rds"
                 , "SAT_Label_2000_Adipocytes.Rds"
)
col.groups1 = c("Senescence.GO.score", "Senescence.score", "SenMayo")
file.type1 = "Senescence_SAT_Label_2000/Senescence.Marker_ModuleScore_Label_metadata.rds"
col.groups = c("Annotation", "Subtype")
row.groups = "Dataset"
file.type = "Subtype.metadata.Rds"
export.folder = "ModuleScore"
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Plotting_Figure_Analysis_ModuleScore_Pre_Post.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Plotting_Figure_Analysis_ModuleScore_Pre_Post', 
                                                 '.html')))
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Plotting_Figure_Analysis_ModuleScore_NoTreatment.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Plotting_Figure_Analysis_ModuleScore_NoTreatment', 
                                                 '.html')))
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Plotting_Figure_Analysis_ModuleScore_Pre_Post_Average.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Plotting_Figure_Analysis_ModuleScore_Pre_Post_Average', 
                                                 '.html')))
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Plotting_Figure_Analysis_ModuleScore_NoTreatment_Average.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Plotting_Figure_Analysis_ModuleScore_NoTreatment_Average', 
                                                 '.html')))
rm(seurat.files, col.groups, row.groups, file.type, col.groups1, file.type1, export.folder)



# Average then compare by 1 way ANOVA, 2 way ANOVA,
seurat.files = c("SAT_Label_2000.Rds"
                 , "SAT_Label_2000_FAP.Rds"
                 , "SAT_Label_2000_Endothelial.Rds"
                 , "SAT_Label_2000_Immune.Rds"
                 , "SAT_Label_2000_Adipocytes.Rds"
)
col.groups1 = c("Senescence.GO.score", "Senescence.score", "SenMayo")
file.type1 = "Senescence_SAT_Label_2000/Senescence.Marker_ModuleScore_Label_metadata.rds"
col.groups = c("Annotation", "Subtype")
row.groups = "Dataset"
file.type = "Subtype.metadata.Rds"
export.folder = "ModuleScore"
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Plotting_Figure_Analysis_ModuleScore_Pre_Post_Average.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Plotting_Figure_Analysis_ModuleScore_Pre_Post_Average', 
                                                 '.html')))
rm(seurat.files, col.groups, row.groups, file.type, col.groups1, file.type1, export.folder)
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Plotting_Figure_Analysis_ModuleScore_NoTreatment_Average.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Plotting_Figure_Analysis_ModuleScore_NoTreatment_Average', 
                                                 '.html')))
rm(seurat.files, col.groups, row.groups, file.type, col.groups1, file.type1, export.folder)


# Add Senescence in metadata 
seurat.file = "SAT_Label_2000.Rds"
file.type1 = "Senescence_SAT_Label_2000/Senescence.Marker_ModuleScore_Label_metadata.rds"
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Metadata_Add_Senescence.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Metadata_Senescence_', stringr::str_split(seurat.file, ".Rds")[[1]][1], '.html')) )# 11∶24∶58 PM 
rm(seurat.file)




# 9. Composition Analysis -----------------------------------------
seurat.files = c("SAT_Label_2000.Rds",
                 "SAT_Label_2000_FAP.Rds", "SAT_Label_2000_Endothelial.Rds",
                 "SAT_Label_2000_Immune.Rds", "SAT_Label_2000_Adipocytes.Rds")
col.groups = c("Annotation", "Subtype")
row.group = "Dataset"
file.type = "Subtype.metadata.Rds"
# need extra file: Samples.df.Rds, _Sample.meta.data.Rds
for (seurat.file in seurat.files) {
  # seurat.file = "SAT_Label_2000.Rds"
  rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Analysis_Composition.rmd"), output_format= "html_document",
                    output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                            paste0('Analysis_Composition_', stringr::str_split(seurat.file, ".Rds")[[1]][1],"_",
                                                   stringr::str_split(file.type, ".Rds")[[1]][1], '.html')))
  # print(seurat.file)
}
rm(seurat.files, col.groups, row.group, file.type)

# combine all composition
export.folder = "Compositions"
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Combine_Composition.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Combine_Composition', 
                                                 # stringr::str_split(file.type, ".Rds")[[1]][1], 
                                                 '.html')))
rm(export.folder)

# # "Senescence"
# seurat.files = c("SAT_Label_2000.Rds",
#                  "SAT_Label_2000_FAP.Rds", "SAT_Label_2000_Endothelial.Rds",
#                  "SAT_Label_2000_Immune.Rds", "SAT_Label_2000_Adipocytes.Rds")
# # "Annotation", "Subtype"
# col.groups1 = c("Senescence.GO")
# file.type1 = "Senescence_SAT_Label_2000/Senescence_ModuleScore_metadata.Rds"
# col.groups = c("Annotation", "Subtype")
# row.group = "Dataset"
# file.type = "Subtype.metadata.Rds"
# for (seurat.file in seurat.files) {
#   # seurat.file = "SAT_Label_2000.Rds"
#   rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Analysis_Composition_Level2.rmd"), output_format= "html_document",
#                     output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
#                                             paste0('Analysis_Composition_', stringr::str_split(seurat.file, ".Rds")[[1]][1],"_",
#                                                    stringr::str_split(file.type, ".Rds")[[1]][1], "_",
#                                                    col.groups1, '.html')))
#   # print(seurat.file)
# }
# rm(seurat.file, seurat.files, col.groups, row.group, file.type, col.groups1, file.type1)

# seurat.files = c("SAT_Label_2000.Rds",
#                  "SAT_Label_2000_FAP.Rds", "SAT_Label_2000_Endothelial.Rds",
#                  "SAT_Label_2000_Immune.Rds", "SAT_Label_2000_Adipocytes.Rds")
# # "Annotation", "Subtype"
# file.type1 = "Senescence.name.metadata.Rds"; col.groups1 = "sene.name"
# col.groups = c("Annotation", "Subtype")
# row.group = "Dataset"
# file.type = "Subtype.metadata.Rds"
# for (seurat.file in seurat.files) {
#   # seurat.file = "SAT_Label_2000.Rds"
#   rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Analysis_Composition_Level2.rmd"), output_format= "html_document",
#                     output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
#                                             paste0('Analysis_Composition_', stringr::str_split(seurat.file, ".Rds")[[1]][1], "_",
#                                                    stringr::str_split(file.type, ".Rds")[[1]][1], "_",
#                                                    col.groups1, '.html')))
#   # print(seurat.file)
# }
# rm(seurat.file, seurat.files, col.groups, row.group, file.type, col.groups1, file.type1)





# 10. Plotting and stats analysis of all variables by Datasets -------------
## 1. plotting by group -------------
export.folder = "Compositions"
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Plotting_Figure_Composition_NoTreatment_Level2.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Plotting_Figure_Composition_Group_Level2_', 
                                                 # stringr::str_split(file.type, ".Rds")[[1]][1], 
                                                 '.html')))
rm(export.folder)

export.folder = "Compositions"
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Plotting_Figure_Composition_Pre_Post_Level2.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Plotting_Figure_Composition_Pre_Post_Level2_', 
                                                 # stringr::str_split(file.type, ".Rds")[[1]][1], 
                                                 '.html')))
rm(export.folder)

## 2. plot bar graph -------
file.variables = "SAT_Label_2000_Subtype_Composition.combine.Rds";
Sample.meta.file = "SAT_Sample.meta.data.Rds";
Select_Group = "Compare_Group1"; compare.group = "BA_Group"; y.lab = "Percentage of Cells (%)";
export.folder = "Variables"
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Plotting_Analysis_Variables_NoTreatment.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Plotting_Analysis_Variables_NoTreatment_', 
                                                 stringr::str_split(file.variables, ".Rds")[[1]][1], 
                                                 '.html')))
rm(file.variables, Sample.meta.file, Select_Group, compare.group, y.lab, export.folder)

file.variables = "SAT_Label_2000_Senescence_ModuleScoreAve.df.Rds";
Sample.meta.file = "SAT_Sample.meta.data.Rds";
Select_Group = "Compare_Group1"; compare.group = "BA_Group"; y.lab = "Senescence Score";
export.folder = "Variables"
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Plotting_Analysis_Variables_NoTreatment.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Plotting_Analysis_Variables_NoTreatment_', 
                                                 stringr::str_split(file.variables, ".Rds")[[1]][1], 
                                                 '.html')))
rm(file.variables, Sample.meta.file, Select_Group, compare.group, y.lab, export.folder)


file.variables = "SAT_Label_2000_Subtype_Composition.combine.Rds";
Sample.meta.file = "SAT_Sample.meta.data.Rds";
Select_Group = "Compare_Group2"; 
compare.group1 = "Group"; compare.group2 = "Pre_Post"; 
y.lab = "Percentage of Cells (%)";
export.folder = "Variables"
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Plotting_Analysis_Variables_SGLT2i.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Plotting_Analysis_Variables_SGLT2i_', 
                                                 stringr::str_split(file.variables, ".Rds")[[1]][1], 
                                                 '.html')))
rm(file.variables, Sample.meta.file, Select_Group, compare.group1, compare.group2, y.lab, export.folder)


file.variables = "SAT_Label_2000_Senescence_ModuleScoreAve.df.Rds";
Sample.meta.file = "SAT_Sample.meta.data.Rds";
Select_Group = "Compare_Group2"; 
compare.group1 = "Group"; compare.group2 = "Pre_Post"; 
y.lab = "Senescence Score";
export.folder = "Variables"
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Plotting_Analysis_Variables_SGLT2i.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Plotting_Analysis_Variables_SGLT2i_', 
                                                 stringr::str_split(file.variables, ".Rds")[[1]][1], 
                                                 '.html')))
rm(file.variables, Sample.meta.file, Select_Group, compare.group1, compare.group2, y.lab, export.folder)


## 3. ML in all variables -----------
export.folder = "ML.Variables";
Select_Group = "Select_Group"; 
Sample.meta.file = "SAT_Sample.meta.data.Rds";
file.variables = "SAT_Label_2000_Senescence_ModuleScoreAve.df.Rds";
file.variables2 = "SAT_Label_2000_Subtype_Composition.combine.Rds";
charac.variables = c("Gender", "Ethnicity", "Group", "Pre_Post")
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Main_08_ML_Fat_Variables.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_08_ML_Fat_Variables_', 
                                                 Select_Group, 
                                                 '.html')))
rm(export.folder, Select_Group, Sample.meta.file, file.variables, file.variables2, charac.variables )

export.folder = "ML.Variables";
Select_Group = "Compare_Group2"; 
Sample.meta.file = "SAT_Sample.meta.data.Rds";
file.variables = "SAT_Label_2000_Senescence_ModuleScoreAve.df.Rds";
file.variables2 = "SAT_Label_2000_Subtype_Composition.combine.Rds";
charac.variables = c("Gender", "Ethnicity", "Group", "Pre_Post")
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Main_08_ML_Fat_Variables.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_08_ML_Fat_Variables_', 
                                                 Select_Group, 
                                                 '.html')))
rm(export.folder, Select_Group, Sample.meta.file, file.variables, file.variables2, charac.variables)


export.folder = "ML.Variables";
Select_Group = c("Compare_Group1", "Compare_Group3"); 
Sample.meta.file = "SAT_Sample.meta.data.Rds";
file.variables = "SAT_Label_2000_Senescence_ModuleScoreAve.df.Rds";
file.variables2 = "SAT_Label_2000_Subtype_Composition.combine.Rds";
charac.variables = c("Gender", "Ethnicity", "Group", "Pre_Post")
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Main_08_ML_Fat_Variables.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_08_ML_Fat_Variables_', 
                                                 paste(Select_Group, collapse = "_"),  
                                                 '.html')))
rm(export.folder, Select_Group, Sample.meta.file, file.variables, file.variables2, charac.variables )


export.folder = "ML.Variables";
Select_Group = "BA_Group"; Select_Group.names =  c("Middle_Lean", "Middle_Overweight", "Older_Lean", "Older_Overweight") ; # Healthy aging group
analysis.group = "HealthyGroups"
Sample.meta.file = "SAT_Sample.meta.data.Rds";
file.variables = "SAT_Label_2000_Senescence_ModuleScoreAve.df.Rds";
file.variables2 = "SAT_Label_2000_Subtype_Composition.combine.Rds";
charac.variables = c("Gender", "Ethnicity", "Group", "Pre_Post")
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Main_08_ML_Fat_Variables.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_08_ML_Fat_Variables_', 
                                                 analysis.group, 
                                                 '.html')))
rm(export.folder, Select_Group, Select_Group.names, Sample.meta.file, file.variables, file.variables2, charac.variables )



export.folder = "ML.Variables";
Select_Group = "BA_Group"; Select_Group.names =  c("Middle_Lean", "Middle_Overweight", "Middle_Obese", 
                                                   "Older_PreD_Lean", "Older_PreD_Overweight", 
                                                   "Older_Lean", "Older_Overweight", "Older_PreD_Obese") ; 
analysis.group = "NoDiabetesGroups"
Sample.meta.file = "SAT_Sample.meta.data.Rds";
file.variables = "SAT_Label_2000_Senescence_ModuleScoreAve.df.Rds";
file.variables2 = "SAT_Label_2000_Subtype_Composition.combine.Rds";
charac.variables = c("Gender", "Ethnicity", "Group", "Pre_Post")
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Main_08_ML_Fat_Variables.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_08_ML_Fat_Variables_', 
                                                 analysis.group, 
                                                 '.html')))
rm(export.folder, Select_Group, Select_Group.names, Sample.meta.file, file.variables, file.variables2, charac.variables , analysis.group)



## ML in all variables Random Forest -----------
# export.folder = "Compositions"
# rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Random Forest.R"), 
#                   output_format= "html_document",
#                   output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
#                                           paste0('RandomForest_', 
#                                                  # stringr::str_split(file.type, ".Rds")[[1]][1], 
#                                                  '.html')))
# rm(export.folder)
# 
# export.folder = "Compositions"
# rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Plotting_RandomForest.R"), 
#                   output_format= "html_document",
#                   output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
#                                           paste0('Plotting_RandomForest_', 
#                                                  # stringr::str_split(file.type, ".Rds")[[1]][1], 
#                                                  '.html')))
# rm(export.folder)


## 4. Scatter plotting of selected variables -----------
export.folder = "ML.Variables";
analysis.file <- c( "ML.Variables/Select_Group_Basic.variables_Clinical.ML.RData"  )
x.factor = "Age"; y.factor = "Age"; Select_Group = "Select_Group"; 

Select_Group = "BA_Group"; Select_Group.names =  c("Middle_Lean", "Middle_Overweight", "Older_Lean", "Older_Overweight") ; # Healthy aging group
x.factor = "Age"; Compare.group = c("Age", "BMI", "A1c_Value","Gender", "Ethnicity"); # "Insulin", , "FastingGlucose"
analysis.group = "HealthyGroups"
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Plotting_Figure_Composition_Regression_Level2.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Plotting_Figure_Composition_Regression_Level2_', 
                                                 # stringr::str_split(file.type, ".Rds")[[1]][1], 
                                                 '.html')))
rm(export.folder, x.factor, Compare.group )



## 5.Regression of all variables by patient from ML ----------------------
export.folder = "Compositions"
Select_Group = "BA_Group"; Select_Group.names =  c("Middle_Lean", "Middle_Overweight", "Older_Lean", "Older_Overweight") ; # Healthy aging group
x.factor = "Age"; Compare.group = c("Age", "BMI", "A1c_Value","Gender", "Ethnicity"); # "Insulin", , "FastingGlucose"
analysis.group = "HealthyGroups"
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Plotting_Figure_Composition_Regression_Level2.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Plotting_Figure_Composition_Regression_Level2_', 
                                                 analysis.group, 
                                                 '.html')))
rm(export.folder,Select_Group, Select_Group.names, x.factor, Compare.group , analysis.group )

export.folder = "Compositions"
Select_Group = "Treatment_Group"; Select_Group.names =  c("Notreatment", "NUT_Pre", "SGLT2i_Pre") ; # All no treatment groups
x.factor = "Age"; Compare.group = c("Age", "BMI", "A1c_Value","Gender", "Ethnicity"); # "Insulin", , "FastingGlucose"
analysis.group = "SGLT2iGroups"
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Plotting_Figure_Composition_Regression_Level2.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Plotting_Figure_Composition_Regression_Level2_', 
                                                 analysis.group, 
                                                 '.html')))
rm(export.folder,Select_Group, Select_Group.names, x.factor, Compare.group , analysis.group )


# export.folder = "Compositions"
# Select_Group = "BA_Group"; Select_Group.names =  c("Older_Lean", "Older_Overweight", "Older_PreD_Obese", "Diabetes") ; # BMI group
# x.factor = "Age"; Compare.group = c("Age", "BMI", "A1c_Value","Gender", "Ethnicity"); # "Insulin", , "FastingGlucose"
# analysis.group = "BMIGroups"
# rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Plotting_Figure_Composition_Regression_Level2.R"), 
#                   output_format= "html_document",
#                   output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
#                                           paste0('Plotting_Figure_Composition_Regression_Level2_', 
#                                                  analysis.group, 
#                                                  '.html')))
# rm(export.folder, x.factor, Compare.group )


export.folder = "Compositions"
Select_Group = "BA_Group"; Select_Group.names =  c("Middle_Lean", "Middle_Overweight", "Middle_Obese", 
                                                   "Older_PreD_Lean", "Older_PreD_Overweight", 
                                                   "Older_Lean", "Older_Overweight", "Older_PreD_Obese") ; 
x.factor = "Age"; Compare.group = c("Age", "BMI", "A1c_Value","Gender", "Ethnicity"); # "Insulin", , "FastingGlucose"
analysis.group = "NoDiabetesGroups"
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Plotting_Figure_Composition_Regression_Level2.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Plotting_Figure_Composition_Regression_Level2_', 
                                                 analysis.group, 
                                                 '.html')))
rm(export.folder,Select_Group, Select_Group.names, x.factor, Compare.group , analysis.group )





# 9. Adipogenesis Trajectory analysis -------------------------------------------------------------------------------------------------
## step 1 -----
seurat.file = c("SAT_Label_2000.Rds")
file.type = "Subtype.metadata.Rds";Select_Group = "Select_Group"; # separate by group
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Main_07_Adipogenesis_Fat_1.Trajectory_Step0.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_07_Adipogenesis_Fat_1.Trajectory_Step1', 
                                                 '.html')))
rm(seurat.file, file.type, Select_Group)


# seurat.file = c("SAT_Label_2000.Rds")
# file.type = "Subtype.metadata.Rds";Select_Group = "Select_Group"; # separate by group
# rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Main_07_Adipogenesis_Fat_1.Trajectory_Step0_All.R"), 
#                   output_format= "html_document",
#                   output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
#                                           paste0('Main_07_Adipogenesis_Fat_1.Trajectory_Step0_All', 
#                                                  '.html')))
# rm(seurat.file, file.type, Select_Group)

seurat.file = c("SAT_Label_2000.Rds")
file.type = "Subtype.metadata.Rds";Select_Group = "Select_Group"; # separate by group
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Main_07_Adipogenesis_Fat_1.Trajectory_Step0_All_CPA1.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_07_Adipogenesis_Fat_1.Trajectory_Step0_All_CPA1', 
                                                 '.html')))
rm(seurat.file, file.type, Select_Group)

seurat.file = c("SAT_Label_2000.Rds")
file.type = "Subtype.metadata.Rds";Select_Group = "Select_Group"; # separate by group
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Main_07_Adipogenesis_Fat_1.Trajectory_Step1_CPA1_10000.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_07_Adipogenesis_Fat_1.Trajectory_Step1_CPA1_10000', 
                                                 '.html')))
rm(seurat.file, file.type, Select_Group)

# seurat.file = c("SAT_Label_2000.Rds")
# file.type = "Subtype.metadata.Rds";Select_Group = "Select_Group"; # separate by group
# rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Main_07_Adipogenesis_Fat_1.Trajectory_Step1_5000.R"), 
#                   output_format= "html_document",
#                   output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
#                                           paste0('Main_07_Adipogenesis_Fat_1.Trajectory_Step1_5000', 
#                                                  '.html')))
# rm(seurat.file, file.type, Select_Group)

seurat.file = c("SAT_Label_2000.Rds")
file.type = "Subtype.metadata.Rds";Select_Group = "Select_Group"; # separate by group
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Main_07_Adipogenesis_Fat_1.Trajectory_Step1_2000.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_07_Adipogenesis_Fat_1.Trajectory_Step1_2000', 
                                                 '.html')))
rm(seurat.file, file.type, Select_Group)


name = "3000";
seurat.file = c("SAT_Label_2000.Rds")
file.type = "Subtype.metadata.Rds";Select_Group = "Select_Group"; # separate by group
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Main_07_Adipogenesis_Fat_1.Trajectory_Step1.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_07_Adipogenesis_Fat_1.Trajectory_Step1_', name, 
                                                 '.html')))
rm(seurat.file, file.type, Select_Group, name)

name = "3000";
seurat.file = c("SAT_Label_2000.Rds")
file.type = "Subtype.metadata.Rds";Select_Group = "Select_Group"; # separate by group
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Main_07_Adipogenesis_Fat_1.Trajectory_Step1_AFE.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_07_Adipogenesis_Fat_1.Trajectory_Step1_AFE.', name, 
                                                 '.html')))
rm(seurat.file, file.type, Select_Group, name)

## step 2 ------
groups = c("CPA1")
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Main_07_Adipogenesis_Fat_1.Trajectory_Step2.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_07_Adipogenesis_Fat_1.Trajectory_Step2_', paste0(groups, collapse = "_"),
                                                 '.html')))
rm(groups)

groups = c("CPA1.10000")
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Main_07_Adipogenesis_Fat_1.Trajectory_Step2.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_07_Adipogenesis_Fat_1.Trajectory_Step2_', paste0(groups, collapse = "_"),
                                                 '.html')))
rm(groups)

groups = c("2000")
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Main_07_Adipogenesis_Fat_1.Trajectory_Step2.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_07_Adipogenesis_Fat_1.Trajectory_Step2_', paste0(groups, collapse = "_"),
                                                 '.html')))
rm(groups)

groups = c("3000")
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Main_07_Adipogenesis_Fat_1.Trajectory_Step2.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_07_Adipogenesis_Fat_1.Trajectory_Step2_', paste0(groups, collapse = "_"),
                                                 '.html')))
rm(groups)


groups = c("NoTreatment")
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Main_07_Adipogenesis_Fat_1.Trajectory_Step2.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_07_Adipogenesis_Fat_1.Trajectory_Step2_', paste0(groups, collapse = "_"),
                                                 '.html')))
rm(groups)

groups = c("SGLT2i")
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Main_07_Adipogenesis_Fat_1.Trajectory_Step2.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_07_Adipogenesis_Fat_1.Trajectory_Step2_', paste0(groups, collapse = "_"),
                                                 '.html')))
rm(groups)

## plotting --------------------
groups = c("CPA1");col.groups = c("Subtype"); 
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Plotting_Figures_Adipogenesis_Fat_1.Trajectory.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Plotting_Figures_Adipogenesis_Fat_1.Trajectory_', paste0(groups, collapse = "_"),
                                                 "_",col.groups,
                                                 '.html')))
rm(groups)

groups = c("CPA1");col.groups = c("Annotation"); 
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Plotting_Figures_Adipogenesis_Fat_1.Trajectory.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Plotting_Figures_Adipogenesis_Fat_1.Trajectory_', paste0(groups, collapse = "_"),
                                                 "_",col.groups,
                                                 '.html')))
rm(groups)

groups = c("CPA1.10000");col.groups = c("Subtype"); 
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Plotting_Figures_Adipogenesis_Fat_1.Trajectory.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Plotting_Figures_Adipogenesis_Fat_1.Trajectory_', paste0(groups, collapse = "_"),
                                                 "_",col.groups,
                                                 '.html')))
rm(groups)
groups = c("CPA1.10000");col.groups = c("Annotation"); 
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Plotting_Figures_Adipogenesis_Fat_1.Trajectory.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Plotting_Figures_Adipogenesis_Fat_1.Trajectory_', paste0(groups, collapse = "_"),
                                                 "_",col.groups,
                                                 '.html')))
rm(groups)


groups = c("2000");col.groups = c("Subtype"); 
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Plotting_Figures_Adipogenesis_Fat_1.Trajectory.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Plotting_Figures_Adipogenesis_Fat_1.Trajectory_', paste0(groups, collapse = "_"),
                                                 "_",col.groups,
                                                 '.html')))
rm(groups)
groups = c("2000");col.groups = c("Annotation"); 
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Plotting_Figures_Adipogenesis_Fat_1.Trajectory.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Plotting_Figures_Adipogenesis_Fat_1.Trajectory_', paste0(groups, collapse = "_"),
                                                 "_",col.groups,
                                                 '.html')))
rm(groups)


groups = c("3000");col.groups = c("Subtype"); 
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Plotting_Figures_Adipogenesis_Fat_1.Trajectory.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Plotting_Figures_Adipogenesis_Fat_1.Trajectory_', paste0(groups, collapse = "_"),
                                                 "_",col.groups,
                                                 '.html')))
rm(groups)
groups = c("3000");col.groups = c("Annotation"); 
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Plotting_Figures_Adipogenesis_Fat_1.Trajectory.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Plotting_Figures_Adipogenesis_Fat_1.Trajectory_', paste0(groups, collapse = "_"),
                                                 "_",col.groups,
                                                 '.html')))
rm(groups)


groups = c("NoTreatment");col.groups = c("Subtype"); 
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Plotting_Figures_Adipogenesis_Fat_1.Trajectory.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Plotting_Figures_Adipogenesis_Fat_1.Trajectory_', paste0(groups, collapse = "_"),
                                                 "_",col.groups,
                                                 '.html')))
rm(groups)
groups = c("NoTreatment");col.groups = c("Annotation"); 
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Plotting_Figures_Adipogenesis_Fat_1.Trajectory.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Plotting_Figures_Adipogenesis_Fat_1.Trajectory_', paste0(groups, collapse = "_"),
                                                 "_",col.groups,
                                                 '.html')))
rm(groups)


groups = c("SGLT2i");col.groups = c("Subtype"); 
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Plotting_Figures_Adipogenesis_Fat_1.Trajectory.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Plotting_Figures_Adipogenesis_Fat_1.Trajectory_', paste0(groups, collapse = "_"),
                                                 "_",col.groups,
                                                 '.html')))
rm(groups)
groups = c("SGLT2i");col.groups = c("Annotation"); 
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Plotting_Figures_Adipogenesis_Fat_1.Trajectory.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Plotting_Figures_Adipogenesis_Fat_1.Trajectory_', paste0(groups, collapse = "_"),
                                                 "_",col.groups,
                                                 '.html')))
rm(groups)



## step 3 -------
groups = c("CPA1");Branch="Branch_8";
groups = c("NoTreatment");Branch="Branch_2";
groups = c("SGLT2i");Branch="Branch_5";
groups = c("2000");Branch="Branch_3";
groups = c("3000");Branch="Branch_3";
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Main_07_Adipogenesis_Fat_1.Trajectory_Step3.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_07_Adipogenesis_Fat_1.Trajectory_Step3_', paste0(groups, collapse = "_"),
                                                 '.html')))
rm(groups)

## Analysis Markers -------
### Generate Association Markers ------
seurat.file = "Differentiation_PHATE_3000.Rds"
path = "Branch_3"
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Analysis_AssociationMarkers_Pt.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Analysis_AssociationMarkers_Pt_', 
                                                 unlist(stringr::str_split(seurat.file, ".Rds"))[1], "_",
                                                 paste0(path),
                                                 '.html')))
rm(seurat.file, path)

### Generate Pattern Markers to groups ------
seurat.file = "Differentiation_PHATE_3000.Rds";path = "Branch_3";
Select_Group = "Compare_Group2"; compare.group = "Pre_Post";
groups = c("Pre", "Post");
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Analysis_PatternTestMarkers_Pt_SGLT2i.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Analysis_PatternTestMarkers_Pt_', 
                                                 unlist(stringr::str_split(seurat.file, ".Rds"))[1], "_",
                                                 paste0(path), "_", compare.group, 
                                                 '.html')))
rm(seurat.file, path, Select_Group, compare.group, groups)

### Generate Trajectory Markers ------
seurat.file = "Differentiation_PHATE_3000.Rds";path = "Branch_3";
Select_Group = "Compare_Group2"; compare.group = "Pre_Post";
groups = c("Pre", "Post");
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Analysis_TrajectoryMarkers_Pt_SGLT2i.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Analysis_TrajectoryMarkers_Pt_SGLT2i_', 
                                                 unlist(stringr::str_split(seurat.file, ".Rds"))[1], "_",
                                                 paste0(path), "_", compare.group, 
                                                 '.html')))
rm(seurat.file, path, Select_Group, compare.group, groups)


seurat.file = "Differentiation_PHATE_3000.Rds";path = "Branch_3";
Select_Group = "Compare_Group1"; compare.group = "BA_Group"; 
groups = c("Older_Lean", "Older_Overweight", "Older_PreD_Obese", "Diabetes");
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Analysis_TrajectoryMarkers_Pt.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Analysis_TrajectoryMarkers_Pt_SGLT2i_', 
                                                 unlist(stringr::str_split(seurat.file, ".Rds"))[1], "_",
                                                 paste0(path), "_", compare.group, 
                                                 '.html')))
rm(seurat.file, path, Select_Group, compare.group, groups)


### Generate Pattern Markers to one ------
seurat.file = "Differentiation_PHATE_3000.Rds";path = "Branch_3";
Select_Group = "Compare_Group1"; compare.group = "BA_Group";
groups = c("Older_Lean", "Older_Overweight", "Older_PreD_Obese", "Diabetes");
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Analysis_PatternTestMarkers_Pt.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Analysis_PatternTestMarkers_Pt_', 
                                                 unlist(stringr::str_split(seurat.file, ".Rds"))[1], "_",
                                                 paste0(path), "_", compare.group, 
                                                 '.html')))
rm(seurat.file, path, Select_Group, compare.group, groups)

seurat.file = "Differentiation_PHATE_3000.Rds";path = "Branch_3";
Select_Group = "Compare_Group3"; compare.group = "Age_Group";
groups = c("Middle", "Older");
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Analysis_PatternTestMarkers_Pt.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Analysis_PatternTestMarkers_Pt_', 
                                                 unlist(stringr::str_split(seurat.file, ".Rds"))[1], "_",
                                                 paste0(path), "_", compare.group, 
                                                 '.html')))
rm(seurat.file, path, Select_Group, compare.group, groups)

# seurat.file = "Differentiation_PHATE_3000.Rds";path = "Branch_3";
# Select_Group = "Compare_Group2"; compare.group = "Treatment_Group";
# groups = c("NUT_Post", "NUT_Pre", "SGLT2i_Post", "SGLT2i_Pre");
# rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Analysis_PatternTestMarkers_Pt.R"), 
#                   output_format= "html_document",
#                   output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
#                                           paste0('Analysis_PatternTestMarkers_Pt_', 
#                                                  unlist(stringr::str_split(seurat.file, ".Rds"))[1], "_",
#                                                  paste0(path), "_", compare.group, 
#                                                  '.html')))
# rm(seurat.file, path, Select_Group, compare.group, groups)




### Combine Markers by group ------
markers.files = c(
  "AssocationMarkers/Differentiation_PHATE_3000_Pt.DEGs.Rds"
); 
sample.file.type = "SAT_Sample.meta.data.Rds"; # sample files
Select_Group = "Select_Group"; 
groups = "Differentiation"; row.by = "Dataset"
rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Main_07_AssociationMarkers_CombinebyGroup.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_07_AssociationMarkers_CombinebyGroup_', 
                                                 groups, "_",
                                                 Select_Group, "_", 
                                                 groups, '.html')))
rm(markers.files, sample.file.type, Select_Group, groups,row.by )

markers.files = c(
  "AssocationMarkers/Differentiation_PHATE_3000_Pt.DEGs.Rds"
); 
sample.file.type = "SAT_Sample.meta.data.Rds"; # sample files
Select_Group = "Compare_Group1"; Compare.G = "BA_Group";
groups = "Differentiation"; row.by = "Dataset"
rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Main_07_AssociationMarkers_CombinebyGroup.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_07_AssociationMarkers_CombinebyGroup_', 
                                                 groups, "_",
                                                 Select_Group, "_", 
                                                 groups, '.html')))
rm(markers.files, sample.file.type, Select_Group, groups, Compare.G,row.by )


markers.files = c(
  "AssocationMarkers/Differentiation_PHATE_3000_Pt.DEGs.Rds"
); 
sample.file.type = "SAT_Sample.meta.data.Rds"; # sample files
Select_Group = "Compare_Group2"; Compare.G = "Treatment_Group";
groups = "Differentiation"; row.by = "Dataset"
rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Main_07_AssociationMarkers_CombinebyGroup.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_07_AssociationMarkers_CombinebyGroup_', 
                                                 groups, "_",
                                                 Select_Group, "_", 
                                                 groups, '.html')))
rm(markers.files, sample.file.type, Select_Group, groups,Compare.G  )

markers.files = c(
  "AssocationMarkers/Differentiation_PHATE_3000_Pt.DEGs.Rds"
); 
sample.file.type = "SAT_Sample.meta.data.Rds"; # sample files
Select_Group = "Compare_Group3"; Compare.G = "Age_Group";
groups = "Differentiation"; row.by = "Dataset"
rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Main_07_AssociationMarkers_CombinebyGroup.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_07_AssociationMarkers_CombinebyGroup_', 
                                                 groups, "_",
                                                 Select_Group, "_", 
                                                 groups, '.html')))
rm(markers.files, sample.file.type, Select_Group, groups,Compare.G,row.by  )



markers.files = c(
  "PatternMarkers_SGLT2i/Differentiation_PHATE_3000_Compare_Group2_Pt.DEGs.Rds"
); 
sample.file.type = "SAT_Sample.meta.data.Rds"; # sample files
Select_Group = "Compare_Group2"; Compare.G = "Group";
groups = "Differentiation"; row.by = "patient_id"
rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Main_07_AssociationMarkers_CombinebyGroup.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_07_AssociationMarkers_CombinebyGroup_', 
                                                 groups, "_",
                                                 Select_Group, "_", 
                                                 groups, '.html')))
rm(markers.files, sample.file.type, Select_Group, groups,Compare.G ,row.by )



### Mean & Plotting -----------------
markers.files = c(
  "AssocationMarkers/Differentiation_PHATE_3000_Pt.Compare_Group3_Age_Group_Combined.RData"
  , "AssocationMarkers/Differentiation_PHATE_3000_Pt.Compare_Group2_Treatment_Group_Combined.RData"
  , "AssocationMarkers/Differentiation_PHATE_3000_Pt.Compare_Group1_BA_Group_Combined.RData"
  , "AssocationMarkers/Differentiation_PHATE_3000_Pt.Select_Group_Combined.RData"
); export.folder <- paste0("AssociationMarkers");
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Mean_Plotting_AssociationMarkers_Dataset.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Mean_Plotting_AssociationMarkers_Dataset.', 
                                                  "Annotation",
                                                 '.html')))
rm(markers.files)

markers.files = c(
  "AssocationMarkers/Differentiation_PHATE_3000_Pt.Select_Group_Age_Group_Combined.RData"
); export.folder <- paste0("AssociationMarkers");
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Mean_Plotting_AssociationMarkers_Dataset.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Mean_Plotting_AssociationMarkers_Dataset.', 
                                                 "Annotation",
                                                 '.html')))
rm(markers.files)

markers.files = c(
  "PatternMarkers_SGLT2i/Differentiation_PHATE_3000_Pt.Compare_Group2_Group_Combined.RData"
); export.folder <- paste0("PatternMarkers_SGLT2i");
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Mean_Plotting_AssociationMarkers_Dataset.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Mean_Plotting_AssociationMarkers_Dataset.', 
                                                 "Annotation",
                                                 '.html')))
rm(markers.files)

### Mean & scatter plotting Markers by Group ---------
markers.files = c(
  "AssocationMarkers/Differentiation_PHATE_3000_Pt.Compare_Group1_BA_Group_Combined.RData"
);
Select_Group = "Compare_Group1"; export.folder <- paste0("AssociationMarkers");
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Mean_Plotting_Markers_level2_1ANOVA.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder),
                                          paste0('Mean_Plotting_Markers_level2_sample_id_',
                                                 Select_Group, "_", "Subtype",
                                                 '.html')))
rm(markers.files, Select_Group)

# 11. Stats Analysis -------------------------------------------------------
compare.list <- list(Group = list(compare.group="Group", compare.group.names= list("C", "D")), # 2 way anova factor 1
                     `Event Name` = list(compare.group="Event Name", compare.group.names= list("Visit 5", "Visit 11")),# 2 way anova factor 2
                     Treatment_Group = list(compare.group="Treatment_Group", 
                                            compare.group.names= list("C_Visit 5", "D_Visit 5","C_Visit 11", "D_Visit 11"))
) 
rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Analysis_Composition_Stats.rmd"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Analysis_Composition_Stats', '.html')))
rm(compare.list)


compare.list <- list( `Event Name` = list(compare.group="Event Name", 
                                          compare.group.names= list("Visit 5", "Visit 11"))# 2 way anova factor 2
) 
rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Analysis_Composition_Paired.rmd"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Analysis_Composition_Paired', '.html')))
rm(compare.list)


compare.list <- list( `Event Name` = list(compare.group="Event Name", 
                                          compare.group.names= list("Visit 5", "Visit 11"),
                                          sub.group = "Group",
                                          sub.group.names= list("C", "D") )# 2 way anova factor 2
) 
rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Analysis_Composition_Paired_subgroup.rmd"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Analysis_Composition_Paired_subgroup', '.html')))
rm(compare.list)


# 12. Cell Markers Analysis and Plotting ---------------------------------------------------------------------------------------------------
# row numbers in metadata.Rds should be same with the whole seurat.data

## 1. level 2 markers: FindMarkers in each sample -------------
seurat.files = c("SAT_Label_2000.Rds")
sample.file.type = "SAT_Sample.meta.data.Rds"; # sample files
# Select_Group = "Compare_Group1";
file.type = "Subtype.metadata.Rds"; 
meta.group = "Annotation"; 
meta.group2.1 = "sample_id"; # subset seurat.data, 
rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Analysis_FindMarkers_seurat_level2.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Analysis_FindMarkers_', stringr::str_split(seurat.files, ".Rds")[[1]][1],"_", 
                                                 paste0(meta.group), '.html')))
rm(seurat.files, sample.file.type, meta.group, file.type, meta.group2.1)


seurat.files = c(
    "SAT_Label_2000_Immune.Rds"
  , "SAT_Label_2000_FAP.Rds"
  , "SAT_Label_2000_Endothelial.Rds"
  , "SAT_Label_2000_Adipocytes.Rds"
)
sample.file.type = "SAT_Sample.meta.data.Rds"; # sample files
# Select_Group = "Compare_Group1";
file.type = "Subtype.metadata.Rds"; 
meta.group = "Subtype"; 
meta.group2.1 = "sample_id"; # subset seurat.data, 
rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Analysis_FindMarkers_seurat_level2.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Analysis_FindMarkers_', "SAT_Label_2000", "_", 
                                                 paste0(meta.group), '.html')))
rm(seurat.files, sample.file.type, meta.group, file.type, meta.group2.1)


## level 2 markers paired -------------
# seurat.files = c("SAT_Label_2000.Rds")
# sample.file.type = "SAT_Sample.meta.data.Rds"; # sample files
# Select_Group = "Compare_Group2";
# file.type = "Subtype.metadata.Rds"; 
# meta.group = "Annotation"; 
# meta.group2.1 = "patient_id"; # subset seurat.data, 
# Compare.G = "Pre_Post";Compare.levels = c("Pre", "Post")
# rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Analysis_FindMarkers_level2_Paired.R"), output_format= "html_document",
#                   output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
#                                           paste0('Analysis_FindMarkers_level2_Paired_', stringr::str_split(seurat.files, ".Rds")[[1]][1],"_", 
#                                                  paste0(meta.group), '.html')))
# rm(seurat.files, sample.file.type, meta.group, file.type, meta.group2.1)



## 2. combine level 2 markers by Group -------------
### cell type markers
markers.files = c(
  "FindMarkers_SAT_Label_2000/Annotation.Level2.Markers/SAT_Label_2000_Annotation.sample_id_Markers.Rds"
); 
sample.file.type = "SAT_Sample.meta.data.Rds"; # sample files
Select_Group = "Select_Group"; 
meta.group = "Annotation"; 
rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Main_07_MarkersLevel2_CombinebyGroup.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_07_MarkersLevel2_CombinebyGroup_', 
                                                 Select_Group, "_", meta.group, '.html')))
rm(markers.files, sample.file.type, meta.group, Select_Group)
markers.files = c(
     "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Adipocytes_Subtype.sample_id_Markers.Rds"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Endothelial_Subtype.sample_id_Markers.Rds"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_FAP_Subtype.sample_id_Markers.Rds"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Immune_Subtype.sample_id_Markers.Rds"
); 
sample.file.type = "SAT_Sample.meta.data.Rds"; # sample files
Select_Group = "Select_Group"; 
meta.group = "Subtype"; 
rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Main_07_MarkersLevel2_CombinebyGroup.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_07_MarkersLevel2_CombinebyGroup_', 
                                                 Select_Group, "_", meta.group, '.html')))
rm(markers.files, sample.file.type, meta.group, Select_Group)


markers.files = c(
  "FindMarkers_SAT_Label_2000/Annotation.Level2.Markers/SAT_Label_2000_Annotation.sample_id_Markers.Rds"
); 
sample.file.type = "SAT_Sample.meta.data.Rds"; # sample files
Select_Group = "Compare_Group1"; 
meta.group = "Annotation"; 
rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Main_07_MarkersLevel2_CombinebyGroup.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_07_MarkersLevel2_CombinebyGroup_', 
                                                 Select_Group, "_", meta.group, '.html')))
rm(markers.files, sample.file.type, meta.group, Select_Group)
markers.files = c(
  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Adipocytes_Subtype.sample_id_Markers.Rds"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Endothelial_Subtype.sample_id_Markers.Rds"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_FAP_Subtype.sample_id_Markers.Rds"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Immune_Subtype.sample_id_Markers.Rds"
); 
sample.file.type = "SAT_Sample.meta.data.Rds"; # sample files
Select_Group = "Compare_Group1"; 
meta.group = "Subtype"; 
rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Main_07_MarkersLevel2_CombinebyGroup.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_07_MarkersLevel2_CombinebyGroup_', 
                                                 Select_Group, "_", meta.group, '.html')))
rm(markers.files, sample.file.type, meta.group, Select_Group)


markers.files = c(
  "FindMarkers_SAT_Label_2000/Annotation.Level2.Markers/SAT_Label_2000_Annotation.sample_id_Markers.Rds"
); 
sample.file.type = "SAT_Sample.meta.data.Rds"; # sample files
Select_Group = "Compare_Group3"; 
meta.group = "Annotation"; 
rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Main_07_MarkersLevel2_CombinebyGroup.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_07_MarkersLevel2_CombinebyGroup_', 
                                                 Select_Group, "_", meta.group, '.html')))
rm(markers.files, sample.file.type, meta.group, Select_Group)
markers.files = c(
  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Adipocytes_Subtype.sample_id_Markers.Rds"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Endothelial_Subtype.sample_id_Markers.Rds"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_FAP_Subtype.sample_id_Markers.Rds"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Immune_Subtype.sample_id_Markers.Rds"
); 
sample.file.type = "SAT_Sample.meta.data.Rds"; # sample files
Select_Group = "Compare_Group3"; 
meta.group = "Subtype"; 
rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Main_07_MarkersLevel2_CombinebyGroup.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_07_MarkersLevel2_CombinebyGroup_', 
                                                 Select_Group, "_", meta.group, '.html')))
rm(markers.files, sample.file.type, meta.group, Select_Group)


## 3. Mean & scatter plotting for level2 Cell Type Markers ---------
markers.files = c(
  "FindMarkers_SAT_Label_2000/Annotation.Level2.Markers/SAT_Label_2000_Annotation.Select_Group_Combined.RData"
); 
Select_Group = "Select_Group"; 
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Mean_Plotting_Markers_level2_sample_id.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Mean_Plotting_Markers_level2_sample_id_', 
                                                 Select_Group, "_", "Annotation",
                                                 '.html')))
rm(markers.files, Select_Group)
markers.files = c(
  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Adipocytes_Subtype.Select_Group_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Endothelial_Subtype.Select_Group_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_FAP_Subtype.Select_Group_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Immune_Subtype.Select_Group_Combined.RData"
); 
Select_Group = "Select_Group"
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Mean_Plotting_Markers_level2_sample_id.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Mean_Plotting_Markers_level2_sample_id_', 
                                                 Select_Group, "_", "Subtype",
                                                 '.html')))
rm(markers.files, Select_Group)


markers.files = c(
  "FindMarkers_SAT_Label_2000/Annotation.Level2.Markers/SAT_Label_2000_Annotation.Compare_Group1_Combined.RData"
); 
Select_Group = "Compare_Group1"; 
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Mean_Plotting_Markers_level2_sample_id.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Mean_Plotting_Markers_level2_sample_id_', 
                                                 Select_Group, "_", "Annotation",
                                                 '.html')))
rm(markers.files, Select_Group)
markers.files = c(
  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Adipocytes_Subtype.Compare_Group1_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Endothelial_Subtype.Compare_Group1_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_FAP_Subtype.Compare_Group1_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Immune_Subtype.Compare_Group1_Combined.RData"
);
Select_Group = "Compare_Group1"
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Mean_Plotting_Markers_level2_sample_id.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder),
                                          paste0('Mean_Plotting_Markers_level2_sample_id_',
                                                 Select_Group, "_", "Subtype",
                                                 '.html')))
rm(markers.files, Select_Group)


markers.files = c(
  "FindMarkers_SAT_Label_2000/Annotation.Level2.Markers/SAT_Label_2000_Annotation.Compare_Group3_Combined.RData"
); 
Select_Group = "Compare_Group3"; 
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Mean_Plotting_Markers_level2_sample_id.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Mean_Plotting_Markers_level2_sample_id_', 
                                                 Select_Group, "_", "Annotation",
                                                 '.html')))
rm(markers.files, Select_Group)
markers.files = c(
  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Adipocytes_Subtype.Compare_Group3_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Endothelial_Subtype.Compare_Group3_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_FAP_Subtype.Compare_Group3_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Immune_Subtype.Compare_Group3_Combined.RData"
);
Select_Group = "Compare_Group3"
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Mean_Plotting_Markers_level2_sample_id.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder),
                                          paste0('Mean_Plotting_Markers_level2_sample_id_',
                                                 Select_Group, "_", "Subtype",
                                                 '.html')))
rm(markers.files, Select_Group)


## 4. VlnPlot plotting for level2 Cell Type Markers ---------
###  "Annotation" ---------
seurat.files =  c("SAT_Label_2000.Rds"); file.type = "Subtype.metadata.Rds"; 
sample.file.type = "SAT_Sample.meta.data.Rds"; # sample files
Select_Group = "Select_Group"; group.x = "Annotation";
markers.files = c(
  "FindMarkers_SAT_Label_2000/Annotation.Level2.Markers/SAT_Label_2000_Annotation.Select_Group_Combined.Mean.Rds"
); 
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Plotting_Markers_level2_VlnPlot.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Plotting_level2_Markers_VlnPlot_', 
                                                 unlist(stringr::str_split(seurat.files, ".Rds"))[1], "_",
                                                 Select_Group, "_", group.x,
                                                 '.html')))
rm(seurat.files , file.type, sample.file.type, markers.files, Select_Group, group.x)
seurat.files =  c("SAT_Label_2000.Rds"); file.type = "Subtype.metadata.Rds"; 
sample.file.type = "SAT_Sample.meta.data.Rds"; # sample files
Select_Group = "Compare_Group1"; group.x = "Annotation";
markers.files = c(
  "FindMarkers_SAT_Label_2000/Annotation.Level2.Markers/SAT_Label_2000_Annotation.Compare_Group1_Combined.Mean.Rds"
); 
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Plotting_Markers_level2_VlnPlot.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Plotting_level2_Markers_VlnPlot_', 
                                                 unlist(stringr::str_split(seurat.files, ".Rds"))[1], "_",
                                                 Select_Group, "_", group.x,
                                                 '.html')))
rm(seurat.files , file.type, sample.file.type, markers.files, Select_Group, group.x)
seurat.files =  c("SAT_Label_2000.Rds"); file.type = "Subtype.metadata.Rds"; 
sample.file.type = "SAT_Sample.meta.data.Rds"; # sample files
Select_Group = "Compare_Group3"; group.x = "Annotation";
markers.files = c(
  "FindMarkers_SAT_Label_2000/Annotation.Level2.Markers/SAT_Label_2000_Annotation.Compare_Group3_Combined.Mean.Rds"
); 
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Plotting_Markers_level2_VlnPlot.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Plotting_level2_Markers_VlnPlot_', 
                                                 unlist(stringr::str_split(seurat.files, ".Rds"))[1], "_",
                                                 Select_Group, "_", group.x,
                                                 '.html')))
rm(seurat.files , file.type, sample.file.type, markers.files, Select_Group, group.x)


###  "Subtype", "Select_Group" ---------
seurat.files =  c("SAT_Label_2000_Immune.Rds"); file.type = "Subtype.metadata.Rds"; 
sample.file.type = "SAT_Sample.meta.data.Rds"; # sample files
Select_Group = "Select_Group"; group.x = "Subtype";
markers.files = c(
  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Immune_Subtype.Select_Group_Combined.Mean.Rds"
); 
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Plotting_Markers_level2_VlnPlot.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Plotting_level2_Markers_VlnPlot_', 
                                                 unlist(stringr::str_split(seurat.files, ".Rds"))[1], "_",
                                                 Select_Group, "_", group.x,
                                                 '.html')))
rm(seurat.files , file.type, sample.file.type, markers.files, Select_Group, group.x)
seurat.files =  c("SAT_Label_2000_FAP.Rds"); file.type = "Subtype.metadata.Rds"; 
sample.file.type = "SAT_Sample.meta.data.Rds"; # sample files
Select_Group = "Select_Group"; group.x = "Subtype";
markers.files = c(
  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_FAP_Subtype.Select_Group_Combined.Mean.Rds"
); 
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Plotting_Markers_level2_VlnPlot.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Plotting_level2_Markers_VlnPlot_', 
                                                 unlist(stringr::str_split(seurat.files, ".Rds"))[1], "_",
                                                 Select_Group, "_", group.x,
                                                 '.html')))
rm(seurat.files , file.type, sample.file.type, markers.files, Select_Group, group.x)
seurat.files =  c("SAT_Label_2000_Endothelial.Rds"); file.type = "Subtype.metadata.Rds"; 
sample.file.type = "SAT_Sample.meta.data.Rds"; # sample files
Select_Group = "Select_Group"; group.x = "Subtype";
markers.files = c(
  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Endothelial_Subtype.Select_Group_Combined.Mean.Rds"
); 
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Plotting_Markers_level2_VlnPlot.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Plotting_level2_Markers_VlnPlot_', 
                                                 unlist(stringr::str_split(seurat.files, ".Rds"))[1], "_",
                                                 Select_Group, "_", group.x,
                                                 '.html')))
rm(seurat.files , file.type, sample.file.type, markers.files, Select_Group, group.x)
seurat.files =  c("SAT_Label_2000_Adipocytes.Rds"); file.type = "Subtype.metadata.Rds"; 
sample.file.type = "SAT_Sample.meta.data.Rds"; # sample files
Select_Group = "Select_Group"; group.x = "Subtype";
markers.files = c(
  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Adipocytes_Subtype.Select_Group_Combined.Mean.Rds"
); 
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Plotting_Markers_level2_VlnPlot.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Plotting_level2_Markers_VlnPlot_', 
                                                 unlist(stringr::str_split(seurat.files, ".Rds"))[1], "_",
                                                 Select_Group, "_", group.x,
                                                 '.html')))
rm(seurat.files , file.type, sample.file.type, markers.files, Select_Group, group.x)

###  "Subtype", "Compare_Group1" ---------
seurat.files =  c("SAT_Label_2000_Immune.Rds"); file.type = "Subtype.metadata.Rds"; 
sample.file.type = "SAT_Sample.meta.data.Rds"; # sample files
Select_Group = "Compare_Group1"; group.x = "Subtype";
markers.files = c(
  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Immune_Subtype.Compare_Group1_Combined.Mean.Rds"
); 
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Plotting_Markers_level2_VlnPlot.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Plotting_level2_Markers_VlnPlot_', 
                                                 unlist(stringr::str_split(seurat.files, ".Rds"))[1], "_",
                                                 Select_Group, "_", group.x,
                                                 '.html')))
rm(seurat.files , file.type, sample.file.type, markers.files, Select_Group, group.x)
seurat.files =  c("SAT_Label_2000_FAP.Rds"); file.type = "Subtype.metadata.Rds"; 
sample.file.type = "SAT_Sample.meta.data.Rds"; # sample files
Select_Group = "Compare_Group1"; group.x = "Subtype";
markers.files = c(
  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_FAP_Subtype.Compare_Group1_Combined.Mean.Rds"
); 
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Plotting_Markers_level2_VlnPlot.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Plotting_level2_Markers_VlnPlot_', 
                                                 unlist(stringr::str_split(seurat.files, ".Rds"))[1], "_",
                                                 Select_Group, "_", group.x,
                                                 '.html')))
rm(seurat.files , file.type, sample.file.type, markers.files, Select_Group, group.x)
seurat.files =  c("SAT_Label_2000_Endothelial.Rds"); file.type = "Subtype.metadata.Rds"; 
sample.file.type = "SAT_Sample.meta.data.Rds"; # sample files
Select_Group = "Compare_Group1"; group.x = "Subtype";
markers.files = c(
  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Endothelial_Subtype.Compare_Group1_Combined.Mean.Rds"
); 
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Plotting_Markers_level2_VlnPlot.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Plotting_level2_Markers_VlnPlot_', 
                                                 unlist(stringr::str_split(seurat.files, ".Rds"))[1], "_",
                                                 Select_Group, "_", group.x,
                                                 '.html')))
rm(seurat.files , file.type, sample.file.type, markers.files, Select_Group, group.x)
seurat.files =  c("SAT_Label_2000_Adipocytes.Rds"); file.type = "Subtype.metadata.Rds"; 
sample.file.type = "SAT_Sample.meta.data.Rds"; # sample files
Select_Group = "Compare_Group1"; group.x = "Subtype";
markers.files = c(
  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Adipocytes_Subtype.Compare_Group1_Combined.Mean.Rds"
); 
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Plotting_Markers_level2_VlnPlot.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Plotting_level2_Markers_VlnPlot_', 
                                                 unlist(stringr::str_split(seurat.files, ".Rds"))[1], "_",
                                                 Select_Group, "_", group.x,
                                                 '.html')))
rm(seurat.files , file.type, sample.file.type, markers.files, Select_Group, group.x)


###  "Subtype", "Compare_Group3" ---------
seurat.files =  c("SAT_Label_2000_Immune.Rds"); file.type = "Subtype.metadata.Rds"; 
sample.file.type = "SAT_Sample.meta.data.Rds"; # sample files
Select_Group = "Compare_Group3"; group.x = "Subtype";
markers.files = c(
  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Immune_Subtype.Compare_Group3_Combined.Mean.Rds"
); 
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Plotting_Markers_level2_VlnPlot.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Plotting_level2_Markers_VlnPlot_', 
                                                 unlist(stringr::str_split(seurat.files, ".Rds"))[1], "_",
                                                 Select_Group, "_", group.x,
                                                 '.html')))
rm(seurat.files , file.type, sample.file.type, markers.files, Select_Group, group.x)
seurat.files =  c("SAT_Label_2000_FAP.Rds"); file.type = "Subtype.metadata.Rds"; 
sample.file.type = "SAT_Sample.meta.data.Rds"; # sample files
Select_Group = "Compare_Group3"; group.x = "Subtype";
markers.files = c(
  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_FAP_Subtype.Compare_Group3_Combined.Mean.Rds"
); 
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Plotting_Markers_level2_VlnPlot.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Plotting_level2_Markers_VlnPlot_', 
                                                 unlist(stringr::str_split(seurat.files, ".Rds"))[1], "_",
                                                 Select_Group, "_", group.x,
                                                 '.html')))
rm(seurat.files , file.type, sample.file.type, markers.files, Select_Group, group.x)
seurat.files =  c("SAT_Label_2000_Endothelial.Rds"); file.type = "Subtype.metadata.Rds"; 
sample.file.type = "SAT_Sample.meta.data.Rds"; # sample files
Select_Group = "Compare_Group3"; group.x = "Subtype";
markers.files = c(
  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Endothelial_Subtype.Compare_Group3_Combined.Mean.Rds"
); 
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Plotting_Markers_level2_VlnPlot.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Plotting_level2_Markers_VlnPlot_', 
                                                 unlist(stringr::str_split(seurat.files, ".Rds"))[1], "_",
                                                 Select_Group, "_", group.x,
                                                 '.html')))
rm(seurat.files , file.type, sample.file.type, markers.files, Select_Group, group.x)
seurat.files =  c("SAT_Label_2000_Adipocytes.Rds"); file.type = "Subtype.metadata.Rds"; 
sample.file.type = "SAT_Sample.meta.data.Rds"; # sample files
Select_Group = "Compare_Group3"; group.x = "Subtype";
markers.files = c(
  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Adipocytes_Subtype.Compare_Group3_Combined.Mean.Rds"
); 
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Plotting_Markers_level2_VlnPlot.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Plotting_level2_Markers_VlnPlot_', 
                                                 unlist(stringr::str_split(seurat.files, ".Rds"))[1], "_",
                                                 Select_Group, "_", group.x,
                                                 '.html')))
rm(seurat.files , file.type, sample.file.type, markers.files, Select_Group, group.x)



## 5. Markers to GProfiler Pathway Analysis and Plotting  ---------
### Select_Group
setwd("/media/jianie/DATA")
markers.files = c(
  "FindMarkers_SAT_Label_2000/Annotation.Level2.Markers/SAT_Label_2000_Annotation.Select_Group_Combined.Mean.Rds"
  ,"FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Endothelial_Subtype.Select_Group_Combined.Mean.Rds"
  ,"FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Adipocytes_Subtype.Select_Group_Combined.Mean.Rds"
  ,"FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_FAP_Subtype.Select_Group_Combined.Mean.Rds"
  ,"FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Immune_Subtype.Select_Group_Combined.Mean.Rds"
); Select_Group = "Select_Group";organism = "hsapiens";
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "CellMarkerstoPathway_Analysis_Plotting_GProfiler.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('CellMarkers_GProfiler_', 
                                                 "SAT_Label_2000", "_", Select_Group, 
                                                 '.html')))
rm(markers.files, organism, Select_Group)


### Compare_Group1
markers.files = c(
  "FindMarkers_SAT_Label_2000/Annotation.Level2.Markers/SAT_Label_2000_Annotation.Compare_Group1_Combined.Mean.Rds"
  ,"FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Endothelial_Subtype.Compare_Group1_Combined.Mean.Rds"
  ,"FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Adipocytes_Subtype.Compare_Group1_Combined.Mean.Rds"
  ,"FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_FAP_Subtype.Compare_Group1_Combined.Mean.Rds"
  ,"FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Immune_Subtype.Compare_Group1_Combined.Mean.Rds"
); Select_Group = "Compare_Group1";organism = "hsapiens";
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "CellMarkerstoPathway_Analysis_Plotting_GProfiler.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('CellMarkers_GProfiler_', 
                                                 "SAT_Label_2000", "_", Select_Group, 
                                                 '.html')))
rm(markers.files, organism, Select_Group)

### Compare_Group3
markers.files = c(
  "FindMarkers_SAT_Label_2000/Annotation.Level2.Markers/SAT_Label_2000_Annotation.Compare_Group3_Combined.Mean.Rds"
  ,"FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Endothelial_Subtype.Compare_Group3_Combined.Mean.Rds"
  ,"FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Adipocytes_Subtype.Compare_Group3_Combined.Mean.Rds"
  ,"FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_FAP_Subtype.Compare_Group3_Combined.Mean.Rds"
  ,"FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Immune_Subtype.Compare_Group3_Combined.Mean.Rds"
); Select_Group = "Compare_Group3";organism = "hsapiens";
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "CellMarkerstoPathway_Analysis_Plotting_GProfiler.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('CellMarkers_GProfiler_', 
                                                 "SAT_Label_2000", "_", Select_Group, 
                                                 '.html')))
rm(markers.files, organism, Select_Group)


## 6. Markers to KEGG Pathway Analysis and Plotting  ---------
markers.files = c(
  "FindMarkers_SAT_Label_2000/Annotation.Level2.Markers/SAT_Label_2000_Annotation.Select_Group_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Adipocytes_Subtype.Select_Group_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Endothelial_Subtype.Select_Group_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_FAP_Subtype.Select_Group_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Immune_Subtype.Select_Group_Combined.RData"
); 
Select_Group = "Select_Group"; type = "avg_log2FC";
rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Main_07_MarkersLevel2_KEGGPathway_Analysis_as.group.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_07_MarkersLevel2_KEGGPathway_', 
                                                 "SAT_Label_2000", "_", Select_Group, 
                                                 '.html')))
rm(markers.files, type, Select_Group)


markers.files = c(
  "FindMarkers_SAT_Label_2000/Annotation.Level2.Markers/SAT_Label_2000_Annotation.Compare_Group1_Combined.RData"
,   "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Adipocytes_Subtype.Compare_Group1_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Endothelial_Subtype.Compare_Group1_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_FAP_Subtype.Compare_Group1_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Immune_Subtype.Compare_Group1_Combined.RData"
  
); 
Select_Group = "Compare_Group1"; type = "avg_log2FC";
rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Main_07_MarkersLevel2_KEGGPathway_Analysis_as.group.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_07_MarkersLevel2_KEGGPathway_', 
                                                 "SAT_Label_2000", "_", Select_Group, 
                                                 '.html')))
rm(markers.files,  type, Select_Group)

markers.files = c(
  "FindMarkers_SAT_Label_2000/Annotation.Level2.Markers/SAT_Label_2000_Annotation.Compare_Group3_Combined.RData"
  ,   "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Adipocytes_Subtype.Compare_Group3_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Endothelial_Subtype.Compare_Group3_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_FAP_Subtype.Compare_Group3_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Immune_Subtype.Compare_Group3_Combined.RData"
  
); 
Select_Group = "Compare_Group3"; type = "avg_log2FC";
rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Main_07_MarkersLevel2_KEGGPathway_Analysis_as.group.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_07_MarkersLevel2_KEGGPathway_', 
                                                 "SAT_Label_2000", "_", Select_Group, 
                                                 '.html')))
rm(markers.files, type, Select_Group)



## 7. Markers to wiki Pathway Analysis and Plotting_not use  ---------
markers.files = c(
  "FindMarkers_SAT_Label_2000/sene.name.Markers/Combined_SAT_Label_2000_Immune_Senescence.name.metadata.Rds"
  #, "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_FAP_Senescence.name.metadata.Rds"
  #, "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Adipocytes_Senescence.name.metadata.Rds"
  #, "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Endothelial_Senescence.name.metadata.Rds"
); 
organism = "Homo sapiens" # organism = "Mus musculus"
rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Main_07_MarkerstoPathway_Analysis.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Analysis_MarkerstoPathway_', "SAT_Label_2000_Immune",  '.html')))
rm(markers.files, organism)

pathways.files = c(
  "Pathways/Pathway_Combined_SAT_Label_2000_Immune_Senescence.name.metadata.Rds"
)
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Plotting_Figure_Analysis_Pathways.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Plotting_Figure_Analysis_pathways',  '.html')))
rm(pathways.files)



# 13. Assigned Markers plotting -----------------------------------------------
## 1. Heatmap plotting for Assigned Markers  ------------------------
seurat.file = c("SAT_Label_2000.Rds");
Markers_file =c("Fat_Subtype.xlsx"); sheet = "Annotation_Markers";
sample.file = "SAT_Sample.meta.data.Rds";
file.type = "Subtype.metadata.Rds"; 
Select_Group = "Select_Group"; meta.group1 = c("Annotation"); meta.group2 = c("sample_id")
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Plotting_Figure_Markers_GeneAverageHeatmap_level2.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Plotting_Figure_Markers_GeneAverageHeatmap_level2_', 
                                                 stringr::str_split(seurat.file, ".Rds")[[1]][1], "_", 
                                                 meta.group1, "_", meta.group2,  '.html')))
rm(Markers.files, seurat.file, sample.file, file.type, Select_Group, meta.group1, meta.group2)

## 2. VlnPlot plotting for Assigned Markers ------------------------
seurat.files =  c("SAT_Label_2000.Rds"); file.type = "Subtype.metadata.Rds"; 
sample.file.type = "SAT_Sample.meta.data.Rds"; # sample files
Select_Group = "Select_Group"; group.x = "Annotation";
Markers_file =c("Fat_Subtype.xlsx"); sheet = "Annotation_Markers";
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Plotting_Markers_VlnPlot.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Plotting_Markers_VlnPlot_', 
                                                 Select_Group, "_", "Annotation",
                                                 '.html')))
rm(seurat.files , file.type, sample.file.type, Markers_file, Select_Group, group.x, sheet)


# 14. DEG analysis between treatment and plotting -----------------------------------------------
## 2. combine level 2 markers by Group -------------
markers.files = c(
  "FindMarkers_SAT_Label_2000/Annotation.Level2.Markers/SAT_Label_2000_Annotation.sample_id_Markers.Rds"
); 
sample.file.type = "SAT_Sample.meta.data.Rds"; # sample files
Select_Group = "Compare_Group1"; Compare.G = "BA_Group";
meta.group = "Annotation"; 
rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Main_07_MarkersLevel2_CombinebyGroup.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_07_MarkersLevel2_CombinebyGroup_', 
                                                 Select_Group, "_", meta.group, "_", Compare.G, '.html')))
rm(markers.files, sample.file.type, meta.group, Select_Group, Compare.G)

markers.files = c(
     "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Adipocytes_Subtype.sample_id_Markers.Rds"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Endothelial_Subtype.sample_id_Markers.Rds"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_FAP_Subtype.sample_id_Markers.Rds"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Immune_Subtype.sample_id_Markers.Rds"
); 
sample.file.type = "SAT_Sample.meta.data.Rds"; # sample files
Select_Group = "Compare_Group1"; Compare.G = "BA_Group";
meta.group = "Subtype"; 
rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Main_07_MarkersLevel2_CombinebyGroup.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_07_MarkersLevel2_CombinebyGroup_', 
                                                 Select_Group, "_", meta.group, "_", Compare.G,  '.html')))
rm(markers.files, sample.file.type, meta.group, Select_Group, Compare.G)


markers.files = c(
  "FindMarkers_SAT_Label_2000/Annotation.Level2.Markers/SAT_Label_2000_Annotation.sample_id_Markers.Rds"
); 
sample.file.type = "SAT_Sample.meta.data.Rds"; # sample files
Select_Group = "Compare_Group2"; Compare.G = "Treatment_Group";
meta.group = "Annotation"; 
rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Main_07_MarkersLevel2_CombinebyGroup.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_07_MarkersLevel2_CombinebyGroup_', 
                                                 Select_Group, "_", meta.group,"_", Compare.G, '.html')))
rm(markers.files, sample.file.type, meta.group, Select_Group, Compare.G)

markers.files = c(
  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Adipocytes_Subtype.sample_id_Markers.Rds"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Endothelial_Subtype.sample_id_Markers.Rds"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_FAP_Subtype.sample_id_Markers.Rds"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Immune_Subtype.sample_id_Markers.Rds"
); 
sample.file.type = "SAT_Sample.meta.data.Rds"; # sample files
Select_Group = "Compare_Group2"; Compare.G = "Treatment_Group";
meta.group = "Subtype"; 
rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Main_07_MarkersLevel2_CombinebyGroup.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_07_MarkersLevel2_CombinebyGroup_', 
                                                 Select_Group, "_", meta.group, "_", Compare.G, '.html')))
rm(markers.files, sample.file.type, meta.group, Select_Group, Compare.G)



markers.files = c(
  "FindMarkers_SAT_Label_2000/Annotation.Level2.Markers/SAT_Label_2000_Annotation.sample_id_Markers.Rds"
); 
sample.file.type = "SAT_Sample.meta.data.Rds"; # sample files
Select_Group = "Compare_Group2"; Compare.G = "Group";
meta.group = "Annotation"; 
rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Main_07_MarkersLevel2_CombinebyGroup.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_07_MarkersLevel2_CombinebyGroup_', 
                                                 Select_Group, "_", meta.group,"_", Compare.G, '.html')))
rm(markers.files, sample.file.type, meta.group, Select_Group, Compare.G)

markers.files = c(
  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Adipocytes_Subtype.sample_id_Markers.Rds"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Endothelial_Subtype.sample_id_Markers.Rds"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_FAP_Subtype.sample_id_Markers.Rds"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Immune_Subtype.sample_id_Markers.Rds"
); 
sample.file.type = "SAT_Sample.meta.data.Rds"; # sample files
Select_Group = "Compare_Group2"; Compare.G = "Group";
meta.group = "Subtype"; 
rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Main_07_MarkersLevel2_CombinebyGroup.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_07_MarkersLevel2_CombinebyGroup_', 
                                                 Select_Group, "_", meta.group, "_", Compare.G, '.html')))
rm(markers.files, sample.file.type, meta.group, Select_Group, Compare.G)



markers.files = c(
  "FindMarkers_SAT_Label_2000/Annotation.Level2.Markers/SAT_Label_2000_Annotation.sample_id_Markers.Rds"
); 
sample.file.type = "SAT_Sample.meta.data.Rds"; # sample files
Select_Group = "Compare_Group3"; Compare.G = "Age_Group";
meta.group = "Annotation"; 
rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Main_07_MarkersLevel2_CombinebyGroup.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_07_MarkersLevel2_CombinebyGroup_', 
                                                 Select_Group, "_", meta.group,"_", Compare.G, '.html')))
rm(markers.files, sample.file.type, meta.group, Select_Group, Compare.G)

markers.files = c(
  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Adipocytes_Subtype.sample_id_Markers.Rds"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Endothelial_Subtype.sample_id_Markers.Rds"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_FAP_Subtype.sample_id_Markers.Rds"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Immune_Subtype.sample_id_Markers.Rds"
); 
sample.file.type = "SAT_Sample.meta.data.Rds"; # sample files
Select_Group = "Compare_Group3"; Compare.G = "Age_Group";
meta.group = "Subtype"; 
rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Main_07_MarkersLevel2_CombinebyGroup.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_07_MarkersLevel2_CombinebyGroup_', 
                                                 Select_Group, "_", meta.group, "_", Compare.G, '.html')))
rm(markers.files, sample.file.type, meta.group, Select_Group, Compare.G)



## 3. Mean & scatter plotting for level2 Markers by Group ---------
markers.files = c(
  "FindMarkers_SAT_Label_2000/Annotation.Level2.Markers/SAT_Label_2000_Annotation.Compare_Group1_BA_Group_Combined.RData"
); 
Select_Group = "Compare_Group1"; 
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Mean_Plotting_Markers_level2_sample_id.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Mean_Plotting_Markers_level2_sample_id_', 
                                                 Select_Group, "_", "BA_Group", "_",  "Annotation",
                                                 '.html')))
rm(markers.files, Select_Group)
markers.files = c(
  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Adipocytes_Subtype.Compare_Group1_BA_Group_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Endothelial_Subtype.Compare_Group1_BA_Group_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_FAP_Subtype.Compare_Group1_BA_Group_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Immune_Subtype.Compare_Group1_BA_Group_Combined.RData"
);
Select_Group = "Compare_Group1"
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Mean_Plotting_Markers_level2_sample_id.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder),
                                          paste0('Mean_Plotting_Markers_level2_sample_id_',
                                                 Select_Group, "_", "BA_Group", "_", "Subtype",
                                                 '.html')))
rm(markers.files, Select_Group)


markers.files = c(
  "FindMarkers_SAT_Label_2000/Annotation.Level2.Markers/SAT_Label_2000_Annotation.Compare_Group2_Treatment_Group_Combined.RData"
); 
Select_Group = "Compare_Group2"; 
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Mean_Plotting_Markers_level2_sample_id.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Mean_Plotting_Markers_level2_sample_id_', 
                                                 Select_Group,"_", "Treatment_Group", "_", "Annotation",
                                                 '.html')))
rm(markers.files, Select_Group)
markers.files = c(
  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Adipocytes_Subtype.Compare_Group2_Treatment_Group_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Endothelial_Subtype.Compare_Group2_Treatment_Group_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_FAP_Subtype.Compare_Group2_Treatment_Group_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Immune_Subtype.Compare_Group2_Treatment_Group_Combined.RData"
);
Select_Group = "Compare_Group2"
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Mean_Plotting_Markers_level2_sample_id.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder),
                                          paste0('Mean_Plotting_Markers_level2_sample_id_',
                                                 Select_Group, "_", "Treatment_Group", "_", "Subtype",
                                                 '.html')))
rm(markers.files, Select_Group)

markers.files = c(
  "FindMarkers_SAT_Label_2000/Annotation.Level2.Markers/SAT_Label_2000_Annotation.Compare_Group3_Age_Group_Combined.RData"
); 
Select_Group = "Compare_Group3"; 
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Mean_Plotting_Markers_level2_sample_id.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Mean_Plotting_Markers_level2_sample_id_', 
                                                 Select_Group,"_", "Age_Group", "_", "Annotation",
                                                 '.html')))
rm(markers.files, Select_Group)
markers.files = c(
  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Adipocytes_Subtype.Compare_Group3_Age_Group_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Endothelial_Subtype.Compare_Group3_Age_Group_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_FAP_Subtype.Compare_Group3_Age_Group_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Immune_Subtype.Compare_Group3_Age_Group_Combined.RData"
);
Select_Group = "Compare_Group3";
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Mean_Plotting_Markers_level2_sample_id.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder),
                                          paste0('Mean_Plotting_Markers_level2_sample_id_',
                                                 Select_Group, "_", "Age_Group", "_", "Subtype",
                                                 '.html')))
rm(markers.files, Select_Group)

## 4. 1ANOVA for level2 Markers by Group ---------------------------
markers.files = c(
  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Adipocytes_Subtype.Compare_Group1_BA_Group_Combined.RData"
);
Select_Group = "Compare_Group1";
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Mean_Plotting_Markers_level2_1ANOVA.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder),
                                          paste0('Mean_Plotting_Markers_level2_1ANOVA_',
                                                 Select_Group, "_", "BA_Group", "_","Subtype",
                                                 '.html')))
rm(markers.files, Select_Group)

markers.files = c(
  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Adipocytes_Subtype.Compare_Group1_BA_Group_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Endothelial_Subtype.Compare_Group1_BA_Group_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_FAP_Subtype.Compare_Group1_BA_Group_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Immune_Subtype.Compare_Group1_BA_Group_Combined.RData"
);
Select_Group = "Compare_Group1";
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Mean_Plotting_Markers_level2_1ANOVA.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder),
                                          paste0('Mean_Plotting_Markers_level2_1ANOVA_',
                                                 Select_Group, "_", "BA_Group", "_","Subtype",
                                                 '.html')))
rm(markers.files, Select_Group)

## 5. 1ANOV_All for level2 Markers by Group  -------
markers.files = c(
  "FindMarkers_SAT_Label_2000/Annotation.Level2.Markers/SAT_Label_2000_Annotation.Select_Group_Combined.RData"
); 
sample.file.type = "Samples.df.Rds"; # sample files
Select_Group = "Compare_Group1";  Compare.G = "BA_Group"; compare.groups = c("Older_Lean", "Older_Overweight", "Older_PreD_Obese", "Diabetes");
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Mean_Plotting_Markers_level2_1ANOVA_All.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder),
                                          paste0('Mean_Plotting_Markers_level2_1ANOVA_All_',
                                                 Select_Group, "_", Compare.G, "_", "Subtype",
                                                 '.html')))
rm(markers.files, Select_Group)
markers.files = c(
  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Adipocytes_Subtype.Select_Group_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Endothelial_Subtype.Select_Group_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_FAP_Subtype.Select_Group_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Immune_Subtype.Select_Group_Combined.RData"
); 
sample.file.type = "Samples.df.Rds"; # sample files
Select_Group = "Compare_Group1";  Compare.G = "BA_Group"; compare.groups = c("Older_Lean", "Older_Overweight", "Older_PreD_Obese", "Diabetes");
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Mean_Plotting_Markers_level2_1ANOVA_All.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder),
                                          paste0('Mean_Plotting_Markers_level2_1ANOVA_All_',
                                                 Select_Group, "_", Compare.G, "_", "Subtype",
                                                 '.html')))
rm(markers.files, Select_Group)


markers.files = c(
  "FindMarkers_SAT_Label_2000/Annotation.Level2.Markers/SAT_Label_2000_Annotation.Select_Group_Combined.RData"
); 
sample.file.type = "Samples.df.Rds"; # sample files
Select_Group = "Compare_Group3";  Compare.G = "Age_Group"; compare.groups = c("Middle", "Older");
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Mean_Plotting_Markers_level2_1ANOVA_All.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder),
                                          paste0('Mean_Plotting_Markers_level2_1ANOVA_All_',
                                                 Select_Group, "_", Compare.G, "_", "Subtype",
                                                 '.html')))
rm(markers.files, Select_Group)
markers.files = c(
  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Adipocytes_Subtype.Select_Group_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Endothelial_Subtype.Select_Group_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_FAP_Subtype.Select_Group_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Immune_Subtype.Select_Group_Combined.RData"
);
Select_Group = "Compare_Group3";  Compare.G = "Age_Group"; compare.groups = c("Middle", "Older");
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Mean_Plotting_Markers_level2_1ANOVA_All.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder),
                                          paste0('Mean_Plotting_Markers_level2_1ANOVA_All_',
                                                 Select_Group, "_", Compare.G, "_", "Subtype",
                                                 '.html')))
rm(markers.files, Select_Group)


# markers.files = c(
#   "FindMarkers_SAT_Label_2000/Annotation.Level2.Markers/SAT_Label_2000_Annotation.Select_Group_Combined.RData"
# ); 
# sample.file.type = "Samples.df.Rds"; # sample files
# Select_Group = "Compare_Group2";  Compare.G = "Treatment_Group"; compare.groups = c("NUT_Pre", "NUT_Post", "SGLT2i_Pre", "SGLT2i_Post");
# rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Mean_Plotting_Markers_level2_1ANOVA_All.R"), output_format= "html_document",
#                   output_file = file.path(paste0(Disk, Project.folder, "/", html.folder),
#                                           paste0('Mean_Plotting_Markers_level2_1ANOVA_All_',
#                                                  Select_Group, "_", Compare.G, "_", "Subtype",
#                                                  '.html')))
# rm(markers.files, Select_Group)
# markers.files = c(
#   "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Adipocytes_Subtype.Select_Group_Combined.RData"
#   ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Endothelial_Subtype.Select_Group_Combined.RData"
#   ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_FAP_Subtype.Select_Group_Combined.RData"
#   ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Immune_Subtype.Select_Group_Combined.RData"
# );
# Select_Group = "Compare_Group2";  Compare.G = "Treatment_Group"; compare.groups = c("NUT_Pre", "NUT_Post", "SGLT2i_Pre", "SGLT2i_Post");
# rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Mean_Plotting_Markers_level2_1ANOVA_All.R"), output_format= "html_document",
#                   output_file = file.path(paste0(Disk, Project.folder, "/", html.folder),
#                                           paste0('Mean_Plotting_Markers_level2_1ANOVA_All_',
#                                                  Select_Group, "_", Compare.G, "_", "Subtype",
#                                                  '.html')))
# rm(markers.files, Select_Group)



## 6. 2ANOVA for level2 Markers by Group  ---------------------------
markers.files = c(
  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Adipocytes_Subtype.Compare_Group2_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Endothelial_Subtype.Compare_Group2_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_FAP_Subtype.Compare_Group2_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Immune_Subtype.Compare_Group2_Combined.RData"
);
Select_Group = "Compare_Group2"
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Mean_Plotting_Markers_level2_2ANOVA.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder),
                                          paste0('Mean_Plotting_Markers_level2_2ANOVA_',
                                                 Select_Group, "_", "Subtype",
                                                 '.html')))
rm(markers.files, Select_Group)

## 7. 2ANOVA_All for level2 Markers by Group  ---------------------------
markers.files = c(
  "FindMarkers_SAT_Label_2000/Annotation.Level2.Markers/SAT_Label_2000_Annotation.Select_Group_Combined.RData"
);
sample.file.type = "Samples.df.Rds"; # sample files
Select_Group = "Compare_Group2";  Compare.G = "Treatment_Group"; compare.groups = c("NUT_Pre", "NUT_Post", "SGLT2i_Pre", "SGLT2i_Post");
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Mean_Plotting_Markers_level2_2ANOVA_All.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder),
                                          paste0('Mean_Plotting_Markers_level2_2ANOVA_All_',
                                                 Select_Group, "_", Compare.G, "_", "Subtype",
                                                 '.html')))
rm(markers.files, Select_Group)
markers.files = c(
  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Adipocytes_Subtype.Select_Group_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Endothelial_Subtype.Select_Group_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_FAP_Subtype.Select_Group_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Immune_Subtype.Select_Group_Combined.RData"
);
Select_Group = "Compare_Group2";  Compare.G = "Treatment_Group"; compare.groups = c("NUT_Pre", "NUT_Post", "SGLT2i_Pre", "SGLT2i_Post");
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Mean_Plotting_Markers_level2_2ANOVA_All.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder),
                                          paste0('Mean_Plotting_Markers_level2_2ANOVA_All_',
                                                 Select_Group, "_", Compare.G, "_", "Subtype",
                                                 '.html')))
rm(markers.files, Select_Group)




## 8. Markers to KEGG Pathway Analysis and Plotting  ---------
### in individual group --------
markers.files = c(
  "FindMarkers_SAT_Label_2000/Annotation.Level2.Markers/SAT_Label_2000_Annotation.Compare_Group1_BA_Group_Combined.RData"
  ,   "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Adipocytes_Subtype.Compare_Group1_BA_Group_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Endothelial_Subtype.Compare_Group1_BA_Group_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_FAP_Subtype.Compare_Group1_BA_Group_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Immune_Subtype.Compare_Group1_BA_Group_Combined.RData"
); 
Select_Group = "Compare_Group1"; type = "avg_log2FC"; 
rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Main_07_MarkersLevel2_KEGGPathway_Analysis_as.group.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_07_MarkersLevel2_KEGGPathway_', 
                                                 "SAT_Label_2000", "_", Select_Group, 
                                                 '.html')))
rm(markers.files, Select_Group, type)

markers.files = c(
  "FindMarkers_SAT_Label_2000/Annotation.Level2.Markers/SAT_Label_2000_Annotation.Compare_Group3_Age_Group_Combined.RData"
  ,   "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Adipocytes_Subtype.Compare_Group3_Age_Group_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Endothelial_Subtype.Compare_Group3_Age_Group_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_FAP_Subtype.Compare_Group3_Age_Group_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Immune_Subtype.Compare_Group3_Age_Group_Combined.RData"
); 
Select_Group = "Compare_Group3"; type = "avg_log2FC";
rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Main_07_MarkersLevel2_KEGGPathway_Analysis_as.group.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_07_MarkersLevel2_KEGGPathway_', 
                                                 "SAT_Label_2000", "_", Select_Group, 
                                                 '.html')))
rm(markers.files, Select_Group, type)


### compare subgroup ---------------
markers.files = c(
  "FindMarkers_SAT_Label_2000/Annotation.Level2.Markers/SAT_Label_2000_Annotation.Compare_Group1_BA_Group_Combined.RData"
  ,   "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Adipocytes_Subtype.Compare_Group1_BA_Group_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Endothelial_Subtype.Compare_Group1_BA_Group_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_FAP_Subtype.Compare_Group1_BA_Group_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Immune_Subtype.Compare_Group1_BA_Group_Combined.RData"
); 
Select_Group = "Compare_Group1"; type = "avg_log2FC"; compare.groups = c("Older_Lean", "Older_Overweight", "Older_PreD_Obese", "Diabetes");
rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Main_07_MarkersLevel2_KEGGPathway_Analysis_as.group.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_07_MarkersLevel2_KEGGPathway_', 
                                                 "SAT_Label_2000", "_", Select_Group, 
                                                 '.html')))
rm(markers.files, Select_Group, type, compare.groups)

markers.files = c(
  "FindMarkers_SAT_Label_2000/Annotation.Level2.Markers/SAT_Label_2000_Annotation.Compare_Group3_Age_Group_Combined.RData"
  ,   "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Adipocytes_Subtype.Compare_Group3_Age_Group_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Endothelial_Subtype.Compare_Group3_Age_Group_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_FAP_Subtype.Compare_Group3_Age_Group_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Immune_Subtype.Compare_Group3_Age_Group_Combined.RData"
); 
Select_Group = "Compare_Group3"; type = "avg_log2FC";compare.groups = c("Middle", "Older");
rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Main_07_MarkersLevel2_KEGGPathway_Analysis_as.group.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_07_MarkersLevel2_KEGGPathway_', 
                                                 "SAT_Label_2000", "_", Select_Group, 
                                                 '.html')))
rm(markers.files, Select_Group, type, compare.groups)


### paired group -----------------
markers.files = c(
  "FindMarkers_SAT_Label_2000/Annotation.Level2.Markers/SAT_Label_2000_Annotation.Compare_Group2_Treatment_Group_Combined.RData"
  ,   "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Adipocytes_Subtype.Compare_Group2_Treatment_Group_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Endothelial_Subtype.Compare_Group2_Treatment_Group_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_FAP_Subtype.Compare_Group2_Treatment_Group_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Immune_Subtype.Compare_Group2_Treatment_Group_Combined.RData"
); 
Select_Group = "Compare_Group2"; type = "avg_log2FC"; 
compare.groups = c( "NUT_Pre",  "NUT_Post", "SGLT2i_Pre", "SGLT2i_Post");
compare.groups1 = c( "Pre",  "Post"); # paired group
compare.groups2 = c( "NUT",  "SGLT2i");
rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Main_07_MarkersLevel2_KEGGPathway_Analysis_paired.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_07_MarkersLevel2_KEGGPathway_', 
                                                 "SAT_Label_2000", "_", Select_Group, 
                                                 '.html')))
rm(markers.files, Select_Group, type, compare.groups, compare.groups1, compare.groups2)

## 9. Markers to GProfiler Pathway Analysis and Plotting  ---------
### Compare_Group1
markers.files = c(
  "FindMarkers_SAT_Label_2000/Annotation.Level2.Markers/SAT_Label_2000_Annotation.Compare_Group1_BA_Group_Combined.Mean.Rds"
  , "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Endothelial_Subtype.Compare_Group1_BA_Group_Combined.Mean.Rds"
  ,"FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Adipocytes_Subtype.Compare_Group1_BA_Group_Combined.Mean.Rds"
  , "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_FAP_Subtype.Compare_Group1_BA_Group_Combined.Mean.Rds"
  , "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Immune_Subtype.Compare_Group1_BA_Group_Combined.Mean.Rds"
); Select_Group = "Compare_Group1";organism = "hsapiens";
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "CellMarkerstoPathway_Analysis_Plotting_GProfiler.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('CellMarkers_GProfiler_', 
                                                 "SAT_Label_2000", "_", Select_Group, 
                                                 '.html')))
rm(markers.files, organism, Select_Group)


## 10. DEGs by zlm Regression --------
args = c("RegressionMarkers",  # output_folder
         "SAT_Label_2000_Immune.Rds", # data_name
         "F",  # str_n_gene
         "Human.adipose", # analyte
         NA,  # keep_sex
         "F", # flag_M_O 
         "Compare_Group3" # tissue, or select data set
)
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Analysis_DEGs_zlm.Regression.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Analysis_DEGs_zlm.Regression_', 
                                                 unlist(strsplit(args[2], "[.Rds]"))[1], 
                                                 "_", args[7], '.html')))
rm(args)

args = c("RegressionMarkers",  # output_folder
         "SAT_Label_2000_Immune.Rds", # data_name
         "T",  # str_n_gene
         "Human.adipose", # analyte
         NA,  # keep_sex
         "F", # flag_M_O 
         "Compare_Group3" # tissue, or select data set
)
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Analysis_DEGs_zlm.Regression.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Analysis_DEGs_zlm.Regression_', 
                                                 unlist(strsplit(args[2], "[.Rds]"))[1], 
                                                 "_", args[7], '.html')))
rm(args)

args = c("RegressionMarkers",  # output_folder
         "SAT_Label_2000_Immune.Rds", # data_name
         "F",  # str_n_gene
         "Human.adipose", # analyte
         NA,  # keep_sex
         "F", # flag_M_O 
         "Compare_Group3" # tissue, or select data set,
         , "Compare_Group1" # tissue, or select data set,
         , "BMI+A1c_Value"  # others
); sample.file.type = "SAT_Sample.meta.data.Rds"; # sample files
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Analysis_DEGs_zlm.Regression_Multiple.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Analysis_DEGs_zlm.Regression_Multiple_', 
                                                 unlist(strsplit(args[2], "[.Rds]"))[1], 
                                                 "_", args[7], '.html')))
rm(args)


args = c("RegressionMarkers",  # output_folder
         "SAT_Label_2000_FAP.Rds", # data_name
         "F",  # str_n_gene
         "Human.adipose", # analyte
         NA,  # keep_sex
         "F", # flag_M_O 
         "Compare_Group3" # tissue, or select data set
)
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Analysis_DEGs_zlm.Regression.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Analysis_DEGs_zlm.Regression_', 
                                                 unlist(strsplit(args[2], "[.Rds]"))[1], 
                                                 "_", args[7], '.html')))
rm(args)

args = c("RegressionMarkers",  # output_folder
         "SAT_Label_2000_Endothelial.Rds", # data_name
         "F",  # str_n_gene
         "Human.adipose", # analyte
         NA,  # keep_sex
         "F", # flag_M_O 
         "Compare_Group3" # tissue, or select data set
)
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Analysis_DEGs_zlm.Regression.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Analysis_DEGs_zlm.Regression_', 
                                                 unlist(strsplit(args[2], "[.Rds]"))[1], 
                                                 "_", args[7], '.html')))
rm(args)

args = c("RegressionMarkers",  # output_folder
         "SAT_Label_2000_Adipocytes.Rds", # data_name
         "F",  # str_n_gene
         "Human.adipose", # analyte
         NA,  # keep_sex
         "F", # flag_M_O 
         "Compare_Group3" # tissue, or select data set
)
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Analysis_DEGs_zlm.Regression.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Analysis_DEGs_zlm.Regression_', 
                                                 unlist(strsplit(args[2], "[.Rds]"))[1], 
                                                 "_", args[7], '.html')))
rm(args)



## level 4 markers: based on each patient ----------------------
seurat.files = c(
  "SAT_Label_2000_Immune.Rds"
  , "SAT_Label_2000_Endothelial.Rds"
  , "SAT_Label_2000_FAP.Rds"
  , "SAT_Label_2000_Adipocytes.Rds"
)
Select_Group = "Compare_Group2"
sample.file.type = "SAT_Sample.meta.data.Rds"; # sample files
meta.group1 = "patient_id" # subset seurat.data
meta.group2.1 = "Group"; # subset seurat.data,
meta.group2.2 = "Pre_Post" # compare group, levels
file.type = "Subtype.metadata.Rds"; # with Annotation.file.name # subset seurat.data
meta.group = "Subtype" # subset seurat.data
rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Analysis_FindMarkers_seurat_level4.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder),
                                          paste0('Analysis_FindMarkers_level4_',
                                                 stringr::str_split(seurat.files[1], "_")[[1]][1],"_",
                                                 stringr::str_split(file.type, ".")[[1]][1],'.html')))
rm(seurat.files, sample.file.type, meta.group1, meta.group2.1, meta.group2.2, file.type, meta.group)






markers.files = c(
  #  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Immune_DCs_B_Cells_Senescence.name.metadata.level4.Rds"
  # , "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Immune_PVM_Senescence.name.metadata.level4.Rds"
  # , "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Immune_NKT_Senescence.name.metadata.level4.Rds"
  # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Immune_CEM1_Senescence.name.metadata.level4.Rds"
  # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Immune_Mast_cells_Senescence.name.metadata.level4.Rds"
  # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Immune_LAM_Senescence.name.metadata.level4.Rds"
  # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Immune_CEM2_Senescence.name.metadata.level4.Rds"
  # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Immune_P_LAM_Senescence.name.metadata.level4.Rds"
  # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_FAP_hFAP6_Senescence.name.metadata.level4.Rds"
  # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_FAP_hFAP3_Senescence.name.metadata.level4.Rds"
  # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_FAP_Pre_adip1_Senescence.name.metadata.level4.Rds"
  # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_FAP_CPA_Senescence.name.metadata.level4.Rds"
  # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_FAP_hFAP2_Senescence.name.metadata.level4.Rds"
  # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_FAP_hFAP4_Senescence.name.metadata.level4.Rds"
  # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_FAP_hFAP7_Senescence.name.metadata.level4.Rds"
  # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_FAP_Pre_adip2_Senescence.name.metadata.level4.Rds"
  # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_FAP_hFAP1_Senescence.name.metadata.level4.Rds"
  # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_FAP_Inflammatory_hFAP_Senescence.name.metadata.level4.Rds"
  # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_FAP_hFAP5_Senescence.name.metadata.level4.Rds"
  # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Adipocytes_hAd07_Senescence.name.metadata.level4.Rds"
  # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Adipocytes_hAd03_Senescence.name.metadata.level4.Rds"
  # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Adipocytes_hAd09_Senescence.name.metadata.level4.Rds"
  # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Adipocytes_hAd04_Senescence.name.metadata.level4.Rds"
  # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Adipocytes_hFAP_Ad2_Senescence.name.metadata.level4.Rds"
  # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Adipocytes_hAd05_Senescence.name.metadata.level4.Rds"
  # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Adipocytes_hAd02_Senescence.name.metadata.level4.Rds"
  # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Adipocytes_hFAP_Ad1_Senescence.name.metadata.level4.Rds"
  # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Adipocytes_hAd10_Senescence.name.metadata.level4.Rds"
  # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Adipocytes_hFAP_Ad3_Senescence.name.metadata.level4.Rds"
  # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Adipocytes_hAd01_Senescence.name.metadata.level4.Rds"
     "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Adipocytes_hAd_NK_Senescence.name.metadata.level4.Rds"
  ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Adipocytes_hAd06_Senescence.name.metadata.level4.Rds"
  ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Adipocytes_hAd08_Senescence.name.metadata.level4.Rds"
  ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Adipocytes_hAd_Macrophage_Senescence.name.metadata.level4.Rds"
); 
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Plotting_Figure_Markers_level4.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Plotting_Figure_Markers_', "Senescence",'.html')))
rm(markers.files)



# 14. Markers to pathways analysis and plotting -----------------------------------------------

## level4
# markers.files = c(
#   #    "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Immune_DCs_B_Cells_Senescence.name.metadata.level4.combined.Rds"
#   # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Immune_PVM_Senescence.name.metadata.level4.combined.Rds"
#   # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Immune_NKT_Senescence.name.metadata.level4.combined.Rds"
#   # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Immune_CEM1_Senescence.name.metadata.level4.combined.Rds"
#   # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Immune_Mast_cells_Senescence.name.metadata.level4.combined.Rds"
#   # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Immune_LAM_Senescence.name.metadata.level4.combined.Rds"
#   # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Immune_CEM2_Senescence.name.metadata.level4.combined.Rds"
#   # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Immune_P_LAM_Senescence.name.metadata.level4.combined.Rds"
#   # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_FAP_hFAP6_Senescence.name.metadata.level4.combined.Rds"
#   # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_FAP_hFAP3_Senescence.name.metadata.level4.combined.Rds"
#   # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_FAP_Pre_adip1_Senescence.name.metadata.level4.combined.Rds"
#   # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_FAP_CPA_Senescence.name.metadata.level4.combined.Rds"
#   # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_FAP_hFAP2_Senescence.name.metadata.level4.combined.Rds"
#   # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_FAP_hFAP4_Senescence.name.metadata.level4.combined.Rds"
#   # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_FAP_hFAP7_Senescence.name.metadata.level4.combined.Rds"
#   # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_FAP_Pre_adip2_Senescence.name.metadata.level4.combined.Rds"
#   # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_FAP_hFAP1_Senescence.name.metadata.level4.combined.Rds"
#   # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_FAP_Inflammatory_hFAP_Senescence.name.metadata.level4.combined.Rds"
#   # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_FAP_hFAP5_Senescence.name.metadata.level4.combined.Rds"
#   # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Adipocytes_hAd07_Senescence.name.metadata.level4.combined.Rds"
#   # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Adipocytes_hAd03_Senescence.name.metadata.level4.combined.Rds"
#   # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Adipocytes_hAd09_Senescence.name.metadata.level4.combined.Rds"
#   # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Adipocytes_hAd04_Senescence.name.metadata.level4.combined.Rds"
#   # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Adipocytes_hFAP_Ad2_Senescence.name.metadata.level4.combined.Rds"
#   # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Adipocytes_hAd05_Senescence.name.metadata.level4.combined.Rds"
#   # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Adipocytes_hAd02_Senescence.name.metadata.level4.combined.Rds"
#   # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Adipocytes_hFAP_Ad1_Senescence.name.metadata.level4.combined.Rds"
#   # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Adipocytes_hAd10_Senescence.name.metadata.level4.combined.Rds"
#   # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Adipocytes_hFAP_Ad3_Senescence.name.metadata.level4.combined.Rds"
#   # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Adipocytes_hAd01_Senescence.name.metadata.level4.combined.Rds"
#   # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Adipocytes_hAd_NK_Senescence.name.metadata.level4.combined.Rds"
#   "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Adipocytes_hAd06_Senescence.name.metadata.level4.combined.Rds"
#   ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Adipocytes_hAd08_Senescence.name.metadata.level4.combined.Rds"
#   ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Adipocytes_hAd_Macrophage_Senescence.name.metadata.level4.combined.Rds"
# ); 
# organism = "Homo sapiens" # organism = "Mus musculus"
# type = "avg_log2FC."
# rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Main_07_MarkersLevel4_toPathway_Analysis.R"), output_format= "html_document",
#                   output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
#                                           paste0('Analysis_MarkersLevel4_toPathway_', "SAT_Label_2000_Immune", '.html')))
# rm(markers.files, organism, type)


pathways.files = c(
  # ,  "Pathways/Pathway_avg_log2FC.SAT_Label_2000_Immune_DCs_B_Cells_Senescence.name.metadata.level4.combined.Rds"
  # , "Pathways/Pathway_avg_log2FC.SAT_Label_2000_Immune_NKT_Senescence.name.metadata.level4.combined.Rds"
  # , "Pathways/Pathway_avg_log2FC.SAT_Label_2000_Immune_PVM_Senescence.name.metadata.level4.combined.Rds"
  # ,   "Pathways/Pathway_avg_log2FC.SAT_Label_2000_Immune_CEM2_Senescence.name.metadata.level4.combined.Rds"
  # ,  "Pathways/Pathway_avg_log2FC.SAT_Label_2000_Immune_P_LAM_Senescence.name.metadata.level4.combined.Rds"
  # ,  "Pathways/Pathway_avg_log2FC.SAT_Label_2000_Immune_CEM1_Senescence.name.metadata.level4.combined.Rds"
  # ,  "Pathways/Pathway_avg_log2FC.SAT_Label_2000_Immune_Mast_cells_Senescence.name.metadata.level4.combined.Rds"
  # ,  "Pathways/Pathway_avg_log2FC.SAT_Label_2000_Immune_LAM_Senescence.name.metadata.level4.combined.Rds"
  # ,  "Pathways/Pathway_avg_log2FC.SAT_Label_2000_FAP_hFAP6_Senescence.name.metadata.level4.combined.Rds"
  # ,  "Pathways/Pathway_avg_log2FC.SAT_Label_2000_FAP_hFAP3_Senescence.name.metadata.level4.combined.Rds"
  # ,  "Pathways/Pathway_avg_log2FC.SAT_Label_2000_FAP_Pre_adip1_Senescence.name.metadata.level4.combined.Rds"
  # ,  "Pathways/Pathway_avg_log2FC.SAT_Label_2000_FAP_CPA_Senescence.name.metadata.level4.combined.Rds"
  # ,  "Pathways/Pathway_avg_log2FC.SAT_Label_2000_FAP_hFAP2_Senescence.name.metadata.level4.combined.Rds"
  # ,  "Pathways/Pathway_avg_log2FC.SAT_Label_2000_FAP_hFAP4_Senescence.name.metadata.level4.combined.Rds"
  # ,  "Pathways/Pathway_avg_log2FC.SAT_Label_2000_FAP_hFAP7_Senescence.name.metadata.level4.combined.Rds"
  # ,  "Pathways/Pathway_avg_log2FC.SAT_Label_2000_FAP_Pre_adip2_Senescence.name.metadata.level4.combined.Rds"
  # ,  "Pathways/Pathway_avg_log2FC.SAT_Label_2000_FAP_hFAP1_Senescence.name.metadata.level4.combined.Rds"
  # ,  "Pathways/Pathway_avg_log2FC.SAT_Label_2000_FAP_Inflammatory_hFAP_Senescence.name.metadata.level4.combined.Rds"
  # ,  "Pathways/Pathway_avg_log2FC.SAT_Label_2000_FAP_hFAP5_Senescence.name.metadata.level4.combined.Rds"
  # ,  "Pathways/Pathway_avg_log2FC.SAT_Label_2000_Adipocytes_hAd07_Senescence.name.metadata.level4.combined.Rds"
  # ,  "Pathways/Pathway_avg_log2FC.SAT_Label_2000_Adipocytes_hAd03_Senescence.name.metadata.level4.combined.Rds"
  # ,  "Pathways/Pathway_avg_log2FC.SAT_Label_2000_Adipocytes_hAd09_Senescence.name.metadata.level4.combined.Rds"
  # ,  "Pathways/Pathway_avg_log2FC.SAT_Label_2000_Adipocytes_hAd04_Senescence.name.metadata.level4.combined.Rds"
  # ,  "Pathways/Pathway_avg_log2FC.SAT_Label_2000_Adipocytes_hFAP_Ad2_Senescence.name.metadata.level4.combined.Rds"
  # ,  "Pathways/Pathway_avg_log2FC.SAT_Label_2000_Adipocytes_hAd05_Senescence.name.metadata.level4.combined.Rds"
  # ,  "Pathways/Pathway_avg_log2FC.SAT_Label_2000_Adipocytes_hAd02_Senescence.name.metadata.level4.combined.Rds"
  # ,  "Pathways/Pathway_avg_log2FC.SAT_Label_2000_Adipocytes_hFAP_Ad1_Senescence.name.metadata.level4.combined.Rds"
  # ,  "Pathways/Pathway_avg_log2FC.SAT_Label_2000_Adipocytes_hAd10_Senescence.name.metadata.level4.combined.Rds"
  # ,  "Pathways/Pathway_avg_log2FC.SAT_Label_2000_Adipocytes_hFAP_Ad3_Senescence.name.metadata.level4.combined.Rds"
  # ,  "Pathways/Pathway_avg_log2FC.SAT_Label_2000_Adipocytes_hAd01_Senescence.name.metadata.level4.combined.Rds"
  # ,  "Pathways/Pathway_avg_log2FC.SAT_Label_2000_Adipocytes_hAd_NK_Senescence.name.metadata.level4.combined.Rds"
  "Pathways/Pathway_avg_log2FC.SAT_Label_2000_Adipocytes_hAd06_Senescence.name.metadata.level4.combined.Rds"
  ,  "Pathways/Pathway_avg_log2FC.SAT_Label_2000_Adipocytes_hAd08_Senescence.name.metadata.level4.combined.Rds"
  ,  "Pathways/Pathway_avg_log2FC.SAT_Label_2000_Adipocytes_hAd_Macrophage_Senescence.name.metadata.level4.combined.Rds"
)
rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Plotting_Figure_Analysis_Pathways.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Plotting_Figure_Analysis_pathways_Level4', "SAT_Label_2000_Immune",  '.html')))
rm(pathways.files)


## KEGG pathways -------

markers.files = c(
  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Adipocytes_Subtype.Compare_Group3_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Endothelial_Subtype.Compare_Group3_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_FAP_Subtype.Compare_Group3_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Immune_Subtype.Compare_Group3_Combined.RData"
); 


markers.files = c(
  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Adipocytes_Subtype.Compare_Group2_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Endothelial_Subtype.Compare_Group2_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_FAP_Subtype.Compare_Group2_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Immune_Subtype.Compare_Group2_Combined.RData"
); 


markers.files = c(
  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Adipocytes_Subtype.Compare_Group1_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Endothelial_Subtype.Compare_Group1_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_FAP_Subtype.Compare_Group1_Combined.RData"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level2.Markers/SAT_Label_2000_Immune_Subtype.Compare_Group1_Combined.RData"
); 
Select_Group = "Compare_Group1"; type = "avg_log2FC";
rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Main_07_MarkersLevel2_KEGGPathway_Analysis_as.group.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_07_MarkersLevel2_KEGGPathway_', 
                                                 Select_Group, "_", 
                                                 "SAT_Label_2000_Subtype",
                                                 '.html')))
rm(markers.files, organism, type)


markers.files = c(
     "FindMarkers_SAT_Label_2000/Subtype.Level4.Markers/SAT_Label_2000_Adipocytes_Subtype.metadata.level4.combined.Rds"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level4.Markers/SAT_Label_2000_Endothelial_Subtype.metadata.level4.combined.Rds"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level4.Markers/SAT_Label_2000_FAP_Subtype.metadata.level4.combined.Rds"
  ,  "FindMarkers_SAT_Label_2000/Subtype.Level4.Markers/SAT_Label_2000_Immune_Subtype.metadata.level4.combined.Rds"
); 
type = "avg_log2FC"
rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Main_07_MarkersLevel4_KEGGPathway_Analysis_as.group.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Main_07_MarkersLevel4_KEGGPathway_', 
                                                 "SAT_Label_2000_Immune", 
                                                 '.html')))
rm(markers.files, organism, type)



# 15. Plotting umap -------------------------------------------------------------------------------
seurat.files = c("SAT_Label_2000.Rds")
file.type = "Subtype.metadata.Rds"; meta.group = c("Annotation"); Select_Group = "Select_Group";
groups.list <- list("1" = list(subset.group = c("Compare_Group1"), split.group = c("BA_Group")),
                    "2" = list(subset.group = c("Compare_Group2"), split.group = c("Treatment_Group")))
for (seurat.file in seurat.files) {
  rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Plotting_Figure_Level1_umap.R"), output_format= "html_document",
                    output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                            paste0('Plotting_Figure_Level1_umap_', stringr::str_split(seurat.file, ".Rds")[[1]][1],"_", 
                                                   paste0(meta.group), '.html')))
}
rm(seurat.files, meta.group, file.type, groups.list, Select_Group)
seurat.files = c("SAT_Label_2000.Rds")
file.type = "Subtype.metadata.Rds"; meta.group = c("Annotation"); Select_Group = "Select_Group";
groups.list <- list("3" = list(subset.group = c("Compare_Group3"), split.group = c("Age_Group")) )
for (seurat.file in seurat.files) {
  rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Plotting_Figure_Level1_umap.R"), output_format= "html_document",
                    output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                            paste0('Plotting_Figure_Level1_umap_', stringr::str_split(seurat.file, ".Rds")[[1]][1],"_", 
                                                   paste0(meta.group), '.html')))
}
rm(seurat.files, meta.group, file.type, groups.list, Select_Group)



seurat.files = c("SAT_Label_2000.Rds")
file.type = "Subtype.metadata.Rds"; meta.group = c("Subtype"); Select_Group = "Select_Group";
groups.list <- list("1" = list(subset.group = c("Compare_Group1"), split.group = c("BA_Group")),
                    "2" = list(subset.group = c("Compare_Group2"), split.group = c("Treatment_Group")) )
for (seurat.file in seurat.files) {
  rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Plotting_Figure_Level1_umap.R"), output_format= "html_document",
                    output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                            paste0('Plotting_Figure_Level1_umap_', stringr::str_split(seurat.file, ".Rds")[[1]][1],"_", 
                                                   paste0(meta.group), '.html')))
}
rm(seurat.files, meta.group, file.type, groups.list, Select_Group)
seurat.files = c("SAT_Label_2000.Rds")
file.type = "Subtype.metadata.Rds"; meta.group = c("Subtype"); Select_Group = "Select_Group";
groups.list <- list("3" = list(subset.group = c("Compare_Group3"), split.group = c("Age_Group")) )
for (seurat.file in seurat.files) {
  rmarkdown::render(input= paste0(Disk, Project.folder, "/", codes.folder, "/", "Plotting_Figure_Level1_umap.R"), output_format= "html_document",
                    output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                            paste0('Plotting_Figure_Level1_umap_', stringr::str_split(seurat.file, ".Rds")[[1]][1],"_", 
                                                   paste0(meta.group), '.html')))
}
rm(seurat.files, meta.group, file.type, groups.list, Select_Group)



seurat.files = c("SAT_Label_2000_FAP.Rds", "SAT_Label_2000_Endothelial.Rds",
                 "SAT_Label_2000_Immune.Rds", "SAT_Label_2000_Adipocytes.Rds")
file.type = "Subtype.metadata.Rds"; meta.group = "Subtype"; Select_Group = "Select_Group";
groups.list <- list("1" = list(subset.group = c("Compare_Group1"), split.group = c("BA_Group")),
                    "2" = list(subset.group = c("Compare_Group2"), split.group = c("Treatment_Group")) ,
                    "3" = list(subset.group = c("Compare_Group3"), split.group = c("Age_Group")) )
for (seurat.file in seurat.files) {
  rmarkdown::render(input= paste0(Disk,Project.folder, "/", codes.folder, "/", "Plotting_Figure_Level1_umap.R"), output_format= "html_document",
                    output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                            paste0('Plotting_Figure_Level1_umap_', stringr::str_split(seurat.file, ".Rds")[[1]][1],"_", 
                                                   paste0(meta.group), '.html')))
}
rm(seurat.files, meta.group, file.type, groups.list, Select_Group)

seurat.files = c(
  "SAT_Label_2000.Rds"
)
file.type = "Senescence.name.metadata.Rds"; meta.group = "sene.name"; Select_Group = "Select_Group";
groups.list <- list("1" = list(subset.group = c("Compare_Group1"), split.group = c("BA_Group")),
                    "2" = list(subset.group = c("Compare_Group2"), split.group = c("Treatment_Group"))  )
for (seurat.file in seurat.files) {
  rmarkdown::render(input= paste0(Disk,Project.folder, "/", codes.folder, "/", "Plotting_Figure_Level1_umap.R"), output_format= "html_document",
                    output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                            paste0('Plotting_Figure_Level1_umap_', stringr::str_split(seurat.file, ".Rds")[[1]][1],"_", 
                                                   paste0(meta.group), '.html')))
}
rm(seurat.files, meta.group, file.type)
seurat.files = c(
  "SAT_Label_2000.Rds"
)
file.type = "Senescence.name.metadata.Rds"; meta.group = "sene.name"; Select_Group = "Select_Group";
groups.list <- list("3" = list(subset.group = c("Compare_Group3"), split.group = c("Age_Group")) )
for (seurat.file in seurat.files) {
  rmarkdown::render(input= paste0(Disk,Project.folder, "/", codes.folder, "/", "Plotting_Figure_Level1_umap.R"), output_format= "html_document",
                    output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                            paste0('Plotting_Figure_Level1_umap_', stringr::str_split(seurat.file, ".Rds")[[1]][1],"_", 
                                                   paste0(meta.group), '.html')))
}
rm(seurat.files, meta.group, file.type)
seurat.files = c(
  "SAT_Label_2000_FAP.Rds", "SAT_Label_2000_Endothelial.Rds",
  "SAT_Label_2000_Immune.Rds", "SAT_Label_2000_Adipocytes.Rds"
)
file.type = "Senescence.name.metadata.Rds"; meta.group = "sene.name"; Select_Group = "Select_Group";
groups.list <- list("1" = list(subset.group = c("Compare_Group1"), split.group = c("BA_Group")),
                    "2" = list(subset.group = c("Compare_Group2"), split.group = c("Treatment_Group")) ,
                    "3" = list(subset.group = c("Compare_Group3"), split.group = c("Age_Group")) )
for (seurat.file in seurat.files) {
  rmarkdown::render(input= paste0(Disk,Project.folder, "/", codes.folder, "/", "Plotting_Figure_Level1_umap.R"), output_format= "html_document",
                    output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                            paste0('Plotting_Figure_Level1_umap_', stringr::str_split(seurat.file, ".Rds")[[1]][1],"_", 
                                                   paste0(meta.group), '.html')))
}
rm(seurat.files, meta.group, file.type)



# seurat.files = c(
#   "SAT_Label_2000_Immune_PVM.Rds"
#   , "SAT_Label_2000_Immune_NKT.Rds"
#   , "SAT_Label_2000_Immune_CEM1.Rds"
#   , "SAT_Label_2000_Immune_DCs_B_Cells.Rds"
#   , "SAT_Label_2000_Immune_Mast_cells.Rds"
#   , "SAT_Label_2000_Immune_LAM.Rds"
#   , "SAT_Label_2000_Immune_CEM2.Rds"
#   , "SAT_Label_2000_Immune_P_LAM.Rds"
#   
#   , "SAT_Label_2000_FAP_hFAP6.Rds"
#   , "SAT_Label_2000_FAP_hFAP3.Rds"
#   , "SAT_Label_2000_FAP_Pre_adip1.Rds"
#   , "SAT_Label_2000_FAP_CPA.Rds"
#   , "SAT_Label_2000_FAP_hFAP2.Rds"
#   , "SAT_Label_2000_FAP_hFAP4.Rds"
#   , "SAT_Label_2000_FAP_hFAP7.Rds"
#   , "SAT_Label_2000_FAP_Pre_adip2.Rds"
#   , "SAT_Label_2000_FAP_hFAP1.Rds"
#   , "SAT_Label_2000_FAP_Inflammatory_hFAP.Rds"
#   , "SAT_Label_2000_FAP_hFAP5.Rds"
#   
#   , "SAT_Label_2000_Adipocytes_hAd07.Rds"
#   , "SAT_Label_2000_Adipocytes_hAd03.Rds"
#   , "SAT_Label_2000_Adipocytes_hAd09.Rds"
#   , "SAT_Label_2000_Adipocytes_hAd04.Rds"
#   , "SAT_Label_2000_Adipocytes_hFAP_Ad2.Rds"
#   , "SAT_Label_2000_Adipocytes_hAd05.Rds"
#   , "SAT_Label_2000_Adipocytes_hAd02.Rds"
#   , "SAT_Label_2000_Adipocytes_hFAP_Ad1.Rds"
#   , "SAT_Label_2000_Adipocytes_hAd10.Rds"
#   , "SAT_Label_2000_Adipocytes_hFAP_Ad3.Rds"
#   , "SAT_Label_2000_Adipocytes_hAd01.Rds"
#   , "SAT_Label_2000_Adipocytes_hAd_NK.Rds"
#   , "SAT_Label_2000_Adipocytes_hAd06.Rds"
#   , "SAT_Label_2000_Adipocytes_hAd08.Rds"
#   , "SAT_Label_2000_Adipocytes_hAd_Macrophage.Rds"
#   , "SAT_Label_2000_Adipocytes_hFAP_Ad4.Rds"
#   
#   , "SAT_Label_2000_Endothelial_MSC.Rds"
#   , "SAT_Label_2000_Endothelial_EC.Rds"
#   , "SAT_Label_2000_Endothelial_MC.Rds"
#   , "SAT_Label_2000_Endothelial_Inflammatory_MC.Rds"
#   , "SAT_Label_2000_Endothelial_VEC.Rds"
#   , "SAT_Label_2000_Endothelial_SMC.Rds"
#   , "SAT_Label_2000_Endothelial_Fibromyocyte.Rds"
#   , "SAT_Label_2000_Endothelial_Fibroblast.Rds"
#   , "SAT_Label_2000_Endothelial_Inflammatory_EC.Rds"
#   , "SAT_Label_2000_Endothelial_Inflammatory_SMC.Rds"
# )
# sub.folder <- "SAT_Label_2000_Subtype"
# file.type = "Senescence.name.metadata.Rds"; meta.group = "sene.name"; Select_Group = "Select_Group";
# groups.list <- list("1" = list(subset.group = c("Compare_Group1"), split.group = c("BA_Group")),
#                     "2" = list(subset.group = c("Compare_Group2"), split.group = c("Treatment_Group")) ,
#                     "3" = list(subset.group = c("Compare_Group3"), split.group = c("Age_Group")) )
# for (seurat.file in seurat.files) {
#   rmarkdown::render(input= paste0(Disk,Project.folder, "/", codes.folder, "/", "Plotting_Figure_Level1_umap.R"), output_format= "html_document",
#                     output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
#                                             paste0('Plotting_Figure_Level1_umap_', stringr::str_split(seurat.file, ".Rds")[[1]][1],"_", 
#                                                    paste0(meta.group), '.html')))
# }
# rm(seurat.files, meta.group, file.type, sub.folder)


# 16. CellChat anaylsis -------------------
seurat.files = c(
  # "SAT_Label_2000.Rds",
  "SAT_Label_2000_Immune.Rds"
  # , "SAT_Label_2000_FAP.Rds"
  # , "SAT_Label_2000_Adipocytes.Rds"
  # , "SAT_Label_2000_Endothelial.Rds"
)
sample.file.type = "SAT_Sample.meta.data.Rds"; # sample files
meta.group2.1 = "Group"; # subset seurat.data, 
meta.group2.2 = "Pre_Post" # compare group, levels
file.type = "Subtype.metadata.Rds"; # with Annotation.file.name # subset seurat.data
meta.group = "Subtype" # subset seurat.data
for (seurat.file in seurat.files) {
  rmarkdown::render(input= paste0(Disk,Project.folder, "/", codes.folder, "/", "Plotting_Figure_Level1_umap.R"), output_format= "html_document",
                    output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                            paste0('Plotting_Figure_Level1_umap_', stringr::str_split(seurat.file, ".Rds")[[1]][1],"_", 
                                                   paste0(meta.group), '.html')))
}
rm(seurat.files, meta.group, file.type)

rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Analysis_FindMarkers_seurat_level4.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/", html.folder), 
                                          paste0('Analysis_FindMarkers_level4_', 
                                                 stringr::str_split(seurat.files[1], "_")[[1]][1],"_", 
                                                 meta.group,'.html')))
rm(seurat.files, sample.file.type, meta.group1, meta.group2.1, meta.group2.2, file.type, meta.group)

# 17. Features plotting -----------
seurat.files = c(
  "SAT_Label_2000.Rds",
  "SAT_Label_2000_FAP.Rds", "SAT_Label_2000_Endothelial.Rds",
  "SAT_Label_2000_Immune.Rds", "SAT_Label_2000_Adipocytes.Rds"
)
features = c("CDKN2A", "CDKN1A") # "CDKN2A"-p16, "CDKN1A"-p21
split.groups = c("BA_Group");select.idents = c("Older_Lean", "Older_Overweight", "Older_PreD_Obese", "Diabetes"); file.type = "Subtype.metadata.Rds"
rmarkdown::render(input= paste0(Disk, Project.folder, "/","Plotting_seurat_Features.R"), output_format= "html_document",
                  output_file = file.path(
                    paste0(Disk, Project.folder, "/", html.folder), 
                    paste0('Plotting_seurat_Features', 
                           ifelse(length(features) >2, paste0(paste(features[1:2], collapse = "."), "_", length(features) ),
                                  paste(features, collapse = ".") ), '.html')) )# 11∶24∶58 PM 
rm(seurat.files, features, split.groups, select.idents, file.type)

seurat.files = c(
  "SAT_Label_2000.Rds",
  "SAT_Label_2000_FAP.Rds", "SAT_Label_2000_Endothelial.Rds",
  "SAT_Label_2000_Immune.Rds", "SAT_Label_2000_Adipocytes.Rds"
)
features = c("ADIPOQ", "LEP") # "ADIPOQ"-Adiponectin, "LEP"- leptin
split.groups = c("BA_Group");select.idents = c("Middle", "Older_Lean", "Older_Overweight", "Older_PreD_Obese", "Diabetes"); file.type = "Subtype.metadata.Rds"
rmarkdown::render(input= paste0(Disk, Project.folder, "/","Plotting_seurat_Features.R"), output_format= "html_document",
                  output_file = file.path(
                    paste0(Disk, Project.folder, "/", html.folder), 
                    paste0('Plotting_seurat_Features', 
                           ifelse(length(features) >2, paste0(paste(features[1:2], collapse = "."), "_", length(features) ),
                                  paste(features, collapse = ".") ), '.html')) )# 11∶24∶58 PM 
rm(seurat.files, features, split.groups, select.idents, file.type)


seurat.files = c(
  "SAT_Label_2000.Rds",
  "SAT_Label_2000_FAP.Rds", "SAT_Label_2000_Endothelial.Rds",
  "SAT_Label_2000_Immune.Rds", "SAT_Label_2000_Adipocytes.Rds"
)
features = c("GLP1R", "INSR") # "GLP1R"-glucagon-like peptide 1 receptor, "INSR"- Insulin Receptor
split.groups = c("BA_Group");select.idents = c("Middle", "Older_Lean", "Older_Overweight", "Older_PreD_Obese", "Diabetes"); file.type = "Subtype.metadata.Rds"
rmarkdown::render(input= paste0(Disk, Project.folder, "/","Plotting_seurat_Features.R"), output_format= "html_document",
                  output_file = file.path(
                    paste0(Disk, Project.folder, "/", html.folder), 
                    paste0('Plotting_seurat_Features', 
                           ifelse(length(features) >2, paste0(paste(features[1:2], collapse = "."), "_", length(features) ),
                                  paste(features, collapse = ".") ), '.html')) )# 11∶24∶58 PM 
rm(seurat.files, features, split.groups, select.idents, file.type)

seurat.files = c(
  "SAT_Label_2000.Rds",
  "SAT_Label_2000_FAP.Rds", "SAT_Label_2000_Endothelial.Rds",
  "SAT_Label_2000_Immune.Rds", "SAT_Label_2000_Adipocytes.Rds"
)
features = c("PDGFRA", "SCA", "CD81") # "GLP1R"-glucagon-like peptide 1 receptor, "INSR"- Insulin Receptor
split.groups = c("BA_Group");select.idents = c("Middle", "Older_Lean", "Older_Overweight", "Older_PreD_Obese", "Diabetes"); file.type = "Subtype.metadata.Rds"
rmarkdown::render(input= paste0(Disk, Project.folder, "/","Plotting_seurat_Features.R"), output_format= "html_document",
                  output_file = file.path(
                    paste0(Disk, Project.folder, "/", html.folder), 
                    paste0('Plotting_seurat_Features', 
                           ifelse(length(features) >2, paste0(paste(features[1:2], collapse = "."), "_", length(features) ),
                                  paste(features, collapse = ".") ), '.html')) )# 11∶24∶58 PM 
rm(seurat.files, features, split.groups, select.idents, file.type)

seurat.files = c(
  "SAT_Label_2000.Rds"
)
features = c("PDGFRA","PLIN4", "LIPE","PPARG", "DCN","COL1A1","VWF", "FLT1","PTPRB","CD247", "SKAP1", "ITGAL") #
file.type = "Subtype.metadata.Rds"
rmarkdown::render(input= paste0(Disk, Project.folder, "/","Plotting_seurat_Features.R"), output_format= "html_document",
                  output_file = file.path(
                    paste0(Disk, Project.folder, "/", html.folder), 
                    paste0('Plotting_seurat_Features', 
                           ifelse(length(features) >2, paste0(paste(features[1:2], collapse = "."), "_", length(features) ),
                                  paste(features, collapse = ".") ), '.html')) )# 11∶24∶58 PM 
rm(seurat.files, features, if(exists("split.groups")) {return(split.groups, select.idents)}, file.type)


seurat.files = c(
  "SAT_Label_2000.Rds"
)
features = c("ADIPOQ","VWF", "FGF14","PTPRC") #
file.type = "Subtype.metadata.Rds"; Select_Group = "Select_Group"; 
rmarkdown::render(input= paste0(Disk, Project.folder, "/","Plotting_seurat_Features.R"), output_format= "html_document",
                  output_file = file.path(
                    paste0(Disk, Project.folder, "/", html.folder), 
                    paste0('Plotting_seurat_Features', 
                           ifelse(length(features) >2, paste0(paste(features[1:2], collapse = "."), "_", length(features) ),
                                  paste(features, collapse = ".") ), '.html')) )# 11∶24∶58 PM 
rm(seurat.files, features,  file.type, Select_Group )


seurat.files = c(
  "SAT_Label_2000.Rds"
)
features = c("PLIN5","PLIN2","PLIN1","PLIN3","PLIN4") #
split.groups = c("BA_Group");select.idents = c("Older_Lean", "Older_Overweight", "Older_PreD_Obese", "Diabetes"); 
file.type = "Subtype.metadata.Rds"; Select_Group = "Select_Group"; 
rmarkdown::render(input= paste0(Disk, Project.folder, "/","Plotting_seurat_Features.R"), output_format= "html_document",
                  output_file = file.path(
                    paste0(Disk, Project.folder, "/", html.folder), 
                    paste0('Plotting_seurat_Features', 
                           ifelse(length(features) >2, paste0(paste(features[1:2], collapse = "."), "_", length(features) ),
                                  paste(features, collapse = ".") ), '.html')) )# 11∶24∶58 PM 
rm(seurat.files, features,  file.type, Select_Group )