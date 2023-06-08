# 1. Create Parameters file --------------------------------
#' create `Project Parameters.R` files for each project
##' counts foler should be `ProjectName_aggr`
rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memory and report the memory usage.
memory.limit(size = 1e+13)
Disk <- c("/media/jianie/Extreme SSD1/") 
# Disk <- c("D:/")
Project.folder <- c("2022-09-01 STARR_SGLT2 Combine")


# 2. prerun ----------------------------------------------------------------------------------------------------------
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Step01_prerun_Creat_matrices.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), paste0('Step01_prerun_Creat_matrices.html'))) # 
#' generate: combine `project name_mat.rds`, feature.names.rds,  only need run once.

rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Step02_prerun_QC_Aggr_all.genes.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), paste0('Step02_prerun_QC_Aggr_all.genes.html'))) # 
#' generate: individual `sample number_mat.Rds`, combined project `dataset.group.Rds`, combined project `mat.Rds`(list), in the counts folder. only need run once.

# 3. QC ----------------------------------------------------------------------------------------------------------
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Step03_QC_Aggr_all.genes_step1.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), paste0('Step03_QC_Aggr_all.genes_step1.html'))) # 
#' QC, filter, and generate: `Ambient.Rds`, `QC files_Genes.RData`, `QC files_Cells.RData`, `Final.sce.list.Rds`, .only need run once.

rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Step03_QC_Aggr_all.genes_step2.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), paste0('Step03_QC_Aggr_all.genes_step2.html'))) # 
#' QC, filter, and generate: `sce.list folder`, `sce.list2 folder`, 
#' `Final.sce.list_part1.Rds`, `Final.sce.list_part2.Rds`, `Final.sce.list_part3.Rds`. only need run once.

# if samples <= 36
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Step04_QC_Aggr_all.genes_step3.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), paste0('Step04_QC_Aggr_all.genes_step3.html'))) # 
#' QC, filter, and generate: `rescaled.Rds`, `seurat.list folder`. only need run once.

# if samples > 36
# run  Step04_QC_Aggr_all.genes_step3.collab.R on colab:https://colab.research.google.com/#create=true&language=r
#' QC, filter, and generate:  `rescaled.Rds` or `rescaled_part1.Rds`, `rescaled_part2.Rds`, `rescaled_part3.Rds`,
#' `seurat.list.combine`: `seurat.list.rds`, or `seurat.list_part1.Rds`, `seurat.list_part2.Rds`,`seurat.list_part3.Rds`,`seurat.list_part4.Rds`. only need run once.
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Step04_QC_Aggr_all.genes_step4_after.collab.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), paste0('Step04_QC_Aggr_all.genes_step4_after.collab.html'))) # this step DONE in Lunix PC
#' seperate to `seurat.list folder`


rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Step05_Combine_seurat_data.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), paste0('Step05_Combine_seurat_data.html'))) # this step DONE in Lunix PC
#' `SAT folder`: combined `SAT1.Rds`, combined `SAT2.Rds`, or combined `SAT.Rds` (depending on sample number).only need run once.


rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Step06_Combine_seurat_data.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), paste0('Step06_Combine_seurat_data.html'))) # this step DONE in Lunix PC
#' `Seurat folder`: combined `SAT.Rds` (if sample number > 18 & <= 36).only need run once. 
#'run Step06_Combine_seurat_data.colab.R on colab:https://colab.research.google.com/#create=true&language=r
#' 25 min upload, 35 min running, 35 min download. 
#' save a copy of SAT.rds file (optional)

rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Step07_QC_Aggr_all.genes_step3.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), paste0('Step07_QC_Aggr_all.genes_step3.html'))) # 
#' import `SAT.Rds` for Seurat analysis. generate: `SAT.basic.meta.data.Rds`, 
#' `SAT_basic.meta.data.Rds`
#' `SAT_Label_1000.Rds`, `SAT_1000.meta.data.Rds`.
#' `SAT_Label_2000.Rds`, `SAT_2000.meta.data.Rds`.
#' `SAT_Label_3000.Rds`, `SAT_3000.meta.data.Rds`.

# GCCRI
rmarkdown::render(input= paste0("Step07_QC_Aggr_all.genes_step3_GCCRI.R"), output_format= "html_document",
                  output_file = file.path(paste0('html.result/Step07_QC_Aggr_all.genes_step3.html'))) # 
#' import `SAT.Rds` for Seurat analysis. generate: `SAT.basic.meta.data.Rds`, 
#' `SAT_Label_1000.Rds`, `SAT_1000.meta.data.Rds`.
#' `SAT_Label_2000.Rds`, `SAT_2000.meta.data.Rds`.
#' `SAT_Label_3000.Rds`, `SAT_3000.meta.data.Rds`.


#' plot QC
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Step08_QC_Aggr_all.genes_step4.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), paste0('Step08_QC_Aggr_all.genes_step4.html'))) # 

## Add sample metadata -------------------------------------
rmarkdown::render(input= paste0(Disk, Project.folder, "/","Metadata_Add_SGLT2_Samples.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                          paste0('Metadata_Add_SGLT2_Samples', '.html')) )# 11∶24∶58 PM 

## plot samples
export.folder = "Samples"
rmarkdown::render(input= paste0(Disk, Project.folder, "/", "Plotting_SGLT2_Samples.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                          paste0('Plotting_SGLT2_Samples', '.html')))
rm(export.folder)



seurat.files = c("SAT_Label_1000.Rds", "SAT_Label_2000.Rds", "SAT_Label_3000.Rds")
rmarkdown::render(input= paste0(Disk, Project.folder, "/","Metadata_Add_SGLT2_Samples_seurat_Plotting.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                          paste0('Metadata_Add_SGLT2_Samples_seurat_Plotting', '.html')) )# 11∶24∶58 PM 
rm(seurat.files)



# 4. Test Resolution --------------------------------------------------------------
## test all resolutions
seurat.file = "SAT_Label_1000.Rds"
res.range <- c(0.005, 0.01); res.by = 0.001
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_00_Step1_Test_Resolution.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                          paste0('Main_00_Step1_Test_Resolution', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 "_",  res.range[1], "_", res.range[2], '.html'))) # 

seurat.file = "SAT_Label_2000.Rds"
res.range <- c(0.005, 0.01); res.by = 0.001
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_00_Step1_Test_Resolution.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                          paste0('Main_00_Step1_Test_Resolution', stringr::str_split(seurat.file, ".Rds")[[1]][1],
                                                 "_",  res.range[1], "_", res.range[2], '.html'))) # 


# 5. cell type analysis -------------------------------------------------------------------------------------------------
# Overall: Annotation
#' check Resolution output files for best resolution and variable gene number
## compare clusters and select best resolution
seurat.file = "SAT_Label_2000.Rds"
res = "_0.005_0.01"  
select.res = "originalexp_snn_res.0.006"
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_01_Human_Fat_Annotation_Compare.rmd"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                          paste0('Main_01_Human_Fat_Annotation_Compare_', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 stringr::str_split(select.res, "originalexp_snn_res")[[1]][2],  '.html')) ) # 11∶17∶45 PM
rm(seurat.file, res, select.res)

## apply final resolution to seurat data
seurat.file = "SAT_Label_2000.Rds"
res = "_0.005_0.01"
select.res = "originalexp_snn_res.0.006"
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_00_Step2_Final_Resolution2.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                          paste0('Main_00_Step2_Final_Resolution_', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 stringr::str_split(select.res, "originalexp_snn_res")[[1]][2], '.html')) ) # 
rm(seurat.file, res, select.res)

## Add annotation metadata --------------------------------------------------------------------------
seurat.file = "SAT_Label_2000.Rds"
# need extra files: Fat_Subtype.xlsx, _Final_Res.metadata.Rds
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Metadata_Add_Annotation.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                          paste0('Metadata_Annotation_', stringr::str_split(seurat.file, ".Rds")[[1]][1], '.html')) )# 11∶24∶58 PM 
rm(seurat.file)




# 6. cell subtype analysis -------------------------------------------------------------------------------------------------
## FAP ------------------------------------------------------
## create seurat data
seurat.file = "SAT_Label_2000.Rds"
celltype.file =  "FAP.Rds"
rmarkdown::render(input= paste0(Disk, Project.folder, "/","Main_04_FAP_Fat_1_Step1.rmd"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                          paste0('Main_04_FAP_Fat_1_Step1_', stringr::str_split(seurat.file, ".Rds")[[1]][1], '.html'))) # 11∶32∶10 PM 

## test all resolutions
seurat.file = "SAT_Label_2000_FAP.Rds"
res.range <- c(0.02, 0.16); res.by = 0.02
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_00_Step1_Test_Resolution.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                          paste0('Main_00_Step1_Test_Resolution_', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 "_",  res.range[1], "_", res.range[2], '.html'))) # 
rm(seurat.file, res.range, res.by)

seurat.file = "SAT_Label_2000_FAP.Rds"
res.range <- c(0.1, 0.5); res.by = 0.05
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_00_Step1_Test_Resolution.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                          paste0('Main_00_Step1_Test_Resolution_', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 "_",  res.range[1], "_", res.range[2], '.html'))) # 
rm(seurat.file, res.range, res.by)

## compare clusters and select best resolution
seurat.file = "SAT_Label_2000_FAP.Rds"
res = "_0.02_0.16"  
select.res = "originalexp_snn_res.0.12"
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_04_Human_FAP_Subtype_Compare.rmd"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                          paste0('Main_04_Human_FAP_Subtype_Compare_', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 stringr::str_split(select.res, "originalexp_snn_res")[[1]][2], '.html'))) # 11∶36∶24 PM 
rm(seurat.file, res, select.res)

seurat.file = "SAT_Label_2000_FAP.Rds"
res = "_0.1_0.5"  
select.res = "originalexp_snn_res.0.2"
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_04_Human_FAP_Subtype_Compare.rmd"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                          paste0('Main_04_Human_FAP_Subtype_Compare_', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 stringr::str_split(select.res, "originalexp_snn_res")[[1]][2], '.html'))) # 11∶36∶24 PM 
rm(seurat.file, res, select.res)

seurat.file = "SAT_Label_2000_FAP.Rds"
res = "_0.1_0.5"  
select.res = "originalexp_snn_res.0.5"
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_04_Human_FAP_Subtype_Compare.rmd"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                          paste0('Main_04_Human_FAP_Subtype_Compare_', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 stringr::str_split(select.res, "originalexp_snn_res")[[1]][2], '.html'))) # 11∶36∶24 PM 
rm(seurat.file, res, select.res)

seurat.file = "SAT_Label_2000_FAP.Rds"
res = "_0.1_0.5"  
select.res = "originalexp_snn_res.0.45"
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_04_Human_FAP_Subtype_Compare.rmd"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                          paste0('Main_04_Human_FAP_Subtype_Compare_', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 stringr::str_split(select.res, "originalexp_snn_res")[[1]][2], '.html'))) # 11∶36∶24 PM 
rm(seurat.file, res, select.res)

## apply final resolution to seurat data
seurat.file = "SAT_Label_2000_FAP.Rds"
res = "_0.1_0.5"
select.res = "originalexp_snn_res.0.45"
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_00_Step2_Final_Resolution2.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                          paste0('Main_00_Step2_Final_Resolution_', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 stringr::str_split(select.res, "originalexp_snn_res")[[1]][2],  '.html')) ) # 
rm(seurat.file, res, select.res)

## Adipocytes ------------------------------------------------------------
## create seurat data
seurat.file = "SAT_Label_2000.Rds"
celltype.file =  "Adipocytes.Rds"
rmarkdown::render(input= paste0(Disk, Project.folder, "/","Main_05_Adipocytes_Fat_1_Step1.rmd"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                          paste0('Main_05_Adipocytes_Fat_1_Step1_', stringr::str_split(seurat.file, ".Rds")[[1]][1], '.html'))) # 11∶59∶03 PM

## test all resolutions
seurat.file = "SAT_Label_2000_Adipocytes.Rds"
res.range <- c(0.01, 0.05); res.by = 0.005
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_00_Step1_Test_Resolution.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                          paste0('Main_00_Step1_Test_Resolution', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 "_",  res.range[1], "_", res.range[2], '.html'))) # 
rm(seurat.file, res.range, res.by)

seurat.file = "SAT_Label_2000_Adipocytes.Rds"
res.range <- c(0.05, 0.5); res.by = 0.05
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_00_Step1_Test_Resolution.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                          paste0('Main_00_Step1_Test_Resolution', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 "_",  res.range[1], "_", res.range[2], '.html'))) # 
rm(seurat.file, res.range, res.by)

seurat.file = "SAT_Label_2000_Adipocytes.Rds"
res.range <- c(0.5, 1); res.by = 0.1
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_00_Step1_Test_Resolution.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                          paste0('Main_00_Step1_Test_Resolution', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 "_",  res.range[1], "_", res.range[2], '.html'))) # 
rm(seurat.file, res.range, res.by)

# compare clusters and select best resolution
seurat.file = "SAT_Label_2000_Adipocytes.Rds"
res = "_0.05_0.5"   
select.res = "originalexp_snn_res.0.05"
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_05_Human_Adipocytes_Subtype_Compare.rmd"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"),
                                          paste0('Main_05_Human_Adipocytes_Subtype_Compare_', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 stringr::str_split(select.res, "originalexp_snn_res")[[1]][2], '.html'))) # 12∶07∶55 AM
rm(seurat.file, res, select.res)

seurat.file = "SAT_Label_2000_Adipocytes.Rds"
res = "_0.05_0.5" 
select.res = "originalexp_snn_res.0.4"
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_05_Human_Adipocytes_Subtype_Compare.rmd"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"),
                                          paste0('Main_05_Human_Adipocytes_Subtype_Compare_', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 stringr::str_split(select.res, "originalexp_snn_res")[[1]][2], '.html'))) # 12∶07∶55 AM
rm(seurat.file, res, select.res)

seurat.file = "SAT_Label_2000_Adipocytes.Rds"
res = "_0.05_0.5" 
select.res = "originalexp_snn_res.0.5"
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_05_Human_Adipocytes_Subtype_Compare.rmd"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"),
                                          paste0('Main_05_Human_Adipocytes_Subtype_Compare_', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 stringr::str_split(select.res, "originalexp_snn_res")[[1]][2], '.html'))) # 12∶07∶55 AM
rm(seurat.file, res, select.res)

## apply final resolution to seurat data
seurat.file = "SAT_Label_2000_Adipocytes.Rds"
res = "_0.05_0.5"
select.res = "originalexp_snn_res.0.4"
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_00_Step2_Final_Resolution2.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                          paste0('Main_00_Step2_Final_Resolution_', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 stringr::str_split(select.res, "originalexp_snn_res")[[1]][2],  '.html')) ) # 
rm(seurat.file, res, select.res)


## Endothelial --------------------------------------------------------------------------------------
## create seurat data
seurat.file = "SAT_Label_2000.Rds"
celltype.file =  "Endothelial.Rds"
rmarkdown::render(input= paste0(Disk, Project.folder, "/","Main_03_Endothelial_Fat_1_Step1.rmd"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                          paste0('Main_03_Endothelial_Fat_1_Step1_', stringr::str_split(seurat.file, ".Rds")[[1]][1], '.html'))) # 12∶36∶31 AM

## test all resolutions
seurat.file = "SAT_Label_2000_Endothelial.Rds"
res = "_0.01_0.1" 
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_00_Step1_Test_Resolution.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                          paste0('Main_00_Step1_Test_Resolution', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 "_",  res.range[1], "_", res.range[2], '.html'))) # 
rm(seurat.file, res.range, res.by)

seurat.file = "SAT_Label_2000_Endothelial.Rds"
res = "_0.1_1"
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_00_Step1_Test_Resolution.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                          paste0('Main_00_Step1_Test_Resolution', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 "_",  res.range[1], "_", res.range[2], '.html'))) # 
rm(seurat.file, res.range, res.by)

# compare clusters and select best resolution
seurat.file = "SAT_Label_2000_Endothelial.Rds"
res = "_0.01_0.1"
select.res = "originalexp_snn_res.0.01"
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_03_Human_Endothelial_Subtype_Compare.rmd"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"),
                                          paste0('Main_03_Human_Endothelial_Subtype_Compare_', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 stringr::str_split(select.res, "originalexp_snn_res")[[1]][2], '.html'))) #  12∶38∶06 AM
rm(seurat.file, res, select.res)

seurat.file = "SAT_Label_2000_Endothelial.Rds"
res = "_0.01_0.1"
select.res = "originalexp_snn_res.0.02"
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_03_Human_Endothelial_Subtype_Compare.rmd"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"),
                                          paste0('Main_03_Human_Endothelial_Subtype_Compare_', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 stringr::str_split(select.res, "originalexp_snn_res")[[1]][2], '.html'))) #  12∶38∶06 AM
rm(seurat.file, res, select.res)

seurat.file = "SAT_Label_2000_Endothelial.Rds"
res = "_0.1_1"
select.res = "originalexp_snn_res.0.4"
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_03_Human_Endothelial_Subtype_Compare.rmd"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"),
                                          paste0('Main_03_Human_Endothelial_Subtype_Compare_', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 stringr::str_split(select.res, "originalexp_snn_res")[[1]][2], '.html'))) #  12∶38∶06 AM
rm(seurat.file, res, select.res)

seurat.file = "SAT_Label_2000_Endothelial.Rds"
res = "_0.1_1"
select.res = "originalexp_snn_res.0.1"
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_03_Human_Endothelial_Subtype_Compare.rmd"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"),
                                          paste0('Main_03_Human_Endothelial_Subtype_Compare_', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 stringr::str_split(select.res, "originalexp_snn_res")[[1]][2], '.html'))) #  12∶38∶06 AM
rm(seurat.file, res, select.res)

## apply final resolution to seurat data
seurat.file = "SAT_Label_2000_Endothelial.Rds"
res = "_0.01_0.1"
select.res = "originalexp_snn_res.0.02"
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_00_Step2_Final_Resolution2.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                          paste0('Main_00_Step2_Final_Resolution_', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 stringr::str_split(select.res, "originalexp_snn_res")[[1]][2], '.html')) ) # 
rm(seurat.file, res, select.res)


## Myeloid ---------------------------------------------------------------------------------------------
## create seurat data
seurat.file = "SAT_Label_2000.Rds"
celltype.file =  "Myeloid.Rds"
rmarkdown::render(input= paste0(Disk, Project.folder, "/","Main_02_Myeloid_Fat_1_Step1.rmd"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                          paste0('Main_02_Myeloid_Fat_1_Step1_', stringr::str_split(seurat.file, ".Rds")[[1]][1], '.html'))) # 12∶22∶38 AM

## test all resolutions
seurat.file = "SAT_Label_2000_Myeloid.Rds"
res.range <- c(0.03, 0.1); res.by = 0.01
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_00_Step1_Test_Resolution.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                          paste0('Main_00_Step1_Test_Resolution', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 "_",  res.range[1], "_", res.range[2], '.html'))) # 
rm(seurat.file, res.range, res.by)

seurat.file = "SAT_Label_2000_Myeloid.Rds"
res.range <- c(0.01, 0.1); res.by = 0.01
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_00_Step1_Test_Resolution.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                          paste0('Main_00_Step1_Test_Resolution', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 "_",  res.range[1], "_", res.range[2], '.html'))) # 
rm(seurat.file, res.range, res.by)

## compare clusters and select best resolution
seurat.file = "SAT_Label_2000_Myeloid.Rds"
res = "_0.01_0.1"
select.res = "originalexp_snn_res.0.02"
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_02_Human_Immune_Subtype_Compare.rmd"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                          paste0('Main_02_Human_Immune_Subtype_Compare_', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 stringr::str_split(select.res, "originalexp_snn_res")[[1]][2], '.html'))) # 12∶25∶47 AM
rm(seurat.file, res, select.res)

seurat.file = "SAT_Label_2000_Myeloid.Rds"
res = "_0.01_0.1"
select.res = "originalexp_snn_res.0.1"
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_02_Human_Immune_Subtype_Compare.rmd"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                          paste0('Main_02_Human_Immune_Subtype_Compare_', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 stringr::str_split(select.res, "originalexp_snn_res")[[1]][2], '.html'))) # 12∶25∶47 AM
rm(seurat.file, res, select.res)


## apply final resolution to seurat data
seurat.file = "SAT_Label_2000_Myeloid.Rds"
res = "_0.01_0.1"
select.res = "originalexp_snn_res.0.02"
rmarkdown::render(input= paste0(Disk, Project.folder, "/","Main_00_Step2_Final_Resolution2.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                          paste0('Main_00_Step2_Final_Resolution_', stringr::str_split(seurat.file, ".Rds")[[1]][1], 
                                                 stringr::str_split(select.res, "originalexp_snn_res")[[1]][2], '.html')) ) # 
rm(seurat.file, res, select.res)



## FAP and Adipocytes combine --------------------------------------------------------
# rmarkdown::render(input= paste0(Disk, Project.folder, "/","Main_05_FAP_Adipocytes_Fat_1_Step1.rmd"), output_format= "html_document") 
# rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_05_FAP_Adipocytes_Fat_1_Step2_for.annotated.rmd"), output_format= "html_document") 


# 7. Add subtype metadata  ------------------------------------------------------------------------
seurat.files = c("SAT_Label_2000_FAP.Rds", "SAT_Label_2000_Endothelial.Rds",
                 "SAT_Label_2000_Myeloid.Rds", "SAT_Label_2000_Adipocytes.Rds")
# need extra files: Fat_Subtype.xlsx, Annotation.metadata.Rds, Final_Res.metadata.Rds,
for (seurat.file in seurat.files) {
  # seurat.file = "SAT_Label_2000_FAP.Rds"
  rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Metadata_Add_Subtype.R"), output_format= "html_document",
                    output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                            paste0('Metadata_Subtype_', stringr::str_split(seurat.file, ".Rds")[[1]][1], '.html')))
}


## combine
seurat.file = "SAT_Label_2000.Rds"
rmarkdown::render(input= paste0(Disk, Project.folder, "/", "Metadata_Combine_Subtype.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                          paste0('Metadata_Combine_Subtype_', stringr::str_split(seurat.file, ".Rds")[[1]][1], '.html')))
rm(seurat.file)



## subset seurat ----
seurat.files = c(
  "SAT_Label_2000_Myeloid.Rds"
  , "SAT_Label_2000_FAP.Rds"
  , "SAT_Label_2000_Adipocytes.Rds"
  , "SAT_Label_2000_Endothelial.Rds"
)
file.type = "Subtype.metadata.Rds"; # with Annotation.file.name # subset seurat.data
subset.group = "Subtype" # subset seurat.data
rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Subset_seurat.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                          paste0('Subset_seurat_', 
                                                 stringr::str_split(seurat.files[1], "_")[[1]][1],"_", 
                                                 stringr::str_split(file.type, ".")[[1]][1],'.html')))
rm(seurat.files, file.type, subset.group)


# 8. cell Senescence analysis -------------------------------------------------------------------------------------------------
seurat.file = "SAT_Label_2000.Rds"
file.type = "Subtype.metadata.Rds"
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_06_Senescence_GeneModule_Compare.rmd"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                          paste0('Main_06_Senescence_', stringr::str_split(seurat.file, ".Rds")[[1]][1], '.html')))
rm(seurat.file)


# compare by Event_Name: V5, V11.
seurat.files = c("SAT_Label_2000.Rds"
                 , "SAT_Label_2000_FAP.Rds"
                 , "SAT_Label_2000_Endothelial.Rds"
                 , "SAT_Label_2000_Myeloid.Rds"
                 , "SAT_Label_2000_Adipocytes.Rds"
)
col.groups1 = c("Senescence.GO.score", "Senescence.score", "SenMayo")
file.type1 = "Senescence_SAT_Label_2000/Senescence.Marker_ModuleScore_Label_metadata.rds"
col.groups = c("Annotation", "Subtype")
row.groups = "Dataset"
file.type = "Subtype.metadata.Rds"
export.folder = "ModuleScore"
rmarkdown::render(input= paste0(Disk, Project.folder, "/", "Plotting_Figure_Analysis_ModuleScore.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                          paste0('Plotting_Figure_Analysis_ModuleScore_', 
                                                 '.html')))
rm(seurat.files, col.groups, row.groups, file.type, col.groups1, file.type1, export.folder)


# compare by Group: C & D
seurat.files = c("SAT_Label_2000.Rds"
                 , "SAT_Label_2000_FAP.Rds"
                 , "SAT_Label_2000_Endothelial.Rds"
                 , "SAT_Label_2000_Myeloid.Rds"
                 , "SAT_Label_2000_Adipocytes.Rds"
)
col.groups1 = c("Senescence.GO.score", "Senescence.score", "SenMayo")
file.type1 = "Senescence_SAT_Label_2000/Senescence.Marker_ModuleScore_Label_metadata.rds"
col.groups = c("Annotation", "Subtype")
row.groups = "Dataset"
file.type = "Subtype.metadata.Rds"
export.folder = "ModuleScore"
rmarkdown::render(input= paste0(Disk, Project.folder, "/", "Plotting_Figure_Analysis_ModuleScore_plot2.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                          paste0('Plotting_Figure_Analysis_ModuleScore_', 
                                                 '.html')))
rm(seurat.files, col.groups, row.groups, file.type, col.groups1, file.type1, export.folder)


# Average then compare by 1 way ANOVA, 2 way ANOVA,
seurat.files = c("SAT_Label_2000.Rds"
                 , "SAT_Label_2000_FAP.Rds"
                 , "SAT_Label_2000_Endothelial.Rds"
                 , "SAT_Label_2000_Myeloid.Rds"
                 , "SAT_Label_2000_Adipocytes.Rds"
)
col.groups1 = c("Senescence.GO.score", "Senescence.score", "SenMayo")
file.type1 = "Senescence_SAT_Label_2000/Senescence.Marker_ModuleScore_Label_metadata.rds"
col.groups = c("Annotation", "Subtype")
row.groups = "Dataset"
file.type = "Subtype.metadata.Rds"
export.folder = "ModuleScore"
rmarkdown::render(input= paste0(Disk, Project.folder, "/", "Plotting_Figure_Analysis_ModuleScore_Average.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                          paste0('Plotting_Figure_Analysis_ModuleScore_', 
                                                 '.html')))
rm(seurat.files, col.groups, row.groups, file.type, col.groups1, file.type1, export.folder)


# Add Senescence in metadata 
seurat.file = "SAT_Label_2000.Rds"
file.type1 = "Senescence_SAT_Label_2000/Senescence.Marker_ModuleScore_Label_metadata.rds"
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Metadata_Add_Senescence.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                          paste0('Metadata_Senescence_', stringr::str_split(seurat.file, ".Rds")[[1]][1], '.html')) )# 11∶24∶58 PM 
rm(seurat.file)



# 9. cell Trajetory analysis -------------------------------------------------------------------------------------------------

# 10. Composition Analysis -----------------------------------------
seurat.files = c("SAT_Label_2000.Rds",
                 "SAT_Label_2000_FAP.Rds", "SAT_Label_2000_Endothelial.Rds",
                 "SAT_Label_2000_Myeloid.Rds", "SAT_Label_2000_Adipocytes.Rds")
# "Annotation", "Subtype"
col.groups = c("Annotation", "Subtype")
row.group = "Dataset"
file.type = "Subtype.metadata.Rds"
# need extra file: Samples.df.Rds, _Sample.meta.data.Rds
for (seurat.file in seurat.files) {
  # seurat.file = "SAT_Label_2000.Rds"
  rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Analysis_Composition.rmd"), output_format= "html_document",
                    output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                            paste0('Analysis_Composition_', stringr::str_split(seurat.file, ".Rds")[[1]][1],"_",
                                                   stringr::str_split(file.type, ".Rds")[[1]][1], '.html')))
  # print(seurat.file)
}
rm(seurat.files, col.groups, row.group, file.type)


# "Senescence"
seurat.files = c("SAT_Label_2000.Rds",
                 "SAT_Label_2000_FAP.Rds", "SAT_Label_2000_Endothelial.Rds",
                 "SAT_Label_2000_Myeloid.Rds", "SAT_Label_2000_Adipocytes.Rds")
# "Annotation", "Subtype"
col.groups1 = c("Senescence.GO")
file.type1 = "Senescence_SAT_Label_2000/Senescence_ModuleScore_metadata.Rds"
col.groups = c("Annotation", "Subtype")
row.group = "Dataset"
file.type = "Subtype.metadata.Rds"
for (seurat.file in seurat.files) {
  # seurat.file = "SAT_Label_2000.Rds"
  rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Analysis_Composition_Level2.rmd"), output_format= "html_document",
                    output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                            paste0('Analysis_Composition_', stringr::str_split(seurat.file, ".Rds")[[1]][1],"_",
                                                   stringr::str_split(file.type, ".Rds")[[1]][1], "_",
                                                   col.groups1, '.html')))
  # print(seurat.file)
}
rm(seurat.file, seurat.files, col.groups, row.group, file.type, col.groups1, file.type1)

seurat.files = c("SAT_Label_2000.Rds",
                 "SAT_Label_2000_FAP.Rds", "SAT_Label_2000_Endothelial.Rds",
                 "SAT_Label_2000_Myeloid.Rds", "SAT_Label_2000_Adipocytes.Rds")
# "Annotation", "Subtype"
file.type1 = "Senescence.name.metadata.Rds"; col.groups1 = "sene.name"
col.groups = c("Annotation", "Subtype")
row.group = "Dataset"
file.type = "Subtype.metadata.Rds"
for (seurat.file in seurat.files) {
  # seurat.file = "SAT_Label_2000.Rds"
  rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Analysis_Composition_Level2.rmd"), output_format= "html_document",
                    output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                            paste0('Analysis_Composition_', stringr::str_split(seurat.file, ".Rds")[[1]][1], "_",
                                                   stringr::str_split(file.type, ".Rds")[[1]][1], "_",
                                                   col.groups1, '.html')))
  # print(seurat.file)
}
rm(seurat.file, seurat.files, col.groups, row.group, file.type, col.groups1, file.type1)



## plotting
export.folder = "Compositions"
file.type = "Subtype.metadata.Rds"
rmarkdown::render(input= paste0(Disk, Project.folder, "/", "Plotting_Figure_Composition_Event.name.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                          paste0('Plotting_Figure_Composition_Event.name_', 
                                                 stringr::str_split(file.type, ".Rds")[[1]][1], '.html')))
rm(export.folder, file.type)

export.folder = "Compositions"
rmarkdown::render(input= paste0(Disk, Project.folder, "/", "Plotting_Figure_Composition_Group_Level2.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                          paste0('Plotting_Figure_Composition_Group_Level2_', 
                                                 # stringr::str_split(file.type, ".Rds")[[1]][1], 
                                                 '.html')))
rm(export.folder)

export.folder = "Compositions"
rmarkdown::render(input= paste0(Disk, Project.folder, "/", "Plotting_Figure_Composition_Event.name_Level2.R"), 
                  output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                          paste0('Plotting_Figure_Composition_Event.name_Level2_', 
                                                 # stringr::str_split(file.type, ".Rds")[[1]][1], 
                                                 '.html')))
rm(export.folder)

# 11. Stats Analysis -------------------------------------------------------
compare.list <- list(Group = list(compare.group="Group", compare.group.names= list("C", "D")), # 2 way anova factor 1
                     `Event Name` = list(compare.group="Event Name", compare.group.names= list("Visit 5", "Visit 11")),# 2 way anova factor 2
                     Treatment_Group = list(compare.group="Treatment_Group", 
                                            compare.group.names= list("C_Visit 5", "D_Visit 5","C_Visit 11", "D_Visit 11"))
) 
rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Analysis_Composition_Stats.rmd"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                          paste0('Analysis_Composition_Stats', '.html')))
rm(compare.list)


compare.list <- list( `Event Name` = list(compare.group="Event Name", 
                                          compare.group.names= list("Visit 5", "Visit 11"))# 2 way anova factor 2
) 
rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Analysis_Composition_Paired.rmd"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                          paste0('Analysis_Composition_Paired', '.html')))
rm(compare.list)


compare.list <- list( `Event Name` = list(compare.group="Event Name", 
                                          compare.group.names= list("Visit 5", "Visit 11"),
                                          sub.group = "Group",
                                          sub.group.names= list("C", "D") )# 2 way anova factor 2
) 
rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Analysis_Composition_Paired_subgroup.rmd"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                          paste0('Analysis_Composition_Paired_subgroup', '.html')))
rm(compare.list)


# 12. cell marker analysis ---------------------------------------------------------------------------------------------------
# row numbers in metadata.Rds should be same with the whole seurat.data
## level 1 markers
### cell type markers
seurat.files = c("SAT_Label_2000.Rds")
file.type = "Subtype.metadata.Rds"; meta.group = "Annotation"
rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Analysis_FindMarkers_seurat.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                          paste0('Analysis_FindMarkers_', stringr::str_split(seurat.files, ".Rds")[[1]][1],"_", 
                                                 paste0(meta.group), '.html')))
rm(seurat.files, meta.group, file.type)


seurat.files = c("SAT_Label_2000.Rds",
                 "SAT_Label_2000_FAP.Rds", "SAT_Label_2000_Endothelial.Rds",
                 "SAT_Label_2000_Myeloid.Rds", "SAT_Label_2000_Adipocytes.Rds")
file.type = "Subtype.metadata.Rds"; meta.group = "Subtype"
rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Analysis_FindMarkers_seurat.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                          paste0('Analysis_FindMarkers_', stringr::str_split(seurat.files[1], ".Rds")[[1]][1],"_", 
                                                 paste0(meta.group), '.html')))
rm(seurat.files, meta.group, file.type)

### Senescence markers
seurat.files = c(
  # "SAT_Label_2000.Rds",
    "SAT_Label_2000_Myeloid.Rds"
  , "SAT_Label_2000_FAP.Rds"
  , "SAT_Label_2000_Adipocytes.Rds"
  , "SAT_Label_2000_Endothelial.Rds"
  )
file.type = "Senescence.name.metadata.Rds"; meta.group = "sene.name"
rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Analysis_FindMarkers_seurat.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                          paste0('Analysis_FindMarkers_', stringr::str_split(seurat.files[1], ".Rds")[[1]][1],"_", 
                                                 paste0(meta.group), '.html')))
rm(seurat.files, meta.group, file.type)



## level 3 markers
seurat.files = c(
  # "SAT_Label_2000.Rds",
  "SAT_Label_2000_Myeloid.Rds"
  , "SAT_Label_2000_FAP.Rds"
  , "SAT_Label_2000_Adipocytes.Rds"
  , "SAT_Label_2000_Endothelial.Rds"
)
sample.file.type = "SAT_Sample.meta.data.Rds"; # sample files
meta.group2.1 = "Group"; # subset seurat.data, 
meta.group2.2 = "Event.Name" # compare group, levels
file.type = "Subtype.metadata.Rds"; # with Annotation.file.name # subset seurat.data
meta.group = "Subtype" # subset seurat.data
rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Analysis_FindMarkers_seurat_level3.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                          paste0('Analysis_FindMarkers_level3_', 
                                                 stringr::str_split(seurat.files[1], "_")[[1]][1],"_", 
                                                 stringr::str_split(file.type, ".")[[1]][1],'.html')))
rm(seurat.files, sample.file.type, meta.group2.1, meta.group2.2, file.type, meta.group)


seurat.files = c(
  # "SAT_Label_2000.Rds",
  "SAT_Label_2000_Myeloid.Rds"
  # , "SAT_Label_2000_FAP.Rds"
  # , "SAT_Label_2000_Adipocytes.Rds"
  # , "SAT_Label_2000_Endothelial.Rds"
)
sample.file.type = "SAT_Sample.meta.data.Rds"; # sample files
meta.group2.1 = "Group"; # subset seurat.data, 
meta.group2.2 = "Event.Name" # compare group, levels
file.type = "Senescence.name.metadata.Rds"; # with Annotation.file.name # subset seurat.data
meta.group = "sene.name" # subset seurat.data
rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Analysis_FindMarkers_seurat_level3.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                          paste0('Analysis_FindMarkers_level3_', 
                                                 stringr::str_split(seurat.files[1], "_")[[1]][1],"_", 
                                                 stringr::str_split(file.type, ".")[[1]][1],'.html')))
rm(seurat.files, sample.file.type, meta.group2.1, meta.group2.2, file.type, meta.group)


## level4 markers: based on each patient
seurat.files = c(
    "SAT_Label_2000_Myeloid_PVM.Rds"
  , "SAT_Label_2000_Myeloid_NKT.Rds"
  , "SAT_Label_2000_Myeloid_CEM1.Rds"
  , "SAT_Label_2000_Myeloid_DCs_B_Cells.Rds"
  , "SAT_Label_2000_Myeloid_Mast_cells.Rds"
  , "SAT_Label_2000_Myeloid_LAM.Rds"
  , "SAT_Label_2000_Myeloid_CEM2.Rds"
  , "SAT_Label_2000_Myeloid_P_LAM.Rds"

  , "SAT_Label_2000_FAP_hFAP6.Rds"
  , "SAT_Label_2000_FAP_hFAP3.Rds"
  , "SAT_Label_2000_FAP_Pre_adip1.Rds"
  , "SAT_Label_2000_FAP_CPA.Rds"
  , "SAT_Label_2000_FAP_hFAP2.Rds"
  , "SAT_Label_2000_FAP_hFAP4.Rds"
  , "SAT_Label_2000_FAP_hFAP7.Rds"
  , "SAT_Label_2000_FAP_Pre_adip2.Rds"
  , "SAT_Label_2000_FAP_hFAP1.Rds"
  , "SAT_Label_2000_FAP_Inflammatory_hFAP.Rds"
  , "SAT_Label_2000_FAP_hFAP5.Rds"

  , "SAT_Label_2000_Adipocytes_hAd07.Rds"
  , "SAT_Label_2000_Adipocytes_hAd03.Rds"
  , "SAT_Label_2000_Adipocytes_hAd09.Rds"
  , "SAT_Label_2000_Adipocytes_hAd04.Rds"
  , "SAT_Label_2000_Adipocytes_hFAP_Ad2.Rds"
  , "SAT_Label_2000_Adipocytes_hAd05.Rds"
  , "SAT_Label_2000_Adipocytes_hAd02.Rds"
  , "SAT_Label_2000_Adipocytes_hFAP_Ad1.Rds"
  , "SAT_Label_2000_Adipocytes_hAd10.Rds"
  , "SAT_Label_2000_Adipocytes_hFAP_Ad3.Rds"
  , "SAT_Label_2000_Adipocytes_hAd01.Rds"
  , "SAT_Label_2000_Adipocytes_hAd_NK.Rds"
  , "SAT_Label_2000_Adipocytes_hAd06.Rds"
  , "SAT_Label_2000_Adipocytes_hAd08.Rds"
  , "SAT_Label_2000_Adipocytes_hAd_Macrophage.Rds"
  , "SAT_Label_2000_Adipocytes_hFAP_Ad4.Rds"

  , "SAT_Label_2000_Endothelial_MSC.Rds"
  , "SAT_Label_2000_Endothelial_EC.Rds"
  , "SAT_Label_2000_Endothelial_MC.Rds"
  , "SAT_Label_2000_Endothelial_Inflammatory_MC.Rds"
  , "SAT_Label_2000_Endothelial_VEC.Rds"
  , "SAT_Label_2000_Endothelial_SMC.Rds"
  , "SAT_Label_2000_Endothelial_Fibromyocyte.Rds"
  , "SAT_Label_2000_Endothelial_Fibroblast.Rds"
  , "SAT_Label_2000_Endothelial_Inflammatory_EC.Rds"
  , "SAT_Label_2000_Endothelial_Inflammatory_SMC.Rds"
)
sub.folder <- "SAT_Label_2000_Subtype"
sample.file.type = "SAT_Sample.meta.data.Rds"; # sample files
meta.group1 = "patient_id" # subset seurat.data
meta.group2.1 = "Group"; # subset seurat.data, 
meta.group2.2 = "Event.Name" # compare group, levels
file.type = "Senescence.name.metadata.Rds"; # with Annotation.file.name # subset seurat.data
meta.group = "sene.name" # subset seurat.data
rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Analysis_FindMarkers_seurat_level4.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                          paste0('Analysis_FindMarkers_level4_', 
                                                 stringr::str_split(seurat.files[1], "_")[[1]][1],"_", 
                                                 stringr::str_split(file.type, ".")[[1]][1],'.html')))
rm(seurat.files, sample.file.type, meta.group1, meta.group2.1, meta.group2.2, file.type, meta.group)

seurat.files = c(
  #  "SAT_Label_2000.Rds",
  "SAT_Label_2000_Myeloid.Rds"
  , "SAT_Label_2000_FAP.Rds"
  , "SAT_Label_2000_Adipocytes.Rds"
  , "SAT_Label_2000_Endothelial.Rds"
)
sample.file.type = "SAT_Sample.meta.data.Rds"; # sample files
meta.group1 = "patient_id" # subset seurat.data
meta.group2.1 = "Group"; # subset seurat.data, 
meta.group2.2 = "Event.Name" # compare group, levels
file.type = "Senescence.name.metadata.Rds"; # with Annotation.file.name # subset seurat.data
meta.group = "sene.name" # subset seurat.data
rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Analysis_FindMarkers_seurat_level4.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                          paste0('Analysis_FindMarkers_level4_', 
                                                 stringr::str_split(seurat.files[1], "_")[[1]][1],"_", 
                                                 stringr::str_split(file.type, ".")[[1]][1],'.html')))
rm(seurat.files, sample.file.type, meta.group1, meta.group2.1, meta.group2.2, file.type, meta.group)


seurat.files = c(
 #  "SAT_Label_2000.Rds",
  "SAT_Label_2000_Myeloid.Rds"
  # , "SAT_Label_2000_FAP.Rds"
  # , "SAT_Label_2000_Adipocytes.Rds"
  # , "SAT_Label_2000_Endothelial.Rds"
)
sub.folder <- "SAT_Label_2000_Subtype"
sample.file.type = "SAT_Sample.meta.data.Rds"; # sample files
meta.group1 = "patient_id" # subset seurat.data
meta.group2.1 = "Group"; # subset seurat.data, 
meta.group2.2 = "Event.Name" # compare group, levels
file.type = "Subtype.metadata.Rds"; # with Annotation.file.name # subset seurat.data
meta.group = "Subtype" # subset seurat.data
rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Analysis_FindMarkers_seurat_level4.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                          paste0('Analysis_FindMarkers_level4_', 
                                                 stringr::str_split(seurat.files[1], "_")[[1]][1],"_", 
                                                 stringr::str_split(file.type, ".")[[1]][1],'.html')))
rm(seurat.files, sample.file.type, meta.group1, meta.group2.1, meta.group2.2, file.type, meta.group)


seurat.files = c(
    "SAT_Label_2000.Rds"
)
sub.folder <- "SAT_Label_2000_Subtype"
sample.file.type = "SAT_Sample.meta.data.Rds"; # sample files
meta.group1 = "patient_id" # subset seurat.data
meta.group2.1 = "Group"; # subset seurat.data, 
meta.group2.2 = "Event.Name" # compare group, levels
file.type = "Subtype.metadata.Rds"; # with Annotation.file.name # subset seurat.data
meta.group = "Annotation" # subset seurat.data
rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Analysis_FindMarkers_seurat_level4.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                          paste0('Analysis_FindMarkers_level4_', 
                                                 stringr::str_split(seurat.files[1], "_")[[1]][1],"_", 
                                                 meta.group,'.html')))
rm(seurat.files, sample.file.type, meta.group1, meta.group2.1, meta.group2.2, file.type, meta.group)

# 13. Markers plotting -----------------------------------------------
# scatter plotting
markers.files = c(
    "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Myeloid_Senescence.name.metadata.Rds"
  #, "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_FAP_Senescence.name.metadata.Rds"
  #, "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Adipocytes_Senescence.name.metadata.Rds"
  #, "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Endothelial_Senescence.name.metadata.Rds"
); 
rmarkdown::render(input= paste0(Disk, Project.folder,  "/", "Plotting_Figure_Markers_level3.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder), 
                                          paste0('Plotting_Figure_Markers_', "Senescence",'.html')))
rm(markers.files)


markers.files = c(
 #  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Myeloid_DCs_B_Cells_Senescence.name.metadata.level4.Rds"
 # , "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Myeloid_PVM_Senescence.name.metadata.level4.Rds"
 # , "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Myeloid_NKT_Senescence.name.metadata.level4.Rds"
 # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Myeloid_CEM1_Senescence.name.metadata.level4.Rds"
 # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Myeloid_Mast_cells_Senescence.name.metadata.level4.Rds"
 # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Myeloid_LAM_Senescence.name.metadata.level4.Rds"
 # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Myeloid_CEM2_Senescence.name.metadata.level4.Rds"
 # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Myeloid_P_LAM_Senescence.name.metadata.level4.Rds"
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
rmarkdown::render(input= paste0(Disk, Project.folder,  "/", "Plotting_Figure_Markers_level4.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder), 
                                          paste0('Plotting_Figure_Markers_', "Senescence",'.html')))
rm(markers.files)


# 14. Plotting umap -------------------------------------------------------------------------------
seurat.files = c("SAT_Label_2000.Rds")
file.type = "Subtype.metadata.Rds"; meta.group = c("Annotation"); split.group = "Treatment_Group"
for (seurat.file in seurat.files) {
  rmarkdown::render(input= paste0(Disk, Project.folder, "/", "Plotting_Figure_Level1_umap.R"), output_format= "html_document",
                    output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                            paste0('Plotting_Figure_Level1_umap_', stringr::str_split(seurat.file, ".Rds")[[1]][1],"_", 
                                                   paste0(meta.group), '.html')))
}
rm(seurat.files, meta.group, file.type)



seurat.files = c("SAT_Label_2000_FAP.Rds", "SAT_Label_2000_Endothelial.Rds",
                 "SAT_Label_2000_Myeloid.Rds", "SAT_Label_2000_Adipocytes.Rds")
file.type = "Subtype.metadata.Rds"; meta.group = "Subtype"; split.group = "Treatment_Group"
for (seurat.file in seurat.files) {
  rmarkdown::render(input= paste0(Disk,Project.folder, "/", "Plotting_Figure_Level1_umap.R"), output_format= "html_document",
                    output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                            paste0('Plotting_Figure_Level1_umap_', stringr::str_split(seurat.file, ".Rds")[[1]][1],"_", 
                                                   paste0(meta.group), '.html')))
}
rm(seurat.files, meta.group, file.type)

seurat.files = c(
 # "SAT_Label_2000.Rds",
                 "SAT_Label_2000_FAP.Rds", "SAT_Label_2000_Endothelial.Rds",
                 "SAT_Label_2000_Myeloid.Rds", "SAT_Label_2000_Adipocytes.Rds"
                 )
file.type = "Senescence.name.metadata.Rds"; meta.group = "sene.name"; split.group = "Treatment_Group"
for (seurat.file in seurat.files) {
  rmarkdown::render(input= paste0(Disk,Project.folder, "/", "Plotting_Figure_Level1_umap.R"), output_format= "html_document",
                    output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                            paste0('Plotting_Figure_Level1_umap_', stringr::str_split(seurat.file, ".Rds")[[1]][1],"_", 
                                                   paste0(meta.group), '.html')))
}
rm(seurat.files, meta.group, file.type)


seurat.files = c(
  "SAT_Label_2000_Myeloid_PVM.Rds"
  , "SAT_Label_2000_Myeloid_NKT.Rds"
  , "SAT_Label_2000_Myeloid_CEM1.Rds"
  , "SAT_Label_2000_Myeloid_DCs_B_Cells.Rds"
  , "SAT_Label_2000_Myeloid_Mast_cells.Rds"
  , "SAT_Label_2000_Myeloid_LAM.Rds"
  , "SAT_Label_2000_Myeloid_CEM2.Rds"
  , "SAT_Label_2000_Myeloid_P_LAM.Rds"
  
  , "SAT_Label_2000_FAP_hFAP6.Rds"
  , "SAT_Label_2000_FAP_hFAP3.Rds"
  , "SAT_Label_2000_FAP_Pre_adip1.Rds"
  , "SAT_Label_2000_FAP_CPA.Rds"
  , "SAT_Label_2000_FAP_hFAP2.Rds"
  , "SAT_Label_2000_FAP_hFAP4.Rds"
  , "SAT_Label_2000_FAP_hFAP7.Rds"
  , "SAT_Label_2000_FAP_Pre_adip2.Rds"
  , "SAT_Label_2000_FAP_hFAP1.Rds"
  , "SAT_Label_2000_FAP_Inflammatory_hFAP.Rds"
  , "SAT_Label_2000_FAP_hFAP5.Rds"
  
  , "SAT_Label_2000_Adipocytes_hAd07.Rds"
  , "SAT_Label_2000_Adipocytes_hAd03.Rds"
  , "SAT_Label_2000_Adipocytes_hAd09.Rds"
  , "SAT_Label_2000_Adipocytes_hAd04.Rds"
  , "SAT_Label_2000_Adipocytes_hFAP_Ad2.Rds"
  , "SAT_Label_2000_Adipocytes_hAd05.Rds"
  , "SAT_Label_2000_Adipocytes_hAd02.Rds"
  , "SAT_Label_2000_Adipocytes_hFAP_Ad1.Rds"
  , "SAT_Label_2000_Adipocytes_hAd10.Rds"
  , "SAT_Label_2000_Adipocytes_hFAP_Ad3.Rds"
  , "SAT_Label_2000_Adipocytes_hAd01.Rds"
  , "SAT_Label_2000_Adipocytes_hAd_NK.Rds"
  , "SAT_Label_2000_Adipocytes_hAd06.Rds"
  , "SAT_Label_2000_Adipocytes_hAd08.Rds"
  , "SAT_Label_2000_Adipocytes_hAd_Macrophage.Rds"
  , "SAT_Label_2000_Adipocytes_hFAP_Ad4.Rds"
  
  , "SAT_Label_2000_Endothelial_MSC.Rds"
  , "SAT_Label_2000_Endothelial_EC.Rds"
  , "SAT_Label_2000_Endothelial_MC.Rds"
  , "SAT_Label_2000_Endothelial_Inflammatory_MC.Rds"
  , "SAT_Label_2000_Endothelial_VEC.Rds"
  , "SAT_Label_2000_Endothelial_SMC.Rds"
  , "SAT_Label_2000_Endothelial_Fibromyocyte.Rds"
  , "SAT_Label_2000_Endothelial_Fibroblast.Rds"
  , "SAT_Label_2000_Endothelial_Inflammatory_EC.Rds"
  , "SAT_Label_2000_Endothelial_Inflammatory_SMC.Rds"
)
sub.folder <- "SAT_Label_2000_Subtype"
file.type = "Senescence.name.metadata.Rds"; meta.group = "sene.name"; split.group = "Treatment_Group"
for (seurat.file in seurat.files) {
  rmarkdown::render(input= paste0(Disk,Project.folder, "/", "Plotting_Figure_Level1_umap.R"), output_format= "html_document",
                    output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                            paste0('Plotting_Figure_Level1_umap_', stringr::str_split(seurat.file, ".Rds")[[1]][1],"_", 
                                                   paste0(meta.group), '.html')))
}
rm(seurat.files, meta.group, file.type, sub.folder)



# 15. Markers to pathways analysis and plotting -----------------------------------------------
markers.files = c(
    "FindMarkers_SAT_Label_2000/sene.name.Markers/Combined_SAT_Label_2000_Myeloid_Senescence.name.metadata.Rds"
  #, "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_FAP_Senescence.name.metadata.Rds"
  #, "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Adipocytes_Senescence.name.metadata.Rds"
  #, "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Endothelial_Senescence.name.metadata.Rds"
); 
organism = "Homo sapiens" # organism = "Mus musculus"
rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Main_07_MarkerstoPathway_Analysis.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                          paste0('Analysis_MarkerstoPathway_', "SAT_Label_2000_Myeloid",  '.html')))
rm(markers.files, organism)


pathways.files = c(
  "Pathways/Pathway_Combined_SAT_Label_2000_Myeloid_Senescence.name.metadata.Rds"
 )
rmarkdown::render(input= paste0(Disk, Project.folder, "/", "Plotting_Figure_Analysis_Pathways.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                          paste0('Plotting_Figure_Analysis_pathways',  '.html')))
rm(pathways.files)


## level4
markers.files = c(
  #    "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Myeloid_DCs_B_Cells_Senescence.name.metadata.level4.ombined.Rds"
  # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Myeloid_PVM_Senescence.name.metadata.level4.combined.Rds"
  # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Myeloid_NKT_Senescence.name.metadata.level4.combined.Rds"
  # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Myeloid_CEM1_Senescence.name.metadata.level4.combined.Rds"
  # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Myeloid_Mast_cells_Senescence.name.metadata.level4.combined.Rds"
  # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Myeloid_LAM_Senescence.name.metadata.level4.combined.Rds"
  # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Myeloid_CEM2_Senescence.name.metadata.level4.combined.Rds"
  # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Myeloid_P_LAM_Senescence.name.metadata.level4.combined.Rds"
  # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_FAP_hFAP6_Senescence.name.metadata.level4.combined.Rds"
  # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_FAP_hFAP3_Senescence.name.metadata.level4.combined.Rds"
  # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_FAP_Pre_adip1_Senescence.name.metadata.level4.combined.Rds"
  # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_FAP_CPA_Senescence.name.metadata.level4.combined.Rds"
  # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_FAP_hFAP2_Senescence.name.metadata.level4.combined.Rds"
  # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_FAP_hFAP4_Senescence.name.metadata.level4.combined.Rds"
  # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_FAP_hFAP7_Senescence.name.metadata.level4.combined.Rds"
  # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_FAP_Pre_adip2_Senescence.name.metadata.level4.combined.Rds"
  # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_FAP_hFAP1_Senescence.name.metadata.level4.combined.Rds"
  # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_FAP_Inflammatory_hFAP_Senescence.name.metadata.level4.combined.Rds"
  # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_FAP_hFAP5_Senescence.name.metadata.level4.combined.Rds"
  # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Adipocytes_hAd07_Senescence.name.metadata.level4.combined.Rds"
  # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Adipocytes_hAd03_Senescence.name.metadata.level4.combined.Rds"
  # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Adipocytes_hAd09_Senescence.name.metadata.level4.combined.Rds"
  # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Adipocytes_hAd04_Senescence.name.metadata.level4.combined.Rds"
  # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Adipocytes_hFAP_Ad2_Senescence.name.metadata.level4.combined.Rds"
  # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Adipocytes_hAd05_Senescence.name.metadata.level4.combined.Rds"
  # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Adipocytes_hAd02_Senescence.name.metadata.level4.combined.Rds"
  # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Adipocytes_hFAP_Ad1_Senescence.name.metadata.level4.combined.Rds"
  # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Adipocytes_hAd10_Senescence.name.metadata.level4.combined.Rds"
  # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Adipocytes_hFAP_Ad3_Senescence.name.metadata.level4.combined.Rds"
  # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Adipocytes_hAd01_Senescence.name.metadata.level4.combined.Rds"
  # ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Adipocytes_hAd_NK_Senescence.name.metadata.level4.combined.Rds"
    "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Adipocytes_hAd06_Senescence.name.metadata.level4.combined.Rds"
  ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Adipocytes_hAd08_Senescence.name.metadata.level4.combined.Rds"
  ,  "FindMarkers_SAT_Label_2000/sene.name.Markers/SAT_Label_2000_Adipocytes_hAd_Macrophage_Senescence.name.metadata.level4.combined.Rds"
  ); 
organism = "Homo sapiens" # organism = "Mus musculus"
type = "avg_log2FC."
rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Main_07_MarkersLevel4_toPathway_Analysis.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                          paste0('Analysis_MarkersLevel4_toPathway_', "SAT_Label_2000_Myeloid", '.html')))
rm(markers.files, organism, type)

pathways.files = c(
  # ,  "Pathways/Pathway_avg_log2FC.SAT_Label_2000_Myeloid_DCs_B_Cells_Senescence.name.metadata.level4.combined.Rds"
  # , "Pathways/Pathway_avg_log2FC.SAT_Label_2000_Myeloid_NKT_Senescence.name.metadata.level4.combined.Rds"
  # , "Pathways/Pathway_avg_log2FC.SAT_Label_2000_Myeloid_PVM_Senescence.name.metadata.level4.combined.Rds"
  # ,   "Pathways/Pathway_avg_log2FC.SAT_Label_2000_Myeloid_CEM2_Senescence.name.metadata.level4.combined.Rds"
  # ,  "Pathways/Pathway_avg_log2FC.SAT_Label_2000_Myeloid_P_LAM_Senescence.name.metadata.level4.combined.Rds"
  # ,  "Pathways/Pathway_avg_log2FC.SAT_Label_2000_Myeloid_CEM1_Senescence.name.metadata.level4.combined.Rds"
  # ,  "Pathways/Pathway_avg_log2FC.SAT_Label_2000_Myeloid_Mast_cells_Senescence.name.metadata.level4.combined.Rds"
  # ,  "Pathways/Pathway_avg_log2FC.SAT_Label_2000_Myeloid_LAM_Senescence.name.metadata.level4.combined.Rds"
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
rmarkdown::render(input= paste0(Disk, Project.folder, "/", "Plotting_Figure_Analysis_Pathways.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                          paste0('Plotting_Figure_Analysis_pathways_Level4', "SAT_Label_2000_Myeloid",  '.html')))
rm(pathways.files)


# 16. CellChat anaylsis -------------------
seurat.files = c(
  # "SAT_Label_2000.Rds",
  "SAT_Label_2000_Myeloid.Rds"
  # , "SAT_Label_2000_FAP.Rds"
  # , "SAT_Label_2000_Adipocytes.Rds"
  # , "SAT_Label_2000_Endothelial.Rds"
)
sample.file.type = "SAT_Sample.meta.data.Rds"; # sample files
meta.group2.1 = "Group"; # subset seurat.data, 
meta.group2.2 = "Event.Name" # compare group, levels
file.type = "Subtype.metadata.Rds"; # with Annotation.file.name # subset seurat.data
meta.group = "Subtype" # subset seurat.data
for (seurat.file in seurat.files) {
  rmarkdown::render(input= paste0(Disk,Project.folder, "/", "Plotting_Figure_Level1_umap.R"), output_format= "html_document",
                    output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                            paste0('Plotting_Figure_Level1_umap_', stringr::str_split(seurat.file, ".Rds")[[1]][1],"_", 
                                                   paste0(meta.group), '.html')))
}
rm(seurat.files, meta.group, file.type)

rmarkdown::render(input= paste0(Disk,"00_Functions_Refs", "/FunctionCodes/", "Analysis_FindMarkers_seurat_level4.R"), output_format= "html_document",
                  output_file = file.path(paste0(Disk, Project.folder, "/html.result/"), 
                                          paste0('Analysis_FindMarkers_level4_', 
                                                 stringr::str_split(seurat.files[1], "_")[[1]][1],"_", 
                                                 meta.group,'.html')))
rm(seurat.files, sample.file.type, meta.group1, meta.group2.1, meta.group2.2, file.type, meta.group)
