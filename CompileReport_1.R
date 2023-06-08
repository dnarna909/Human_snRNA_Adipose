# /media/jianie/Extreme SSD1/ 
# D:/
rmarkdown::render(input= paste0("/media/jianie/Extreme SSD1/00_Functions_Refs/Single cell altas of human adipose tissue/2.Calculate Markers.r"), output_format= "html_document") 

#' creat `Project Parameters.R` files for each project
##' counts foler should be `ProjectName_aggr`
#' creat `CompileParameters.R` files for each project
source(paste0("/media/jianie/Extreme SSD1/2022-09-01 STARR_SGLT2 Combine/CompileParameters.R"), local = knitr::knit_global())

# prerun ----------------------------------------------------------------------------------------------------------
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_00_prerun_Creat_matrices.R"), output_format= "html_document") # 
#' generate: combine `project_mat.rds`,feature.names.rds,  only need run once.
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_00_prerun_QC_Aggr_all.genes.R"), output_format= "html_document") # 
#' generate: individual `sample_mat.Rds`, combined `dataset.group.Rds`, combined `mat.Rds`(list), in the counts folder.only need run once.

# QC ----------------------------------------------------------------------------------------------------------
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_00_QC_Aggr_all.genes_step1.R"), output_format= "html_document") # 
#' QC, filter, and generate: `Ambient.Rds`, `QC files_Genes.RData`, `QC files_Cells.RData`, `Final.sce.list.Rds`, `rescaled.Rds`.only need run once.
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_00_QC_Aggr_all.genes_step2.R"), output_format= "html_document") # 
#' generate: individual `seurat.list.rds`.only need run once.
rmarkdown::render(input= paste0(Disk, Project.folder, "/","Main_00_Step1_Combine_seurat_data.R"), output_format= "html_document") # 
#' combined `SAT1.Rds`, combined `SAT2.Rds`, or combined `SAT.Rds` (depending on sample number).only need run once.
rmarkdown::render(input= paste0(Disk, "00_Functions_Refs", "/FunctionCodes/", "Main_00_QC_Aggr_all.genes_step3.R"), output_format= "html_document") # 
#' import `SAT.Rds` for seurat analysis. generate: `meta.data.Rds`, `SAT.Rds`.

#' check `Main_00_Step2_Test_Resolution.R` files for best resolution and variable gene number


# cell type analysis -------------------------------------------------------------------------------------------------
# Overall
# source(paste0("/media/jianie/Extreme SSD1/2022-09-01 STARR_SGLT2 Combine/CompileParameters.R"), local = knitr::knit_global())
# rmarkdown::render(input= paste0(Disk, Project.folder, "/","Main_01_Overall_Fat_1_Step1.R"), output_format= "html_document") # 10∶26∶54 PM - 11∶14∶30 PM 
source(paste0("/media/jianie/Extreme SSD1/2022-09-01 STARR_SGLT2 Combine/CompileParameters.R"), local = knitr::knit_global())
rmarkdown::render(input= paste0(Disk, Project.folder, "/","Main_01_Overall_Fat_1_Step2_for.annotated.rmd"), output_format= "html_document") # 11∶17∶45 PM

# Overall
source(paste0("/media/jianie/Extreme SSD1/2022-09-01 STARR_SGLT2 Combine/CompileParameters.R"), local = knitr::knit_global()) 
rmarkdown::render(input= paste0(Disk, Project.folder, "/","Main_01_Overall_Fat_1_Step3_name.R"), output_format= "html_document")# 11∶24∶58 PM 

# cell subtype analysis -------------------------------------------------------------------------------------------------
# FAP
source(paste0("/media/jianie/Extreme SSD1/2022-09-01 STARR_SGLT2 Combine/CompileParameters.R"), local = knitr::knit_global())
rmarkdown::render(input= paste0(Disk, Project.folder, "/","Main_04_FAP_Fat_1_Step1.rmd"), output_format= "html_document") # 11∶32∶10 PM 
source(paste0("/media/jianie/Extreme SSD1/2022-09-01 STARR_SGLT2 Combine/CompileParameters.R"), local = knitr::knit_global())
rmarkdown::render(input= paste0(Disk, Project.folder, "/","Main_04_FAP_Fat_1_Step2_for.annotated.rmd"), output_format= "html_document") # 11∶36∶24 PM 

# Adipocytes
source(paste0("/media/jianie/Extreme SSD1/2022-09-01 STARR_SGLT2 Combine/CompileParameters.R"), local = knitr::knit_global())
rmarkdown::render(input= paste0(Disk, Project.folder, "/","Main_05_Adipocytes_Fat_1_Step1.rmd"), output_format= "html_document") # 11∶59∶03 PM
source(paste0("/media/jianie/Extreme SSD1/2022-09-01 STARR_SGLT2 Combine/CompileParameters.R"), local = knitr::knit_global())
rmarkdown::render(input= paste0(Disk, Project.folder, "/","Main_05_Adipocytes_Fat_1_Step2_for.annotated.rmd"), output_format= "html_document") # 12∶07∶55 AM 

# Myeloid
source(paste0("/media/jianie/Extreme SSD1/2022-09-01 STARR_SGLT2 Combine/CompileParameters.R"), local = knitr::knit_global())
rmarkdown::render(input= paste0(Disk, Project.folder, "/","Main_02_Myeloid_Fat_1_Step1.rmd"), output_format= "html_document") # 12∶22∶38 AM
source(paste0("/media/jianie/Extreme SSD1/2022-09-01 STARR_SGLT2 Combine/CompileParameters.R"), local = knitr::knit_global())
rmarkdown::render(input= paste0(Disk, Project.folder, "/","Main_02_Myeloid_Fat_1_Step2_for.annotated.rmd"), output_format= "html_document") # 12∶25∶47 AM

# Endothelial
source(paste0("/media/jianie/Extreme SSD1/2022-09-01 STARR_SGLT2 Combine/CompileParameters.R"), local = knitr::knit_global())
rmarkdown::render(input= paste0(Disk, Project.folder, "/","Main_03_Endothelial_Fat_1_Step1.rmd"), output_format= "html_document") # 12∶36∶31 AM
source(paste0("/media/jianie/Extreme SSD1/2022-09-01 STARR_SGLT2 Combine/CompileParameters.R"), local = knitr::knit_global())
rmarkdown::render(input= paste0(Disk, Project.folder, "/","Main_03_Endothelial_Fat_1_Step2_for.annotated.rmd"), output_format= "html_document") #  12∶38∶06 AM

# # FAP and Adipocytes combine
# source(paste0("/media/jianie/Extreme SSD1/2022-09-01 STARR_SGLT2 Combine/CompileParameters.R"), local = knitr::knit_global())
# rmarkdown::render(input= paste0(Disk, Project.folder, "/","Main_05_FAP_Adipocytes_Fat_1_Step1.rmd"), output_format= "html_document") 
# source(paste0("/media/jianie/Extreme SSD1/2022-09-01 STARR_SGLT2 Combine/CompileParameters.R"), local = knitr::knit_global())
# rmarkdown::render(input= paste0(Disk, Project.folder, "/","Main_05_FAP_Adipocytes_Fat_1_Step2_for.annotated.rmd"), output_format= "html_document") 


# add name to cell subtype  ------------------------------------------------------------------------
# FAP
source(paste0("/media/jianie/Extreme SSD1/2022-09-01 STARR_SGLT2 Combine/CompileParameters.R"), local = knitr::knit_global())
rmarkdown::render(input= paste0(Disk, Project.folder, "/","Main_04_FAP_Fat_1_Step3_name.R"), output_format= "html_document") # 11∶45∶46 PM
# Adipocytes
source(paste0("/media/jianie/Extreme SSD1/2022-09-01 STARR_SGLT2 Combine/CompileParameters.R"), local = knitr::knit_global()) 
rmarkdown::render(input= paste0(Disk, Project.folder, "/","Main_05_Adipocytes_Fat_1_Step3_name.R"), output_format= "html_document") # 12∶18∶28 AM 
# Myeloid
source(paste0("/media/jianie/Extreme SSD1/2022-09-01 STARR_SGLT2 Combine/CompileParameters.R"), local = knitr::knit_global())
rmarkdown::render(input= paste0(Disk, Project.folder, "/","Main_02_Myeloid_Fat_1_Step3_name.R"), output_format= "html_document") # 12∶33∶55 AM
# Endothelial
source(paste0("/media/jianie/Extreme SSD1/2022-09-01 STARR_SGLT2 Combine/CompileParameters.R"), local = knitr::knit_global())
# rmarkdown::render(input= paste0(Disk, Project.folder, "/","Main_03_Endothelial_Fat_1_Step3_name.R"), output_format= "html_document")
# # FAP and Adipocytes combine
# source(paste0("/media/jianie/Extreme SSD1/2022-09-01 STARR_SGLT2 Combine/CompileParameters.R"), local = knitr::knit_global()) 
# rmarkdown::render(input= paste0(Disk, Project.folder, "/","Main_05_FAP_Adipocytes_Fat_1_Step3_name.R"), output_format= "html_document") 


source(paste0("/media/jianie/Extreme SSD1/2022-09-01 STARR_SGLT2 Combine/CompileParameters.R"), local = knitr::knit_global())
rmarkdown::render(input= paste0(Disk, Project.folder, "/","Main_04_FAP_Fat_1_Step3_name.R"), output_format= "html_document") # 11∶45∶46 PM
source(paste0("/media/jianie/Extreme SSD1/2022-09-01 STARR_SGLT2 Combine/CompileParameters.R"), local = knitr::knit_global()) 
rmarkdown::render(input= paste0(Disk, Project.folder, "/","Main_05_Adipocytes_Fat_1_Step3_name.R"), output_format= "html_document") # 12∶18∶28 AM 
source(paste0("/media/jianie/Extreme SSD1/2022-09-01 STARR_SGLT2 Combine/CompileParameters.R"), local = knitr::knit_global())
rmarkdown::render(input= paste0(Disk, Project.folder, "/","Main_02_Myeloid_Fat_1_Step3_name.R"), output_format= "html_document") # 12∶33∶55 AM
