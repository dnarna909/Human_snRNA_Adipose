
# set parameters for QC
# Disk <- c("/media/jianie/Extreme SSD1/") 
# Disk <- c("D:/")
# counts.folder <- c("01_STARR")
counts.folder <- c("02_SGLT2")
counts.subfolder <- c("STARR-SGLT2-1_aggr_lowest_16", "STARR-SGLT2-2_aggr_lowest_17_32", "SGLT2_Final_aggr")
# filefolder.aggr = stringr::str_split(counts.subfolder[1], "_")[[1]][1] # substr(counts.subfolder, 0, 7)  #don't change #  check matrices.All file
# Project.folder <- c("2022-09-01 STARR_SGLT2 Combine")
Rds.folder <- c("Rds files_Aggr_all")
figures.folder <- c("Figures")
species.ensembl <- "hsapiens" # Human genes, changes in line 99
# species.ensembl <- "mmusculus" # Mouse genes, changes in line 99
mt.genes <- "^Mt-" # Human genes
# mt.genes <- "^mt-" # Mouse genes


# set parameters for Aggregation Information
# aggr file:
sample.info.disk1 <- "/media/jianie/DATA/"
# sample.info.disk1 <- "C:/Users/niej/Documents/"
sample.info.path1 <- "UTHSC_cellranger_codes_files/AggrFiles/"
sample.file.name1 <- c("STARR.SGLT2.1_aggr_lowest_16.rds", "STARR.SGLT2.2_aggr_lowest_17_32.rds", "SGLT2_Final.rds") # from aggregate # need to add "Dataset", "List", "sample_id

# samples information by All Sample Information
sample.info.disk2 <- "/media/jianie/DATA/"
# sample.info.disk2 <- "C:/Users/niej/Documents/"
sample.info.path2 <- "UTHSC_cellranger_codes_files/RedCap_SubjectsInfo/"
sample.file.name2 <- c("STARR.all.subject.info.rds", "SGLT2.all.subject.info.rds")  # contain "Dataset" (Sample01)

# # samples information by researcher
# sample.info.disk3 <- Disk
# sample.info.path3 <- paste0(Project.folder, "/")
# sample.file.name3 <- "MarmosetIslet.csv" # contain "Dataset" (Sample01)

# Load libraries 
# packageList <- scan( paste0(Disk, "00_Functions_Refs/SingleCellPreprocessPackageList.txt") , character(), quote = "")
# for(package in packageList){
#   if(!require(package,character.only = TRUE)){
#     install.packages(package);require(package,character.only = TRUE);}
# }
# # library("velocyto.R") # linux -specific
