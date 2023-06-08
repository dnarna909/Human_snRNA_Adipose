

sce.list <- readRDS("/media/jianie/Extreme SSD1/2022-09-01 STARR_SGLT2 Combine/Rds files_Aggr_all/Final.sce.list.Rds")
# saveRDS(sce.list[1:16], file = paste0(Disk, Project.folder, "/",  "Rds files_Aggr_all", "/Final.sce.list_part1.Rds"))
saveRDS(sce.list[17:32], file = paste0(Disk, Project.folder, "/",  "Rds files_Aggr_all", "/Final.sce.list_part2.Rds"))
saveRDS(sce.list[33:55], file = paste0(Disk, Project.folder, "/",  "Rds files_Aggr_all", "/Final.sce.list_part3.Rds"))

for (ii in names(sce.list)){
  saveRDS(sce.list[[ii]], file = paste0(Disk, Project.folder, "/",  "Rds files_Aggr_all/", "sce.list","/", ii,".Rds"))
  print(ii)
}
