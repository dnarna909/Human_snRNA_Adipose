
Result <- readRDS("/media/jianie/Extreme SSD/2022-09-01 STARR_SGLT2 Combine/Rds files_Aggr_all/FindMarkers_SAT_Label_2000/Subtype.Level4.Markers/SAT_Label_2000_Endothelial_Subtype.metadata.level4.combined.Rds")
cell.id = "EC4"
cell.id = "Inflammatory_EC"
tt <- Result[[paste0(cell.id,"_SGLT2i_Post:Pre")]][[paste0(cell.id,"_avg_log2FC.SGLT2i:NUT")]]

tt <- Result[[paste0(cell.id,"_SGLT2i_Post:Pre")]][[paste0(cell.id,"_avg_log2FC.SGLT2i:NUT")]] %>% dplyr::filter(Direction == "Up")

library("AnnotationDbi")
library("org.Hs.eg.db")
columns(org.Hs.eg.db)

tt$ENSEMBL = mapIds(org.Hs.eg.db,
                    keys=row.names(tt), 
                    column="ENSEMBL" ,
                    keytype="SYMBOL",
                    multiVals="first")
tt$name =   mapIds(org.Hs.eg.db,
                    keys=row.names(tt), 
                    column="GENENAME",
                    keytype="SYMBOL",
                    multiVals="first")
tt$entrez = mapIds(org.Hs.eg.db,
                   keys=row.names(tt), 
                   column="ENTREZID",
                   keytype="SYMBOL",
                   multiVals="first")
head(tt, 10)


library(pathview)
library(gage)
library(gageData)
data(kegg.sets.hs)
data(sigmet.idx.hs)

kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]
head(kegg.sets.hs, 3)

foldchanges = tt$mean.x
names(foldchanges) = tt$entrez
head(foldchanges)

# Get the results
keggres = gage(foldchanges, gsets=kegg.sets.hs, same.dir=TRUE)
head(keggres$greater)

pathways = data.frame(id=rownames(keggres$greater), keggres$greater)
head(pathways)
library(dplyr)
pathways.sig <- pathways %>% dplyr::filter(p.val <= 0.05)


sim.cpd.data=sim.mol.data(mol.type="cpd", nmol=3000) # compound nodes
data(cpd.simtypes)

export.folder0 <- paste0("Pathways")
dir.create(file.path(paste0(Disk, Project.folder, "/", figures.folder, "/"), export.folder0), showWarnings = FALSE)
dir <- paste0(Disk, Project.folder, "/", figures.folder, "/", export.folder0, "/")

setwd(dir)
for (i in 1:nrow(pathways.sig)){
  pathview(gene.data = foldchanges , pathway.id = stringr::str_split(rownames(pathways.sig[i, ]), " ")[[1]][1],
           out.suffix = cell.id,  kegg.native = T)# native original KEGG graph, png
  pathview(gene.data = foldchanges , pathway.id = stringr::str_split(rownames(pathways.sig[i, ]), " ")[[1]][1],
           out.suffix = cell.id,  kegg.native = F) # Graphviz engine, pdf
  pathview(gene.data = foldchanges , pathway.id = stringr::str_split(rownames(pathways.sig[i, ]), " ")[[1]][1],
           out.suffix = cell.id,  kegg.native = F, same.layer = F) # Graphviz engine, pdf, page 1 is the main graph, page 2 is the legen
  pathview(gene.data = foldchanges , pathway.id = stringr::str_split(rownames(pathways.sig[i, ]), " ")[[1]][1],
           out.suffix = cell.id,  kegg.native = F, split.group = T) # # Graphviz engine, pdf, split the node groups into individual detached nodes
  pathview(gene.data = foldchanges , pathway.id = stringr::str_split(rownames(pathways.sig[i, ]), " ")[[1]][1],
           out.suffix = cell.id,  kegg.native = F, split.group = T, expand.node = T) # Graphviz engine, pdf, expand the multiple-gene nodes into individual genes
  pv.out <- pathview(gene.data = foldchanges , pathway.id = stringr::str_split(rownames(pathways.sig[i, ]), " ")[[1]][1]
           , cpd.data = sim.cpd.data,
           out.suffix = cell.id,  keys.align = "y", kegg.native = T) # native original KEGG graph, png, Compound and gene data
}



# --------------------
tt <- Result[["EC4_SGLT2i_Post:Pre"]][["EC4_SGLT2i_Post:Pre_data"]]
tt$ENSEMBL = mapIds(org.Hs.eg.db,
                    keys=tt$rowname, 
                    column="ENSEMBL" ,
                    keytype="SYMBOL",
                    multiVals="first")
tt$name =   mapIds(org.Hs.eg.db,
                   keys=tt$rowname, 
                   column="GENENAME",
                   keytype="SYMBOL",
                   multiVals="first")
tt$entrez = mapIds(org.Hs.eg.db,
                   keys=tt$rowname, 
                   column="ENTREZID",
                   keytype="SYMBOL",
                   multiVals="first")
head(tt, 10)
tt <- tt[!is.na(tt$entrez), ] 
grep('SGLT2i', colnames(tt), ignore.case =TRUE)
dd1 <- tt[, c(grep('entrez', colnames(tt), ignore.case =TRUE), grep('avg_log2FC', colnames(tt), ignore.case =TRUE))]

tt <- Result[["EC4_SGLT2i_Post:Pre"]][["EC4_NUT_Post:Pre_data"]]
tt$ENSEMBL = mapIds(org.Hs.eg.db,
                    keys=tt$rowname, 
                    column="ENSEMBL" ,
                    keytype="SYMBOL",
                    multiVals="first")
tt$name =   mapIds(org.Hs.eg.db,
                   keys=tt$rowname, 
                   column="GENENAME",
                   keytype="SYMBOL",
                   multiVals="first")
tt$entrez = mapIds(org.Hs.eg.db,
                   keys=tt$rowname, 
                   column="ENTREZID",
                   keytype="SYMBOL",
                   multiVals="first")
head(tt, 10)
tt <- tt[!is.na(tt$entrez), ] 
dd2 <- tt[, c(grep('entrez', colnames(tt), ignore.case =TRUE), grep('avg_log2FC', colnames(tt), ignore.case =TRUE))] 

dd <- dd1 %>% full_join(dd2, by = "entrez")%>% tibble::column_to_rownames(var = "entrez") %>% as.matrix()
keggres = gage(dd, gsets=kegg.sets.hs, same.dir=TRUE, ref =  grep('NUT', colnames(dd), ignore.case =TRUE), compare='1ongroup' )# grep('mean.y', colnames(Result[[nn]][[nn1]]), ignore.case =TRUE)
# keggres = gage(dd, gsets=kegg.sets.hs, same.dir=TRUE, ref =  grep('NUT', colnames(dd), ignore.case =TRUE), compare='as.group' )# grep('mean.y', colnames(Result[[nn]][[nn1]]), ignore.case =TRUE)
head(keggres$greater)
pathways = data.frame(id=rownames(keggres$greater), keggres$greater)
head(pathways)
pathview(gene.data = dd , pathway.id = stringr::str_split(rownames(pathways[6, ]), " ")[[1]][1],
         out.suffix = "test",  kegg.native = T,  multi.state = T, same.layer = T)# native original KEGG graph, png


df <- data.frame(SGLT = rowSums(dd[, grep('SGLT2i', colnames(dd), ignore.case =TRUE)] )/ length(grep('SGLT2i', colnames(dd), ignore.case =TRUE)),
NUT = rowSums(dd[, grep('NUT', colnames(dd), ignore.case =TRUE)] )/length(grep('NUT', colnames(dd), ignore.case =TRUE)) )

df <- data.frame(SGLT = rowMeans(dd[, grep('SGLT2i', colnames(dd), ignore.case =TRUE)] ),
                 NUT = rowMeans(dd[, grep('NUT', colnames(dd), ignore.case =TRUE)] ))
pathview(gene.data = df , pathway.id = stringr::str_split(rownames(pathways[6, ]), " ")[[1]][1],
         out.suffix = "test",  kegg.native = T,  multi.state = T, same.layer = T)# native original KEGG graph, png
