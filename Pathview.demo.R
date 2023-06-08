# https://pathview.r-forge.r-project.org/pathview.pdf


library(pathview)

## load example data, and check
data(gse16873.d)
gse16873.d[1:4,1:4]
# DCIS_1      DCIS_2      DCIS_3      DCIS_4
# 10000 -0.3076448 -0.14722769 -0.02378481 -0.07056193
# 10001  0.4158680 -0.33477259 -0.51313691 -0.16653712
# 10002  0.1985493  0.03789588  0.34186534 -0.08527420
# 10003 -0.2315530 -0.09659311 -0.10472728 -0.04801404
#  
## load example pathways, and check
data(demo.paths)
demo.paths
# $sel.paths
# # [1] "04110" "00620" "00640"
# 
# $kpos1
# # [1] "topright"    "bottomright" "topright"   
# 
# $kpos2
# # [1] "topright"    "topright"    "bottomright"
# 
# $spos
# # [1] "bottomleft" "bottomleft" "topright"  
# 
# $offs
# # [1] -1.0 -1.0 -0.8


## plot expression data of first sample
## on the first pathway (="04110")
i <- 1

## create actual plot. png will be saved in working directory
pv.out <- pathview(gene.data = gse16873.d[, 1],    ## plot expression data from first sample
                   pathway.id = demo.paths$sel.paths[i],    ## plot on the first ("04110") of the 3 demo pathways
                   species = "hsa",    ## hsa = KEGG identifier for humanm
                   out.suffix = "gse16873",    ## add this suffix to name of png file
                   kegg.native = TRUE)    ## overlay expression data on KEGG-like pathway (and do not create a graphviz object) 
# 'select()' returned 1:1 mapping between keys and columns
# Info: Working in directory E:/test
# Info: Writing image file hsa04110.gse16873.png

pv.out <- pathview(gene.data = gse16873.d[, 1], pathway.id = demo.paths$sel.paths[i],
                    species = "hsa", out.suffix = "gse16873", kegg.native = F,
                    sign.pos = demo.paths$spos[i])

str(pv.out)
# List of 2
# $ plot.data.gene:'data.frame': 92 obs. of  10 variables:
#   ..$ kegg.names: chr [1:92] "1029" "51343" "4171" "4998" ...
# ..$ labels    : chr [1:92] "CDKN2A" "FZR1" "MCM2" "ORC1" ...
# ..$ all.mapped: chr [1:92] "1029" "51343" "4171,4172,4173,4174,4175,4176" "4998,4999,5000,5001,23594,23595" ...
# ..$ type      : chr [1:92] "gene" "gene" "gene" "gene" ...
# ..$ x         : num [1:92] 532 919 553 494 919 919 188 432 123 77 ...
# ..$ y         : num [1:92] 124 536 556 556 297 519 519 191 704 687 ...
# ..$ width     : num [1:92] 46 46 46 46 46 46 46 46 46 46 ...
# ..$ height    : num [1:92] 17 17 17 17 17 17 17 17 17 17 ...
# ..$ mol.data  : num [1:92] 0.129 -0.404 -0.42 0.986 0.936 ...
# ..$ mol.col   : chr [1:92] "#BEBEBE" "#5FDF5F" "#5FDF5F" "#FF0000" ...
# $ plot.data.cpd : NULL

head(pv.out$plot.data.gene)
# kegg.names labels                                  all.mapped type   x   y
# 4       1029 CDKN2A                                        1029 gene 532 124
# 5      51343   FZR1                                       51343 gene 919 536
# 6       4171   MCM2               4171,4172,4173,4174,4175,4176 gene 553 556
# 7       4998   ORC1             4998,4999,5000,5001,23594,23595 gene 494 556
# 8        996  CDC27 996,8697,8881,10393,25847,25906,29882,51433 gene 919 297
# 9        996  CDC27 996,8697,8881,10393,25847,25906,29882,51433 gene 919 519
# width height   mol.data mol.col
# 4    46     17  0.1291987 #BEBEBE
# 5    46     17 -0.4043256 #5FDF5F
# 6    46     17 -0.4202181 #5FDF5F
# 7    46     17  0.9864873 #FF0000
# 8    46     17  0.9363018 #FF0000
# 9    46     17  0.9363018 #FF0000

sim.cpd.data=sim.mol.data(mol.type="cpd", nmol=3000) # compound nodes
data(cpd.simtypes)
i=3
pv.out <- pathview(gene.data = gse16873.d[, 1], cpd.data = sim.cpd.data,
                    pathway.id = demo.paths$sel.paths[i], species = "hsa", out.suffix = "gse16873.cpd",
                    keys.align = "y", kegg.native = T, key.pos = demo.paths$kpos1[i])

pv.out <- pathview(gene.data = gse16873.d[, 1], cpd.data = sim.cpd.data,
                    pathway.id = demo.paths$sel.paths[i], species = "hsa", out.suffix = "gse16873.cpd",
                    keys.align = "y", kegg.native = F, key.pos = demo.paths$kpos2[i],
                    sign.pos = demo.paths$spos[i], cpd.lab.offset = demo.paths$offs[i])

sim.cpd.data2 = matrix(sample(sim.cpd.data, 18000,
                               replace = T), ncol = 6)
rownames(sim.cpd.data2) = names(sim.cpd.data)
 colnames(sim.cpd.data2) = paste("exp", 1:6, sep = "")
 
pv.out <- pathview(gene.data = gse16873.d[, 1:3],
                   cpd.data = sim.cpd.data2[, 1:2], pathway.id = demo.paths$sel.paths[i],
                    species = "hsa", out.suffix = "gse16873.cpd.3-2s", keys.align = "y",
                    kegg.native = T, match.data = F, multi.state = T, same.layer = T)

## version information
packageVersion("pathview")
# [1] ‘1.38.0’

sessionInfo()