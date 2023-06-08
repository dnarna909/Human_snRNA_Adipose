
colnames(Result.all[[1]])
types = c("avg_log2FC", "pct.dif")
for (ff in names(Result.all)){
  for (type in types) {
    Result.all[[ff]]  <- Result.all[[ff]] %>% 
      rename( !!paste0(type, ".cl" ) := (!!sym(paste0(ff, "_",type, ".cl" ) )))
  }
}
colnames(Result.all[[1]])
save(Result.all,
     Result.anova ,
     Result.t.test ,
     Result.t.test2 ,
     Result.summary , Result.summary3 ,
     Result ,
     Result.df_long , 
     sample.meta,
     file = paste0("/media/jianie/Extreme SSD/2022-09-01 STARR_SGLT2 Combine/Rds files_Aggr_all/FindMarkers_SAT_Label_2000/Annotation.Level2.Markers/SAT_Label_2000_Annotation_Compare_Group1_BA_Group.Mean.ANOVA.RData")
)




Result.all <- lapply(Result.all, function(x) {
  x <- x  %>%
    arrange((!!sym(paste0("anova.", types[2], ".pval" ) ))) %>%
    mutate(  (!!sym(paste0("Rank.ID.", types[2]))) := 1:n()) %>%
    arrange((!!sym(paste0("anova.", types[1], ".pval" ) ))) %>%
    mutate( (!!sym(paste0("Rank.ID.", types[1])))  := 1:n()) %>%
    mutate( (!!sym(paste0("mlog10.", types[1], ".pval")))  :=  -log10( (!!sym(paste0("anova.", types[1], ".pval" ) )) ) ) %>%
    mutate( (!!sym(paste0("mlog10.", types[2], ".pval")))  :=  -log10( (!!sym(paste0("anova.", types[2], ".pval" ) )) ) ) 
  
  x$mlog10.pval <- x[[paste0("mlog10.", types[1], ".pval")]]+x[[paste0("mlog10.", types[2], ".pval")]]
  x <- x %>%
    arrange(-(!!sym(paste0("mlog10.pval" ) ))) %>%
    mutate(  (!!sym(paste0("Rank.ID1"))) := 1:n())
})
colnames(Result.all[[1]])


colnames( Result.all[[1]])
for (ff in names(Result.all)){
  Result.all[[ff]] <- Result.all[[ff]] %>% 
    full_join( Result.summary3[[ff]][, c (grep("Direction_avg_log2FC_mean", colnames(Result.summary3[[ff]])), 
                                          grep("SetDiff.value_avg_log2FC_mean", colnames(Result.summary3[[ff]])),
                                          grep("rowname", colnames(Result.summary3[[ff]])))], by = "rowname") 
  Result.all[[ff]][["Direction"]] <- Result.all[[ff]][, c (grep("Direction_avg_log2FC_mean", colnames(Result.all[[ff]])) )]
  Result.all[[ff]][["SetDiff.value"]] <- Result.all[[ff]][, c (grep("SetDiff.value", colnames(Result.all[[ff]])) )]
  Result.all[[ff]] <- Result.all[[ff]] %>% 
    mutate(Rank.ID.SetDiff = abs(SetDiff.value))
}
rm(ff)  
colnames( Result.all[[1]])
# Result.all <- lapply(Result.all, function(x) {
#   x$Rank.ID <- x[[paste0("Rank.ID.", types[2])]]+x[[paste0("Rank.ID.", types[1])]]
#   x <- x %>%
#     arrange(Rank.ID) 
# })
Result.all <- lapply(Result.all, function(x) {
  x <- x %>%
    arrange(-(!!sym(paste0("Rank.ID.SetDiff" ) ))) %>%
    mutate(  (!!sym(paste0("Rank.ID2"))) := 1:n()) %>%
    mutate(  (!!sym(paste0("Rank.ID"))) :=  (!!sym(paste0("Rank.ID1"))) + (!!sym(paste0("Rank.ID2"))) ) %>%
    arrange( (!!sym(paste0("Rank.ID"))) ) 
})
colnames( Result.all[[1]])

#'  Gene list functional enrichment analysis with gost
#'  Enrichment analysis  
#'  Data sources and versions
# Available data sources and their abbreviations are:
#   
#   Gene Ontology (GO or by branch GO:MF, GO:BP, GO:CC)
# KEGG (KEGG)
# Reactome (REAC)
# WikiPathways (WP)
# TRANSFAC (TF)
# miRTarBase (MIRNA)
# Human Protein Atlas (HPA)
# CORUM (CORUM)
# Human phenotype ontology (HP)
# https://cran.r-project.org/web/packages/gprofiler2/vignettes/gprofiler2.html
library(gprofiler2)


genes.df <- list()
genes <- list()
for (dd in c("Up", "Down", "No")) {
  genes.df[[dd]] <-  Result.all[[1]] %>%
    dplyr::filter(Direction ==  dd,# "Down", # # "Down", "No
                  anova.avg_log2FC.pval < 0.05 & anova.pct.dif.pval < 0.05 ) %>%
    arrange(Rank.ID ) %>%
    dplyr::select(one_of(c("rowname", "Rank.ID1")))
  genes[[dd]] <- genes.df[[dd]][["Rank.ID1"]]
  names(genes[[dd]]) = genes.df[[dd]][["rowname"]]
  genes[[dd]] = sort(genes[[dd]], decreasing =  FALSE)
  head(genes[[dd]])
}
rm(dd)
genes <- genes[sapply(genes, length)>0]
head(genes[[1]])




genes.df <-  Result.all[[1]] %>%
  dplyr::filter(
    anova.avg_log2FC.pval < 0.05 & anova.pct.dif.pval < 0.05 ) %>%
  arrange(Rank.ID ) %>%
  dplyr::select(one_of(c("rowname", "Rank.ID1")))
genes <- genes.df[["Rank.ID1"]]
names(genes) = genes.df[["rowname"]]
genes = sort(genes, decreasing =  FALSE)
head(genes)


gostres <- gost(query = genes, # a (named) list of gene identifiers
                organism = "hsapiens", 
                ordered_query = TRUE, # input genes are decreasingly ordered
                multi_query = FALSE, # multi_query = TRUE, # input queries are grouped according to term IDs for better comparison
                significant = T, 
                exclude_iea = FALSE, # exclude the electronic GO annotations
                measure_underrepresentation = FALSE, # measure under-representation instead of over-representation
                evcodes = TRUE, # includes the evidence codes to the results
                user_threshold = 0.05, 
                # as_short_link = TRUE, # query results can also be gathered into a short-link to the g:Profiler web tool.
                correction_method = "g_SCS", #' synonyms g_SCS and analytical, Bonferroni correction (correction_method = "bonferroni") or FDR (correction_method = "fdr")
                # domain_scope = "annotated", custom_bg = NULL, #' only the genes with at least one annotation are considered to be part of the full domain
                # domain_scope = "known", #' then all the genes of the given organism are considered to be part of the domain.
                domain_scope = "custom", custom_bg = (Result.all[[1]][["rowname"]]) , #' gost provides the means to define a custom background as a (mixed) list of gene identifiers with the parameter custom_bg
                # domain_scope = "custom_annotated",  custom_bg = sample(Result.all[[1]][["rowname"]], length(unlist(genes))) ,#' which will use the set of genes that are annotated in the data source and are also included in the user provided background list.
               
                numeric_ns = "", sources = NULL, as_short_link = FALSE)
names(gostres);head(gostres$result, 3);names(gostres$meta)
for (rr in 1:nrow(gostres$result)) {
  gostres$result[rr, "genes.matched"] <- paste(names(genes[as.vector(genes) %in% as.numeric(unlist(stringr::str_split(gostres[["result"]]$intersection[rr], ",")))]), collapse = ",")
}

gostplot(gostres, capped = TRUE, interactive = TRUE)
p <- gostplot(gostres, capped = FALSE, interactive = FALSE)
p

pp <- publish_gostplot(p, highlight_terms = c(gostres[["result"]]  %>% dplyr::filter(source== "KEGG", p_value <0.05) %>% top_n(-8, p_value) %>% pull(term_id)), 
                       width = NA, height = NA, filename = NULL )
pp

pd <- publish_gosttable(gostres, highlight_terms = gostres[["result"]]  %>% dplyr::filter(source== "KEGG", p_value <0.05) %>% top_n(-10, p_value), 
                  use_colors = TRUE, 
                  show_columns = c("source", "term_name", "term_size", "intersection_size"),
                  filename = NULL)

if (nrow(gostres[["result"]]) >0 ) {
  # p.adjust or pvalue
  df.all <- gostres[["result"]] %>% dplyr::filter( p_value < 0.05
                                        # ,Description != "Cytoplasmic ribosomal proteins"
  ) %>% 
    rename(Direction = query, Description = term_name) %>% 
    mutate(y = ifelse(Direction == "Up", -log10(p_value), 
                      ifelse(Direction == "Down", log10(p_value), NA)),
           Direction = factor(Direction, levels = c("Up", "Down"))
    )
  df.split <-split(df.all, df.all[, "source"])
  df.split <- df.split[sapply(df.split, nrow)>0]
  
  for (dd in names(df.split)){
    # dd = names(df.split[1]);dd
    df <- df.split[[dd]]
    ## plot ggboxplot for all data
    p <- ggbarplot(df %>% top_n(-50, p_value), x = "Description", y = "y", 
                   fill = "Direction",               # change fill color by cyl
                   color = "white",            # Set bar border colors to white
                   palette = c("red3", "skyblue3"),            # jco journal color palett. see ?ggpar
                   width = 1,
                   #main = paste0("Top ", "significant pathways"), 
                   #font.main = c(14, "bold", "black"), 
                   xlab = FALSE,
                   font.x = c(14, "bold", "black"),
                   ylab = "-log10(Pvalue) for Induced / log10(Pvalue) for Repressed", # "-log10(AdjPvalue) for Induced / log10(AdjPvalue) for Repressed", 
                   # font.y = c(0, "bold", "white"), # error:Error in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  : polygon edge not found
                   # facet.by = "clusters", 
                   scales = "free",  ncol = 1,
                   #panel.labs = list(clusters = c("","","","","")),
                   
                   #sort.val = "asc",          # Sort the value in ascending order
                   #sort.by.groups = TRUE,     # Don't sort inside each group
                   
                   legend.title = "Group",
                   font.legend = c(20, "plain", "black"),
                   legend ="top",
                   # x.text.angle = 90,           # Rotate vertically x axis texts
                   rotate = TRUE,
                   #orientation = "horizontal",
                   ggtheme = theme_minimal() )  +
      geom_hline(yintercept = -log10(0.05), color = "black", linetype="dashed", linewidth=0.5) + 
      geom_hline(yintercept = log10(0.05), color = "black", linetype="dashed", linewidth=0.5) +
      border() +
      theme(axis.text.y = element_text(size = 15, colour = "black", face = "bold"),
            strip.background = element_blank(),
            strip.text.x = element_blank(),
            panel.spacing.y = unit(0, "lines")) 
    p
    png(filename = paste0(dir,  sub("_[^_]+$", "", stringr::str_split(pathways.file, "/")[[1]][2]), "_",
                          stringr::str_split(nn, ":")[[1]][1], ".", "bar.plot.png"), 
        width = 18, 
        height= nrow(df) * 0.22 + 1 , res = 300, units = "in")
    print(p)
    dev.off()
  }

}



genes.convert <- gconvert(query = Result.all[[1]][["rowname"]], organism = "hsapiens", 
                          target="ENTREZGENE_ACC" , mthreshold = Inf, filter_na = TRUE) # "ENSG"
select.col <- setdiff(grep("_avg_log2FC_mean", colnames(Result.summary3[[1]])),grep("Direction", colnames(Result.summary3[[1]])) )
select.col <- setdiff(select.col, grep("SetDiff.value", colnames(Result.summary3[[1]])) );select.col
colnames(Result.summary3[[1]])[select.col]

dd <- genes.convert %>% 
  as.data.frame()%>% 
  distinct(input, .keep_all = TRUE) %>% 
  distinct(target, .keep_all = TRUE) %>% 
  left_join(Result.summary3[[1]][, c(select.col,
                                     grep("rowname", colnames(Result.summary3[[1]])))], by = c("input"="rowname")) %>% 
  dplyr::select(one_of(c(colnames(Result.summary3[[1]])[select.col], "target")))  %>% 
  tibble::column_to_rownames(var = "target") %>% 
  as.matrix()
colnames(dd) 
# colnames(dd) <- c("exp_1", "exp_2", "exp_3", "exp_4")
keggres = gage(dd, gsets=kegg.sets.hs, same.dir=TRUE, ref = grep('Older_Lean', colnames(dd), ignore.case =TRUE) )
pathview(gene.data = dd, pathway.id = ids[5],   species = "hsa",
         kegg.native = T,
         same.layer = T
         )
kegg.ids <- gostres[["result"]]  %>% dplyr::filter(source== "KEGG", p_value <0.05) %>% top_n(-20, p_value) %>% pull(term_id)
ids <- stringr::str_split_fixed(kegg.ids, ":", 2)[,2];ids

pv.out <- pathview(gene.data = dd, pathway.id = "04911",   species = "hsa",
         map.symbol=F
)
str(pv.out)
hist(dd)
pathview(gene.data = dd[,3], pathway.id = "01100" )
dd <- dd +0.5
pathview(gene.data = dd, pathway.id = "04261",   species = "hsa",
         kegg.native = T,
         multi.state = TRUE, match.data = FALSE, same.layer = T
)
gse16873.d =gse16873.d+2
data(gse16873.d); dim(gse16873.d)
pathview(gene.data = gse16873.d[,1], pathway.id = "04911" )
pathview(gene.data = gse16873.d, pathway.id = "04911",   species = "hsa",
         kegg.native = T,
         multi.state = TRUE, same.layer = T
)
