### DISCLAIMER
#   Some of the algorithms are non-deterministic making the results slightly different from run to run.
# Many of the algorithms are frequently updated with new releases making the results slightly different from version to version.
# Depending on your system, this code below may not produce an exact 1:1 copy of the results. 
# This adversely affects clustering, filtering etc. if the code is naively copy-pasted.
# For better reproducibility, please download quality-filtered Rds objects and work on downstream analyses. 
# <br>
#   change "D:/" to "/media/jianie/Extreme SSD1/" # if in LUNIX computer

### Load libraries -------------------------------------------------------------------------------------------------------------
#rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.
memory.limit(size = 1e+13)

# source(paste0("/media/jianie/Extreme SSD1/2022-09-01 STARR Combine/Project Parameters.R"), local = knitr::knit_global())
source(paste0(Disk, Project.folder, "/", codes.folder, "/", "Project Parameters.R"), local = knitr::knit_global())

# # Prepare for KEGGPathways analysis
library(pathview)
library(gage)
library(gageData)

# KEGG
data(kegg.sets.hs) # ist of 229 elements, a single KEGG pathway
data(sigmet.idx.hs) # a index numbers of sinaling and metabolic pathways in kegg.set.gs
kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs] # "cleaner" gene sets of sinaling and metabolic pathways only.
head(kegg.sets.hs, 3)
# data(kegg.gs)
#' kegg.gs only include the subset of canonical signaling and metabolic pathways from KEGG pathway database, and 
#' 

data(kegg.gs.dise)
#' kegg.gs.dise is the subset of disease pathways. 
#' 
#' And it is recommended to do KEGG pathway analysis with either kegg.gs or kegg.gs.dise seperately (rather than combined altogether) for better defined results. 
#' 

# GO
data(go.sets.hs) # list of 17202 elements, a single Gene Ontology term
data(go.gs) # list of 1000 elements, a single Gene Ontology term
data(go.subs.hs) # BP, CC, MF
gobpsets = go.sets.hs[go.subs.hs$BP]
goccsets = go.sets.hs[go.subs.hs$CC]
gomfsets = go.sets.hs[go.subs.hs$MF]
#' go.gs derived from human GO database only includes 1000 gene sets due to size limit. 
#' 
#' For full go.gs or gene sets data for other species, we may always use the gageData package. 
#' 
#' If you don’t find the gene set data for your target species, you may use kegg.gsets or go.gsets to compile pathway gene set data any time needed 
#' as long as it is one of the 3000 KEGG species or a species with gene annotation package supplied in Bioconductor or the user. 
#' You may want to save the new gene set data for later use. 
#' 
#' An another advantage of using kegg.gsets is that you get the most updated pathway gene set data as it is retrieved from KEGG in real time.

# BioCarta
data(carta.hs) # a named list of 259 elements

# prepare folder
dir1 <- paste0(Disk, Project.folder, "/", Rds.folder, "/");dir1

export.folder <- paste0("Pathways")
dir.create(file.path(paste0(Disk, Project.folder, "/", Rds.folder, "/"), export.folder), showWarnings = FALSE)
dir1.1 <- paste0(Disk, Project.folder, "/", Rds.folder, "/", export.folder, "/");dir1.1

dir.create(file.path(paste0(Disk, Project.folder, "/", figures.folder, "/"), export.folder), showWarnings = FALSE)
dir2 <- paste0(Disk, Project.folder, "/", figures.folder, "/", export.folder, "/");dir2


for (markers.file in markers.files) {
  # markers.file = markers.files[1]; markers.file
  
  Result <- readRDS(paste0(Disk, Project.folder, "/", Rds.folder, "/", markers.file))  
  names(Result); 
  
  # prepare folder
  export.folder0 <- stringr::str_split(markers.file, "/")[[1]] ;export.folder0
  dir1.2 <- dir1
  for (ee in 1:(length(export.folder0) -1) ) {
    dir.create(file.path(dir1, export.folder0[ee]), showWarnings = FALSE)
    dir1.2 <- paste0(dir1.2, export.folder0[ee], "/")
  }
  dir1.2
  
  Compare.variable = stringr::str_split(stringr::str_split(markers.file, "_Combined.Mean.Rds")[[1]][1], paste0(Select_Group, "_"))[[1]][2];
  if (is.na(Compare.variable)){Compare.variable = ""} else {Compare.variable = Compare.variable};Compare.variable
  dir.create(file.path(paste0(Disk, Project.folder, "/", figures.folder, "/", export.folder, "/"), 
                       paste0(stringr::str_split(unlist(stringr::str_split(markers.file, "[.]"))[3], "_")[[1]][4], ".", Select_Group)  ), showWarnings = FALSE)
  dir2.1 <- paste0(Disk, Project.folder, "/", figures.folder, "/", export.folder, "/",  
                   stringr::str_split(unlist(stringr::str_split(markers.file, "[.]"))[3], "_")[[1]][4], ".", Select_Group, "/"); dir2.1
  
  gostres_list <- list()
  gostres_plot_list  <- list()
  gostres_plot_list_combine <- list()
  gostres_plot_list_combine.plot <- list()
  
  Result.Path.list <- list()
  Result.Data.list <- list()
  Result.Path.list.final <- list()
  Result.Path.list.final.sig <- list()
  Result.Data.list.df <- list()
  Result.Path.list.final.plot <- list()
  
  library(gprofiler2)
  genes.convert <- gconvert(query = Result[[1]][[1]][["rowname"]], organism = organism, 
                            target="ENTREZGENE_ACC" , mthreshold = Inf, filter_na = TRUE) # "ENSG"
  
  # calculate results ---------------------------
  for (ll in names(Result) ) {
    # ll = names(Result)[1];ll
    
    # Gene list functional enrichment analysis with gost --------------------------------------------------
    for (nn in names(Result[[ll]]) ) {
      # nn = names(Result[[ll]])[1];nn
      if (nrow(Result[[ll]][[nn]]) > 0 ) {
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
        
        # generate gene file
        genes.df <- list()
        genes <- list()
        
        genes.df[["Up"]] <-  Result[[ll]][[nn]] %>% 
          dplyr::filter(Direction == "Up") %>%
          arrange(Rank.ID ) %>%
          dplyr::select(one_of(c("rowname", "Rank.ID")))
        genes[["Up"]] <- genes.df[["Up"]][["Rank.ID"]]
        names(genes[["Up"]]) = genes.df[["Up"]][["rowname"]]
        genes[["Up"]] = sort(genes[["Up"]], decreasing =  FALSE)
        head(genes[["Up"]])
        
        genes.df[["Down"]] <-  Result[[ll]][[nn]] %>% 
          dplyr::filter(Direction == "Down") %>%
          arrange(Rank.ID ) %>%
          dplyr::select(one_of(c("rowname", "Rank.ID.neg")))
        genes[["Down"]] <- genes.df[["Down"]][["Rank.ID.neg"]]
        names(genes[["Down"]]) = genes.df[["Down"]][["rowname"]]
        genes[["Down"]] = sort(genes[["Down"]], decreasing =  FALSE)
        head(genes[["Down"]])
        
        genes <- genes[sapply(genes, length)>0]
        if(length(genes) != 0){
          head(genes[[1]])
          
          # combine up and down genes
          # genes.df <-   Result[[ll]][[nn]] %>%
          #   dplyr::select(one_of(c("rowname", "Rank.ID.neg")))
          # genes <- genes.df[["Rank.ID.neg"]]
          # names(genes) = genes.df[["rowname"]]
          # genes = sort(genes, decreasing =  FALSE)
          # head(genes)
          
          # create file from gost 
          gostres <- gost(query = genes, # a (named) list of gene identifiers
                          organism = organism, 
                          ordered_query = TRUE, # input genes are decreasingly ordered
                          multi_query = FALSE, # multi_query = TRUE, # input queries are grouped according to term IDs for better comparison
                          significant = F, 
                          exclude_iea = FALSE, # exclude the electronic GO annotations
                          measure_underrepresentation = FALSE, # measure under-representation instead of over-representation
                          evcodes = TRUE, # includes the evidence codes to the results
                          user_threshold = 0.05, 
                          # as_short_link = TRUE, # query results can also be gathered into a short-link to the g:Profiler web tool.
                          correction_method = "g_SCS", #' synonyms g_SCS and analytical, Bonferroni correction (correction_method = "bonferroni") or FDR (correction_method = "fdr")
                          # domain_scope = "annotated", custom_bg = NULL, #' only the genes with at least one annotation are considered to be part of the full domain
                          # domain_scope = "known", #' then all the genes of the given organism are considered to be part of the domain.
                          domain_scope = "custom", custom_bg = (Result[[ll]][[nn]][["rowname"]]) , #' gost provides the means to define a custom background as a (mixed) list of gene identifiers with the parameter custom_bg
                          # domain_scope = "custom_annotated",  custom_bg = sample(Result.all[[1]][["rowname"]], length(unlist(genes))) ,#' which will use the set of genes that are annotated in the data source and are also included in the user provided background list.
                          
                          numeric_ns = "", sources = NULL, as_short_link = FALSE)
          if (!is.null(gostres)){
            names(gostres);head(gostres$result, 3);names(gostres$meta)
            genes.all <- c(genes[["Up"]], genes[["Down"]])
            for (rr in 1:nrow(gostres$result)) {
              gostres$result[rr, "genes.matched"] <- paste(names(genes.all[as.vector(genes.all) %in% as.numeric(unlist(stringr::str_split(gostres[["result"]]$intersection[rr], ",")))]), collapse = ",")
            }
            gostres_list[[ll]][[nn]] <- gostres
            
            # plotting
            gostplot(gostres, capped = TRUE, interactive = TRUE)
            if( exists("p1") ){rm(p1)}
            p1 <- gostplot(gostres, capped = FALSE, interactive = FALSE)
            p1
            
            # select sig: p.adjust or pvalue
            df.all <- gostres[["result"]] %>% 
              dplyr::filter( p_value < 0.05,
                             !term_id %in% c("KEGG:00000", "WP:000000", "REAC:0000000", "HP:0000001", "CORUM:0000000")
                             # ,Description != "Cytoplasmic ribosomal proteins"
              ) %>% 
              rename(Direction = query, Description = term_name) %>% 
              mutate(y = ifelse(Direction == "Up", -log10(p_value), 
                                ifelse(Direction == "Down", log10(p_value), NA)),
                     Direction = factor(Direction, levels = c("Up", "Down"))
              )
            if (nrow(df.all) > 0 ) {
              df.split <-split(df.all, df.all[, "source"])
              df.split <- df.split[sapply(df.split, nrow)>0]
              for (dd in names(df.split)){
                # dd = names(df.split)[7];dd
                df.split[[dd]] <- df.split[[dd]]%>%
                  mutate(path.group = dd, 
                         compared.groups = nn);dim(df)
              }
              for (dd in names(df.split)){
                gostres_plot_list[[ll]][[dd]][[nn]] <- df.split[[dd]]
              }
              for (dd in names(df.split)){
                # dd = names(df.split)[7];dd
                df <- df.split[[dd]]; dim(df)
                
                if(length(Result[[ll]]) == 1){
                  plot.pathways <- c(df %>% top_n(-5, p_value) %>% pull(term_id));plot.pathways
                  # plot.pathways <- plot.pathways[!plot.pathways %in% c("KEGG:00000", "WP:000000", "REAC:0000000", "HP:0000001")];plot.pathways
                  if( exists("pp") ){rm(pp)};
                  pp <- publish_gostplot(gostplot(gostres, capped = FALSE, interactive = FALSE), highlight_terms = plot.pathways, 
                                         width = NA, height = NA, filename = NULL )
                  pp
                  png(filename = paste0(dir2.1, 
                                        sub("_[^_]+$", "", stringr::str_split(markers.file, "/")[[1]][3]), "_",
                                        paste0(ll, ".", nn, ".", make.names(dd) ), "_l2_gprofiler.gostplot.png"), 
                      width =10, 
                      height= 10 , res = 300, units = "in")
                  print(pp)
                  dev.off()
                  
                  
                  plot.pathways <- c(df %>% top_n(-20, p_value) %>% pull(term_id));plot.pathways
                  # plot.pathways <- plot.pathways[!plot.pathways %in% c("KEGG:00000", "WP:000000", "REAC:0000000", "HP:0000001")];plot.pathways 
                  if( exists("pd") ){rm(pd)}
                  pd <- publish_gosttable(gostres, highlight_terms = plot.pathways, 
                                          use_colors = TRUE, 
                                          show_columns = c("source", "term_name", "term_size", "intersection_size", "p_value"),
                                          filename = NULL)
                  pd
                  png(filename = paste0(dir2.1, 
                                        sub("_[^_]+$", "", stringr::str_split(markers.file, "/")[[1]][3]), "_",
                                        paste0(ll, ".", nn, ".", make.names(dd) ), "_l2_gprofiler.gosttable.png"), 
                      width =20, 
                      height= 10 , res = 300, units = "in")
                  print(pd)
                  dev.off()
                }
              }
            }
          }
        }
      }
    }
    
    # calculate pathway from gage --------------
    sum.Path.list <- list();sum.data.list <- list(); sum.Path.list.final <- list(); sum.Path.list.final.sig <- list();sum.data.list.df<-data.frame()
    for (nn in names(Result[[ll]]) ) {
      # nn = names(Result[[ll]])[1];nn
      if (nrow(Result[[ll]][[nn]]) > 0 ) {
        select.col <- grep("_avg_log2FC", colnames(Result[[ll]][[nn]]))
        colnames(Result[[ll]][[nn]])[select.col]
        dd <- genes.convert %>% 
          as.data.frame()%>% 
          distinct(input, .keep_all = TRUE) %>% 
          distinct(target, .keep_all = TRUE) %>% 
          inner_join(Result[[ll]][[nn]]%>% 
                       dplyr::filter(Direction %in% c("Up", "Down")) %>% 
                       dplyr::select(one_of(colnames(Result[[ll]][[nn]])[c(select.col,
                                                                           grep("rowname", colnames(Result[[ll]][[nn]])))])) , by = c("input"="rowname")) %>% 
          dplyr::select(one_of(c(colnames(Result[[ll]][[nn]])[select.col], "target")))  %>% 
          tibble::column_to_rownames(var = "target") %>% 
          as.matrix()
        colnames(dd) ;head(dd)
        keggres = gage(dd, gsets=kegg.sets.hs, same.dir=TRUE , compare='as.group' )
        gobpres = gage(dd, gsets=gobpsets, same.dir=TRUE , compare='as.group' )
        goccres = gage(dd, gsets=goccsets, same.dir=TRUE , compare='as.group' )
        gomfres = gage(dd, gsets=gomfsets, same.dir=TRUE , compare='as.group' )
        cartares = gage(dd, gsets=carta.hs, same.dir=TRUE , compare='as.group' )
        kegg.diseres = gage(dd, gsets=kegg.gs.dise, same.dir=TRUE , compare='as.group' )
        
        sum.Path.list[["BioCarta"]][[nn]] <- cartares
        sum.Path.list[["GO.MF"]][[nn]] <- gomfres
        sum.Path.list[["GO.CC"]][[nn]] <- goccres
        sum.Path.list[["GO.BP"]][[nn]] <- gobpres
        sum.Path.list[["KEGG"]][[nn]] <- keggres
        sum.Path.list[["KEGG.dise"]][[nn]] <- kegg.diseres 
        sum.data.list[[nn]] <- dd
      }
    }
    sum.Path.list.final <- lapply(sum.Path.list, function(x){
      lapply(x, function(xx){
        for (gg in names(xx)){
          # gg = names(xx)[1];gg
          colnames( xx[[gg]] ) <- paste0(gg, "_" , colnames( xx[[gg]]) ) 
          xx[[gg]] <-  xx[[gg]] %>% as.data.frame() %>%  tibble::rownames_to_column(var = "id")
        }
        x <- Reduce(full_join, xx)
      })
    })
    sum.Path.list.final <- lapply(sum.Path.list.final, function(x){
      for (gg in names(x)){
        x[[gg]]$compared.groups <- gg
      }
      x <- do.call(rbind, x)
      x <- x %>% 
        dplyr::filter(complete.cases(.)) %>% 
        mutate(
          path.name = stringr::str_split_fixed(id, " ", 2)[, 2],
          hsa.name = stringr::str_split_fixed(id, " ", 2)[, 1] # substr(id, start=1, stop=8)
        )%>% 
        mutate(y.greater =  -log10((!!sym( paste0("greater", "_" , "p.val") ))),
               y.less =  -log10((!!sym( paste0("less", "_" , "p.val") ))) ) %>% 
        mutate(diff.y = y.greater - y.less)%>% 
        mutate(Direction = ifelse(diff.y > 0, "Up", 
                                  ifelse(diff.y < 0, "Down", "NS")) ) %>% 
        mutate(y = ifelse(Direction == "Up", -log10((!!sym( paste0("greater", "_" , "p.val") ))  ), 
                          ifelse(Direction == "Down", log10((!!sym( paste0("less", "_" , "p.val") ))  ), NA)),
               Direction = factor(Direction, levels = c("Up", "Down", "NS")) )
    })
    sum.Path.list.final <- sum.Path.list.final[sapply(sum.Path.list.final, nrow)>0] # Remove empty zero length rows in list with R
    for ( path.name1 in names(sum.Path.list.final) ) {
      sum.Path.list.final[[path.name1]][["path.name"]] <- path.name1
    }
    sum.Path.list.final.sig <- lapply(sum.Path.list.final, function(x){
      x <- x %>% 
        group_by(compared.groups) %>% 
        dplyr::filter((!!sym( paste0("greater", "_" , "p.val") )) < 0.05 | (!!sym( paste0("less", "_" , "p.val") )) < 0.05 ) 
    })
    
    # calculate average of each group
    for (gg in names(sum.data.list)){
      sum.data.list[[gg]] <- as.data.frame(sum.data.list[[gg]]) # rowname: entrez
      sum.data.list[[gg]]$compared.groups <- gg
    }
    head(sum.data.list[[1]])
    
    sum.data.list2 <- lapply(sum.data.list, function(x){
      x[[paste0("Mean_", unique(x$compared.groups), "_group")]] <-  rowMeans(dplyr::select(x, contains(unique(x$compared.groups))), na.rm = TRUE)
      x[[paste0("SD_", unique(x$compared.groups), "_group")]] <-  rowSds(as.matrix(dplyr::select(x, contains(unique(x$compared.groups))), na.rm = TRUE))
      x[[paste0("SEM_", unique(x$compared.groups), "_group")]] <-   rowSds(as.matrix(dplyr::select(x, contains(unique(x$compared.groups))), na.rm = TRUE))/sqrt(length(grep(unique(x$compared.groups), colnames(x), ignore.case =TRUE)))
      x[["entrez"]] <- rownames(x)  
      x <- x[, -c(grep("avg_log2FC", colnames(x), ignore.case =TRUE), grep("compared.groups", colnames(x), ignore.case =TRUE))] 
    })
    head(sum.data.list2[[1]])
    sum.data.list.df <- Reduce(full_join, sum.data.list2); 
    rownames(sum.data.list.df) <- NULL
    sum.data.list.df <-  sum.data.list.df %>%  tibble::column_to_rownames(var = "entrez")
    head(sum.data.list.df)
    
    if(length(sum.data.list) >0 ){ Result.Data.list[[ll]] <- sum.data.list}
    if(nrow(sum.data.list.df) >0 ){ Result.Data.list.df[[ll]] <- sum.data.list.df}
    if(length(sum.Path.list) >0 ){ Result.Path.list[[ll]] <- sum.Path.list}
    if(length(sum.Path.list.final) >0 ){ Result.Path.list.final[[ll]]  <- sum.Path.list.final}
    if(length(sum.Path.list.final.sig) >0 ){ Result.Path.list.final.sig[[ll]]  <- sum.Path.list.final.sig}
    
  }
  
  gostres_plot_list_combine <- lapply(gostres_plot_list, function(x){
    lapply(x, function(xx){
      for (gg in names(xx)){
        # gg = names(xx)[1];gg
        xx[[gg]]$compared.groups <- gg
      }
      xx <- do.call(rbind, xx)
    })
  })  
  gostres_plot_list_combine <- lapply(gostres_plot_list_combine, function(x){
    lapply(x, function(xx){
      xx <- xx %>%  
        mutate(compared.groups = factor(compared.groups, levels = rev(names(Result[[ll]]))),
               hsa.name = stringr::str_split_fixed(term_id, ":", 2)[,2])
      df.sum2 <- xx %>%
        mutate(Direction = factor(Direction, levels = c("Down", "Up" )))%>%
        group_by(compared.groups, term_id, Description) %>%
        arrange_at(c("Direction")) %>% 
        reframe(diff.y = sum(y)) # second one - first one
      xx <- xx %>%
        left_join(df.sum2, by = c("compared.groups","term_id", "Description"))
    })
  })  
  save(gostres_list,  
       gostres_plot_list ,  
       gostres_plot_list_combine ,  
       
       Result.Path.list,  
       Result.Data.list ,  
       Result.Path.list.final ,  
       Result.Path.list.final.sig ,  
       Result.Data.list.df ,  
       
       file = paste0(Disk, Project.folder, "/", Rds.folder, "/", stringr::str_split(markers.file, ".Rds")[[1]][1],
                     ".GProfiler.Pathways.RData") 
  )
  
  # plotting --------------------------------------
  for (ll in names(gostres_plot_list_combine) ) {
    # ll = names(gostres_plot_list_combine)[1];ll
    for (nn in names(gostres_plot_list_combine[[ll]]) ) {
      # nn = names(gostres_plot_list_combine[[ll]])[7];nn
      if (nrow(gostres_plot_list_combine[[ll]][[nn]]) > 0 ) {
        pathways <- gostres_plot_list_combine[[ll]][[nn]]
      }
    }  
    for (nn in names(Result.Path.list.final[[ll]]) ) {
      # nn = names(gostres_plot_list_combine[[ll]])[7];nn
      if (nrow(Result.Path.list.final[[ll]][[nn]]) > 0 ) {
        ## plot ggboxplot for all gost results -----------
        if( exists("p1") ){ rm(p1) }
        df <- data.frame()
        colnames(pathways)
        
        # get rid of duplicated rows
        df.sum <- pathways %>%
          group_by(intersection, compared.groups) %>% 
          dplyr::filter(!(duplicated(intersection))) %>% # |duplicated(intersection, fromLast = TRUE)
          ungroup()
        #table(df.sum$intersection)
        
        if (length(names(Result[[ll]])) > 1){
          df.sum <- df.sum %>%
            group_by(term_id, Description ) %>%
            dplyr::filter(n_distinct((!!sym("compared.groups"))) > 1) 
          # table(df.sum$compared.groups, df.sum$Description)
          if (nrow(df.sum) >0){
            df.sum2 <-df.sum %>%
              group_by(term_id, Description) %>%
              reframe(diff.y.path = diff(range(diff.y)) )%>% 
              arrange(-diff.y.path) %>%  
              mutate(rank = 1:n())
            # select.rows <- df.sum %>% dplyr::filter(dense_rank(mean.y) <= 10 | dense_rank(dplyr::desc(mean.y)) <= 10);
            select.rows <- df.sum %>% 
              left_join(df.sum2, by = c("term_id", "Description"))  %>% 
              dplyr::filter(abs(diff.y.path) >= 0.5) %>%  
              dplyr::filter(abs(y) >= -log10(0.05)) %>%
              arrange(-y) %>%  
              mutate(rank.y = 1:n())  
            df <- pathways %>%
              left_join(df.sum2, by = c("term_id", "Description"))  %>%
              mutate(select.plot = ifelse(term_id %in% c(select.rows%>% ungroup() %>% 
                                                           dplyr::filter( dense_rank(dplyr::desc(diff.y.path)) <= 20) %>% 
                                                           pull(term_id)), "Yes", NA),
                     select = ifelse(term_id %in% c(select.rows %>% pull(term_id)), "Yes", NA)
              ) %>%
              arrange(rank) 
            df <-  df %>%    
              mutate( Description = factor(Description, levels = unique(df$Description)),
                      compared.groups = factor(compared.groups, levels = rev(names(Result[[ll]]))))
            table(df$select); table(df$select.plot) 
          }
        } else {
          select.rows <- pathways %>% 
            dplyr::filter(abs(diff.y) >= -log10(0.05)) %>%  
            arrange(-abs(diff.y)) %>%  
            mutate(rank = 1:n())
          if (nrow(select.rows) >0){
            df <-  pathways %>%  
              mutate(select.plot = ifelse(term_id %in% c(select.rows %>%
                                                           dplyr::filter(dense_rank(diff.y) <= ceiling(20/length(names(Result[[ll]])) ) | dense_rank(dplyr::desc(diff.y)) <= ceiling(20/length(names(Result[[ll]])) ))%>% 
                                                           pull(term_id)), "Yes", NA),
                     select = ifelse(term_id %in% c(select.rows%>% pull(term_id)), "Yes", NA) )  %>%  
              arrange(-abs(diff.y)) %>%  
              mutate(rank = 1:n())
            # dplyr::filter(select == "Yes") 
            df <-  df %>%  
              mutate( Description = factor(Description, levels = unique(df$Description)),
                      compared.groups = factor(compared.groups, levels = rev(names(Result[[ll]]))))
          }
        }
        if (nrow(df) >0){
          gostres_plot_list_combine.plot[[ll]][[nn]] <- df
          df <- df  %>% dplyr::filter(select.plot == "Yes") 
          p1 <- ggbarplot(df %>% dplyr::filter(!grepl("disease", Description)), x = "compared.groups", y = "y", 
                          fill = "Direction",               # change fill color by cyl
                          color = "white",            # Set bar border colors to white
                          palette = c("red3", "skyblue3"),            # jco journal color palett. see ?ggpar
                          width = 1,
                          main = paste0("Top ", "significant ", nn, " pathways"," in ",  ll, ", compare to ",  Compare.variable , ", by ", "gost ", "function" ), 
                          #font.main = c(14, "bold", "black"), 
                          xlab = FALSE,
                          font.x = c(14, "bold", "black"),
                          ylab = "-log10(Pvalue) for Induced / log10(Pvalue) for Repressed", # "-log10(AdjPvalue) for Induced / log10(AdjPvalue) for Repressed", 
                          # font.y = c(0, "bold", "white"), # error:Error in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  : polygon edge not found
                          facet.by = c( "Description", "path.group"), # "compared.groups", # panel.labs = list(compared.groups = unique(df$compared.groups)),
                          short.panel.labs = TRUE,  strip.position = "bottom",# ncol = 2,
                          panel.labs.font.y = list(face = "bold", color = "black", size = 15, angle = 0, hjust = 0, vjust = 0.5),
                          scales = "free",  # "fixed", # 
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
            theme(axis.text.y = element_text(size = 12, colour = "black", face = "bold"),
                  # strip.background = element_blank(),
                  strip.background = element_rect(colour = "black", fill = "gray87"),
                  # strip.text.x = element_blank(),
                  strip.text.x = element_text(size = 15, colour = "black", face = "bold"),# , angle = 45, hjust = 1, vjust = 0.7
                  panel.spacing.y = unit(0, "lines")) + scale_x_discrete(position = "bottom") ;p1
          png(filename = paste0(dir2.1,  
                                sub("_[^_]+$", "", stringr::str_split(markers.file, "/")[[1]][3]), "_",
                                paste0(ll, ".", nn), "_l2_gprofiler.ggbarplot.png"),  
              width = 13, 
              height= nrow(df) * 0.23 + 2.5 , res = 300, units = "in")
          print(p1)
          dev.off()
        }
        
        ## plot pathview for KEGG from gost  -----------------------------
        if( nn %in% c("KEGG", "KEGG.dise") ) {
          df_list <- list()
          if (length(names(Result[[ll]])) > 1){
            df_list[["Sig"]] <- df %>%
              ungroup() %>% as.data.frame() %>%
              dplyr::filter(!(duplicated(Description))) %>% 
              dplyr::filter(select == "Yes") %>%
              top_n( -10, rank)
          } else {
            df1 <- df %>%
              dplyr::filter(Direction == "Up"   ) %>%
              ungroup() %>% as.data.frame() %>%
              dplyr::filter(select == "Yes") %>%
              top_n( -5, rank)
            df2 <-df %>%
              dplyr::filter(Direction == "Down"   ) %>%
              ungroup() %>% as.data.frame() %>%
              dplyr::filter(select == "Yes") %>%
              top_n( 5, rank)
            if(nrow(df1) > 0) {df_list[["Up"]] <- df1}
            if(nrow(df2) > 0) {df_list[["Down"]] <- df2}
          }
          for (vv in names(df_list)) {
            # vv= names(df_list)[1]; vv
            
            for (dd.row in 1:nrow(df_list[[vv]])){
              # dd.row = (1:nrow(df_list[[vv]]))[1];dd.row
              dd <- Result.Data.list.df[[ll]] %>%
                dplyr::select(contains("Mean"))
              setwd(dir2.1)
              head(df_list[[vv]])
              pathview(gene.data = as.matrix(dd) , pathway.id = as.character(df_list[[vv]][dd.row, "hsa.name"]), # stringr::str_split(rownames(df[i, ]), " ")[[1]][1],
                       out.suffix = paste0(sub("_[^_]+$", "", stringr::str_split(markers.file, "/")[[1]][3]), "_",
                                           paste0(ll, ".", nn, ".", vv),  "_l2_gost.ori" ), kegg.native = T, multi.state = T )# native original KEGG graph, png   # "_",Compare.Group1,
              
              pathview(gene.data = as.matrix(dd) , pathway.id = as.character(df_list[[vv]][dd.row, "hsa.name"]), # stringr::str_split(rownames(df[i, ]), " ")[[1]][1],
                       out.suffix = paste0(sub("_[^_]+$", "", stringr::str_split(markers.file, "/")[[1]][3]), "_",
                                           paste0(ll, ".", nn, ".", vv),  "_l2_gost.Graphviz" ),   kegg.native = F, same.layer = F, multi.state = T) # Graphviz engine, pdf, page 1 is the main graph, page 2 is the legend
              setwd("~/")
            }
          }
        }
        
        ## plot ggboxplot for gage results -----------
        if( exists("p1") ){ rm(p1) }
        pathways <- Result.Path.list.final[[ll]][[nn]]
        df <- data.frame()
        colnames(pathways)
        
        # get rid of duplicated rows
        df.sum <- pathways %>%
          group_by(id, compared.groups) %>% 
          dplyr::filter(!(duplicated(id))) %>% # |duplicated(intersection, fromLast = TRUE)
          ungroup()
        #table(df.sum$intersection)
        
        if (length(names(Result[[ll]])) > 1){
          df.sum <- df.sum %>%
            group_by(id, path.name ) %>%
            dplyr::filter(n_distinct((!!sym("compared.groups"))) > 1) 
          # table(df.sum$compared.groups, df.sum$Description)
          if (nrow(df.sum) >0){
            df.sum2 <-df.sum %>%
              group_by(id, path.name) %>%
              reframe(diff.y.path = diff(range(diff.y)) )%>% 
              arrange(-diff.y.path) %>%  
              mutate(rank = 1:n())
            # select.rows <- df.sum %>% dplyr::filter(dense_rank(mean.y) <= 10 | dense_rank(dplyr::desc(mean.y)) <= 10);
            select.rows <- df.sum %>% 
              left_join(df.sum2, by = c("id", "path.name"))  %>% 
              dplyr::filter(abs(diff.y.path) >= 0.5) %>%  
              dplyr::filter(abs(y) >= -log10(0.05)) %>%
              arrange(-y) # %>%  
              # mutate(rank.y = 1:n())  
            df <- pathways %>%
              left_join(df.sum2, by = c("id", "path.name"))  %>%
              mutate(select.plot = ifelse(id %in% c(select.rows%>% ungroup() %>% 
                                                      dplyr::filter( dense_rank(dplyr::desc(diff.y.path)) <= 20) %>% 
                                                      pull(id)), "Yes", NA),
                     select = ifelse(id %in% c(select.rows %>% pull(id)), "Yes", NA)
              ) %>%
              arrange(rank) 
            df <-  df %>%    
              mutate( path.name = factor(path.name, levels = unique(df$path.name)),
                      compared.groups = factor(compared.groups, levels = rev(names(Result[[ll]]))))
            table(df$select); table(df$select.plot) 
          }
        } else {
          select.rows <- pathways %>% 
            dplyr::filter(abs(diff.y) >= -log10(0.05)) %>%  
            arrange(-abs(diff.y)) %>%  
            mutate(rank = 1:n())
          if (nrow(select.rows) >0){
            df <-  pathways %>%  
              mutate(select.plot = ifelse(term_id %in% c(select.rows %>%
                                                           dplyr::filter(dense_rank(diff.y) <= ceiling(20/length(names(Result[[ll]])) ) | dense_rank(dplyr::desc(diff.y)) <= ceiling(20/length(names(Result[[ll]])) ))%>% 
                                                           pull(term_id)), "Yes", NA),
                     select = ifelse(term_id %in% c(select.rows%>% pull(term_id)), "Yes", NA) )  %>%  
              arrange(-abs(diff.y)) %>%  
              mutate(rank = 1:n())
            # dplyr::filter(select == "Yes") 
            df <-  df %>%  
              mutate( Description = factor(Description, levels = unique(df$Description)),
                      compared.groups = factor(compared.groups, levels = rev(names(Result[[ll]]))))
          }
        }
        if (nrow(df  %>% dplyr::filter(select.plot == "Yes")) >0){
          Result.Path.list.final.plot[[ll]][[nn]] <- df
          df <- df  %>% dplyr::filter(select.plot == "Yes") 
          p1 <- ggbarplot(df , x = "compared.groups", y = "y", 
                          fill = "Direction",               # change fill color by cyl
                          color = "white",            # Set bar border colors to white
                          palette = c("red3", "skyblue3"),            # jco journal color palett. see ?ggpar
                          width = 1,
                          main = paste0("Top ", "significant ", nn, " pathways"," in ",  ll, ", compare to ",  Compare.variable , ", by ", "gost ", "function" ), 
                          #font.main = c(14, "bold", "black"), 
                          xlab = FALSE,
                          font.x = c(14, "bold", "black"),
                          ylab = "-log10(Pvalue) for Induced / log10(Pvalue) for Repressed", # "-log10(AdjPvalue) for Induced / log10(AdjPvalue) for Repressed", 
                          # font.y = c(0, "bold", "white"), # error:Error in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  : polygon edge not found
                          facet.by = c( "path.name", "path.group"), # "compared.groups", # panel.labs = list(compared.groups = unique(df$compared.groups)),
                          short.panel.labs = TRUE,  strip.position = "bottom",# ncol = 2,
                          panel.labs.font.y = list(face = "bold", color = "black", size = 15, angle = 0, hjust = 0, vjust = 0.5),
                          scales = "free",  # "fixed", # 
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
            theme(axis.text.y = element_text(size = 12, colour = "black", face = "bold"),
                  # strip.background = element_blank(),
                  strip.background = element_rect(colour = "black", fill = "gray87"),
                  # strip.text.x = element_blank(),
                  strip.text.x = element_text(size = 15, colour = "black", face = "bold"),# , angle = 45, hjust = 1, vjust = 0.7
                  panel.spacing.y = unit(0, "lines")) + scale_x_discrete(position = "bottom") ;p1
          png(filename = paste0(dir2.1,  
                                sub("_[^_]+$", "", stringr::str_split(markers.file, "/")[[1]][3]), "_",
                                paste0(ll, ".", nn ), "_l2_gage.ggbarplot.png"),  
              width = 13, 
              height= nrow(df) * 0.23 + 2.5 , res = 300, units = "in")
          print(p1)
          dev.off()
        }
        
        ## plot pathview for KEGG from gage  -----------------------------
        if( nn %in% c("KEGG", "KEGG.dise") ) {
          df_list <- list()
          if (length(names(Result[[ll]])) > 1){
            df_list[["Sig"]] <- df %>%
              ungroup() %>% as.data.frame() %>%
              dplyr::filter(!(duplicated(path.name))) %>% 
              dplyr::filter(select == "Yes") %>%
              top_n( -10, rank)
          } else {
            df1 <- df %>%
              dplyr::filter(Direction == "Up"   ) %>%
              ungroup() %>% as.data.frame() %>%
              dplyr::filter(select == "Yes") %>%
              top_n( -5, rank)
            df2 <-df %>%
              dplyr::filter(Direction == "Down"   ) %>%
              ungroup() %>% as.data.frame() %>%
              dplyr::filter(select == "Yes") %>%
              top_n( 5, rank)
            if(nrow(df1) > 0) {df_list[["Up"]] <- df1}
            if(nrow(df2) > 0) {df_list[["Down"]] <- df2}
          }
          for (vv in names(df_list)) {
            # vv= names(df_list)[1]; vv
            
            for (dd.row in 1:nrow(df_list[[vv]])){
              # dd.row = (1:nrow(df_list[[vv]]))[1];dd.row
              dd <- Result.Data.list.df[[ll]] %>%
                dplyr::select(contains("Mean"))
              setwd(dir2.1)
              head(df_list[[vv]])
              pathview(gene.data = as.matrix(dd) , pathway.id = as.character(df_list[[vv]][dd.row, "hsa.name"]), # stringr::str_split(rownames(df[i, ]), " ")[[1]][1],
                       out.suffix = paste0(sub("_[^_]+$", "", stringr::str_split(markers.file, "/")[[1]][3]), "_",
                                           paste0(ll, ".", nn, ".", vv),  "_l2_gage.ori" ), kegg.native = T, multi.state = T )# native original KEGG graph, png   # "_",Compare.Group1,
              
              pathview(gene.data = as.matrix(dd) , pathway.id = as.character(df_list[[vv]][dd.row, "hsa.name"]), # stringr::str_split(rownames(df[i, ]), " ")[[1]][1],
                       out.suffix = paste0(sub("_[^_]+$", "", stringr::str_split(markers.file, "/")[[1]][3]), "_",
                                           paste0(ll, ".", nn, ".", vv),  "_l2_gage.Graphviz" ),   kegg.native = F, same.layer = F, multi.state = T) # Graphviz engine, pdf, page 1 is the main graph, page 2 is the legend
              setwd("~/")
            }
          }
        }
      }
    }
  }
  save(gostres_list,  
       gostres_plot_list ,  
       gostres_plot_list_combine ,  
       gostres_plot_list_combine.plot ,  
       
       Result.Path.list,  
       Result.Data.list ,  
       Result.Path.list.final ,  
       Result.Path.list.final.sig ,  
       Result.Data.list.df ,  
       Result.Path.list.final.plot,
       
       file = paste0(Disk, Project.folder, "/", Rds.folder, "/", stringr::str_split(markers.file, ".Rds")[[1]][1],
                     ".GProfiler.Pathways.RData") 
  )
  print(markers.file)
}

rm(bods, carta.hs, dd, df, df.all, df.split, genes, genes.convert, genes.df, 
   go.gs, go.sets.hs, go.subs.hs, gobpsets, goccsets, gomfsets, 
   gostres, gostres_list, kegg.gs.dise, kegg.sets.hs, p, pd, pp, Result,
   dir1, dir1.1, dir1.2, dir2, dir2.1, ee, export.folder, export.folder0, ids,
   kegg.ids, ll, nn, plot.pathways, rr, sigmet.idx.hs, select.col, gene.all,  kegg.list,
   cpd.names, kegg.met, keggres, pathways, pathways.sig, rn.list, gg
)
gc()
sessionInfo()