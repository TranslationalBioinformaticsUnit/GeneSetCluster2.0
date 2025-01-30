### Functions --------------------------------------------------------------------
# These functions are not meant to be invoked directly by the user.
# See the PlotPathwayCluster method instead.
obtainDFmetadata <- function(df, mat_cor) {
  ngroups <- unique(df$Groups)
  info <- vector("list", length(ngroups))

  for (i in 1:length(ngroups)) {
    info[[i]] <- df[which(df$Groups==ngroups[i]), "Pathways"]
  }

  df_metadata = as.data.frame(matrix(0, nrow = nrow(mat_cor), ncol = length(ngroups)))
  rownames(df_metadata) = colnames(mat_cor)
  colnames(df_metadata) = ngroups

  for (i in 1:nrow(mat_cor)) {
    for (j in 1:length(ngroups)) {
      if (rownames(mat_cor)[i] %in% info[[j]]) {
        df_metadata[i, j] <- 1
      }
    }
  }

  return(df_metadata)
}


obtainDefCluster <- function(mat_using){

  candidates_clustering <- c()
  definitive_clusters <- list()
  j <- 1
  num_sim <- 0.6
  for (i in 1:(ncol(mat_using)-1)) {
    if (i < j-1 | names(mat_using[i,])[j] %in% candidates_clustering) {  ## Avoid repeating same cluster
      next
    }

    if (i == ncol(mat_using)-1) {
      if (mat_using[i,i+1] >= num_sim) {
        candidates_clustering <- c(names(mat_using[i,])[i], names(mat_using[i,])[i+1])
        definitive_clusters <- append(definitive_clusters, list(candidates_clustering))
      }
    } else {

      if (mat_using[i,i+1] >= num_sim) {
        candidates_clustering <- c(names(mat_using[i,])[i])

        for (j in (i+1):(ncol(mat_using))) {

          if (j==ncol(mat_using)) {
            candidates_clustering <- c(candidates_clustering, names(mat_using[i,])[j])
            definitive_clusters <- append(definitive_clusters, list(candidates_clustering))
            break
          }

          if(mat_using[i,j] >= num_sim) {
            candidates_clustering <- c(unique(candidates_clustering), names(mat_using[i,])[j])
            if (j == (ncol(mat_using)-2)) {
              append(definitive_clusters, list(candidates_clustering))
            }

          } else if (mat_using[i+1,j] >= num_sim){ # The lower this value is, the larger group of clusters is obtained
            candidates_clustering <- c(unique(candidates_clustering), names(mat_using[i,])[j])
            if (j == (ncol(mat_using)-2)) {
              append(definitive_clusters, list(candidates_clustering))
            }

          } else {
            definitive_clusters <- append(definitive_clusters, list(candidates_clustering))
            break
          }
        }
      }
    }
  }

  return(definitive_clusters)
}



# This function has been adapted from simplifyEnrichment package (Gu, Z. & Hubschmann, D. 2022).
#' @import simplifyEnrichment
#' @import GetoptLong
keywordEnrichment <- function(term_id, tdm, min_bg = 5, min_term = 2) {
  GO_EXCLUDE_WORDS <- c("via", "protein", "factor", "side", "type", "specific", "biological", "process")

  tdm2 <- tdm[slam::row_sums(tdm) >= 5, ]
  l <- colnames(tdm2) %in% term_id

  n <- nrow(tdm2)
  n_term <- numeric(n)
  n_bg <- numeric(n)
  p <- numeric(n)
  for(i in seq_len(n)) {
    if(interactive() && simplifyEnrichment::se_opt$verbose) {
      if(i %% 100 == 0 || i == n) {
        message(strrep("\r", 100), appendLF = FALSE)
        message(GetoptLong::qq("performing keyword enrichment, @{i}/@{n}..."), appendLF = FALSE)
      }
    }
    v <- as.vector(tdm2[i, ])
    s11 <- sum(v & l)
    if(s11 < 2) {
      next
    }
    s12 <- sum(!v & l)
    s21 <- sum(v & !l)
    s22 <- sum(!v & !l)

    n_term[i] <- s11
    n_bg[i] <- s11 + s21

    p[i] <- fisher.test(cbind(c(s11, s21), c(s12, s22)), alternative = "greater")$p.value

  }
  if(interactive() && simplifyEnrichment::se_opt$verbose) {
    message("")
  }

  df <- data.frame(keyword = rownames(tdm2), n_term = n_term, n_bg = n_bg, p = p)
  df <- df[df$n_term >= 2, , drop = FALSE]
  df$padj <- p.adjust(df$p)
  df[order(df$padj, df$p), , drop = FALSE]

  df <- df[which(!df$keyword %in% GO_EXCLUDE_WORDS),]

  return(df)
}


obtainGOdescription <- function(go_ids) {
  tryCatch(
    {
      go_description <- AnnotationDbi::select(GO.db, keys=go_ids, columns="TERM")
    },
    error=function(e) {
      message(paste0(go_ids, " not found."))
      go_description <- "Not found"

    }
  )

  return(go_description)
}


#' @import simplifyEnrichment
#' @import ggwordcloud
wordcloud_generation <- function(go_id_list, min_stat=5) {
  env_tdm_GO <- readRDS(system.file("extdata", "tdm_GO.rds", package = "simplifyEnrichment"))
  df <- keywordEnrichment(go_id_list, env_tdm_GO)
  df <- df[df$p <= min_stat, , drop = FALSE]
  df <- data.frame(df[, 1], -log10(df$p))
  if (nrow(df) > 10) {
    df <- df[1:10,]
  }

  wordplot <- ggplot(df, aes(label = df[,1], size = df[,2])) +
                             # color = factor(sample.int(10, nrow(df), replace = TRUE)))) + #to add colour
                    geom_text_wordcloud(shape="square") +
                    theme_minimal() +
                    labs(title = "")+
                    theme(
                      plot.title = element_blank(),
                      legend.title = element_blank(),
                      legend.background = element_blank(),
                      legend.box.spacing = unit(0, "mm")
                    )

  return(wordplot)
}

obtainOrg <- function(Object) {
  mm <- c("MUS MUSCULUS" ,"MM", "ORG.MM.EG.DB")
  hs <- c("HOMO SAPIENS", "HS", "ORG.HS.EG.DB")
  if (toupper(Object@metadata$organism[1]) %in% hs) {
    usingOrg <- org.Hs.eg.db
  } else if (toupper(Object@metadata$organism[1]) %in% mm) {
    usingOrg <- org.Mm.eg.db
  } else {
    message("Not recognised organism. Will use Homo sapiens data base as default.")
    usingOrg <- org.Hs.eg.db
  }
  return(usingOrg)
}

#' @import AnnotationDbi
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @import clusterProfiler
#' @import stringr
#' @export
ORAperCluster <- function(Object, doORA, clusters, clus, completeRes = FALSE) {
  options <- c("SYMBOL", "ENTREZID", "ENSEMBLID", "ENTREZ", "ENSEMBL")
  entrezOptions <- c("ENTREZID", "ENTREZ")
  ensemblOptions <- c("ENSEMBL", "ENSEMBLID")

  # internal check of parameters
  if (!toupper(Object@metadata$structure[1]) %in% options | doORA == FALSE) {
    if (!toupper(Object@metadata$structure[1]) %in% options) {
      message("Genes are not included as genes Symbol, EntrezID or ENSEMBLID.\nORA of the genes involved in the cluster pathway will NOT be performed.")
    }
    keywords <- rep(".", times=clusters)
    return(keywords)
  }

  #organism
  usingOrg <- obtainOrg(Object)

  #obtain genes EntrezID for ORA
  info_df <- data.frame(Object@Data[[1]]$Pathways, Object@Data[[1]]$Molecules)
  # keywords <- vector()
  message(paste0("Performing ORA for ", clusters, " clusters and choosing the most overrepresented pathway..."))
  keywords <- list()

  for (i in unique(clus)) {
    paths_of_clust <- names(clus[clus == i])
    genes2check <- unique(unlist(str_split(paste0(info_df[which(info_df[,1] %in% paths_of_clust),2], collapse = ","), Object@metadata$seperator[1])))

    if (toupper(Object@metadata$structure[1]) == "SYMBOL") {
      if (table(genes2check==toupper(genes2check))[1] < length(genes2check)/2 | table(genes2check==tolower(genes2check))[1] < length(genes2check)/2 | all(genes2check==toupper(genes2check)) | all(genes2check==tolower(genes2check))){ # To correct capital letters of gene symbol
        genes2check <- tolower(genes2check)
        genes2check <- paste(toupper(substr(tolower(genes2check), 1, 1)), substr(tolower(genes2check), 2, nchar(tolower(genes2check))), sep = "")
      }
      entrezID <- AnnotationDbi::select(usingOrg, keys=genes2check, columns='ENTREZID', keytype='SYMBOL')
      # entrezID <- AnnotationDbi::select(usingOrg, keys=toupper(genes2check), columns='ENTREZID', keytype='SYMBOL')
      entrezID <- entrezID$ENTREZID
    } else if (toupper(Object@metadata$structure[1]) %in% entrezOptions) {
      entrezID <- genes2check
    } else if (toupper(Object@metadata$structure[1]) %in% ensemblOptions) {
      entrezID <- AnnotationDbi::select(usingOrg, keys=genes2check, columns='ENTREZID', keytype='ENSEMBL')
      entrezID <- entrezID$ENTREZID
    }

    ora <- clusterProfiler::enrichGO(na.omit(entrezID), OrgDb=usingOrg, keyType="ENTREZID",
                                     ont="BP", pvalueCutoff=1, pAdjustMethod="BH", qvalueCutoff=1)

    if (completeRes == TRUE) {
      keywords[[i]] <- ora
    } else {

      if (nchar(ora@result$Description[1]) > 50) {
        keywords[i] <- paste0(substr(ora@result$Description[1], 1, 40), "... ", ora@result$ID[1])
      } else {
        keywords[i] <- ora@result$Description[1]
      }
    }
  }
  return(keywords)
}

#' @import GO.db
getterm <- function(goid){
  termid <- GO.db::GOTERM[[goid]]
  if(is.null(termid)) {
    return("NA")
  } else {
    return(AnnotationDbi::Term(termid))
  }
}

convertNA <- function(x) {
  na_count <- sum(x == "NA")
  ifelse(x == "NA", paste0("NA.", seq_len(na_count)), x)
}

checkGO <- function(Object) {
  if (substr(Object@Data[[1]]$Pathways[1], 1, 3) == "GO:" & !is.na(is.numeric(substr(Object@Data[[1]]$Pathways[1], 4, 6)))) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

parallelization.message <- function(threads) {
  message(paste0("Parallelizing processes, ", threads, " cores have been selected. "), sep="")

  if (threads > detectCores()) {
    message("More cores have been selected than detected.")
    threads <- detectCores()
    message(paste0("The maximum number of cores will be used:", threads, " cores."))
  }

  return(threads)
}

scale_values <- function(x, ...) {(x - min(x, ...)) / (max(x, ...) - min(x, ...))}

makeSymm <- function(m) {
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  return(m)
}


# guess the separator of a csv file
guessDelim <- function(fileLocation) {
  allDelim <- c(";", "\t", ",")
  delimReturn <- NULL
  for (delim in allDelim)
  {
    file <- read.csv(fileLocation, sep = delim)
    if (ncol(file) > 5)
    {
      delimReturn <- delim
    }
  }
  if (is.null(delimReturn))
  {
    message("The file introduced is not correct. Please check it.")
    stop()
  }

  return(delimReturn)
}


#scaling the matrix
scaleCorMatrix <- function(mat)
  {
  mat[lower.tri(mat, diag = F)] <- NA
  mat_mod <- apply(mat, 1, function(x) scale_values(x, na.rm=T))
  mat_sym <- makeSymm(mat_mod)
  mat_sym[ncol(mat_sym), ncol(mat_sym)] <- 1

  if (anyNA(mat_sym)) { #because of the scaling method when the last values of the diagonal have the maximum value they are setted as NA. We are fixing it.
    mat_sym[which(is.na(mat_sym))] <- 1
  }

  return(mat_sym)
}


#check if introduced tissues names is correct
checkTissues <- function(tissues, tissue_names)
  {
  if (tissues[1] != "ALL")
  {
    if (all(tissues %in% tissue_names)) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  } else {
    return(TRUE)
  }

}


prepare_plotting_df <- function(df) {
  df_tissue_z[is.na(df_tissue_z)] <- 0
  plot_melt_mtx_z <- melt(df_tissue_z)
  colnames(plot_melt_mtx_z)[1:3] <- c("cluster", "Dataset", "value")
  plot_cluster_size_z <- aggregate(value ~ cluster, data = plot_melt_mtx_z, FUN = sum)
  colnames(plot_cluster_size_z) <- c("cluster", "value")

  return(plot_cluster_size_z)
}

# this function filter the genes of GTEx database and select ENTREZID
getGenesTissue <- function(Object, mol.unique.df, dic)
{
  entrezOptions <- c("ENTREZID", "ENTREZ")
  ensemblOptions <- c("ENSEMBL", "ENSEMBLID")
  genes <- rownames(mol.unique.df)

  if (substr(genes[1], 1, 3) == "ENS") {
    dic.selected <- dic[which(dic$GTEx.median.TPM.Name %in% toupper(genes)),]
    genes.selected <- dic.selected$GTEx.median.TPM.Name

    if (length(genes.selected) == 0) {
      stop("The genes used to create the object are not correct. Please make sure they are introduced either as GENE SYMBOL, ENTREZID or ENSEMBLID.")
    }

    touse <- mol.unique.df[which(toupper(genes) %in% dic$GTEx.median.TPM.Name),]

    #label number of genes per cluster
    ngenes <- apply(touse, 2, function(x) length(which(x==1)))
    colnames(touse) <- paste0(names(ngenes), " (", ngenes, ")")

    touse$ENSID <- rownames(touse)
    meltedtouse <- melt(touse, id="ENSID")
    meltedtouse <- meltedtouse[which(meltedtouse$value==1), c("variable", "ENSID")]
    colnames(meltedtouse) <- c("Term", "Gene")

  } else if (toupper(Object@metadata$structure[1]) %in% entrezOptions) {
    usingOrg <- obtainOrg(Object)
    resSymbol <- AnnotationDbi::select(usingOrg, keys=genes, columns='SYMBOL', keytype='ENTREZID')

    dic.selected <- dic[which(dic$GTEx.median.TPM.Description %in% toupper(resSymbol$SYMBOL)),]
    genes.selected <- dic.selected$GTEx.median.TPM.Name

    if (length(genes.selected) == 0) {
      stop("The genes used to create the object are not correct. Please make sure they are introduced either as GENE SYMBOL, ENTREZID or ENSEMBLID.")
    }

    colnames(dic.selected) <- c("ENSEMBLID", "SYMBOL")
    dicreference <- merge(resSymbol, dic.selected, by="SYMBOL")
    rownames(dicreference) <- dicreference$ENSEMBLID
    touse <- dicreference[genes.selected,]

    mol.unique.df$ENTREZID <- rownames(mol.unique.df)
    touse <- merge(mol.unique.df, touse, by="ENTREZID")
    rownames(touse) <- touse$ENSEMBLID
    touse <- subset(touse, select=-c(ENTREZID, SYMBOL))

    #label number of genes per cluster
    ngenes <- apply(touse, 2, function(x) length(which(x==1)))
    colnames(touse)[1:ncol(touse)-1] <- paste0(names(ngenes), " (", ngenes, ")")[1:ncol(touse)-1]

    meltedtouse <- melt(touse, id="ENSEMBLID")
    meltedtouse <- meltedtouse[which(meltedtouse$value==1), c("variable", "ENSEMBLID")]
    colnames(meltedtouse) <- c("Term", "Gene")

  } else {
    mol.unique.df$GTEx.median.TPM.Description <- toupper(rownames(mol.unique.df))
    dic.selected <- dic[which(dic$GTEx.median.TPM.Description %in% toupper(genes)),]
    genes.selected <- dic.selected$GTEx.median.TPM.Name

    if (length(genes.selected) == 0) {
      stop("The genes used to create the object are not correct. Please make sure they are introduced either as GENE SYMBOL, ENTREZID or ENSEMBLID.")
    }

    touse <- merge(mol.unique.df, dic.selected, by="GTEx.median.TPM.Description")
    touse <- subset(touse, select = -GTEx.median.TPM.Description)

    #label number of genes per cluster
    ngenes <- apply(touse, 2, function(x) length(which(x==1)))[1:ncol(touse)-1]
    newnames <- paste0(names(ngenes), " (", ngenes, ")")
    colnames(touse) <- c(newnames, colnames(touse[ncol(touse)]))

    meltedtouse <- melt(touse, id="GTEx.median.TPM.Name")
    meltedtouse <- meltedtouse[which(meltedtouse$value==1), c("variable", "GTEx.median.TPM.Name")]
    colnames(meltedtouse) <- c("Term", "Gene")

  }

  return(list(meltedtouse, touse))
}

#checks if the cluster independent slot is null
checkIndependent <- function(Object)
{
  if (is.null(Object@cIndependentMethod))
  {
    message("The cluster independent method is not defined. Please firstly run ClusterIndependentGeneSet in order to select the ordering method.")
    stop()
  }
}

# return the number of groups generated for a given minimum of elements
numberGroups <- function(res, nPathway)
{
  return(length(res[which(lapply(res, function(x) length(x)) >= nPathway)]))
}


#gets the optimal number of minimum elements for group
getOptimalNumber <- function(res)
{
  if (numberGroups(res, 4) == 0) {
    # message(paste0("There is no groups with at least 4 pathways. Using the minimum value, 2 pathways per group."))
    optimalNumber <- 2
  } else if (numberGroups(res, 5) == 0)
  {
    optimalNumber <- 4
  } else if (numberGroups(res, 6) == 0)
  {
    optimalNumber <- 5
  } else {
    if (numberGroups(res, 10) > 3)
    {
      optimalNumber <- 10
    } else if (numberGroups(res, 9) > 4 & 9 >= numberGroups(res, 9))
    {
      optimalNumber <- 9
    } else if (numberGroups(res, 8) > 4 & 9 >= numberGroups(res, 8))
    {
      optimalNumber <- 8
    } else if (numberGroups(res, 7) > 4 & 9 >= numberGroups(res, 7))
    {
      optimalNumber <- 7
    } else  if (numberGroups(res, 6) > 4)
    {
      optimalNumber <- 6
    } else {
      optimalNumber <- 4
    }
  }
  
  return(optimalNumber)
}

# get color palette
getPal <- function()
{
    pal.c <-  c(brewer.pal(n = 8, name ="Accent" ),
                brewer.pal(n = 8, name ="Dark2" ),
                brewer.pal(n = 8, name ="Set3"),
                brewer.pal(n = 8, name ="Set1"),
                brewer.pal(n = 8, name ="Set2"),
                brewer.pal(n = 8, name ="Greens"),
                brewer.pal(n = 8, name ="BrBG"))
    return(pal.c)
}


#This function performs the ora per cluster from a given object

PerformORAGOperCluster <- function(Object, uniquePathways, clusterIndependent)
{
  
  go_id_list_new <- obtainGOidList(Object, uniquePathways, clusterIndependent)
  #cluster-independent
  if (clusterIndependent)
  {
    clusters <- length(go_id_list_new)
    clus <- vector()
    for (i in 1:clusters) {
      paths <- c(rep(i, length(go_id_list_new[[i]])))
      names(paths) <- go_id_list_new[[i]]
      clus <- c(clus, paths)
    }
    #cluster-dependent
  } else {
    if (uniquePathways)
    {
      clus <- Object@plot$aka2Unique$Cluster
      names(clus) <- rownames(Object@plot$aka2Unique)
      clusters <- length(go_id_list_new)
    } else {
      clus <- Object@Data[[1]]$cluster
      names(clus) <- Object@Data[[1]]$Pathways
      clusters <- length(go_id_list_new)
    }
  }
  
  # message("\nCalculating ORA per cluster...")
  # keywords_ora <- ORAperCluster(Object, doORA=TRUE, clusters,
  #                               clus, completeRes=TRUE)
  
  message("\nPerforming Semantic enrichment per cluster...")
  if (checkGO(Object) == FALSE) {
    message("No GO terms have been detected in the pathways. The semantic enrichment word cloud will not be generated.")
    wordcloud <- FALSE
    term <- NULL
  } else {
    #adapted from anno_word_cloud_from_GO function
    env_tdm_GO <- readRDS(system.file("extdata", "tdm_GO.rds", package = "simplifyEnrichment"))
    names(go_id_list_new) <- as.character(1:length(go_id_list_new))
    
    # keyword enrichment
    message(paste0("Performing keyword enrichment for "), length(go_id_list_new), " group(s) of pathways.")
    term <- lapply(go_id_list_new, function(x, min_stat=5) {
      df <- keywordEnrichment(x, env_tdm_GO)
      df <- df[df$p <= min_stat, , drop = FALSE]
      data.frame(df[, 1], -log10(df$p))
    })
  }
  
  #return(list(ORA=keywords_ora, GO=term))
  return(list(GO=term))
  
}



obtainGOidList <- function(Object, uniquePathways, clusterIndependent)
{
  #cluster-independent
  if (clusterIndependent)
  {
    if (uniquePathways)
    {
      mat_sym <- scaleCorMatrix(Object@DataPathways.RR)
      #mat_cor <- mat_sym[Object@cIndependentMethod[[2]][[1]],
      #                   Object@cIndependentMethod[[2]][[1]]]
      mat_cor <- mat_sym[Object@cIndependentMethod[[1]][[1]],
                         Object@cIndependentMethod[[1]][[1]]]

      res <- obtainDefCluster(mat_cor)
      #go_id_list <- res[which(lapply(res,
      #                         function(x) length(x))>=Object@plot$aka3IndependentUnique$optimalNumberPathway)]
      
      go_id_list <- res[which(lapply(res,
                                     function(x) length(x))>=getOptimalNumber(res))]

    } else {
      mat_sym <- scaleCorMatrix(Object@Data.RR)
      mat_cor <- mat_sym[Object@cIndependentMethod[[1]][[1]],
                         Object@cIndependentMethod[[1]][[1]]]

      res <- obtainDefCluster(mat_cor)
      #go_id_list <- res[which(lapply(res,
      #                               function(x) length(x))>=Object@plot$aka3IndependentUnique$optimalNumber)]
      go_id_list <- res[which(lapply(res,
                                     function(x) length(x))>=getOptimalNumber(res))]
      go_id_list_new <- go_id_list
    }
    
    
    if (uniquePathways)
    {
      go_id_list_new <- go_id_list
    } else {
      go_id_list_new <- list()
      for (i in 1:length(go_id_list)) {
        go_id_list_new[[i]] <- Object@Data[[1]]$Pathways[Object@Data[[1]]$RR_name %in% go_id_list[[i]]]
      }
    }

  #cluster-dependent
  } else {
    if(uniquePathways == F)
    {
      clus <- Object@Data[[1]]$cluster
      names(clus) <- Object@Data[[1]]$Pathways
    }else{
      clus <- Object@plot$aka2Unique$Cluster
      names(clus) <- rownames(Object@plot$aka2Unique)
    }
    go_id_list <- list()
    for (z in unique(clus))
    {
      go_id_list[[z]] <- names(clus[clus==z])
    }
    go_id_list_new <- go_id_list

  }

  return(go_id_list_new)
}


obtainGREATtable <- function(GREATobject)
{
  if (!class(GREATobject)[1] == "GreatObject")
  {
    stop("The introduced object is not GREAT object.")
  }

  great_1 <- rGREAT::getEnrichmentTable(Trans1_GREAT)
  zgreat_1 <- great_1[great_1$p_adjust<0.05,]
  pathways_1 <- great_1[, "id"]
  genes_1 <- rGREAT::getRegionGeneAssociations(Trans1_GREAT, use_symbols = FALSE)
  genes_1 <- unique(unlist(genes_1$annotated_genes))
  gene_sets1 <- data.frame(id=c(0), Genes=c(0))
  for(i in 1:length(Trans1_GREAT@gene_sets)){
    id <- names(Trans1_GREAT@gene_sets[i])
    genes <- Trans1_GREAT@gene_sets[[i]][which(Trans1_GREAT@gene_sets[[i]] %in% genes_1)]
    gene_sets1[i,] <- c(id, paste(genes, collapse = "/"))
  }

  return(merge(great_1, gene_sets1, by="id"))
}


printTemplate <- function()
{
  message("The TEMPLATE should be a data frame with the following structure...")
  message("The following column names: Pathways, Pvalues, EnrichmentScore, Molecules. Here is an example.")
  dfexample <- data.frame(Pathways = c("GO:0043269","GO:0010273","GO:1990169","GO:0014888","GO:0046323","GO:0071280"),
                          Pvalues = c(5.584191e-06,1.144264e-05,1.144264e-05,1.340369e-05,1.643183e-05,2.235679e-05),
                          EnrichmentScore = c(-0.3426331,0.8733950,0.873395,0.6782518,-0.5851564,0.7640066),
                          Molecules = c("SYT15/SHISA7/PTPN6/KCNJ4/CLIC1/EDNRA/GRM7/VAMP8/ABCB1/JSRP1/SLC8A1",
                                        "MT1X/MT1H/MT1M/MT1G/MT2A/MT1E/MT1F/MT1DP",
                                        "KLF15/TRIM63/ERRFI1/MYOG/MYOZ2/FOXO1/MLIP/MSTN/ATP2A2/CFLAR/PPARG",
                                        "SLC2A10/SLC27A4/HK2/PTPN11/SLC2A1/RAB4B/CLTCL1/SLC2A14/GRB10/RTN2",
                                        "CYP1A1/MT1X/MT1H/MT1M/MT1G/MT2A/MT1E/MT1F/MT1DP",
                                        "KLF15/LEP/IRS2/PRKAG2/SORBS1/PIK3R1/ENPP1/MYC/IRS1/OPN3"))

  print(knitr::kable(dfexample, format = "markdown"))

}


optimalDist <- function(mat_sym, all_methods=FALSE, use_method=NULL) {
  set.seed(123)
  if (is.null(use_method))
  {

    seriation_methods <- c("GW", "GW_average", "GW_complete", "GW_single", "GW_ward", "HC",
                           "HC_average", "HC_complete", "HC_single", "HC_ward", "isoMDS",
                           "MDS", "MDS_angle", "monoMDS", "OLO", "OLO_average", "OLO_complete",
                           "OLO_single", "OLO_ward", "QAP_2SUM", "QAP_BAR", "R2E", "Sammon_mapping",
                           "Spectral", "Spectral_norm", "TSP", "Identity", "SPIN_NH", "SPIN_STS")

    seriation_res <- apply_seriation(mat_sym,                     # apply seriation methods
                                     seriation_methods)
    statistics_res <- extract_statistics(seriation_res, mat_sym)  # obtain statistics
    statistics_res_scaled <- generate_scale(statistics_res)

    # weights for combined 40% BBanded anti-Robinson form criterion (BAR), 40% Hamiltonian path length and 20% total number of gene-set in clusters
    statistics_res_scaled$combined <- 0.4*statistics_res_scaled$BAR+0.4*statistics_res_scaled$Path_length+0.2*(1-statistics_res_scaled$TotalElements)
    use_method <- rownames(statistics_res_scaled[order(statistics_res_scaled$combined, decreasing = FALSE),])[1]
    optimal_seriation <- seriation_res[use_method]
    i_def <- get_order(optimal_seriation[[1]])
    mode <- "automatic_selection"


  } else {
    if (use_method %in% list_seriation_methods('dist'))
    {
      d <- as.dist(1 - mat_sym)
      seriation_res <- seriate(d, use_method)
      i_def <- get_order(seriation_res)
      mode <- "manual_selection"


    } else {
      message("The selected method is not correct or it is not among the possible methods.")
      stop()
    }
  }

  out <- list(i_def, use_method, mode)
  names(out) <- c("Order", "Used method", "Mode")
  return(out)
}


extract_numberClusters <- function(mat_sym, s) {
  order <- get_order(s)
  res <- obtainDefCluster(mat_sym[order,order])
  go_id_list <- res[which(lapply(res, function(x) length(x)) >= 6)]
  cluster_num <- unlist(lapply(go_id_list, length))

  return(cluster_num)
}


apply_seriation <- function(mat_sym, seriation_methods)
{
  d <- as.dist(1- mat_sym)
  o <- pbsapply(seriation_methods, FUN = function(m) seriate(d, m))
  o <- ser_align(o)

  return(o)
}


extract_statistics <- function(o, mat_sym)
{
  d <- as.dist(1- mat_sym)
  message("\nComputing the loss functions...")
  crit <- pbsapply(o, FUN = function(x) criterion(d, x))
  clust_size <- sapply(o, FUN = function(x) extract_numberClusters(mat_sym, x))
  df_include <- data.frame(Ncluster=sapply(clust_size, length), TotalElements=sapply(clust_size, sum),
                           meanCluster=sapply(clust_size, mean))
  df_out <- rbind(crit, t(df_include))

  return(t(df_out))
}


generate_scale <- function(df)
{
  df <- as.data.frame(df)
  scaled_df <- as.data.frame(lapply(df, function(x) {
    min_val <- min(x)
    max_val <- max(x)
    scaled <- (x - min_val) / (max_val - min_val)
    return(scaled)
  }))

  colnames(scaled_df) <- colnames(df)
  rownames(scaled_df) <- rownames(df)

  return(scaled_df)
}



