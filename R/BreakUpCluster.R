#' @title BreakUpCluster
#' @description
#' Split a cluster in a desired number of subclusters.
#'
#' @import methods
#' @import RColorBrewer
#' @importFrom limma strsplit2
#' @importFrom stats kmeans hclust cutree dist
#'
#' @param Object A PathwayObject
#' @param breakup.cluster The cluster to be broken up
#' @param sub.cluster The amount of subclusters to be generated
#' @param uniquePathways A boolean to specify the approach of unique Gene-set or unique Pathways.
#' FALSE if for unique Gene-set and TRUE for unique Pathways.
#'
#' @return PathwayObject
#'
#' @examples
#' GSEA.files <- c(system.file("extdata", "Covid19AI_Healthy_GSEA.csv",
#'                              package = "GeneSetCluster"),
#'                 system.file("extdata", "Covid193Mo_Healthy_GSEA.csv",
#'                              package = "GeneSetCluster"),
#'                 system.file("extdata", "Covid196Mo_Healthy_GSEA.csv",
#'                              package = "GeneSetCluster"))
#'
#' GSEA.Object1 <- LoadGeneSets(file_location = GSEA.files,
#'                               groupnames= c("IA", "Mo3", "Mo6"),
#'                               P.cutoff = 0.05,
#'                               Mol.cutoff = 5,
#'                               Source = "GSEA",
#'                               Great.Background = FALSE,
#'                               type = NA,
#'                               topranks = NA,
#'                               structure = "ENTREZID",
#'                               Organism = "org.Hs.eg.db",
#'                               seperator = "/")
#' \dontrun{
#' GSEA.Object2 <- CombineGeneSets(Object = GSEA.Object1, threads=2)
#'
#' GSEA.Object3 <- ClusterGeneSets(Object = GSEA.Object2,
#'                 clusters = 3,
#'                 method = "kmeans")
#'
#' OptimalGeneSets(Object = GSEA.Object3, uniquePathway = FALSE, cluster = 2,
#'                 method = "silhouette", max_cluster= 24, cluster_method = "kmeans",
#'                 main= "Kmeans for 24 clusters of Cluster 2 subclustering Gene-set")
#'
#' GSEA.Object3.1 <- BreakUpCluster(GSEA.Object3, breakup.cluster=2, sub.cluster=2, uniquePathways=FALSE)
#' }
#'
#' @export
#'
setGeneric(name="BreakUpCluster",
           def=function(Object = Object, breakup.cluster=1, sub.cluster=3, uniquePathways=FALSE)
           {
             standardGeneric("BreakUpCluster")
           }
)

setMethod(f="BreakUpCluster",
          definition = function(Object = Object, breakup.cluster=1, sub.cluster=3, uniquePathways=FALSE)
{

  message("[=========================================================]")
  message("[<<<<            BreakUpCluster START                >>>>>]")
  message("-----------------------------------------------------------")


  if(is.na(Object@metadata$cluster.method[1]))
  {
    message("Make sure youre object has been combined by CombineGeneSets")
    stop()
  }


  order <- Object@metadata$order.group[1]
  method <- Object@metadata$cluster.method[1]


  #############################################
  ###----Kmeans Cluster---------------------###
  #############################################

  if(Object@metadata$cluster.method[1] == "kmeans")
  {

    if (uniquePathways)
    {

      # breakup.cluster<-sub("Cluster_","",breakup.cluster)
      canonical.df <- Object@plot$aka2Unique
      RR <- Object@DataPathways.RR
      Data.RR.clus <- RR[Object@plot$aka2Unique$Cluster==breakup.cluster,Object@plot$aka2Unique$Cluster==breakup.cluster]
      canonical.df$Cluster[Object@plot$aka2Unique$Cluster==breakup.cluster] <- paste0(breakup.cluster, ".",kmeans(x = Data.RR.clus, centers = sub.cluster)$cluster)
      Object@plot$aka2Unique <- canonical.df

    } else {
      canonical.df <- Object@Data[[1]]

      Data.RR.clus <- Object@Data.RR[Object@Data[[1]]$cluster==breakup.cluster,Object@Data[[1]]$cluster==breakup.cluster]
      canonical.df$cluster[Object@Data[[1]]$cluster==breakup.cluster] <- paste0(breakup.cluster, ".",kmeans(x = Data.RR.clus, centers = sub.cluster)$cluster)

      RR<-Object@Data.RR
      #Data.RR.clus <- Object@Data.RR[Object@Data[[1]]$cluster==breakup.cluster,Object@Data[[1]]$cluster==breakup.cluster]
      Data.RR.clus <- RR[Object@Data[[1]]$cluster==breakup.cluster,Object@Data[[1]]$cluster==breakup.cluster]
      canonical.df$cluster[Object@Data[[1]]$cluster==breakup.cluster] <- paste0(breakup.cluster, ".",kmeans(x = Data.RR.clus, centers = sub.cluster)$cluster)
      Object@Data[[1]] <- canonical.df


      #Object@metadata[,"cluster.method"] <- paste0(Object@metadata[,"cluster.method"], "_breakup")
    }
  }

  #############################################
  ###----Hierarchical Cluster---------------###
  #############################################
  if(Object@metadata$cluster.method[1] == "Hierarchical")
  {
    if(uniquePathways){
      breakup.cluster<-sub("Cluster_","",breakup.cluster)
      canonical.df<-Object@plot$aka2Unique
      RR<-Object@DataPathways.RR
      Data.RR.clus <- RR[Object@plot$aka2Unique$Cluster==breakup.cluster,Object@plot$aka2Unique$Cluster==breakup.cluster]
      temp.clx <- cutree(hclust(dist(t(Data.RR.clus)), method = "ward.D2"), k = sub.cluster)
      canonical.df$Cluster[Object@plot$aka2Unique$Cluster==breakup.cluster] <- paste0(breakup.cluster, ".",temp.clx)
      Object@plot$aka2Unique <- canonical.df

    }else{
      canonical.df <- Object@Data[[1]]
      RR<-Object@Data.RR
      #Data.RR.clus <- Object@Data.RR[Object@Data[[1]]$cluster==breakup.cluster,Object@Data[[1]]$cluster==breakup.cluster]
      Data.RR.clus <- RR[Object@Data[[1]]$cluster==breakup.cluster,Object@Data[[1]]$cluster==breakup.cluster]
      temp.clx <- cutree(hclust(dist(t(Data.RR.clus)), method = "ward.D2"), k = sub.cluster)
      canonical.df$cluster[Object@Data[[1]]$cluster==breakup.cluster] <- paste0(breakup.cluster, ".",temp.clx)
      Object@Data[[1]] <- canonical.df
    }


    #Object@metadata[,"cluster.method"] <- paste0(Object@metadata[,"cluster.method"], "_breakup")
  }

  ############################################
  #---------Adding cluster info--------------#
  ############################################
  if(!uniquePathways){
    canonical.df$mean.RR.cl <- NA
    canonical.df$sum.RR.cl <- NA


    for(cl.i in unique(canonical.df$cluster))
    {
      DF.cl <- canonical.df[canonical.df$cluster == cl.i,]
      idx <- DF.cl$RR_name
      rr.x <- Object@Data.RR[idx, idx]
      rr.x <- data.frame(sapply(rr.x, function(x) as.numeric(as.character(x))))
      rr.x <- as.matrix(rr.x)
      DF.cl$mean.RR.cl <- round(mean(rr.x),digits = 4)
      DF.cl$sum.RR.cl <- round(sum(rr.x),digits = 4)
      canonical.df[canonical.df$cluster == cl.i,] <- DF.cl

    }
  }
  ###############################
  #--------Order clusters-------#
  ###############################
  message("Ordering pathway clusters")
  if(!uniquePathways){
    if(order == "group")
    {
      Object@Data <- list(canonical.df[order(canonical.df$Group,canonical.df$cluster,  decreasing = F),])

      Object@metadata[,"order.group"] <- rep(order, times = nrow(Object@metadata))

    }
    if(order == "cluster")
    {
      Object@Data <- list(canonical.df[order(canonical.df$cluster,  decreasing = F),])

      Object@metadata[,"order.group"] <- rep(order, times = nrow(Object@metadata))
    }
    if(order == "mean.RR")
    {
      Object@Data <- list(canonical.df[order(canonical.df$mean.RR.cl,  decreasing = F),])

      Object@metadata[,"order.group"] <- rep(order, times = nrow(Object@metadata))

    }
    if(order == "Sum.RR")
    {
      Object@Data <- list(canonical.df[order(canonical.df$sum.RR.cl,  decreasing = F),])

      Object@metadata[,"order.group"] <- rep(order, times = nrow(Object@metadata))

    }
    Object@metadata[,"cluster.method"] <- rep(method, times = nrow(Object@metadata))


    Object@Data.RR <- Object@Data.RR[Object@Data[[1]]$RR_name,Object@Data[[1]]$RR_name]
  }

  ################################################
  ####-----------Plotting Info----------------####
  ################################################


  pal.c <-  c(brewer.pal(n = 8, name ="Accent" ),
              brewer.pal(n = 8, name ="Dark2" ),
              brewer.pal(n = 8, name ="Set3"),
              brewer.pal(n = 8, name ="Set1"))

  if(uniquePathways){
    #unique pathways
    aka2Unique <- Object@plot$aka2Unique[,-ncol(Object@plot$aka2Unique)]
    if(length(unique(Object@plot$aka2Unique$Cluster)) > 32)
    {
      message("Warning: number of cluster is larger than colours supported")
    }
    groups.col <- brewer.pal(n = 8, name ="Set2" )[1:ncol(aka2Unique)]
    names(groups.col) <- colnames(aka2Unique)
    pal <- pal.c[1:length(unique(Object@plot$aka2Unique$Cluster))]
    names(pal) <- sub("Cluster_","",unique(Object@plot$aka2Unique$Cluster))


    aka3Unique <- list(Group = groups.col,
                      Cluster= pal)
    names(aka3Unique) <-  c("Group", "Cluster")

    ##################################
    #-----------Return---------------#
    ##################################
    Object@plot$aka3Unique = aka3Unique

  }else{
    aka2 <- matrix(data = NA, nrow = nrow((Object@Data.RR)), ncol = 2)
    colnames(aka2) <- c("Group", "Cluster")
    rownames(aka2) <- Object@Data[[1]]$RR_name

    aka2[,"Group"] <- as.character(Object@Data[[1]]$Groups)
    aka2[,"Cluster"] <- sub("Cluster_","",as.character(Object@Data[[1]]$cluster))

    aka2 <- as.data.frame(aka2)
    if(length(unique(Object@Data[[1]]$cluster)) > 32)
    {
      message("Warning: number of cluster is larger than colours supported")
    }


    groups.col <- brewer.pal(n = 8, name ="Set2" )[1:length(unique(aka2[,"Group"]))]
    names(groups.col) <- unique(aka2[,"Group"])
    pal <- pal.c[1:length(unique(aka2[,"Cluster"]))]
    names(pal) <- sub("Cluster_","",unique(aka2[,"Cluster"]))


    aka3 <- list(Group = groups.col,
                Cluster= pal)
    names(aka3) <-  c("Group", "Cluster")

    ##################################
    #-----------Return---------------#
    ##################################

    Object@plot$aka2 = aka2
    Object@plot$aka3 = aka3
  }

  if (uniquePathways)
  {
    message("\nCalculating ORA and Semantic enrichment per cluster for Unique Gene-sets...")
    functional_annotationUnique <- PerformORAGOperCluster(Object, uniquePathways = TRUE, clusterIndependent = FALSE)
    functional_annotation <- Object@functionalAnnot$Geneset

  } else {
    message("\nCalculating ORA and Semantic enrichment per cluster for Unique Pathways...")
    functional_annotation <- PerformORAGOperCluster(Object, uniquePathways = FALSE, clusterIndependent = FALSE)
    functional_annotationUnique <- Object@functionalAnnot$Pathway
  }


  Object@functionalAnnot <- list(Geneset=functional_annotation,
                                 Pathway=functional_annotationUnique)


  message("-----------------------------------------------------------")
  message("[<<<<<             BreakUpCluster END               >>>>>>]")
  message("[=========================================================]")
  message("[You may want to process HighlightGeneSets next.           ]")
  message("[You may want to plot the results using PlotGeneSets next. ]")


  return(Object)
})



