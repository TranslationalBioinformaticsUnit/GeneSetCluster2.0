#' @title OptimalGeneSets
#'
#' @description
#' There are many ways to determine what is the best separation of the clusters.
#' The function OptimalGeneSets can help with determining what is the best separation.
#' There are 3 different statistics the function can plot that help with this; Gap, Elbow and Silhouette.
#' Gap: Compares the total within intra-cluster variation for different values of k with their expected values under null reference distribution of the data
#' Elbow: For each k, calculate the total within cluster sum of squares
#' Silhouette: Determines how wel each object lies within its cluster. Higher the better.
#'
#' @import cluster
#' @import factoextra
#' @import ggplot2
#'
#' @param Object A PathwayObject.
#' @param uniquePathways A boolean to specify unique Gene-set or unique Pathways.
#' FALSE if for unique Gene-set and TRUE for unique Pathways.
#' unique Pathway if it is TRUE.
#' @param cluster the cluster for sub-clustering.
#' If NA calculates the optimal for all the data.
#' @param method Which method to determing optimal number of clusters. gap, elbow or silhouette.
#' @param max_cluster Max number of clusters to test
#' @param cluster_method kmeans or hcut. Which clustering method is used
#' @param main A string to be used as title in the plot
#'
#' @return plot
#'
#' @export
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
#' OptimalGeneSets(Object = GSEA.Object2, uniquePathways = FALSE, cluster=NULL,
#'                 method = "silhouette", max_cluster= 24, cluster_method = "kmeans",
#'                 main= "Kmeans for 24 clusters")
#'
#' OptimalGeneSets(Object = GSEA.Object2, uniquePathways = TRUE, cluster=NULL,
#'                 method = "silhouette", max_cluster= 24, cluster_method = "kmeans",
#'                 main= "Kmeans for 24 clusters")
#'
#' GSEA.Object3 <- ClusterGeneSets(Object = GSEA.Object2,
#'                 clusters = 3,
#'                 method = "kmeans")
#'
#' OptimalGeneSets(Object = GSEA.Object3, uniquePathways = TRUE, cluster= 1,
#'                 method = "silhouette", max_cluster= 24, cluster_method = "kmeans",
#'                 main= "Kmeans for 24 clusters for Cluster 1")
#' }
#'
setGeneric(name="OptimalGeneSets",
           def=function(Object, uniquePathways=FALSE, cluster=NULL,
                        method = "silhouette", max_cluster=24,
                        cluster_method="kmeans", main = "")
           {
             standardGeneric("OptimalGeneSets")
           }
)


setMethod(f="OptimalGeneSets",
          definition = function(Object, uniquePathways=FALSE, cluster=NULL,
                                method = "silhouette", max_cluster=24,
                                cluster_method="kmeans", main = "")
{

  message("[=========================================================]")
  message("[<<<<            OptimalGeneSets START               >>>>>]")
  message("-----------------------------------------------------------")
  if(method %in% c("Elbow", "elbow")){method <- "elbow"}
  if(method %in% c("Silhouette", "silhouette")){method <- "silhouette"}
  if(method %in% c("Gap", "gap")){method <- "gap"}
  if(cluster_method %in% c("HierArchical", "hierarchical", "Hierarchical", "hclust")){cluster_method <- "hcut"}


  message("Finding optimal number of clusters")
  message(paste("Clustering method = ", cluster_method))
  message(paste("Optimizing method = ", method))


  if (uniquePathways == TRUE) {
    if (is.null(cluster)) {
      RRData <- Object@DataPathways.RR
    } else {
      if (any(cluster %in% Object@plot$aka2Unique$Cluster))
      {
        RRData <- Object@DataPathways.RR[which(Object@plot$aka2Unique$Cluster == cluster)]
      } else {
        stop("The introduced cluster is not correct.")
      }
    }
  } else {
    if (is.null(cluster)) {
      RRData <- Object@Data.RR
    } else {
      if (any(cluster %in% Object@plot$aka2$Cluster))
      {
        RRData <- Object@Data.RR[which(Object@plot$aka2$Cluster == cluster)]
      } else {
        stop("The introduced cluster is not correct.")
      }
    }
  }

  ######################################
  ##---------elbow--------------------##
  ######################################

  if(method == "elbow")
  {
    plot.x <-   fviz_nbclust(RRData, get(cluster_method), method = "wss", k.max = max_cluster) +ggtitle(paste(main,": The Elbow method",sep=""))
  }
  ######################################
  ##---------gap----------------------##
  ######################################

  if(method == "gap")
  {
    gap_stat <- clusGap(RRData, get(cluster_method), K.max = max_cluster, B = 50)
    plot.x <- fviz_gap_stat(gap_stat) + ggtitle(paste(main,": Gap Statistics",sep=""))
  }

  ######################################
  ##---------silhoutte----------------##
  ######################################

  if(method == "silhouette")
  {
    plot.x <-  fviz_nbclust(x = RRData, get(cluster_method), method = "silhouette", k.max = max_cluster) + ggtitle(paste(main,": The Silhouette Plot",sep=""))
  }

  message("-----------------------------------------------------------")
  message("[<<<<<             OptimalGeneSets END              >>>>>>]")
  message("[=========================================================]")
  message("[You may want to process ClusterGeneSets next.            ]")

  return(plot(plot.x))

})
