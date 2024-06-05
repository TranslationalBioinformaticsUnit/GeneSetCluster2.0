OptimalGeneSets_2 <- function(object, method, max_cluster, cluster_method, main, uniquePathways)
{

  message("[=========================================================]")
  message("[<<<<            OptimalGeneSets_2 START               >>>>>]")
  message("-----------------------------------------------------------")
  if(method %in% c("Elbow", "elbow")){method <- "elbow"}
  if(method %in% c("Silhouette", "silhouette")){method <- "silhouette"}
  if(method %in% c("Gap", "gap")){method <- "gap"}
  if(cluster_method %in% c("HierArchical", "hierarchical", "Hierarchical", "hclust")){cluster_method <- "hcut"}

  message("Finding optimal number of clusters")
  message(paste("Clustering method = ", cluster_method))
  message(paste("Optimizing method = ", method))



  ######################################
  ##---------elbow--------------------##
  ######################################
  if(uniquePathways){
    RR<-object@DataPathways.RR
  }else{
    RR<-object@Data.RR
  }

  if(method == "elbow")
  {
    plot.x <-   fviz_nbclust(RR, get(cluster_method), method = "wss", k.max = max_cluster) +ggtitle(paste(main,": The Elbow method",sep=""))
  }
  ######################################
  ##---------gap----------------------##
  ######################################

  if(method == "gap")
  {
    gap_stat <- clusGap(RR, get(cluster_method), K.max = max_cluster, B = 50)
    plot.x <-    fviz_gap_stat(gap_stat) + ggtitle(paste(main,": Gap Statistics",sep=""))
  }

  ######################################
  ##---------silhoutte----------------##
  ######################################

  if(method == "silhouette")
  {
    plot.x <-  fviz_nbclust(x = RR,get(cluster_method), method = "silhouette", k.max = max_cluster) + ggtitle(paste(main,": The Silhouette Plot",sep=""))
  }

  message("-----------------------------------------------------------")
  message("[<<<<<             OptimalGeneSets_2 END              >>>>>>]")
  message("[=========================================================]")
  message("[You may want to process ClusterGeneSets next.            ]")

  return(plot.x)

}
