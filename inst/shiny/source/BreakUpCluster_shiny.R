#' BreakUpCluster_shiny
#'
#' @import stats
#' @import methods
#' @import limma
#'
#' @importFrom stats dist hclust kmeans
#'
#' @param Object A PathwayObject
#' @param breakup.cluster An integer with the number of the cluster to be broken up
#' @param sub.cluster The amount of subclusters to be generated
#' @param uniquePathways Boolean to merge or not unique pathways
#' @return PathwayObject
#' @export
#'
#' @examples
#' Great.files <- c(system.file("extdata", "MM10.GREAT.KO.uGvsMac.bed.tsv",
#'                              package = "GeneSetCluster"),
#' system.file("extdata", "MM10.GREAT.KO.uGvsMac.bed_BCKGRND.tsv", package = "GeneSetCluster"),
#' system.file("extdata", "MM10.GREAT.WT.uGvsMac.bed.tsv", package = "GeneSetCluster"),
#' system.file("extdata", "MM10.GREAT.WT.uGvsMac.bed_BCKGRND.tsv", package = "GeneSetCluster"))
#' Great.files.bckgrnd <- Great.files[grepl("BCKGRND", Great.files)]
#'
#'
#' Great.bckgnrd.Object1 <- LoadGeneSets(file_location = Great.files.bckgrnd,
#'                                       groupnames= c("KO", "WT"),
#'                                       P.cutoff = 0.05,
#'                                       Mol.cutoff = 5,
#'                                       Source = "Great",
#'                                       Great.Background = TRUE,
#'                                       type = "Canonical_Pathways",
#'                                     topranks = 20,
#'                                    structure = "SYMBOL",
#'                                    Organism = "org.Mm.eg.db",
#'                                    seperator = ",")
#' man.Great.Object1 <- ManageGeneSets(Object = Great.bckgnrd.Object1,
#'                                    keep.type =c("Disease Ontology",
#'                                    "GO Biological Process" ),
#'                                    exclude.type="")
#' man.Great.Object2 <- CombineGeneSets(Object = man.Great.Object1)
#' man.Great.Object3 <- ClusterGeneSets(Object = man.Great.Object2,
#'                                      clusters = 5,
#'                                      method = "kmeans")
#'                                      
#' man.Great.Object3 <- (Object = man.Great.Object3, breakup.cluster = 3, sub.cluster=3)                                     
#'                                      

BreakUpCluster_shiny <- function(Object = Object, breakup.cluster = 6, sub.cluster=3, uniquePathways=F)
{
  
  message("[=========================================================]")
  message("[<<<<            BreakUpCluster START                >>>>>]")
  message("-----------------------------------------------------------")
  
  require(RColorBrewer)
  
  
  warning("Currently only working for Kmeans and Hierarchical")
  if(is.na(Object@metadata$cluster.method[1]))
  {
    message("Make sure youre object has been combined by CombineGeneSets")
    stop()
  }  
  canonical.df <- Object@Data[[1]]
  
  
  order <- Object@metadata$order.group[1]
  method <- Object@metadata$cluster.method[1]
  
  
  #############################################
  ###----Kmeans Cluster--------------------###
  #############################################
  
  if(Object@metadata$cluster.method[1] == "kmeans")
  {
    breakup.cluster2<-sub("Cluster_","",breakup.cluster)
    if(uniquePathways){
      canonical.df<-Object@plot$aka2Unique
      RR<-Object@DataPathways.RR
      Data.RR.clus <- RR[Object@plot$aka2Unique$Cluster==breakup.cluster2,Object@plot$aka2Unique$Cluster==breakup.cluster2]
      temp.clx<-kmeans(x = Data.RR.clus, centers = sub.cluster)$cluster
      canonical.df$Cluster[Object@plot$aka2Unique$Cluster==breakup.cluster2] <- paste0(breakup.cluster2, ".",temp.clx)
      Object@plot$aka2Unique<-canonical.df
    }else{
      canonical.df<-Object@plot$aka2
      RR<-Object@Data.RR
      Data.RR.clus <- RR[Object@Data[[1]]$cluster==breakup.cluster,Object@Data[[1]]$cluster==breakup.cluster]
      temp.clx<-kmeans(x = Data.RR.clus, centers = sub.cluster)$cluster
      canonical.df$Cluster[Object@plot$aka2$Cluster==breakup.cluster2] <- paste0(breakup.cluster2, ".",temp.clx)
      Object@Data[[1]]$cluster[Object@Data[[1]]$cluster==breakup.cluster] <- paste0(breakup.cluster, ".",temp.clx)
      Object@plot$aka2<-canonical.df
      
    }
    
    #Data.RR.clus <- Object@Data.RR[Object@Data[[1]]$cluster==breakup.cluster,Object@Data[[1]]$cluster==breakup.cluster]
    #canonical.df$cluster[Object@Data[[1]]$cluster==breakup.cluster] <- paste0(breakup.cluster, ".",kmeans(x = Data.RR.clus, centers = sub.cluster)$cluster)
    
    #Object@metadata[,"cluster.method"] <- paste0(Object@metadata[,"cluster.method"], "_breakup")
  }
  
  #############################################
  ###----Hierarchical Cluster--------------------###
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
      Object@plot$aka2Unique<-canonical.df
      
    }else{
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
      #Object@Data <- list(canonical.df[order(canonical.df$Cluster,  decreasing = F),])
      
      #Object@metadata[,"order.group"] <- rep(order, times = nrow(Object@metadata))
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
  }else{
    #unique pathways = TRUE and ordering by cluster Object@Data
    
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
    aka2Unique<-Object@plot$aka2Unique[,-ncol(Object@plot$aka2Unique)]
    if(length(unique(Object@plot$aka2Unique$Cluster)) > 32)
    {
      message("Warning: number of cluster is larger than colours supported")
    }
    groups.col <- brewer.pal(n = 8, name ="Set2" )[1:ncol(aka2Unique)]
    names(groups.col) <- colnames(aka2Unique)
    pal <- pal.c[1:length(unique(Object@plot$aka2Unique$Cluster))]
    names(pal) <- sub("Cluster_","",unique(Object@plot$aka2Unique$Cluster))
    
    
    aka3Unique = list(Group = groups.col,
                Cluster= pal)
    names(aka3Unique) <-  c("Group", "Cluster")
    
    ##################################
    #-----------Return---------------#
    ##################################
    Object@plot$aka3Unique = aka3Unique
    
  }else{
    
    
    #aka2 <- matrix(data = NA, nrow = nrow((Object@Data.RR)), ncol = 2)
    #colnames(aka2) <- c("Group", "Cluster")
    #rownames(aka2) <- Object@Data[[1]]$RR_name
    
    #aka2[,"Group"] <- as.character(Object@Data[[1]]$Groups)
    #aka2[,"Cluster"] <- sub("Cluster_","",as.character(Object@Data[[1]]$cluster))
    
    #aka2 <- as.data.frame(aka2)
    aka2 <- as.data.frame(Object@plot$aka2)
    if(length(unique(Object@Data[[1]]$cluster)) > 32)
    {
      message("Warning: number of cluster is larger than colours supported")
    }
    
    
    groups.col <- brewer.pal(n = 8, name ="Set2" )[1:length(unique(aka2[,"Group"]))]
    names(groups.col) <- unique(aka2[,"Group"])
    pal <- pal.c[1:length(unique(aka2[,"Cluster"]))]
    names(pal) <- sub("Cluster_","",unique(aka2[,"Cluster"]))
    
    aka3 = list(Group = groups.col,
                    Cluster= pal)
    names(aka3) <-  c("Group", "Cluster")

    ##################################
    #-----------Return---------------#
    ##################################

    Object@plot$aka2 = aka2
    Object@plot$aka3 = aka3
  }
  
  message("-----------------------------------------------------------")
  message("[<<<<<             BreakUpCluster END               >>>>>>]")
  message("[=========================================================]")
  message("[You may want to process HighlightGeneSets next.           ]")
  message("[You may want to plot the results using PlotGeneSets next. ]")
  
  
  return(Object)
}



