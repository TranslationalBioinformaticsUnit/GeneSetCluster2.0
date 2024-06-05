#' @title PlotGeneNetworks
#'
#' @description
#' Plots a heatmap with the distances per cluster
#'
#' @import factoextra
#' @import network
#' @import sna
#' @importFrom GGally ggnet2
#'
#' @param Object A PathwayObject.
#' @param uniquePathways A boolean to specify the approach of unique Gene-set or unique Pathways.
#' FALSE if for unique Gene-set and TRUE for unique Pathways.
#' @param clusterIndependent A boolean to specify the clustering approach.
#' TRUE for seriation-based clustering approach and FALSE for classic approach.
#' @param labels Boolean. TRUE to labels, or a vector with labels to be shown
#' @param RRmin filter for different RR cutoff for what is an edge. Minimum value 0 and Maximum value 1.
#'
#' @return plot
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
#' GSEA.Object3 <- ClusterGeneSets(Object = GSEA.Object2,
#'                 clusters = 3,
#'                 method = "kmeans")
#' GSEA.Object3 <- ClusterIndependentGeneSet(Object = GSEA.Object3)
#'
#' PlotGeneNetworks(GSEA.Object3, uniquePathways=FALSE, clusterIndependent=FALSE)
#' PlotGeneNetworks(GSEA.Object3, uniquePathways=FALSE, clusterIndependent=TRUE)
#' }
#' @export
#'
setGeneric(name="PlotGeneNetworks",
           def=function(Object, uniquePathways=FALSE, clusterIndependent=FALSE,
                        labels =FALSE, RRmin = 0.1)
           {
             standardGeneric("PlotGeneNetworks")
           }
)

setMethod(f="PlotGeneNetworks",
          signature="PathwayObject",
          definition=function(Object, uniquePathways=FALSE, clusterIndependent=FALSE,
                              labels =FALSE, RRmin = 0.1)
          {

            if (uniquePathways)
            {
              RR.data <- scaleCorMatrix(Object@DataPathways.RR)

              if (clusterIndependent)
              {
                clusters <- Object@plot$aka2IndependentUnique[colnames(RR.data),]
                y <- Object@plot$aka3IndependentUnique$Cluster[order(names(Object@plot$aka3IndependentUnique$Cluster))]
                y["0"] <- "grey"


              } else {
                clusters <- Object@plot$aka2Unique[colnames(RR.data),]
                y <- Object@plot$aka3Unique$Cluster

              }

            } else {
                RR.data <- scaleCorMatrix(Object@Data.RR)

              if (clusterIndependent)
              {
                clusters <- Object@plot$aka2Independent[colnames(RR.data),]
                y <- Object@plot$aka3Independent$Cluster
                y["0"] <- "grey"


              } else {
                clusters <- Object@plot$aka2[colnames(RR.data),]
                y <- Object@plot$aka3$Cluster
              }
            }


            rownames(RR.data) <- paste("Cluster", clusters$Cluster,sep="")
            colnames(RR.data) <- paste("Cluster", clusters$Cluster,sep="")

            for(i in 1:nrow(RR.data))
            {
              RR.data[i,RR.data[i,] <= RRmin] <- 0
            }

            RR.net = network::network(RR.data,  directed = TRUE)

            # Cluster affiliation
            x <- as.factor(paste("Cluster", clusters$Cluster, sep=""))
            RR.net %v% "Cluster" = as.character(x)

            # color palette
            names(y) <- paste("Cluster", names(y), sep="")


            if(sum(labels %in% T) == 1)
            {
              labels = Object@Data[[1]]$Pathways
            }
            if(length(labels) > 1 & length(labels)<= ncol(RR.data))
            {
              labels.2 <- rep("", times =ncol(RR.data))
              labels.2[Object@Data[[1]]$Pathways %in% labels] <- labels
              labels <- labels.2
            }


            # network plot
            ggnet2(RR.net, color = "Cluster", palette = y, alpha = 0.75, size = 4,
                   edge.alpha = 0.5, label = labels)
})
