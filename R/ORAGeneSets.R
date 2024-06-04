#' @title ORAGeneSets
#'
#' @description
#' Obtains the results of the ORAs performed by the genes that make up each cluster.
#'
#' @param Object A PathwayObject.
#' @param uniquePathway A boolean to specify the approach of unique Gene-set or unique Pathways.
#' FALSE if for unique Gene-set and TRUE for unique Pathways.
#' @param clusterIndependent A boolean to specify the clustering approach.
#' TRUE for seriation-based clustering approach and FALSE for classic approach.
#'
#' @return list with the over represented gene sets for all the clusters.
#'
#' @examples
#'
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
#' GSEA.Object2 <- CombineGeneSets(Object=GSEA.Object1, threads=2)
#' GSEA.Object3 <- ClusterIndependentGeneSet(Object=GSEA.Object2)
#' }
#' @export
#'
setGeneric(name="ORAGeneSets",
           def=function(Object, uniquePathway=FALSE, clusterIndependent=FALSE)
           {
             standardGeneric("ORAGeneSets")
           }
)


setMethod(f="ORAGeneSets",
          signature = "PathwayObject",
          definition = function(Object, uniquePathway=FALSE, clusterIndependent=FALSE)
          {

          message("[=========================================================]")
          message("[<<<<               ORAGeneSets START                >>>>>]")
          message("-----------------------------------------------------------")

          if (clusterIndependent)
          {
            checkIndependent(Object)
            if (uniquePathway)
            {
              ora_list <- Object@functionalAnnotIndependent$Pathway$ORA
              names(ora_list) <- paste0("Cluster_", 1:length(Object@functionalAnnotIndependent$Pathway$ORA))

            } else {
              ora_list <- Object@functionalAnnotIndependent$Geneset$ORA
              names(ora_list) <- paste0("Cluster_", 1:length(Object@functionalAnnotIndependent$Geneset$ORA))

            }
          } else {
            checkDependent(Object)
            if (uniquePathway)
            {
              ora_list <- Object@functionalAnnot$Pathway$ORA
              names(ora_list) <- paste0("Cluster_", 1:length(Object@functionalAnnot$Pathway$ORA))

            } else {
              ora_list <- Object@functionalAnnot$Geneset$ORA
              names(ora_list) <- paste0("Cluster_", 1:length(Object@functionalAnnot$Geneset$ORA))

            }
          }

          message("-----------------------------------------------------------")
          message("[<<<<<             ORAGeneSets END              >>>>>>]")
          message("[=========================================================]")
          message("[You may want to process HighlightGeneSets next.            ]")

          return(ora_list)
})
