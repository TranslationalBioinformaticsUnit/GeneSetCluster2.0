#' @title ClusterIndependentGeneSet
#'
#' @description
#' Performs seriation-based clustering of Gene-sets and pathways based on
#' previously calculated distance similarity. The function applies 29 different
#' seriation algorithms and selects the one with the highest score to represent
#' the ordering of the distance correlation matrix. Depending on the  minimum
#' number of Gene-sets estimated to form a cluster, the clusters will be defined.
#' The minimum number of Gene-sets is estimated automatically depending on the number
#' of Gene-sets and the number of clusters generated. But a custom minimum number
#' of Gene-sets can be selected using the method SetPathway.
#'
#' @import seriation
#' @import pbapply
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#'
#' @param Object PathwayObject.
#' @param use_method the method to use for ordering the correlation matrix.
#' @param nPathways a numeric or character. Optimal as default.
#'
#' @return PathwayObject
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
#' GSEA.Object2 <- CombineGeneSets(Object = GSEA.Object1, threads=2)
#' GSEA.Object3 <- ClusterIndependentGeneSet(Object = GSEA.Object2)
#' }
#'
#' @export
#'
setGeneric(name="ClusterIndependentGeneSet",
           def=function(Object, use_method = NULL, nPathways = "optimal")
           {
             standardGeneric("ClusterIndependentGeneSet")
           }
)

setMethod(f="ClusterIndependentGeneSet",
          signature = "PathwayObject",
          definition = function(Object, use_method = NULL, nPathways = "optimal")
          {
            message("Performing the estimation for all the pathways...\n")
            mat_sym <- scaleCorMatrix(Object@Data.RR)
            # select the best ordering method for the correlation matrix
            ordergo <- optimalDist(mat_sym, use_method = use_method)

            message("\nPerforming the estimation for unique pathways...\n")
            mat_symunique <- scaleCorMatrix(Object@DataPathways.RR)
            # select the best ordering method for the correlation matrix
            ordergounique <- optimalDist(mat_symunique, use_method = use_method)

            Object@cIndependentMethod <- list(ordergo, ordergounique)
            message("\n")
            Object <- SetPathway(Object, nPathways = nPathways)


            return(Object)
          }
)
