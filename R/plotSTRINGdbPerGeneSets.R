#' @title plotSTRINGdbPerGeneSets
#'
#' @description
#' Plot the protein-protein interaction predcited using the STRINGdb plugin.
#' Before GetSTRINGdbPerGeneSets must be run.
#'
#'
#' @import STRINGdb
#'
#' @param Object a pathway object
#' @param plot.cluster The cluster for which the string image should be plotted.
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
#' \donttest{
#' GSEA.Object2 <- CombineGeneSets(Object = GSEA.Object1, threads=2)
#' GSEA.Object3 <- ClusterGeneSets(Object = GSEA.Object2,
#'                 clusters = 3,
#'                 method = "kmeans")
#' GSEA.Object3 <- ClusterIndependentGeneSet(Object = GSEA.Object3)
#'
#' GSEA.Object4 <- GetSTRINGdbPerGeneSets(Object = GSEA.Object3, uniquePathways=FALSE,
#'                                        clusterIndependent=FALSE)
#' GSEA.Object4 <- GetSTRINGdbPerGeneSets(Object = GSEA.Object4, uniquePathways=FALSE,
#'                                        clusterIndependent=TRUE)
#' plotSTRINGdbPerGeneSets(Object = GSEA.Object4, plot.cluster = 1)
#' plotSTRINGdbPerGeneSets(Object = GSEA.Object4, plot.cluster = 2)
#' plotSTRINGdbPerGeneSets(Object = GSEA.Object4, plot.cluster = 3)
#' }
#'
#' @export
#'

setGeneric(name="plotSTRINGdbPerGeneSets",
           def=function(Object, plot.cluster=1)
           {
             standardGeneric("plotSTRINGdbPerGeneSets")
           }
)

setMethod(f="plotSTRINGdbPerGeneSets",
          signature = "PathwayObject",
          definition = function(Object, plot.cluster=1)
          {

          if (is.null( Object@stringdbInfo))
          {
            message("There is no information processed of STRINGdb.\nPlease first run GetSTRINGdbPerGeneSets.")
            stop()
          }

          StringObject <- Object@stringdbInfo[[1]]

          if (plot.cluster > length(StringObject))
          {
            message(paste0("There are only ", length(StringObject), " clusters for the parameters that have been selected."))
            stop()
          }


          graphics::plot(1:(dim(StringObject[[plot.cluster]]$img)[2]), type = "n", xaxt = "n", yaxt = "n",
               xlab = "", ylab = "",
               ylim = c(1, dim(StringObject[[plot.cluster]]$img)[1]),
               xlim = c(1, (dim(StringObject[[plot.cluster]]$img)[2])), asp = 1)

          graphics::mtext(StringObject[[plot.cluster]]$input, cex = 0.7)

          graphics::mtext(StringObject[[plot.cluster]]$link, cex = 0.7,side = 1)

          graphics::rasterImage(StringObject[[plot.cluster]]$img, 1, 1,
                               (dim(StringObject[[plot.cluster]]$img)[2]),
                                dim(StringObject[[plot.cluster]]$img)[1])
})
