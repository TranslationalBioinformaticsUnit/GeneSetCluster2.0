#' @title HighlightGeneSets
#' @description
#' Adds a highlight score if the Gene-Set overlaps with a gene subset which is supplied by the user.
#'
#' @param Object A PathwayObject.
#' @param highligt.genes A vector with genes from the subset the user is interested in. e.g. a list of ROS genes.
#' @param name The name of the subset which will be added to the score calculated.
#' @param clusterIndependent A boolean to specify the clustering approach.
#' TRUE for seriation-based clustering approach and FALSE for classic approach.
#'
#' @return PathwayObject
#'
#' @examples
#' IPA.files <- c(system.file("extdata",
#'                            "MM10.IPA.KO.uGvsMac.Canonical_pathways.xls",
#'                             package = "GeneSetCluster"),
#'              system.file("extdata",
#'                             "MM10.IPA.WT.uGvsMac.Canonical_pathways.xls",
#'                              package = "GeneSetCluster"),
#'              system.file("extdata",
#'                              "MM10.IPA.KO.uGvsMac.Functional_annotations.xls",
#'                              package = "GeneSetCluster"),
#'              system.file("extdata",
#'                              "MM10.IPA.WT.uGvsMac.Functional_annotations.xls",
#'                              package = "GeneSetCluster"))
#' canonical.files <- IPA.files[grep("Canonical", IPA.files)]
#'
#' IPA.object1 <- LoadGeneSets(file_location = canonical.files,
#'                          groupnames= c("KO", "WT"),
#'                          P.cutoff = 1.3,
#'                          Mol.cutoff = 5,
#'                          Source = "IPA",
#'                          type = "Canonical_Pathways",
#'                          structure = "SYMBOL",
#'                          seperator = ",",
#'                          Organism = "org.Mm.eg.db")
#' \dontrun{
#' IPA.object2 <- CombineGeneSets(Object = IPA.object1)
#'
#' IPA.object3 <- ClusterGeneSets(Object = IPA.object2,
#'                               clusters = 6,
#'                               method = "kmeans")
#' system.file("data", "Redox.genes", package = "GeneSetCluster")
#' IPA.object3.highlight <- HighlightGeneSets(Object = IPA.object3,
#'                                           highligt.genes = Redox.genes,
#'                                           name = "Ros", clusterIndependent = FALSE)
#' }
#' @export
#'
setGeneric(name="HighlightGeneSets",
           def=function(Object, highligt.genes, name = "",
                        clusterIndependent = FALSE)
           {
             standardGeneric("HighlightGeneSets")
           }
)

setMethod(f="HighlightGeneSets",
          definition=function(Object, highligt.genes, name = "",
                              clusterIndependent = FALSE)
{
  message("[=========================================================]")
  message("[<<<<<             HighlightGeneSets start          >>>>>>]")
  message("[=========================================================]")

  if(is.na(unique(Object@metadata[,"cluster.method"])))
  {
    message("Make sure youre object has been clustered by ClusterGeneSets")
    stop()
  }

  message("make sure that the structure of the highlight.genes is the same as the data")
  #Make sure that the highlight genes have the same structure as the object
  if(unique(as.character(Object@metadata[,"structure"])) == "SYMBOL")
  {
    message("transforming all highligt.genes to upper case, make sure this doesnt change the data")
    message(paste( "raw data has " , length(unique(highligt.genes))," highligt.genes", sep=""))
    highligt.genes <- toupper(highligt.genes)
    message(paste( "Transformed data has " , length(unique(highligt.genes))," highligt.genes", sep=""))
  }


  ################################################
  #-----------Seperate per cluster---------------#
  ################################################

  message("Calculating overlap between pathways and highlight genes")


  if (clusterIndependent == FALSE)
  {
    checkDependent(Object)
    data <- vector()
    for(cl.i in unique(Object@Data[[1]]$cluster))
    {
      DF.cl <- Object@Data[[1]][Object@Data[[1]]$cluster == cl.i,]

      highligt.score <- vector()
      for(genes.cl.i in 1:nrow(DF.cl))
      {
        genes.x <- as.vector(strsplit2(as.character(DF.cl[genes.cl.i,]$Molecules), split=Object@metadata[1,"seperator"]))
        y <- sum(highligt.genes %in% genes.x)
        highligt.score[genes.cl.i] <- y/length(genes.x)
      }
      DF.cl$Highlight <- highligt.score
      DF.cl$Highlight.mean <- mean(highligt.score)
      data <- rbind(data, DF.cl)
    }
    message("Highlight mean values per cluster:")
    for (i in 1:max(data$cluster))
    {
      message(paste0("\tCluster ", i, ": ", data[which(data$cluster==i), "Highlight.mean"][1]))

    }
  } else {
    checkIndependent(Object)
    data <- vector()
    tmp.df <- Object@Data[[1]]
    rownames(tmp.df) <- tmp.df$RR_name
    tmp.df <- tmp.df[rownames(Object@plot$aka2Independent),]
    tmp.df$cluster_seriation <- as.numeric(Object@plot$aka2Independent$Cluster)

    for(cl.i in 1:max(tmp.df$cluster_seriation))
    {
      DF.cl <- tmp.df[tmp.df$cluster_seriation == cl.i,]

      highligt.score <- vector()
      for(genes.cl.i in 1:nrow(DF.cl))
      {
        genes.x <- as.vector(strsplit2(as.character(DF.cl[genes.cl.i,]$Molecules), split=Object@metadata[1,"seperator"]))
        y <- sum(highligt.genes %in% genes.x)
        highligt.score[genes.cl.i] <- y/length(genes.x)
      }
      DF.cl$Highlight_seriationBased <- highligt.score
      DF.cl$Highlight.mean_seriationBased <- mean(highligt.score)
      data <- rbind(data, DF.cl)

    }

    message("Highlight mean values per cluster:")
    for (i in 1:max(data$cluster_seriation))
    {
      message(paste0("\tCluster ", i, ": ", data[which(data$cluster_seriation==i), "Highlight.mean_seriationBased"][1]))
    }
}

  Object@Data <- list(data)

  message("-----------------------------------------------------------")
  message("[<<<<<             HighlightGeneSets END            >>>>>>]")
  message("[=========================================================]")
  return(Object)
})
