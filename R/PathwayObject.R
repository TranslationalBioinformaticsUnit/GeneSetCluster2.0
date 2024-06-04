setClassUnion("DForNULL", c("data.frame", "NULL", "matrix"))
setClassUnion("LISTorNULL", c("list", "NULL"))


#' @title PathwayObject
#' @aliases PathwayObject-class
#' @description Stores the all the information for GeneSetCluster. Essentially
#' pathways IDs, list of the genes, distance matrix and functional annotation
#' of the clusters.
#' @slot Data list. All the introduced data
#' @slot PData data.frame.
#' @slot metadata data.frame. Metdata of the object.
#' @slot plot list. Information for plotting.
#' @slot Data.RR data.frame. Distance matrix for unique Gene-sets.
#' @slot DataPathways.RR data.frame. Distance matrix for unique Pathways.
#' @slot cIndependentMethod list. Optimal seriation method for the data.
#' @slot dfTissue list. Information about tissue enrichment in traditional clustering.
#' @slot dfTissueIndependent list. Information about tissue enrichment in seriation-based clustering.
#' @slot functionalAnnot list. Functional annotation about the groups generated in traditional clustering.
#' @slot functionalAnnotIndependent list. Functional annotation about the groups generated in seriation-based clustering.
#' @slot stringdbInfo list. Information of the STRING database and the parameters used.
#'
#' @return PathwayObject
#' @exportClass PathwayObject
#' @export

PathwayObject <- setClass(
  Class="PathwayObject",
  slots=list(
    Data="list",
    PData="data.frame",
    metadata="data.frame",
    plot="list",
    Data.RR="data.frame",
    DataPathways.RR="DForNULL",
    cIndependentMethod="LISTorNULL",
    dfTissue="LISTorNULL",
    dfTissueIndependent="LISTorNULL",
    functionalAnnot="LISTorNULL",
    functionalAnnotIndependent="LISTorNULL",
    stringdbInfo="LISTorNULL"),
  prototype=list(
    Data=NULL,
    PData=NULL,
    metadata = NULL,
    plot=NULL,
    Data.RR=NULL,
    DataPathways.RR=NULL,
    cIndependentMethod=NULL,
    dfTissue=NULL,
    dfTissueIndependent=NULL,
    functionalAnnot=NULL,
    functionalAnnotIndependent=NULL,
    stringdbInfo=NULL)
)


#' Show method for a \code{\link{PathwayObject-class}} object
#' @title show method for PathwayObject
#' @param object \code{\link{PathwayObject-class}} object
#' @return summarize information about the object
setMethod("show", "PathwayObject", function(object) {

  cat(
    "The PathwayObject contians:",
    "\nTotal number of pathways:", length(object@Data[[1]]$Pathways),
    "\nGroups:", paste0(unique(object@Data[[1]]$Groups), sep=" "),
    "\nDistances calculated:", ifelse(nrow(object@Data.RR)==0, "NOT processed", "YES"),
    "\nTraditional clustering performed:", ifelse(nrow(object@plot$aka2)==0, "NOT processed", paste0(max(object@plot$aka2$Cluster), " clusters in unique Gene-Sets and ", max(object@plot$aka2Unique$Cluster),  " clusters in unique Pathways")),
    "\nSeriation-based clustering performed:", ifelse(is.null(object@cIndependentMethod), "NOT processed", paste0(max(object@plot$aka2Independent$Cluster), " clusters in unique Gene-Sets and ", max(object@plot$aka2IndependentUnique$Cluster),  " clusters in unique Pathways")),
    "\nTissue enrichment for traditional clustering performed:", ifelse(nrow(object@dfTissue)==0, "NOT processed", "YES"),
    "\nTissue enrichment for seriation-based clustering performed:", ifelse(is.null(object@dfTissueIndependent), "NOT processed", "YES"),
    "\nSTRINGdb performed:", ifelse(is.null(object@stringdbInfo), "NOT processed", "YES"),
    "\n"
  )
}
)
