#' ClusterIndependentGeneSet_shiny
#'
#' Choose the most appropiate ordering method for the DataPathways.RR correlation matrix.
#' And perform the functional annotation per cluster.
#' @import simplifyEnrichment
#' @import slam
#' @import GetoptLong
#' @import seriation
#'
#' @param Object a pathway object
#' @param iterations number of iterations to perform per method
#' @param use_method the method to use for ordering the correlation matrix
#' @param nPathways a numeric or character. Optimal as default.
#' @param uniquePathways boolean for unique pathways
#'
#' @return PathwayObject
#'
#' @export
#'
setGeneric(name="ClusterIndependentGeneSet_shiny",
           def=function(Object, iterations = 3, use_method = NULL, all_methods = FALSE, nPathways = "optimal", uniquePathways = FALSE)
           {
             standardGeneric("ClusterIndependentGeneSet_shiny")
           }
)

#' ClusterIndependentGeneSet_shiny
#'
#' @param Object a pathway object
#' @param iterations number of iterations to perform per method
#' @param use_method the method to use for ordering the correlation matrix
#' @param nPathways a numeric or character. Optimal as default.
#'
#' @return PathwayObject
#'
#' @export
#'

setMethod(f="ClusterIndependentGeneSet_shiny",
          signature = "PathwayObject",
          definition = function(Object, iterations = 3, use_method = NULL, all_methods = FALSE, nPathways = "optimal", uniquePathways = FALSE)
          {
            
            set.seed(123)
            if(uniquePathways){
              message("\nPerforming the estimation for UNIQUE pathways...\n")
              mat_symunique <- scaleCorMatrix(Object@DataPathways.RR)
              # select the best ordering method for the correlation matrix
              ordergo <- optimalDist(mat_symunique, all_methods = F, use_method = use_method)
            }else{
              message("Performing the estimation for ALL the pathways...\n")
              mat_sym <- scaleCorMatrix(Object@Data.RR)
              # select the best ordering method for the correlation matrix
              ordergo <- optimalDist(mat_sym, all_methods = F, use_method = use_method)
            }
            
            #Object@cIndependentMethod <- list(ordergo, ordergounique)
            Object@cIndependentMethod <- list(ordergo)

            Object <- SetPathway_shiny(Object, nPathways = nPathways, uniquePathways = uniquePathways)


            return(Object)
          }
)
