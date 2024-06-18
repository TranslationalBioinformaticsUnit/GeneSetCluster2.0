
setGeneric(name="ClusterIndependentGeneSet_shiny",
           def=function(Object, iterations = 3, use_method = NULL, all_methods = FALSE, nPathways = "optimal", uniquePathways = FALSE)
           {
             standardGeneric("ClusterIndependentGeneSet_shiny")
           }
)


setMethod(f="ClusterIndependentGeneSet_shiny",
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
