
setGeneric(name="ShowPathwayCluster_shiny",
           def=function(Object, cluster = "ALL", nPathway = 4, uniquePathways = FALSE)
           {
             standardGeneric("ShowPathwayCluster_shiny")
           }
)



setMethod(f="ShowPathwayCluster_shiny",
          definition = function(Object, cluster="ALL", nPathway = 4, uniquePathways = FALSE)
{

  message("[=============================================================]")
  message("[<<<<               ShowPathwayCluster START                >>>>>]")
  message("---------------------------------------------------------------")
  if (is.null(Object@cIndependentMethod))
  {
    message("The cluster independent methos is not defines. Please firstly run ClusterIndependentGeneSet in order to select the ordering method.")
    stop()
  }

  message("Make sure you are using the same nPathway value as in PlotPathwayCluster.")

  tryCatch(
    {
      if (toupper(cluster) != "ALL") {
        as.numeric(cluster)
      }
    },
    error=function(e) {
      message("The argument included for cluster parameter is not valid. Instead showing all the clusters.")
      cluster <- "ALL"
    }
  )

  #scaling the matrix
  if(uniquePathways){
    mat <- Object@DataPathways.RR
    nPathway <-Object@plot$aka3IndependentUnique$optimalNumberPathway
  }else{
    mat <- Object@Data.RR
    nPathway <-Object@plot$aka3Independent$optimalNumberPathway
  }
  mat[lower.tri(mat, diag = F)] <- NA
  mat_mod <- apply(mat, 1, function(x) scale_values(x, na.rm=T))
  mat_sym <- makeSymm(mat_mod)
  mat_sym[ncol(mat_sym), ncol(mat_sym)] <- 1

  if (anyNA(mat_sym)) {
    mat_sym[which(is.na(mat_sym))] <- 1
  }

  mat_cor <- mat_sym[Object@cIndependentMethod[[1]][[1]], Object@cIndependentMethod[[1]][[1]]]


  res <- obtainDefCluster(mat_cor)

  if (length(res[which(lapply(res, function(x) length(x))>=nPathway)]) == 0) {
    message(paste0("There is no groups with at least ", nPathway, " pathways. Using the minimum value, 2 pathways per group."))
    go_id_list <- res[which(lapply(res, function(x) length(x))>=2)]
  } else {
    go_id_list <- res[which(lapply(res, function(x) length(x))>=nPathway)]
  }

  if (((cluster > length(go_id_list) | cluster <= 0)) & toupper(cluster) != "ALL")  {
    paste0(message("There is no cluster ", cluster, ". Instead showing all the clusters."))
    cluster = "ALL"
  }

  #if (checkGO(Object)) {
  #  GOres <- lapply(go_id_list, obtainGOdescription)
  #} else {
  #  GOres <- go_id_list
  #}

  GOres <- go_id_list

  if (toupper(cluster) == "ALL") {
    finalRes <- GOres
  } else {
    finalRes <- GOres[[cluster]]
  }

  message("[=============================================================]")
  message("[<<<<               ShowPathwayCluster END                >>>>>]")
  message("---------------------------------------------------------------")

  return(finalRes)
}
)

