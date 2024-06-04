#' @title SetPathway
#'
#' @description
#' Set the minimum number of Gene-sets to be considered as group
#' in seriation-method clustering. Then the clusters will be definded and will
#' annotated.
#'
#' @import slam
#' @import GetoptLong
#' @import simplifyEnrichment
#' @import grid
#' @import RColorBrewer
#' @importFrom stringr str_split
#'
#' @param Object a PathwayObject.
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
#' \donttest{
#' GSEA.Object2 <- CombineGeneSets(Object = GSEA.Object1, threads=2)
#' GSEA.Object3 <- ClusterIndependentGeneSet(Object = GSEA.Object2)
#' GSEA.Object3.1 <- SetPathway(GSEA.Object3, nPathways=12)
#' }
#'
#' @export
#'

setGeneric(name="SetPathway",
           def=function(Object, nPathways="optimal")
           {
             standardGeneric("SetPathway")
           }
)

setMethod(f="SetPathway",
          signature="PathwayObject",
          definition=function(Object, nPathways="optimal")
          {
            message("[=============================================================]")
            message("[<<<<               SetPathway START                >>>>>]")
            message("---------------------------------------------------------------")

            #check if cluster independent slot is null
            checkIndependent(Object)

            mat_cor <- scaleCorMatrix(Object@Data.RR)[Object@cIndependentMethod[[1]][[1]],
                                                      Object@cIndependentMethod[[1]][[1]]]

            mat_corPathway <- scaleCorMatrix(Object@DataPathways.RR)[Object@cIndependentMethod[[2]][[1]],
                                                                     Object@cIndependentMethod[[2]][[1]]]

            res <- obtainDefCluster(mat_cor)
            resPathway <- obtainDefCluster(mat_corPathway)

            if (nPathways == "optimal")
            {
            #scaling the matrix
            optimalNumber <- getOptimalNumber(res)
            optimalNumberPathway <- getOptimalNumber(resPathway)

            message(paste0("The optimal number of minimum Gene-sets for Unique Gene-set approach is ", optimalNumber))
            message(paste0("The optimal number of minimum Pathways for Unique Pathway approach is ", optimalNumberPathway))
            } else {
              optimalNumber <- nPathways
              optimalNumberPathway <- nPathways
              message(paste0("Setting ", nPathways, " as minimum number of elements to consider a cluster...\n"))
            }

            go_id_list <- res[which(lapply(res, function(x) length(x)) >= optimalNumber)]
            go_id_listPathway <- resPathway[which(lapply(resPathway, function(x) length(x)) >= optimalNumberPathway)]

            extractedData <- Object@Data[[1]][which(Object@Data[[1]]$RR_name %in% unlist(go_id_list)), ]


            aka2 <- matrix(data = NA, nrow = nrow(Object@Data[[1]]), ncol = 2)
            colnames(aka2) <- c("Group", "Cluster")
            rownames(aka2) <- Object@Data[[1]]$RR_name

            aka2[,"Group"] <- as.character(Object@Data[[1]]$Groups)
            #order matrix
            aka2 <- aka2[rownames(mat_cor),]

            for (i in 1:length(go_id_list)) {
              aka2[which(rownames(aka2) %in% go_id_list[[i]]), "Cluster"] <- i
            }

            aka2 <- as.data.frame(aka2)
            aka2$Cluster[is.na(aka2$Cluster)] <- 0

            if(length(unique(aka2$Cluster)) > 32)
            {
              message("Warning: number of cluster is larger than colours supported")
            }

            pal.c <-  getPal()

            groups.col <- brewer.pal(n = 8, name ="Set2" )[1:length(unique(aka2[,"Group"]))]
            names(groups.col) <- unique(aka2[,"Group"])
            pal <- pal.c[1:length(unique(aka2[,"Cluster"]))]
            names(pal) <- unique(aka2[,"Cluster"])
            pal["0"] <- "white"

            aka3 <- list(Group = groups.col,
                        Cluster= pal,
                        optimalNumberPathway = optimalNumber)

            Object@plot$aka2Independent <- aka2
            Object@plot$aka3Independent <- aka3


            aka2Unique <- obtainDFmetadata(Object@Data[[1]], Object@DataPathways.RR)
            #order matrix
            aka2Unique <- aka2Unique[rownames(mat_corPathway),]
            aka2Unique$Cluster <- 0
            for (i in 1:length(go_id_listPathway))
            {
              aka2Unique[go_id_listPathway[[i]], "Cluster"] <- i
            }

            clusterOrder <- sort(unique(aka2Unique[,"Cluster"]))
            clusterOrder <- c(1:(length(clusterOrder)-1), 0)
            palUnique <- pal.c[1:length(clusterOrder)]
            names(palUnique) <- clusterOrder
            palUnique["0"] <- "white"


            aka3Unique <- list(Group = groups.col,
                               Cluster = palUnique,
                               optimalNumberPathway = optimalNumberPathway)

            Object@plot$aka2IndependentUnique <- aka2Unique
            Object@plot$aka3IndependentUnique <- aka3Unique


            message("\nCalculating ORA and Semantic enrichment per cluster for Unique Gene-sets...")
            functional_annotation <- PerformORAGOperCluster(Object, uniquePathways = FALSE, clusterIndependent = TRUE, minNumber = nPathways)
            message("\nCalculating ORA and Semantic enrichment per cluster for Unique Pathways...")
            functional_annotationUnique <- PerformORAGOperCluster(Object, uniquePathways = TRUE, clusterIndependent = TRUE, minNumber = nPathways)

            Object@functionalAnnotIndependent <- list(Geneset=functional_annotation,
                                                      Pathway=functional_annotationUnique)


            message("--------------------------------------------------------------")
            message("[<<<<               SetPathway ends                >>>>>]")
            message("[============================================================]")

            return(Object)
}
)
