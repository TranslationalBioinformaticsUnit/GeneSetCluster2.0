#' GenesPerGeneSet
#'
#' @importFrom limma strsplit2
#'
#' Extracts a data.frame from the object with every gene in per cluster
#' @param Object a PathwayObject
#'
#' @return dataframe of pathwayobjects
#' @export
#'

setGeneric(name="GenesPerGeneSet",
           def=function(Object, uniquePathways=FALSE, clusterIndependent=FALSE)
           {
             standardGeneric("GenesPerGeneSet")
           }
)

#' GenesPerGeneSet
#'
#' @param Object a PathwayObject
#' @param PathwayObject  a PathwayObject
#'
#' @return dataframe of pathwayobjects
#' @export
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
#'                          seperator = ",")
#'IPA.object2 <- CombineGeneSets(Object = IPA.object1)
#'IPA.object3 <- ClusterGeneSets(Object = IPA.object2,
#'                               clusters = 7,
#'                               method = "kmeans")
#' GenesPerGeneSet(Object =IPA.object3 )
setMethod(f="GenesPerGeneSet",
          signature="PathwayObject",
          definition=function(Object, uniquePathways=FALSE, clusterIndependent=FALSE)
          {
            message("[=============================================================]")
            message("[<<<<               GenesPerGeneSet START                >>>>>]")
            message("---------------------------------------------------------------")

            clus.mol.ls <- list()

            if (uniquePathways)
            {
              if (clusterIndependent)
              {
                message("Obtaining genes from Cluster-independent approach...")
                forIteration <- 1:max(Object@plot$aka2IndependentUnique$Cluster)

                message("Unique Pathways selected.")
                message(paste("Extracting genes for all ", length(forIteration), " clusters", sep=""))

                pathwaysUse <- rownames(Object@plot$aka2IndependentUnique[Object@plot$aka2IndependentUnique$Cluster>0,])
                use_cluster <- Object@plot$aka2IndependentUnique$Cluster[Object@plot$aka2IndependentUnique$Cluster>0]

              } else {
                forIteration <- unique(Object@plot$aka2Unique$Cluster)

                message("Unique Pathways selected.")
                message(paste("Extracting genes for all ", length(forIteration), " clusters", sep=""))
                pathwaysUse <- rownames(Object@plot$aka2Unique)
                use_cluster <- Object@plot$aka2Unique$Cluster
              }

            } else {
              if (clusterIndependent)
              {
                message("Obtaining genes from Cluster-independent approach...")
                forIteration <- 1:max(Object@plot$aka2Independent$Cluster)

                message("Unique Pathways selected.")
                message(paste("Extracting genes for all ", length(forIteration), " clusters", sep=""))
                dfTmp <- Object@Data[[1]]
                rownames(dfTmp) <- Object@Data[[1]]$RR_name
                dfTmp <- dfTmp[rownames(Object@plot$aka2Independent),]
                pathwaysUse <- Object@plot$aka2Independent
                pathwaysUse$Pathways <- dfTmp$RR_name
                pathwaysUse <- pathwaysUse$Pathways[pathwaysUse$Cluster>0]
                use_cluster <- Object@plot$aka2Independent$Cluster[Object@plot$aka2Independent$Cluster>0]

              } else {
                #forIteration <- unique(Object@Data[[1]]$cluster)
                forIteration <- gsub("Cluster_","",unique(Object@Data[[1]]$cluster))

                message("Unique Gene-sets selected.")
                message(paste("Extracting genes for all ", length(forIteration), " clusters", sep=""))
                dfTmp <- Object@Data[[1]]
                rownames(dfTmp) <- Object@Data[[1]]$RR_name
                dfTmp <- dfTmp[rownames(Object@plot$aka2),]
                pathwaysUse <- dfTmp$RR_name
                use_cluster <- Object@plot$aka2$Cluster
              }
            }


            for(clus.i in unique(use_cluster))
            {

              pathways.x <- pathwaysUse[use_cluster == clus.i]
              # clus.x <- Object@Data[[1]][use_cluster == unique(use_cluster)[clus.i],]
              if (uniquePathways)
              {
                clus.x <- Object@Data[[1]][Object@Data[[1]]$Pathways %in% pathways.x,]
                clus.x$Molecules <- as.character(clus.x$unionMolecules)
              } else {
                clus.x <- Object@Data[[1]][Object@Data[[1]]$RR_name %in% pathways.x,]
                clus.x$Molecules <- as.character(clus.x$Molecules)
              }
              clus.mol.ls[[clus.i]] <- unique(as.vector(strsplit2(clus.x$Molecules, split = Object@metadata$seperator[1])))
            }
            # names(clus.mol.ls) <- unique(use_cluster)

            unique.mol <- unique(do.call(what = c, args = clus.mol.ls))

            mol.unique.df <- as.data.frame(matrix(0, nrow = length(unique.mol), ncol = length(clus.mol.ls)))
            rownames(mol.unique.df) <- unique.mol
            colnames(mol.unique.df) <-  paste(unique(use_cluster), sep="")
            
            for(clus.i in forIteration)
            {
              mol.unique.df[clus.mol.ls[[clus.i]],clus.i] <- 1
            }
            
            
            mol.unique.df <- mol.unique.df[, order(as.numeric(colnames(mol.unique.df)), decreasing = F)]
            colnames(mol.unique.df) <-  paste("Cluster_", colnames(mol.unique.df), sep="")

            message("--------------------------------------------------------------")
            message("[<<<<               GenesPerGeneSet ends                >>>>>]")
            message("[============================================================]")

            return(mol.unique.df)
}
)
