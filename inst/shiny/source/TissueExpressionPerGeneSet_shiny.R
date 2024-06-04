#' TissueExpressionPerGeneSet
#'
#' Extracts every gen of every cluster from the object. Using the median gene expression from the GTEx database build a data.frame with every tissue per cluster.
#'
#' @import jsonlite
#' @import httr
#' @import reshape2
#' @import dplyr
#' @import utils
#' @import pbapply
#'
#' @param Object a PathwayObject
#' @param localDatabase data frame with the GTEx database information
#' @param uniquePathways  a boolean
#'
#' @return dataframe of pathwayobjects with tissue expression per cluster
#'

setGeneric(name="TissueExpressionPerGeneSet",
           def=function(Object, tissues="ALL", localDatabase=NULL, dic=NULL, clusterIndependent=FALSE, uniquePathways=F)
           {
             standardGeneric("TissueExpressionPerGeneSet")
           }
)

#' TissueExpressionPerGeneSet
#'
#' @param Object a PathwayObject
#' @param localDatabase data frame with the GTEx database information
#' @param uniquePathways  a boolean
#'
#' @return dataframe of pathwayobjects with tissue expression per cluster
#' @export
#'
#' @examples
#' #' IPA.files <- c(system.file("extdata",
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
#' IPA.object2 <- CombineGeneSets(Object = IPA.object1)
#' IPA.object3 <- ClusterGeneSets(Object = IPA.object2,
#'                               clusters = 7,
#'                               method = "kmeans")
#' IPA.object4 <- TissueExpressionPerGeneSet(Object = IPA.object3, threads = 8)
setMethod(f="TissueExpressionPerGeneSet",
          signature="PathwayObject",
          definition=function(Object, tissues="ALL", localDatabase=NULL, dic=NULL, clusterIndependent=FALSE, uniquePathways=F)
          {
          message("[=========================================================]")
          message("[<<<<         ObtainTissueExpression START           >>>>>]")
          message("-----------------------------------------------------------")

          options <- c("SYMBOL", "ENTREZID", "ENSEMBLID", "ENTREZ", "ENSEMBL")


          if (!toupper(Object@metadata$structure[1]) %in% options) {
            message("Genes are not included as genes Symbol, EntrezID or ENSEMBLID.\nSo TissuePerGeneSet will NOT be performed.")
            message("Please make sure that you are using one of the followings: GENE SYMBOL, ENTREZID or ENSEMBLID.")
            stop()
          }

          # data("dic", package = "GeneSetCluster")
          # data("availableTissues", package = "GeneSetCluster")
          # 
          # if (checkTissues(tissues, tissue_names))
          # {
          #   if (tissues[1] == "ALL")
          #   {
          #     tissues <- tissue_names
          #     message("\nSelected ALL tissues.\n")
          # 
          #   } else {
          #     message("Selected tissues:\n")
          #     show_tissues <- lapply(tissues, function(x) paste0(x, "\n"))
          #     message(show_tissues)
          # 
          #   }
          # } else {
          #   tissue_names <- lapply(tissue_names, function(x) paste0(x, "\n"))
          #   warning("Introduced tissue names are incorrect.\nThese are available tissue names")
          #   message(tissue_names)
          #   stop("\nPlease make sure you introduced tissue(s) name(s) correctly.")
          # }

          if (clusterIndependent)
          {
            checkIndependent(Object)
            mol.unique.df <- GenesPerGeneSet(Object, uniquePathways = uniquePathways, clusterIndependent = TRUE)
            #mol.unique.dfPathway <- GenesPerGeneSet(Object, uniquePathways = TRUE, clusterIndependent = TRUE)
          } else {
            mol.unique.df <- GenesPerGeneSet(Object, uniquePathways = uniquePathways, clusterIndependent = FALSE)
            #mol.unique.dfPathway <- GenesPerGeneSet(Object, uniquePathways = TRUE, clusterIndependent = FALSE)
          }


          genestissue <- getGenesTissue(Object, mol.unique.df, dic)
          meltedtouse <- genestissue[[1]]
          touse <- genestissue[[2]]

          # genestissuePathway <- getGenesTissue(Object, mol.unique.dfPathway, dic)
          # meltedtousePathway <- genestissuePathway[[1]]
          # tousePathway <- genestissuePathway[[2]]

          
          tryCatch({
            #localDatabase_use <- localDatabase[tissues]
            message("Performing GSEA for unique Gene-sets...\n")
            results <- pblapply(localDatabase, function(x) clusterProfiler::GSEA(geneList=x, TERM2GENE=meltedtouse, minGSSize=5, maxGSSize=10000, verbose=F))
            #message("Performing GSEA for unique Pathways...\n")
            #resultsPathway <- pblapply(localDatabase_use, function(x) clusterProfiler::GSEA(geneList=x, TERM2GENE=meltedtousePathway, minGSSize=5, maxGSSize=10000, verbose=F))

          }, error=function(e){
            stop("Seems that the database provided is not correct. Please make sure you are using the correct one.\nYou can download it on https://github.com/TranslationalBioinformaticsUnit/GeneSetCluster")
          })

          tissue.df <- data.frame(matrix(0, nrow = length(results), ncol = length(unique(meltedtouse$Term))))
          #tissue.dfPathway <- data.frame(matrix(0, nrow = length(resultsPathway), ncol = length(unique(meltedtousePathway$Term))))
          rownames(tissue.df) <- names(results)
          colnames(tissue.df) <- colnames(touse)[1:ncol(tissue.df)]
          #rownames(tissue.dfPathway) <- names(resultsPathway)
          #colnames(tissue.dfPathway) <- colnames(tousePathway)[1:ncol(tissue.dfPathway)]

          for (i in 1:length(results))
          {
            tissue.df[i, results[[i]]@result$ID] <- results[[i]]@result$NES
          }

          # for (i in 1:length(resultsPathway))
          # {
          #   tissue.dfPathway[i, resultsPathway[[i]]@result$ID] <- resultsPathway[[i]]@result$NES
          # }


          if (clusterIndependent)
          {
            #Object@dfTissueIndependent <- list(Geneset=tissue.df, Pathway=tissue.dfPathway)
            Object@dfTissueIndependent <- tissue.df
          } else {
            #Object@dfTissue <- list(Geneset=tissue.df, Pathway=tissue.dfPathway)
            Object@dfTissue <- tissue.df
          }

          message("\n")
          message("[=========================================================]")
          message("[<<<<         ObtainTissueExpression END             >>>>>]")
          message("-----------------------------------------------------------")
          message("[You may want to plot the results using PlotTissueExpression next.]")

          return(Object)
        }
)
