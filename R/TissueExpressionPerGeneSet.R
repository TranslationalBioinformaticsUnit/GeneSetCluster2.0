#' @title TissueExpressionPerGeneSet
#'
#' @description
#' Extracts every gen of every cluster from the object. Using the median gene
#' expression from the GTEx database build a data.frame with every tissue per cluster.
#' The list of available tissues:
#' Adipose...Subcutaneous; Adipose...Visceral..Omentum.; Adrenal.Gland; Artery...Aorta;
#' Artery...Coronary; Artery...Tibial; Bladder; Brain...Amygdala; Brain...Anterior.cingulate.cortex..BA24.;
#' Brain...Caudate..basal.ganglia.; Brain...Cerebellar.Hemisphere; Brain...Cerebellum; Brain...Cortex;
#' Brain...Frontal.Cortex..BA9.; Brain...Hippocampus; Brain...Hypothalamus; Brain...Nucleus.accumbens..basal.ganglia.;
#' Brain...Putamen..basal.ganglia.; Brain...Spinal.cord..cervical.c.1.; Brain...Substantia.nigra;
#' Breast...Mammary.Tissue; Cells...Cultured.fibroblasts; Cells...EBV.transformed.lymphocytes;
#' Cervix...Ectocervix; Cervix...Endocervix; Colon...Sigmoid; Colon...Transverse;
#' Esophagus...Gastroesophageal.Junction; Esophagus...Mucosa; Esophagus...Muscularis;
#' Fallopian.Tube; Heart...Atrial.Appendage; Heart...Left.Ventricle; Kidney...Cortex;
#' Kidney...Medulla; Liver; Lung; Minor.Salivary.Gland; Muscle...Skeletal; Nerve...Tibial;
#' Ovary; Pancreas; Pituitary; Prostate; Skin...Not.Sun.Exposed..Suprapubic.;
#' Skin...Sun.Exposed..Lower.leg.; Small.Intestine...Terminal.Ileum; Spleen; Stomach;
#' Testis; Thyroid; Uterus; Vagina; Whole.Blood.
#'
#' @import jsonlite
#' @import httr
#' @import reshape2
#' @import utils
#' @import pbapply
#' @importFrom clusterProfiler GSEA
#'
#' @param Object a PathwayObject.
#' @param tissues a character vector with the name of the tissues to use.
#' @param localDatabase the local GTEx database Object. If is NULL then the database will be obtained through custom API
#' Can be downloaded: https://github.com/TranslationalBioinformaticsUnit/GeneSetCluster
#' @param clusterIndependent A boolean to specify the clustering approach.
#' TRUE for seriation-based clustering approach and FALSE for classic approach.
#'
#' @return dataframe of pathwayobjects with tissue expression per cluster
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
#' \dontrun{
#' GSEA.Object2 <- CombineGeneSets(Object = GSEA.Object1, threads=2)
#' GSEA.Object3 <- ClusterGeneSets(Object = GSEA.Object2,
#'                 clusters = 3,
#'                 method = "kmeans")
#' GSEA.Object3 <- ClusterIndependentGeneSet(Object = GSEA.Object3)
#' GSEA.Object4 <- TissueExpressionPerGeneSet(GSEA.Object3, tissues="ALL", localDatabase=NULL, clusterIndependent=FALSE)
#' GSEA.Object4 <- TissueExpressionPerGeneSet(GSEA.Object4, tissues="ALL", localDatabase=NULL, clusterIndependent=TRUE)
#' }
#' @export
#'

setGeneric(name="TissueExpressionPerGeneSet",
           def=function(Object, tissues="ALL", localDatabase=NULL, clusterIndependent=FALSE)
           {
             standardGeneric("TissueExpressionPerGeneSet")
           }
)


setMethod(f="TissueExpressionPerGeneSet",
          signature="PathwayObject",
          definition=function(Object, tissues="ALL", localDatabase=NULL, clusterIndependent=FALSE)
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

          data("dic", package = "GeneSetCluster")
          data("availableTissues", package = "GeneSetCluster")

          if (checkTissues(tissues, tissue_names))
          {
            if (tissues[1] == "ALL")
            {
              tissues <- tissue_names
              message("\nSelected ALL tissues.\n")

            } else {
              message("Selected tissues:\n")
              show_tissues <- lapply(tissues, function(x) paste0(x, "\n"))
              message(show_tissues)

            }
          } else {
            tissue_names <- lapply(tissue_names, function(x) paste0(x, "\n"))
            message("Introduced tissue names are incorrect.\nThese are available tissue names")
            message(tissue_names)
            stop("\nPlease make sure you introduced tissue(s) name(s) correctly.")
          }

          if (clusterIndependent)
          {
            checkIndependent(Object)
            mol.unique.df <- GenesPerGeneSet(Object, uniquePathways = FALSE, clusterIndependent = TRUE)
            mol.unique.dfPathway <- GenesPerGeneSet(Object, uniquePathways = TRUE, clusterIndependent = TRUE)
          } else {
            mol.unique.df <- GenesPerGeneSet(Object, uniquePathways = FALSE, clusterIndependent = FALSE)
            mol.unique.dfPathway <- GenesPerGeneSet(Object, uniquePathways = TRUE, clusterIndependent = FALSE)
          }


          genestissue <- getGenesTissue(Object, mol.unique.df, dic)
          meltedtouse <- genestissue[[1]]
          touse <- genestissue[[2]]

          genestissuePathway <- getGenesTissue(Object, mol.unique.dfPathway, dic)
          meltedtousePathway <- genestissuePathway[[1]]
          tousePathway <- genestissuePathway[[2]]

          if (is.null(localDatabase)) {
            #GTEx API REST query of median Gene Expression: database -> gtex_v8
            message("Performing GTEx API REST query to obtain the GTEx database.")

            url <- "translationalbioinformatic.pythonanywhere.com/GTExdatabase"
            urlNames <- "translationalbioinformatic.pythonanywhere.com/GTExdatabaseNames"
            tryCatch(
              {
                response <- httr::GET(url)
                GTEx.info <- fromJSON(rawToChar(response$content))

                responseNames <- httr::GET(urlNames)
                GTEx.infoNames <- fromJSON(rawToChar(responseNames$content))

                for (i in 1:length(GTEx.info))
                {
                  names(GTEx.info[[i]]) <- GTEx.infoNames[[i]]
                }

                GTEx.info_use <- GTEx.info[tissues]

                message("Performing GSEA for unique Gene-sets...\n")
                results <- pblapply(GTEx.info_use, function(x) clusterProfiler::GSEA(geneList=x, TERM2GENE=meltedtouse, minGSSize=5, maxGSSize=10000, verbose=F))
                message("Performing GSEA for unique Pathways...\n")
                resultsPathway <- pblapply(GTEx.info_use, function(x) clusterProfiler::GSEA(geneList=x, TERM2GENE=meltedtousePathway, minGSSize=5, maxGSSize=10000, verbose=F))

              }, error=function(e) {
                stop("The API service is not currently working... Use the local database instead, you can download it on: INCLUDE LINk")
              }
            )


          } else {
            tryCatch({
              localDatabase_use <- localDatabase[tissues]
              message("Performing GSEA for unique Gene-sets...\n")
              results <- pblapply(localDatabase_use, function(x) clusterProfiler::GSEA(geneList=x, TERM2GENE=meltedtouse, minGSSize=5, maxGSSize=10000, verbose=F))
              message("Performing GSEA for unique Pathways...\n")
              resultsPathway <- pblapply(localDatabase_use, function(x) clusterProfiler::GSEA(geneList=x, TERM2GENE=meltedtousePathway, minGSSize=5, maxGSSize=10000, verbose=F))


            }, error=function(e){
              stop("Seems that the database provided is not correct. Please make sure you are using the correct one.\nYou can download it on https://github.com/TranslationalBioinformaticsUnit/GeneSetCluster")
            })
          }

          tissue.df <- data.frame(matrix(0, nrow = length(results), ncol = length(unique(meltedtouse$Term))))
          tissue.dfPathway <- data.frame(matrix(0, nrow = length(resultsPathway), ncol = length(unique(meltedtousePathway$Term))))
          rownames(tissue.df) <- names(results)
          colnames(tissue.df) <- colnames(touse)[1:ncol(tissue.df)]
          rownames(tissue.dfPathway) <- names(resultsPathway)
          colnames(tissue.dfPathway) <- colnames(tousePathway)[1:ncol(tissue.dfPathway)]

          for (i in 1:length(results))
          {
            tissue.df[i, results[[i]]@result$ID] <- results[[i]]@result$NES
          }

          for (i in 1:length(resultsPathway))
          {
            tissue.dfPathway[i, resultsPathway[[i]]@result$ID] <- resultsPathway[[i]]@result$NES
          }


          if (clusterIndependent)
          {
            Object@dfTissueIndependent <- list(Geneset=tissue.df, Pathway=tissue.dfPathway)
          } else {
            Object@dfTissue <- list(Geneset=tissue.df, Pathway=tissue.dfPathway)
          }

          message("\n")
          message("[=========================================================]")
          message("[<<<<         ObtainTissueExpression END             >>>>>]")
          message("-----------------------------------------------------------")
          message("[You may want to plot the results using PlotTissueExpression next.]")

          return(Object)
        }
)
