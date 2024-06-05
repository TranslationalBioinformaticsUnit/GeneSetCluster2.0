#' @title GetSTRINGdbPerGeneSets
#'
#' @description
#' Connects to the STRING database in order to obtain predictions of
#' protein-protein interactions.
#'
#' @import STRINGdb
#' @importFrom limma strsplit2
#' @importFrom dplyr distinct
#'
#' @param Object A PathwayObject.
#' @param uniquePathways A boolean to specify the approach of unique Gene-set or unique Pathways.
#' FALSE if for unique Gene-set and TRUE for unique Pathways.
#' @param clusterIndependent A boolean to specify the clustering approach.
#' TRUE for seriation-based clustering approach and FALSE for classic approach.
#' @param plot.input A vector with group names of the different gene set experiments
#' @param threshold A numeric for the threshold of the string analysis
#' @param version Version of StringDB to run (currently 12)
#'
#' @return PathwayObject
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
#'
#' GSEA.Object4 <- GetSTRINGdbPerGeneSets(Object = GSEA.Object3, uniquePathways=FALSE,
#'                                        clusterIndependent=FALSE)
#' GSEA.Object4 <- GetSTRINGdbPerGeneSets(Object = GSEA.Object4, uniquePathways=FALSE,
#'                                        clusterIndependent=TRUE)
#' }
#'
#' @export
#'
setGeneric(name="GetSTRINGdbPerGeneSets",
           def=function(Object, uniquePathways=FALSE, clusterIndependent=FALSE,
                        plot.input="All", threshold = 100, version = "12")
           {
             standardGeneric("GetSTRINGdbPerGeneSets")
           }
)

setMethod(f="GetSTRINGdbPerGeneSets",
          definition = function(Object, uniquePathways=FALSE, clusterIndependent=FALSE,
                                plot.input="All", threshold = 100, version = "12")
          {

          message("[==================================================================]")
          message("[<<<<             GetSTRINGdbPerGeneSets START                >>>>>]")
          message("--------------------------------------------------------------------")

          clus.mol.ls <- list()

          if (clusterIndependent)
          {
            checkIndependent(Object)
            if (uniquePathways)
            {
              df_use <- Object@plot$aka2IndependentUnique
              df_use$Pathways <- rownames(df_use)

              df_unique <- Object@Data[[1]] %>%
                dplyr::distinct(Pathways, .keep_all = TRUE)
              rownames(df_unique) <- df_unique$Pathways

              df_use$Molecules <- df_unique[rownames(df_use), "unionMolecules"]

            } else {

              df_use <- Object@plot$aka2Independent
              df_unique <- Object@Data[[1]]
              rownames(df_unique) <- df_unique$RR_name
              df_use$Molecules <- df_unique[rownames(df_use), "Molecules"]

            }
          } else {
            checkDependent(Object)
            if (uniquePathways)
            {
              df_use <- Object@plot$aka2Unique
              df_use$Pathways <- rownames(df_use)

              df_unique <- Object@Data[[1]] %>%
                distinct(Pathways, .keep_all = TRUE)
              rownames(df_unique) <- df_unique$Pathways

              df_use$Molecules <- df_unique[rownames(df_use), "unionMolecules"]
            } else {
              df_use <- Object@Data[[1]]
              df_use$Cluster <- df_use$cluster
            }

          }


          for (clus.i in 1:max(df_use$Cluster)) {
            clus.x <- df_use[df_use$Cluster == clus.i, ]
            clus.x$Molecules <- as.character(clus.x$Molecules)
            clus.mol.ls[[clus.i]] <- unique(as.vector(limma::strsplit2(clus.x$Molecules,
                                                                       split = Object@metadata$seperator[1])))
          }
          unique.mol <- unique(do.call(what = c, args = clus.mol.ls))
          mol.unique.df <- as.data.frame(matrix(0, nrow = length(unique.mol),
                                                ncol = length(clus.mol.ls)))
          rownames(mol.unique.df) <- unique.mol
          colnames(mol.unique.df) <- paste("Cluster_", 1:length(clus.mol.ls),
                                           sep = "")
          for (clus.i in 1:max(Object@Data[[1]]$cluster)) {
            mol.unique.df[clus.mol.ls[[clus.i]], clus.i] <- 1
          }
          df <- mol.unique.df
          if (Object@metadata$organism[1] == "org.Hs.eg.db") {
            species = 9606
          }
          if (Object@metadata$organism[1] == "org.Mm.eg.db") {
            species = 10090
          }
          string_db <- STRINGdb$new(version = version, species = species,
                                    score_threshold = threshold, input_directory = "")
          df.plot <- list()
          for (df.i in 1:ncol(df)) {

            genes.i <- rownames(df[df[, df.i] == 1, ])


            if (length(genes.i) == 0) {
              message(paste("Cluster ", clus.i, " has no unique genes",
                            sep = ""))
              message("Moving to next cluster")
              (next)()
            }
            genes.i <- as.matrix(genes.i)
            colnames(genes.i) <- "gene"
            mapped.genes.i <- string_db$map(genes.i, "gene", removeUnmappedRows = TRUE)
            Enrichment.genes.i <- string_db$ppi_enrichment(string_ids = mapped.genes.i$STRING_id)

            if (plot.input == "All") {
              plot.input <- length(mapped.genes.i$STRING_id)
            }
            list.i <- list(map = mapped.genes.i,
                           img = string_db$get_png(mapped.genes.i$STRING_id[1:plot.input], payload_id = NULL),
                           input = string_db$get_summary(mapped.genes.i$STRING_id[1:plot.input], required_score = threshold)
                           # , link = string_db$get_link(mapped.genes.i$STRING_id[1:plot.input], payload_id = NULL, required_score = 100)
                           )

            list.i$input <- gsub("proteins", "All_GenesPerCluster",
                                 list.i$input)

            df.plot[[df.i]] <- list.i
          }

          message("-------------------------------------------------------------------")
          message("[<<<<             GetSTRINGdbPerGeneSets ends                >>>>>]")
          message("[=================================================================]")
          message("To view the StringPlots use plotSTRINGdbPerGeneSets")

          metadata <- paste0("Performed with the following parameters:\n",
                             "uniquePathways: ", uniquePathways,
                             "\nclusterIndependent: ", clusterIndependent)

          Object@stringdbInfo <- list(df.plot=df.plot,
                                      metdata=metadata)

          return(Object)
})


