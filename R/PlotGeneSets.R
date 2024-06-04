#' @title PlotGeneSets
#' @description
#' Plots the correlation matrix with the distances per cluster.
#' A functional annotation of the group of pathways based on GO terms descriptions.
#' Either Unique Gene-set or Pathway approach can be selected.
#'
#' @import colorRamp2
#' @importFrom simplifyEnrichment anno_word_cloud
#' @importFrom grDevices colorRampPalette
#' @importFrom ComplexHeatmap rowAnnotation HeatmapAnnotation Heatmap Legend draw ht_opt
#'
#' @param Object a PathwayObject
#' @param uniquePathways A boolean to specify the approach of unique Gene-set or unique Pathways.
#' FALSE if for unique Gene-set and TRUE for unique Pathways.
#' @param doORA a boolean to show the ORA per cluster.
#' @param wordcloud a boolean to show the wordcloud annotation per cluster.
#'
#' @return plot
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
#' \donttest{
#' GSEA.Object2 <- CombineGeneSets(Object = GSEA.Object1, threads=2)
#'
#' OptimalGeneSets(Object = GSEA.Object2, uniquePathways = FALSE, cluster=NULL,
#'                 method = "silhouette", max_cluster= 24, cluster_method = "kmeans",
#'                 main= "Kmeans for 24 clusters")
#'
#' OptimalGeneSets(Object = GSEA.Object2, uniquePathways = TRUE, cluster = NULL,
#'                 method = "silhouette", max_cluster= 24, cluster_method = "kmeans",
#'                 main= "Kmeans for 24 clusters")
#'
#' GSEA.Object3 <- ClusterGeneSets(Object = GSEA.Object2,
#'                 clusters = 3,
#'                 method = "kmeans")
#'
#' PlotGeneSets(GSEA.Object3, uniquePathways = FALSE)
#' PlotGeneSets(GSEA.Object3, uniquePathways = TRUE)
#' }
#'
#' @export
#'
setGeneric(name="PlotGeneSets",
           def=function(Object,
                        uniquePathways = FALSE,
                        doORA = TRUE,
                        wordcloud = TRUE)
           {
             standardGeneric("PlotGeneSets")
           }
)



setMethod(f="PlotGeneSets",
          signature="PathwayObject",
          definition=function(Object,
                              uniquePathways = FALSE,
                              doORA = TRUE,
                              wordcloud = TRUE)
          {

            if (is.null(Object@functionalAnnot))
            {
              message("Please first you should use ClusterGeneSet in ordert to group Gene-sets in clusters...")
              stop()
            }

            #scaling the matrix
            tmp_data <- Object@Data[[1]][order(Object@Data[[1]]$cluster),]
            Data.RR <- scaleCorMatrix(Object@Data.RR)
            Data.RR <- Data.RR[order(Object@Data[[1]]$cluster), order(Object@Data[[1]]$cluster)]

            DataPathway.RR <- scaleCorMatrix(Object@DataPathways.RR)

            if (uniquePathways == T)
            {
              Object@plot$aka2Unique <- Object@plot$aka2Unique[order(Object@plot$aka2Unique$Cluster), ]
              DataPathway.RR <- DataPathway.RR[rownames(Object@plot$aka2Unique), rownames(Object@plot$aka2Unique)]
            }


            if(uniquePathways == F)
            {
              keywords_ora <- Object@functionalAnnot$Geneset$ORA
              term <- Object@functionalAnnot$Geneset$GO

            }else{
              keywords_ora <- Object@functionalAnnot$Pathway$ORA
              term <- Object@functionalAnnot$Pathway$GO
            }

            keywords_ora <- lapply(keywords_ora, function (x)
              if (nchar(x$Description[1]) > 50) {
                paste0(substr(x$Description[1], 1, 40), "... ", x$ID[1])
              } else {
                x$Description[1]
              }
            )


            if (uniquePathways == F)
            {
              clus2 <- tmp_data$cluster
              names(clus2) <- tmp_data$RR_name
              go_id_list2 <- list()
              for (z in unique(clus2))
              {
                go_id_list2[[z]] <- names(clus2[clus2==z])
              }
              names(go_id_list2) <- as.character(1:length(go_id_list2))
              align_to <- go_id_list2
              for (i in 1:length(align_to)){
                for (j in 1:length(align_to[[i]]))
                  align_to[[i]][j] <- as.numeric(which(colnames(Data.RR) == align_to[[i]][j]))
              }
            }else{
              clus2 <- Object@plot$aka2Unique$Cluster
              names(clus2) <- rownames(Object@plot$aka2Unique)
              go_id_list2 <- list()
              for (z in unique(clus2))
              {
                go_id_list2[[z]] <- names(clus2[clus2==z])
              }

              names(go_id_list2) <- as.character(1:length(go_id_list2))
              align_to <- go_id_list2
              for (i in 1:length(align_to)){
                for (j in 1:length(align_to[[i]]))
                  align_to[[i]][j] <- as.numeric(which(colnames(DataPathway.RR) == align_to[[i]][j]))
              }
            }

            align_to <- lapply(align_to, as.numeric)

            annot_label <- rowAnnotation(keywords = anno_word_cloud(align_to, term),
                                         annotation_name_align=T)


            if (wordcloud == TRUE) {
              if(is.null(unlist(term))) {
                message("The pathways not contain GO IDs. Then the Semantic enrichment wordcloud cannot be generated.")
                wordcloud <- FALSE
                annot_label <- NULL
              } else {
                annot_label <- rowAnnotation(keywords = anno_word_cloud(align_to, term),
                                             annotation_name_align=T)
              }
            } else {
              annot_label <- NULL
            }

            if (uniquePathways == FALSE)
            {

              colors <- list(Group = Object@plot$aka3$Group,
                             Cluster = Object@plot$aka3$Cluster)
              color_fun <- colorRamp2(seq(min(Data.RR), max(Data.RR), length = 9),
                                      brewer.pal(9, "Reds"))
              annot_top <- HeatmapAnnotation(df=Object@plot$aka2[order(Object@Data[[1]]$cluster),],
                                             show_legend = T, col = colors)
              row_annot <- rowAnnotation(df=Object@plot$aka2[order(Object@Data[[1]]$cluster),],
                                         show_legend = F, col = colors)
              # ht_opt$message <- FALSE


              plot <- Heatmap(matrix = Data.RR, cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F ,
                              name="Similarity", col=color_fun, top_annotation = annot_top, left_annotation = row_annot, right_annotation = annot_label)


              if (doORA == TRUE) {
                col_ora <- Object@plot$aka3$Cluster[order(names(Object@plot$aka3$Cluster))]
                legend <- suppressMessages(
                                           Legend(labels = keywords_ora,
                                           title = "\nORA", legend_gp = gpar(fill = col_ora),
                                            nr=8,  title_position = "leftcenter-rot", )
                                           )


                mergedplot <- draw(plot, annotation_legend_list=legend, annotation_legend_side = "right")

              } else {
                mergedplot <- plot
              }

            }else{

              df_metadata <- Object@plot$aka2Unique[, -ncol(Object@plot$aka2Unique)]
              df_cluster <- as.data.frame(as.character(Object@plot$aka2Unique$Cluster))
              rownames(df_cluster) <- rownames(Object@plot$aka2Unique)
              colnames(df_cluster) <- "Cluster"

              col_cluster <- list(Cluster = Object@plot$aka3Unique$Cluster)
              col_group <- list()
              for (i in colnames(df_metadata))
              {
                tmp_vector <- c("white", Object@plot$aka3Unique$Group[i])
                names(tmp_vector) <- c(0, 1)
                tmp_vector <- list(i=tmp_vector)
                col_group <- append(col_group, tmp_vector)
              }

              names(col_group) <- names(Object@plot$aka3Unique$Group)

              annot_top <- HeatmapAnnotation(df=df_metadata, show_legend = F, col=col_group, annotation_name_gp = gpar(fontsize = 10))
              row_annot <- rowAnnotation(df=df_cluster, show_legend = T, col = col_cluster)


              # ht_opt$message <- FALSE

              color_fun <- colorRamp2(seq(min(DataPathway.RR), max(DataPathway.RR), length = 9), brewer.pal(9, "Reds"))


              plot <- Heatmap(DataPathway.RR, cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F ,
                              name="Similarity", col=color_fun, top_annotation = annot_top, right_annotation = annot_label, left_annotation = row_annot)


              if (doORA == TRUE) {
                col_ora <- Object@plot$aka3Unique$Cluster[order(names(Object@plot$aka3Unique$Cluster))]
                legend <- suppressMessages(
                                           Legend(labels = keywords_ora,
                                           title = "\nORA", legend_gp = gpar(fill = col_ora),
                                           nr=8,  title_position = "leftcenter-rot")
                                           )

                mergedplot <- draw(plot, annotation_legend_list=legend, annotation_legend_side = "bottom")
              } else {
                mergedplot <- plot
              }
            }

            return(mergedplot)

          }
)

