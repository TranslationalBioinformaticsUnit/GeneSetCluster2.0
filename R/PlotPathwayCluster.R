#' @title PlotPathwayCluster
#'
#' @description
#' Plots the correlation matrix clustered by the seriation-based method.
#' A functional annotation of the group of pathways based on GO terms descriptions.
#' Either Unique Gene-set or Pathway approach can be selected.
#'
#' @import colorRamp2
#' @importFrom simplifyEnrichment anno_word_cloud
#' @importFrom grDevices colorRampPalette
#' @importFrom ComplexHeatmap rowAnnotation HeatmapAnnotation Heatmap Legend draw ht_opt anno_mark
#'
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
#' \donttest{
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
#'
#' GSEA.Object2 <- CombineGeneSets(Object = GSEA.Object1, threads=2)
#' GSEA.Object3 <- ClusterIndependentGeneSet(Object = GSEA.Object2)
#' PlotPathwayCluster(GSEA.Object3, uniquePathways = FALSE)
#' PlotPathwayCluster(GSEA.Object3, uniquePathways = TRUE)
#' }
#'
#' @export
#'
setGeneric(name="PlotPathwayCluster",
           def=function(Object, doORA = TRUE, wordcloud = TRUE, uniquePathways = FALSE)
           {
             standardGeneric("PlotPathwayCluster")
           }
)

setMethod(f="PlotPathwayCluster",
          signature = "PathwayObject",
          definition = function(Object, doORA = TRUE, wordcloud = TRUE, uniquePathways = FALSE)
          {

          #check if cluster independent slot is null
          checkIndependent(Object)
          # ht_opt$message <- FALSE
          if (uniquePathways)
          {
            mat_sym <- scaleCorMatrix(Object@DataPathways.RR)
            mat_cor <- mat_sym[Object@cIndependentMethod[[2]][[1]],
                               Object@cIndependentMethod[[2]][[1]]]

            keywords_ora <- Object@functionalAnnotIndependent$Pathway$ORA

            keywords_ora <- lapply(keywords_ora, function (x)
              if (nchar(x$Description[1]) > 50) {
                paste0(substr(x$Description[1], 1, 40), "... ", x$ID[1])
              } else {
                x$Description[1]
              }
            )

            term <- Object@functionalAnnotIndependent$Pathway$GO

            res <- obtainDefCluster(mat_cor)

            go_id_list <- res[which(lapply(res,
                                           function(x) length(x))>=Object@plot$aka3IndependentUnique$optimalNumberPathway)]

          } else {

            mat_cor <- scaleCorMatrix(Object@Data.RR)[Object@cIndependentMethod[[1]][[1]],
                                                      Object@cIndependentMethod[[1]][[1]]]

            keywords_ora <- Object@functionalAnnotIndependent$Geneset$ORA

            keywords_ora <- lapply(keywords_ora, function (x)
                              if (nchar(x$Description[1]) > 50) {
                                paste0(substr(x$Description[1], 1, 40), "... ", x$ID[1])
                              } else {
                                x$Description[1]
                              }
                            )

            term <- Object@functionalAnnotIndependent$Geneset$GO

            res <- obtainDefCluster(mat_cor)
            go_id_list <- res[which(lapply(res,
                                           function(x) length(x))>=Object@plot$aka3Independent$optimalNumberPathway)]

          }

          names(go_id_list) <- 1:length(go_id_list)

          align_to <- go_id_list

          for (i in 1:length(align_to)){
            for (j in 1:length(align_to[[i]]))
              align_to[[i]][j] <- as.numeric(which(colnames(mat_cor) == align_to[[i]][j]))
          }

          align_to <- lapply(align_to, as.numeric)


          if (wordcloud == TRUE) {
            if (is.null(unlist(term))) {
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

          rowNumbers <- rowAnnotation(foo = anno_mark(at = unlist(lapply(align_to, function(x) round(mean(x)))),
                                            labels = rep(paste0("Cluster ", 1:length(align_to))),
                                            link_width = unit(0,"mm")))

          if (uniquePathways)
          {
            df_metadata <- Object@plot$aka2IndependentUnique[, -ncol(Object@plot$aka2IndependentUnique)]
            df_cluster <- as.data.frame(as.character(Object@plot$aka2IndependentUnique$Cluster))
            rownames(df_cluster) <- rownames(Object@plot$aka2IndependentUnique)
            colnames(df_cluster) <- "Cluster"

            if (10 %in% df_cluster$Cluster) { #case there is 10 to not transform it into 1
              df_cluster$Cluster <- gsub(10,"ten",df_cluster$Cluster)
              df_cluster$Cluster <- gsub(0,"",df_cluster$Cluster)
              df_cluster$Cluster <- gsub("ten",10,df_cluster$Cluster)

            } else {
              df_cluster$Cluster <- gsub(0,"",df_cluster$Cluster)
            }

            col_cluster <- list(Cluster = Object@plot$aka3IndependentUnique$Cluster)
            col_cluster$Cluster <- col_cluster$Cluster[-which(col_cluster$Cluster=="white")]

            col_group <- list()
            for (i in colnames(df_metadata))
            {
              tmp_vector <- c("white", Object@plot$aka3IndependentUnique$Group[i])
              names(tmp_vector) <- c(0, 1)
              tmp_vector <- list(i=tmp_vector)
              col_group <- append(col_group, tmp_vector)
            }

            names(col_group) <- names(Object@plot$aka3IndependentUnique$Group)

            annot_top <- HeatmapAnnotation(df=df_metadata, show_legend = F,
                                           col=col_group, annotation_name_gp = gpar(fontsize = 10))
            row_annot <- rowAnnotation(df=df_cluster, show_legend = T, col = col_cluster, na_col="white")

            # ComplexHeatmap::ht_opt$message <- FALSE
            color_fun <- colorRamp2(seq(min(mat_cor), max(mat_cor), length = 9),
                                    brewer.pal(9, "Reds"))

            plot <- Heatmap(mat_cor, cluster_rows = F, cluster_columns = F,
                            show_row_names = F, show_column_names = F ,
                            name="Similarity", col=color_fun,
                            top_annotation = annot_top, right_annotation = annot_label,
                            left_annotation = row_annot) + rowNumbers

            if (doORA == TRUE) {
              legend <- suppressMessages(
                               Legend(labels = keywords_ora,
                               title = "\nORA", legend_gp =  gpar(fill = col_cluster$Cluster),
                               nr=8,  title_position = "leftcenter-rot")
                               )

              mergedplot <- draw(plot, annotation_legend_list = legend, annotation_legend_side = "bottom")
            } else {
              mergedplot <- plot
            }

          } else {

            col_cluster <- list(Cluster = Object@plot$aka3Independent$Cluster)
            col_cluster$Cluster <- col_cluster$Cluster[-which(col_cluster$Cluster=="white")]
            col_all <- list(Cluster = Object@plot$aka3Independent$Cluster,
                            Group = Object@plot$aka3Independent$Group)

            df_cluster <- as.data.frame(Object@plot$aka2Independent)

            if (10 %in% df_cluster$Cluster) { #case there is 10 to not transform it into 1
              df_cluster$Cluster <- gsub(10,"ten",df_cluster$Cluster)
              df_cluster$Cluster <- gsub(0,"",df_cluster$Cluster)
              df_cluster$Cluster <- gsub("ten",10,df_cluster$Cluster)

            } else {
              df_cluster$Cluster <- gsub(0,"",df_cluster$Cluster)
            }


            rowNumbers <- rowAnnotation(foo = anno_mark(at = unlist(lapply(align_to, function(x) round(mean(x)))),
                                                        labels = rep(paste0("Cluster ", 1:length(align_to))),
                                                        link_width = unit(0,"mm")))

            color_fun <- colorRamp2(seq(min(mat_cor), max(mat_cor), length = 9), brewer.pal(9, "Reds"))
            annot_top <- HeatmapAnnotation(df=df_cluster, show_legend = T,
                                           col = col_all,
                                           na_col="white")
            row_annot <- rowAnnotation(df=df_cluster,
                                       show_legend = F, col = col_all,
                                       na_col="white")


            plot <- Heatmap(matrix = mat_cor, cluster_rows = F, cluster_columns = F,
                            show_row_names = F, show_column_names = F ,
                            name="Similarity", col=color_fun, top_annotation = annot_top,
                            left_annotation = row_annot, right_annotation = annot_label) + rowNumbers


            if (doORA == TRUE) {
              legend <- suppressMessages(
                               Legend(labels = keywords_ora,
                               title = "\nORA", legend_gp = gpar(fill = col_cluster$Cluster),
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

