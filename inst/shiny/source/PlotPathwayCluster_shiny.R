#' PlotPathwayCluster_shiny
#'
#' Plots a correlation matrix showing the correlation between the pathways and a functional annotation of the group of pathways based on GO terms descriptions.
#' @import simplifyEnrichment
#' @import slam
#' @import GetoptLong
#' @import ComplexHeatmap
#' @import grid
#' @import RColorBrewer
#' @import colorRamp2
#' @import seriation
#'
#'
#' @param Object a pathway object
#' @param doORA boolean perform ORA analysis
#' @param wordcloud boolean perform boolean to add wordcloud annotations of each cluster
#' @param wordclouds a list with wordclouds.
#' @param uniquePathways logical to use either unique Gene-set or Pathways approach. If FALSE unique Gene-sets.
#' @param keywords_ora_inde a list with key words ora.
#'
#' @return plot
#'
#' @export
#'
setGeneric(name="PlotPathwayCluster_shiny",
           def=function(Object, doORA = TRUE, wordcloud = TRUE, wordclouds = "", uniquePathways = FALSE, keywords_ora_inde="")
           {
             standardGeneric("PlotPathwayCluster_shiny")
           }
)

#' PlotPathwayCluster_shiny
#'
#' @param Object a pathway object
#' @param doORA boolean to perform ORA analysis
#' @param wordcloud boolean to add wordcloud annotations of each cluster
#' @param wordclouds a list with wordclouds.
#' @param uniquePathways logical to use either unique Gene-set or Pathways approach. If FALSE unique Gene-sets.
#' @param keywords_ora_inde a list with key words ora.
#'
#' @return plot
#'
#' @examples

setMethod(f="PlotPathwayCluster_shiny",
          signature = "PathwayObject",
          definition = function(Object, doORA = TRUE, wordcloud = TRUE, wordclouds = "", uniquePathways = FALSE, keywords_ora_inde="")
          {

          #check if cluster independent slot is null
          checkIndependent(Object)

          if (uniquePathways)
          {
            mat_sym <- scaleCorMatrix(Object@DataPathways.RR)
            mat_cor <- mat_sym[Object@cIndependentMethod[[1]][[1]],
                               Object@cIndependentMethod[[1]][[1]]]

            # keywords_ora <- Object@functionalAnnotIndependent$Pathway$ORA
            # 
            # keywords_ora <- lapply(keywords_ora, function (x)
            #   if (nchar(x@result$Description[1]) > 50) {
            #     paste0(substr(x@result$Description[1], 1, 40), "... ", x@result$ID[1])
            #   } else {
            #     x@result$Description[1]
            #   }
            # )

            #term <- Object@functionalAnnotIndependent$Pathway$GO
            res <- obtainDefCluster(mat_cor)

            go_id_list <- res[which(lapply(res,
                                           function(x) length(x))>=Object@plot$aka3IndependentUnique$optimalNumberPathway)]

          } else {
            mat_sym <- scaleCorMatrix(Object@Data.RR)
            mat_cor <- mat_sym[Object@cIndependentMethod[[1]][[1]],
                               Object@cIndependentMethod[[1]][[1]]]

            # keywords_ora <- Object@functionalAnnotIndependent$Geneset$ORA
            # 
            # keywords_ora <- lapply(keywords_ora, function (x)
            #                   if (nchar(x@result$Description[1]) > 50) {
            #                     paste0(substr(x@result$Description[1], 1, 40), "... ", x@result$ID[1])
            #                   } else {
            #                     x@result$Description[1]
            #                   }
            #                 )
            # 
            # term <- Object@functionalAnnotIndependent$Geneset$GO
            
            res <- obtainDefCluster(mat_cor)

            go_id_list <- res[which(lapply(res,
                                           function(x) length(x))>=Object@plot$aka3Independent$optimalNumberPathway)]

          }


          print(go_id_list)
          names(go_id_list) <- 1:length(go_id_list)

          align_to <- go_id_list

          for (i in 1:length(align_to)){
            for (j in 1:length(align_to[[i]]))
              align_to[[i]][j] <- as.numeric(which(colnames(mat_cor) == align_to[[i]][j]))
          }

          align_to <- lapply(align_to, as.numeric)


          if (wordcloud == TRUE) {
            term<-wordclouds$GO
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

            #df_cluster$Cluster <- gsub(0,"",df_cluster$Cluster)
            
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
            
            names(col_group) <- names(Object@plot$aka3Unique$Group)
            
            annot_top <- HeatmapAnnotation(df=df_metadata, show_legend = F,
                                           col=col_group, annotation_name_gp = gpar(fontsize = 10))
            row_annot <- rowAnnotation(df=df_cluster, show_legend = T, col = col_cluster, na_col="white")
            
            ht_opt$message <- FALSE
            color_fun <- colorRamp2(seq(min(mat_cor), max(mat_cor), length = 9),
                                    brewer.pal(9, "Reds"))
            
            plot <- Heatmap(mat_cor, cluster_rows = F, cluster_columns = F,
                            show_row_names = F, show_column_names = F ,
                            name="Similarity", col=color_fun,
                            top_annotation = annot_top, right_annotation = annot_label,
                            left_annotation = row_annot) + rowNumbers

            if (doORA == TRUE) {
              #legend <- Legend(labels = rep(paste0(1:length(keywords_ora), "- ", keywords_ora)),
              # legend <- Legend(labels = keywords_ora_inde,                
              #                  title = "\nORA", legend_gp = gpar(fill = Object@plot$aka3IndependentUnique$Cluster[!names(Object@plot$aka3IndependentUnique$Cluster)=="0"]),
              #                  nr=8,  title_position = "leftcenter-rot")
              legend <- suppressMessages(
                Legend(labels = keywords_ora_inde,
                       title = "\nORA", legend_gp =  gpar(fill = col_cluster$Cluster),
                       nr=8,  title_position = "leftcenter-rot")
              )

              mergedplot <- draw(plot, annotation_legend_list = legend, annotation_legend_side = "bottom")
            } else {
              mergedplot <- plot
            }

          } else {

              # col_cluster <- list(Cluster = Object@plot$aka3Independent$Cluster)
              # col_cluster$Cluster <- col_cluster$Cluster[-which(col_cluster$Cluster=="white")]
              # col_all <- list(Cluster = Object@plot$aka3Independent$Cluster,
              #                 Group = Object@plot$aka3Independent$Group)
              # 
              # df_cluster <- as.data.frame(Object@plot$aka2Independent)
              # df_cluster$Cluster <- gsub(0,"",df_cluster$Cluster)
            
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
              
              
              color_fun <- colorRamp2(seq(min(mat_cor), max(mat_cor), length = 9), brewer.pal(9, "Reds"))
              annot_top <- HeatmapAnnotation(df=df_cluster, show_legend = T,
                                             col = col_all,
                                             na_col="white")
              row_annot <- rowAnnotation(df=df_cluster,
                                         show_legend = F, col = col_all,
                                         na_col="white")
              
              
              rowNumbers <- rowAnnotation(foo = anno_mark(at = unlist(lapply(align_to, function(x) round(mean(x)))),
                                                          labels = rep(paste0("Cluster ", 1:length(align_to))),
                                                          link_width = unit(0,"mm")))
              
              plot <- Heatmap(matrix = mat_cor, cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F ,
                              name="Similarity", col=color_fun, top_annotation = annot_top, left_annotation = row_annot, right_annotation = annot_label) + rowNumbers
              
            
              if (doORA == TRUE) {
                #legend <- Legend(labels = rep(paste0(1:length(keywords_ora), "- ", keywords_ora)),
                # legend <- Legend(labels = keywords_ora_inde,  
                #                  title = "\nORA", legend_gp = gpar(fill = col_cluster$Cluster),
                #                  nr=8,  title_position = "leftcenter-rot")
                
                legend <- suppressMessages(
                  Legend(labels = keywords_ora_inde,
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

