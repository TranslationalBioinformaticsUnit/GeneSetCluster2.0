setGeneric(name="PlotGeneSets_shiny",
           def=function(Object,
                        uniquePathways = F,
                        doORA = T,
                        wordcloud = T,
                        wordclouds = "",
                        fontsize = 5,
                        legend = T,
                        annotation.mol=F,
                        main="",
                        RR.max = "",
                        cluster.order = "",
                        keywords_ora="")
           {
             standardGeneric("PlotGeneSets_shiny")
           }
)


setMethod(f="PlotGeneSets_shiny",
          definition=function(Object,
                              uniquePathways = F,
                              doORA = T,
                              wordcloud = T,
                              wordclouds = "",
                              fontsize = 5,
                              legend = T,
                              annotation.mol=F,
                              main="",
                              RR.max = "",
                              cluster.order = "",
                              keywords_ora="")
          {
            #scaling the matrix
            Data.RR <- scaleCorMatrix(Object@Data.RR)
            DataPathway.RR <- scaleCorMatrix(Object@DataPathways.RR)

            if(!RR.max == "")
            {
              RR.max <- as.numeric(as.character(RR.max))
              for(rows.i in 1:nrow(Data.RR))
              {
                idx <- Data.RR[rows.i,] > RR.max
                names(idx) <- NULL
                Data.RR[rows.i,idx] <- RR.max
              }
            }
            if(!cluster.order[1]== "")
            {
              if(!length(cluster.order) == max(as.numeric(as.character(Object@plot$aka2$Cluster))))
              {
                message("Make sure youre the order of supplied clusters is the same length as the number of clusters")

              }else{
                order.x <- vector()
                if (uniquePathways == F) {

                  for(order.i in cluster.order)
                  {
                    order.x <- c(order.x, which(order.i == as.numeric(as.character(Object@plot$aka2$Cluster))))
                  }
                  Object@plot$aka2 <- Object@plot$aka2[order.x,]

                  Data.RR <- Data.RR[order.x, order.x]
                }else{
                  for(order.i in cluster.order)
                  {
                    order.x <- c(order.x, which(order.i == as.numeric(as.character(Object@plot$aka2Unique$Cluster))))
                  }
                  Object@plot$aka2Unique <- Object@plot$aka2Unique[order.x,]

                  DataPathway.RR <- DataPathway.RR[order.x, order.x]
                }
              }
            }else{
              if (uniquePathways == T)
              {
                Object@plot$aka2Unique <- Object@plot$aka2Unique[order(Object@plot$aka2Unique$Cluster), ]
                DataPathway.RR <- DataPathway.RR[rownames(Object@plot$aka2Unique), rownames(Object@plot$aka2Unique)]
              }
            }

            #performing ora
            if(uniquePathways == F)
            {
              clus <- Object@Data[[1]]$cluster
              names(clus) <- Object@Data[[1]]$Pathways
              # clus<-gsub("Cluster_", "", clus)
              #keywords_ora <- ORAperCluster(Object, doORA, length(unique(clus)), clus)
            }else{
              clus <- Object@plot$aka2Unique$Cluster
              names(clus) <- rownames(Object@plot$aka2Unique)
              # clus<-gsub("Cluster_", "", clus)
              #keywords_ora <- ORAperCluster(Object, doORA, length(unique(clus)), clus)
            }


            go_id_list <- list()
            for (z in unique(clus))
            {
              go_id_list[[z]] <- names(clus[clus==z])
            }


            if (wordcloud == T)
            {
              #performing keyword erichment to generate wordcloud
              if (checkGO(Object) == FALSE) {
                message("No GO terms have been detected in the pathways. The semantic enrichment word cloud will not be generated.")
                wordcloud = FALSE
              } else {
                #adapted from anno_word_cloud_from_GO function
                env_tdm_GO <- readRDS(system.file("extdata", "tdm_GO.rds", package = "simplifyEnrichment"))
                names(go_id_list) <- as.character(1:length(go_id_list))

                # keyword enrichment
                message(paste0("Performing keyword enrichment for "), length(go_id_list), " group(s) of pathways.")
                # term <- lapply(go_id_list, function(x, min_stat=5) {
                #   df <- keywordEnrichment(x, env_tdm_GO)
                #   df <- df[df$p <= min_stat, , drop = FALSE]
                #   data.frame(df[, 1], -log10(df$p))
                # })
                #term<-PerformORAGOperCluster(Object, uniquePathways=uniquePathways, clusterIndependent=FALSE)
                term<-wordclouds$GO
              }


              if (uniquePathways == F)
              {
                #plot labels
                clus2 <- Object@Data[[1]]$cluster
                names(clus2) <- Object@Data[[1]]$RR_name
                # clus2<-gsub("Cluster_", "", clus2)
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

            } else {
              annot_label <- NULL
            }

            #rowNumbers <- rowAnnotation(foo = anno_mark(at = unlist(lapply(align_to, function(x) round(mean(x)))),
            #                                            labels = rep(paste0("Cluster ", 1:length(align_to))),
            #                                            link_width = unit(0,"mm")))

            if (uniquePathways == F)
            {

              colors <- list(Group = Object@plot$aka3$Group, Cluster = Object@plot$aka3$Cluster)
              color_fun <- colorRamp2(seq(min(Data.RR), max(Data.RR), length = 9), brewer.pal(9, "Reds"))
              annot_top <- HeatmapAnnotation(df=Object@plot$aka2, show_legend = T, col = colors)
              row_annot <- rowAnnotation(df=Object@plot$aka2, show_legend = F, col = colors)

              if(annotation.mol==T)
              {
                plot <- Heatmap(matrix = Data.RR, cluster_rows = F, cluster_columns = F, show_row_names = T, show_column_names = F ,
                                name="Similarity", col=color_fun, top_annotation = annot_top, left_annotation = NULL, right_annotation = annot_label) + rowNumbers

              } else {

                plot <- Heatmap(matrix = Data.RR, cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F ,
                                name="Similarity", col=color_fun, top_annotation = annot_top, left_annotation = row_annot, right_annotation = annot_label)
              }

              if (doORA == TRUE) {
                legend <- Legend(labels = keywords_ora,
                                 title = "\nORA", legend_gp = gpar(fill = Object@plot$aka3$Cluster[order(names(Object@plot$aka3$Cluster))]),
                                 nr=length(Object@plot$aka3$Cluster[order(names(Object@plot$aka3$Cluster))]),  title_position = "leftcenter-rot")

                mergedplot <- draw(plot, annotation_legend_list=legend, annotation_legend_side = "bottom")
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

              # colors <- list(Group = Object@plot$aka3$Group, Cluster = Object@plot$aka3$Cluster)
              # color_fun <- colorRamp2(seq(min(Data.RR), max(Data.RR), length = 9), brewer.pal(9, "Reds"))
              # annot_top <- HeatmapAnnotation(df=Object@plot$aka2, show_legend = T, col = colors)


              ht_opt$message <- FALSE

              color_fun <- colorRamp2(seq(min(DataPathway.RR), max(DataPathway.RR), length = 9), brewer.pal(9, "Reds"))

              if (annotation.mol == T)
              {
                plot <- Heatmap(DataPathway.RR, cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = T ,
                                name="Similarity", col=color_fun, top_annotation = annot_top, right_annotation = annot_label, left_annotation = row_annot) + rowNumbers
              }else{
                plot <- Heatmap(DataPathway.RR, cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F ,
                                name="Similarity", col=color_fun, top_annotation = annot_top, right_annotation = annot_label, left_annotation = row_annot)
              }

              if (doORA == TRUE) {
                #legend <- Legend(labels = rep(paste0(1:length(keywords_ora), "- ", keywords_ora)),
                legend <- Legend(labels = keywords_ora,
                                 title = "\nORA", legend_gp = gpar(fill = Object@plot$aka3Unique$Cluster[order(names(Object@plot$aka3Unique$Cluster))]),
                                 nr=length(Object@plot$aka3Unique$Cluster[order(names(Object@plot$aka3Unique$Cluster))]),  title_position = "leftcenter-rot")

                mergedplot <- draw(plot, annotation_legend_list=legend, annotation_legend_side = "bottom")
              } else {
                mergedplot <- plot
              }
            }

            return(mergedplot)

          }
)

