
setGeneric(name="PlotTissueExpression",
           def=function(Object, all = FALSE, uniquePathways=FALSE, clusterIndependent=FALSE, showZscore=FALSE)
           {
             standardGeneric("PlotTissueExpression")
           }
)

setMethod(f="PlotTissueExpression",
          signature = "PathwayObject",
          definition = function(Object, all=FALSE, uniquePathways=FALSE, clusterIndependent=FALSE, showZscore=FALSE)
  {

    if (clusterIndependent)
    {
      if (is.null(Object@dfTissueIndependent) | anyNA(Object@dfTissueIndependent)) {
        stop("First you have to run TissueExpressionPerGeneSet clusterIndependent=TRUE to obtain tissue expression information.")
      }

      # if (uniquePathways)
      # {
      #   df <- Object@dfTissueIndependent$Pathway
      # } else {
      #   df <- Object@dfTissueIndependent$Geneset
      # }
      df <- Object@dfTissueIndependent

    } else {
      if (is.null(Object@dfTissue) | anyNA(Object@dfTissue)) {
        stop("First you have to run TissueExpressionPerGeneSet clusterIndependent=FALSE to obtain tissue expression information.")
      }
      # if (uniquePathways)
      # {
      #   df <- Object@dfTissue$Pathway
      # } else {
      #   df <- Object@dfTissue$Geneset
      # }
      df <- Object@dfTissue
    }


    # df <- Object@dfTissue
    df_keep <- df
    df$tissue <- rownames(df)
    plot_melt_mtx <- melt(df, id="tissue")
    colnames(plot_melt_mtx)[1:3] <- c("cluster", "Dataset", "value")
    plot_cluster_size <- aggregate(value ~ cluster, data = plot_melt_mtx, FUN = sum)
    colnames(plot_cluster_size) <- c("cluster", "value")

    if (nrow(df) > 1)
    {
      df_tissue_z <- apply(df_keep, 2, function(x) (x-mean(x))/sd(x))
      df_tissue_z[is.na(df_tissue_z)] <- 0
      plot_melt_mtx_z <- melt(df_tissue_z)
      colnames(plot_melt_mtx_z)[1:3] <- c("cluster", "Dataset", "value")
      plot_cluster_size_z <- aggregate(value ~ cluster, data = plot_melt_mtx_z, FUN = sum)
      colnames(plot_cluster_size_z) <- c("cluster", "value")
      if (nrow(df_keep) < 54)
      {
        all <- TRUE
      }

      if (all == FALSE)
      {
        plot_cluster_size <- plot_cluster_size[order(plot_cluster_size[,"value"],
                                                     decreasing = T),]
        keep <- plot_cluster_size$cluster[1:20] # top 20
        plot_cluster_size <- plot_cluster_size[plot_cluster_size$cluster %in% keep,]

        plot_cluster_size_z <- plot_cluster_size_z[plot_cluster_size_z$cluster %in% keep,]
        plot_cluster_size_z <- plot_cluster_size_z[rownames(plot_cluster_size),]
        plot_melt_mtx_z <- plot_melt_mtx_z[plot_melt_mtx_z$cluster %in% keep,]
        plot_melt_mtx <- plot_melt_mtx[plot_melt_mtx$cluster %in% keep,]
      }

      sorted_labels <- paste(sort(levels(plot_cluster_size$cluster), decreasing = F))
      sorted_labels <- paste(sort(plot_cluster_size$cluster, decreasing = F))

      plot_cluster_size$cluster <- factor(plot_cluster_size$cluster,levels = sorted_labels)
      plot_cluster_size_z$cluster <- factor(plot_cluster_size_z$cluster,levels = sorted_labels)
      plot_melt_mtx_z$cluster <- factor(plot_melt_mtx_z$cluster,levels = sorted_labels)
      plot_melt_mtx_z$tot_Ncells <- ave(plot_melt_mtx_z$value, plot_melt_mtx_z$cluster, FUN=sum)
      plot_melt_mtx_z$percent <- plot_melt_mtx_z$value*100/plot_melt_mtx_z$tot_Ncells

      plot_melt_mtx$cluster <- factor(plot_melt_mtx$cluster,levels = sorted_labels)
      plot_melt_mtx$tot_Ncells <- ave(plot_melt_mtx$value, plot_melt_mtx$cluster, FUN=sum)
      plot_melt_mtx$percent <- plot_melt_mtx$value*100/plot_melt_mtx$tot_Ncells

      cluster_size <- plot_cluster_size
      cluster_size_z <- plot_cluster_size_z


    } else {
      if (showZscore)
      {
        message("Not showing Z-score as only one tissue has been used.")
        showZscore <- FALSE
      }

      sorted_labels <- paste(sort(levels(plot_cluster_size$cluster), decreasing = F))
      sorted_labels <- paste(sort(plot_cluster_size$cluster, decreasing = F))

      plot_cluster_size$cluster <- factor(plot_cluster_size$cluster,levels = sorted_labels)

      plot_melt_mtx$cluster <- factor(plot_melt_mtx$cluster,levels = sorted_labels)
      plot_melt_mtx$tot_Ncells <- ave(plot_melt_mtx$value, plot_melt_mtx$cluster, FUN=sum)
      plot_melt_mtx$percent <- plot_melt_mtx$value*100/plot_melt_mtx$tot_Ncells

      cluster_size <- plot_cluster_size
    }


    if(clusterIndependent){
      if(uniquePathways){
        my_pal <- Object@plot$aka3IndependentUnique$Cluster
        my_pal <- my_pal[-which(my_pal=="white")]
      }else{
        my_pal <- Object@plot$aka3Independent$Cluster
        my_pal <- my_pal[-which(my_pal=="white")]
      }

    }else{
      if(uniquePathways){
        my_pal <- Object@plot$aka3Unique$Cluster
      }else{
        my_pal <- Object@plot$aka3$Cluster
      }
    }
    names(my_pal) <- NULL

    #Plot total proportion
    p1 <- ggplot(plot_cluster_size, aes(x = value, y = reorder(cluster, value)))

    #Plot Barplot by cluster percentages
    p4 <- ggplot(plot_melt_mtx, aes(x=reorder(cluster, value), y=percent, fill=Dataset))

    #Plot Barplot by cluster values NES
    p5 <- ggplot(plot_melt_mtx, aes(x=reorder(cluster, value), y=value, fill=Dataset))

    p1 <- p1 + geom_bar(position="dodge", stat="identity",fill = "grey60", width = 0.85) +
      theme_bw() + xlab("Sum of NES values") +
      geom_text(aes(y=cluster,x=1.1,label=round(value, digits = 2),hjust="bottom"), size=4) + ylab("")+
      theme(axis.title = element_text(size = 16), axis.text.x = element_text(size = 12),
            axis.text.y = element_blank(), axis.ticks.y=element_blank())


    p4 <- p4 + geom_bar(position="dodge", stat="identity", width = 0.7) + theme_bw() + coord_flip() +
      scale_fill_manual(values = my_pal) +
      ylab("NES per cluster") + xlab("Tissue") +
      theme(legend.position="top", axis.title = element_text(size = 16),
            axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 14))


    p5 <- p5 + geom_bar(position="dodge", stat="identity", width = 0.7) + theme_bw() + coord_flip() +
      scale_fill_manual(values = my_pal) +
      ylab("NES per cluster") + xlab("Tissue") + ylim(c(0, 3.5)) +
      theme(legend.position="top", axis.title = element_text(size = 16),
            axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 14))


    #Combine plots
    if (showZscore)
    {
      #Plot total z-score proportion
      p2 <- ggplot(plot_cluster_size_z, aes(x = value, y = reorder(cluster, value)))

      #Plot Barplot by cluster z percentages
      p3 <- ggplot(plot_melt_mtx_z, aes(x=reorder(cluster, value), y=percent, fill=Dataset))

      #Plot Barplot by cluster z values
      p6 <- ggplot(plot_melt_mtx_z, aes(x=reorder(cluster, value), y=value, fill=Dataset))

      p2 <- p2 + geom_bar(position="dodge", stat="identity",fill = "grey60", width = 0.85) +
        theme_bw() + xlab("Sum of Z-score") +
        geom_text(aes(y=cluster,x=1.1,label=round(value, digits = 2),hjust="bottom"), size=4) + ylab("")+
        theme(axis.title = element_text(size = 16), axis.text.x = element_text(size = 12),
              axis.text.y = element_blank(), axis.ticks.y=element_blank())

      p3 <- p3 + geom_bar(position="dodge", stat="identity", width = 0.7) + theme_bw() + coord_flip() +
        scale_fill_manual(values = my_pal) +
        ylab("Z-score per cluster") + xlab("Tissue") +
        theme(legend.position="top", axis.title = element_text(size = 16),
              axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 14))

      p6 <- p6 + geom_bar(position="dodge", stat="identity", width = 0.7) + theme_bw() + coord_flip() +
        scale_fill_manual(values = my_pal) +
        ylab("Z-Score per cluster") +  xlab("")+
        theme(legend.position="none", axis.title = element_text(size = 16),
              axis.text.x = element_text(size = 16),
              axis.text.y = element_blank(), axis.ticks.y=element_blank())


      finalPlot <- p5 + p6 + p2 + p1 + plot_layout(widths = c(4, 4, 2, 2))
    } else {
      finalPlot <- p5 + p1 + plot_layout(widths = c(4, 2))
    }


    message("[You may want to plot the cluster tissue expression using PlotClusterTissueExpression next.]")

    return(finalPlot)
}
)
