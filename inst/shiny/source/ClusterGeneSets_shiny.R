ClusterGeneSets_shiny <- function(Object, clusters = 5, method = "kmeans", order = "group", molecular.signature = "All", user_function)
{
  message("[=========================================================]")
  message("[<<<<            ClusterGeneSets START               >>>>>]")
  message("-----------------------------------------------------------")

  if(nrow(Object@Data.RR)==0)
  {
    message("Make sure youre object has been combines by CombineGeneSets")
    stop()
  }

  #Currently supported ways of clustering
  #Kmeans  #kmeans_group  #mclust #mclust_group #Hierarchical #Hierarchical_group
  #Correct spelling mistakes#
  if(method %in% c("kmeans", "Kmeans")){method <- "kmeans"}
  if(method %in% c("kmeans_group","Kmeans_group","Kmeans_Group", "kmeans group","Kmeans group")){method <- "kmeans_group"}
  #
  if(method %in% c("HierArchical", "hierarchical", "Hierarchical")){method <- "Hierarchical"}
  if(method %in% c("Hierarchical_group","hierarchical_group", "Hierarchical group","hierarchical group", "hierarchical_group")){method <- "Hierarchical_group"}

  #####################################################
  if(clusters > nrow(Object@Data[[1]]))
  {
    message("Too many cluster inputted, reducing it to the number of datapoints")
    clusters <- nrow(Object@Data[[1]])
  }

  if((length(clusters) > 1 & !grepl("group", method)))
  {
    message("Too many different cluster options inputted, reducing it to the number of datapoints")
    clusters <- clusters[1]
  }

  if((!length(clusters) == nrow(Object@PData) & grepl("group", method)))
  {
    message("please input a vector with the same length as the number of groups, reducing it to the 1st number")
    clusters <- clusters[1]
  }

  if(sum(method %in% c("kmeans", "kmeans_group", "Hierarchical", "Hierarchical_group")) < 1)
  {
    message("Method is not a supported method of clustering")
    stop()
  }
  if(length(method) > 1)
  {
    message("Too many methods supplied, pick between Kmeans, kmeans_group, mclust, mclust_group")
    stop()
  }

  ####################################
  #------------Combine DF------------#
  ####################################
  canonical.df <-Object@Data[[1]]
  if(molecular.signature == "All")
  {
    message("Using all Gene sets")
    canonical.df <-Object@Data[[1]]
  }
  if(molecular.signature == "Unique")
  {
    message("Using Unique Gene sets")
    if(Object@metadata[,"display"] == "Condensed")
    {
      warning("Display is set to Condensed, this will lead to a loss of information")
      warning("Advised is to set display to Expanded before filtering the unique molecular.signature")
    }
    Object@Data.RR <- Object@Data.RR[!duplicated(canonical.df$Molecules),
                                     !duplicated(canonical.df$Molecules)]

    canonical.df <- canonical.df[!duplicated(canonical.df$Molecules),]
    Object@metadata[,"mol.signature"] <- rep("Unique", times = nrow(metadata))

  }


  ####################################################
  #------------Kmeans of clusters--------------------#
  ####################################################
  if(method == "kmeans")
  {
    message("Running kmeans")

    if(length(clusters) > 1)
    {
      message("Too many clusters supplied, Taking first cluster")
      clusters <- clusters[1]
    }
    canonical.df$cluster <- kmeans(x = Object@Data.RR, centers = clusters)$cluster
    clusterUnique <- kmeans(x = Object@DataPathways.RR, centers = clusters)$cluster
  }


  ##############################################################
  #------------Kmeans of clusters per group--------------------#
  ##############################################################
  if(method == "kmeans_group")
  {
    print("Running kmeans per group seperatly")
    if(!length(clusters) == nrow(Object@PData))
    {
      message("Wrong number of clusters supplied")
      stop()
    }

    temp.cl2 <- vector()
    for(i in 1:length(unique(as.character(Object@Data[[1]]$Groups))))
    {

      idx.ls <- Object@Data[[1]][as.character(Object@Data[[1]]$Groups) == as.character(Object@PData[i,"Groupnames"]),"RR_name"]

      rr.temp <- Object@Data.RR[idx.ls,idx.ls]

      temp.cl1 <-  kmeans(x = rr.temp, centers = clusters[i])$cluster
      if(i >1){temp.cl1 <- temp.cl1 + max(temp.cl2)}
      temp.cl2 <- c(temp.cl2, temp.cl1)
    }
    canonical.df$cluster <- temp.cl2


  }

  #########################################################
  #------------Hierarchical clustering--------------------#
  #########################################################
  if(method == "Hierarchical")
  {
    message("Running Hierarchical")

    if(length(clusters) > 1){message("Too many clusters supplied, Taking first cluster")}
    canonical.df$cluster <- cutree(hclust(dist(t(Object@Data.RR)), method = "ward.D2"), k = clusters)
    clusterUnique <- cutree(hclust(dist(t(Object@DataPathways.RR)), method = "ward.D2"), k = clusters)

  }

  #########################################################
  #------------Hierarchical clustering--------------------#
  #########################################################
  if(method == "Hierarchical_group")
  {
    message("Running Hierarchical per group")

    if(length(clusters) > 1){message("Too many clusters supplied, Taking first cluster")}
    temp.cl2 <- vector()
    for(i in 1:nrow(Object@PData))
    {

      idx.ls <- Object@Data[[1]][as.character(Object@Data[[1]]$Group) == as.character(Object@PData[i,"Groupnames"]),"RR_name"]

      rr.temp <- Object@Data.RR[idx.ls,idx.ls]

      temp.clx<- cutree(hclust(dist(t(rr.temp)), method = "ward.D2"), k = clusters)

      temp.cl1 <- hclust(dist(t(rr.temp)), method = "ward.D2")$order
      if(i >1){temp.cl1 <- temp.cl1 + max(temp.cl2)}
      temp.cl2 <- c(temp.cl2, temp.cl1)
    }
    canonical.df$cluster <- temp.cl2

  }

  #########################################################
  #------------User defined clustering--------------------#
  #########################################################
  if(method == "User_supplied")
  {
    message("Running User defined clustering function")

    canonical.df$cluster <- user_function(Object@Data.RR)
    clusterUnique <- user_function(Object@DataPathways.RR)

  }


  ############################################
  #---------Adding cluster info--------------#
  ############################################

  canonical.df$mean.RR.cl <- NA
  canonical.df$sum.RR.cl <- NA


  for(cl.i in unique(canonical.df$cluster))
  {
    DF.cl <- canonical.df[canonical.df$cluster == cl.i,]
    idx <- DF.cl$RR_name
    rr.x <- Object@Data.RR[idx, idx]
    rr.x <- data.frame(sapply(rr.x, function(x) as.numeric(as.character(x))))
    rr.x <- as.matrix(rr.x)
    DF.cl$mean.RR.cl <- round(mean(rr.x),digits = 4)
    DF.cl$sum.RR.cl <- round(sum(rr.x),digits = 4)

    canonical.df[canonical.df$cluster == cl.i,] <- DF.cl


  }




  ###############################
  #--------Order clusters-------#
  ###############################
  message("Ordering pathway clusters")
  if(order == "group")
  {
    Object@Data <- list(canonical.df[order(canonical.df$Group,canonical.df$cluster,  decreasing = F),])

    Object@metadata[,"order.group"] <- rep(order, times = nrow(Object@metadata))

  }
  if(order == "cluster")
  {
    Object@Data <- list(canonical.df[order(canonical.df$cluster,  decreasing = F),])

    Object@metadata[,"order.group"] <- rep(order, times = nrow(Object@metadata))
  }
  if(order == "mean.RR")
  {
    Object@Data <- list(canonical.df[order(canonical.df$mean.RR.cl,  decreasing = F),])

    Object@metadata[,"order.group"] <- rep(order, times = nrow(Object@metadata))

  }
  if(order == "Sum.RR")
  {
    Object@Data <- list(canonical.df[order(canonical.df$sum.RR.cl,  decreasing = F),])

    Object@metadata[,"order.group"] <- rep(order, times = nrow(Object@metadata))

  }
  Object@metadata[,"cluster.method"] <- rep(method, times = nrow(Object@metadata))
  Object@Data.RR <- Object@Data.RR[Object@Data[[1]]$RR_name,Object@Data[[1]]$RR_name]

  ##########################
  #--------plot Info-------#
  ##########################

  pal.c <-  c(brewer.pal(n = 8, name ="Accent" ),
              brewer.pal(n = 8, name ="Dark2" ),
              brewer.pal(n = 8, name ="Set3"),
              brewer.pal(n = 8, name ="Set1"),
              brewer.pal(n = 8, name ="Set2"),
              brewer.pal(n = 8, name ="Greens"),
              brewer.pal(n = 8, name ="BrBG"))

  if(Object@metadata$display[1] == "Expanded")
  {
    #Display is expanded meaning groups get marked seperatly
    if(length(unique(as.character(Object@Data[[1]]$Type))) > 1)
    {
      aka2 <- matrix(data = NA,
                     nrow = nrow((Object@Data.RR)),
                     ncol = 3 +nrow(Object@metadata))
      colnames(aka2) <- c("Group", "Cluster", "Type", unique(Object@metadata$Groups))
      aka2[,"Type"] <- as.character(Object@Data[[1]]$Type)
      pal.type <- pal.c[1:length(unique(aka2[,"Type"]))]
      names(pal.type) <- unique(aka2[,"Type"])

      for(groups.i in 1:nrow(Object@metadata))
      {
        aka2[,Object@metadata$Groups[groups.i]] <- Object@Data[[1]][,Object@metadata$Groups[groups.i]]

      }


    }else{
      aka2 <- matrix(data = NA, nrow = nrow((Object@Data.RR)), ncol = 2+ nrow(Object@metadata))
      colnames(aka2) <- c("Group", "Cluster", unique(Object@metadata$Groups))

      for(groups.i in 1:nrow(Object@metadata))
      {
        aka2[,Object@metadata$Groups[groups.i]] <- Object@Data[[1]][,Object@metadata$Groups[groups.i]]

      }
    }
  }else{
    if(length(unique(as.character(Object@Data[[1]]$Type))) > 1)
    {
      aka2 <- matrix(data = NA, nrow = nrow((Object@Data.RR)), ncol = 3)
      colnames(aka2) <- c("Group", "Cluster", "Type")
      aka2[,"Type"] <- as.character(Object@Data[[1]]$Type)
      pal.type <- pal.c[1:length(unique(aka2[,"Type"]))]
      names(pal.type) <- unique(aka2[,"Type"])

    }else{
      aka2 <- matrix(data = NA, nrow = nrow((Object@Data.RR)), ncol = 2)
      colnames(aka2) <- c("Group", "Cluster")
    }
  }
  rownames(aka2) <- Object@Data[[1]]$RR_name

  aka2[,"Group"] <- as.character(Object@Data[[1]]$Groups)
  aka2[,"Cluster"] <- as.character(Object@Data[[1]]$cluster)


  if(!method %in% c("kmeans_group","Kmeans_group","Kmeans_Group", "kmeans group","Kmeans group", "Hierarchical_group","hierarchical_group", "Hierarchical group","hierarchical group", "hierarchical_group"))
  {
    aka2Unique <- obtainDFmetadata(Object@Data[[1]], Object@DataPathways.RR)
    aka2Unique$Cluster <- clusterUnique
  }

  aka2 <- as.data.frame(aka2)
  if(length(unique(Object@Data[[1]]$cluster)) > 32)
  {
    message("Warning: number of cluster is larger than colours supported")
  }

  groups.col <- brewer.pal(n = 8, name ="Set2" )[1:length(unique(aka2[,"Group"]))]
  names(groups.col) <- unique(aka2[,"Group"])
  pal <- pal.c[1:length(unique(aka2[,"Cluster"]))]
  names(pal) <- unique(aka2[,"Cluster"])

  if(!method %in% c("kmeans_group","Kmeans_group","Kmeans_Group", "kmeans group","Kmeans group", "Hierarchical_group","hierarchical_group", "Hierarchical group","hierarchical group", "hierarchical_group"))
  {
    palUnique <- pal.c[1:length(unique(aka2Unique[,"Cluster"]))]
    names(palUnique) <- unique(aka2Unique[,"Cluster"])
  }

  if(Object@metadata$display[1] == "Expanded")
  {

    if(length(unique(as.character(Object@Data[[1]]$Type))) > 1)
    {
      aka3 = list(Group = groups.col,
                  Cluster= pal,
                  Type = pal.type)

      for(groups.i in 1:nrow(Object@metadata))
      {
        vector.x <- c("black","red")
        names(vector.x) <- c(0,1)
        aka3[[Object@metadata$Groups[groups.i]]] <-   vector.x
      }

      names(aka3) <-  c("Group", "Cluster", "Type",  unique(Object@metadata$Groups))

      if(!method %in% c("kmeans_group","Kmeans_group","Kmeans_Group", "kmeans group","Kmeans group", "Hierarchical_group","hierarchical_group", "Hierarchical group","hierarchical group", "hierarchical_group"))
      {
        aka3Unique = list(Group = groups.col,
                          Cluster= palUnique)

        for(groups.i in 1:nrow(Object@metadata))
        {
          vector.x <- c("black","red")
          names(vector.x) <- c(0,1)
          aka3Unique[[Object@metadata$Groups[groups.i]]] <-   vector.x
        }

        names(aka3Unique) <-  c("Group", "Cluster",  unique(Object@metadata$Groups))
      }

    }else{
      aka3 = list(Group = groups.col,
                  Cluster= pal)

      for(groups.i in 1:nrow(Object@metadata))
      {
        vector.x <- c("black","red")
        names(vector.x) <- c(0,1)
        aka3[[Object@metadata$Groups[groups.i]]] <-   vector.x
      }

      names(aka3) <-  c("Group", "Cluster",  unique(Object@metadata$Groups))

      if(!method %in% c("kmeans_group","Kmeans_group","Kmeans_Group", "kmeans group","Kmeans group", "Hierarchical_group","hierarchical_group", "Hierarchical group","hierarchical group", "hierarchical_group"))
      {
        aka3Unique = list(Group = groups.col,
                          Cluster= palUnique)

        for(groups.i in 1:nrow(Object@metadata))
        {
          vector.x <- c("black","red")
          names(vector.x) <- c(0,1)
          aka3Unique[[Object@metadata$Groups[groups.i]]] <-   vector.x
        }

        names(aka3Unique) <-  c("Group", "Cluster",  unique(Object@metadata$Groups))
      }
    }
  }else{

    if(length(unique(as.character(Object@Data[[1]]$Type))) > 1)
    {
      aka3 = list(Group = groups.col,
                  Cluster= pal,
                  Type = pal.type)
      names(aka3) <-  c("Group", "Cluster", "Type")

      if(!method %in% c("kmeans_group","Kmeans_group","Kmeans_Group", "kmeans group","Kmeans group", "Hierarchical_group","hierarchical_group", "Hierarchical group","hierarchical group", "hierarchical_group"))
      {
        aka3Unique = list(Group = groups.col,
                          Cluster= palUnique)

        names(aka3Unique) <-  c("Group", "Cluster")
      }

    }else{
      aka3 = list(Group = groups.col,
                  Cluster= pal)
      names(aka3) <-  c("Group", "Cluster")

      if(!method %in% c("kmeans_group","Kmeans_group","Kmeans_Group", "kmeans group","Kmeans group", "Hierarchical_group","hierarchical_group", "Hierarchical group","hierarchical group", "hierarchical_group"))
      {
        aka3Unique = list(Group = groups.col,
                          Cluster= palUnique)

        names(aka3Unique) <-  c("Group", "Cluster")
      }
    }
  }

  ##################################
  #-----------Return---------------#
  ##################################

  Object@plot$aka2 = aka2
  Object@plot$aka2Unique = aka2Unique
  Object@plot$aka3 = aka3
  Object@plot$aka3Unique = aka3Unique



  message("-----------------------------------------------------------")
  message("[<<<<<             ClusterGeneSets END              >>>>>>]")
  message("[=========================================================]")
  message("[You may want to process HighlightGeneSets next.            ]")
  message("[You may want to plot the results using PlotGeneSets next. ]")


  return(Object)
}
