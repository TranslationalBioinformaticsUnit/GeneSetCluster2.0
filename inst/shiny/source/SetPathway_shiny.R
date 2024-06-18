setGeneric(name="SetPathway_shiny",
           def=function(Object, nPathways="optimal", uniquePathways = FALSE)
           {
             standardGeneric("SetPathway_shiny")
           }
)


setMethod(f="SetPathway_shiny",
          definition=function(Object, nPathways="optimal", uniquePathways = FALSE)
          {
            message("[=============================================================]")
            message("[<<<<               SetPathway START                >>>>>]")
            message("---------------------------------------------------------------")

            #check if cluster independent slot is null
            checkIndependent(Object)
            pal.c <-  getPal()
            groups.col <- brewer.pal(n = 8, name ="Set2" )[1:length(unique(Object@plot$aka2[,"Group"]))]
            names(groups.col) <- unique(Object@plot$aka2[,"Group"])

            if (nPathways == "optimal")
            {
              if(uniquePathways){
                mat_corPathway <- scaleCorMatrix(Object@DataPathways.RR)[Object@cIndependentMethod[[1]][[1]],
                                                                         Object@cIndependentMethod[[1]][[1]]]
                resPathway <- obtainDefCluster(mat_corPathway)
                optimalNumberPathway <- getOptimalNumber(resPathway)
                message(paste0("The optimal number of minimum Pathways for Unique Pathway approach is ", optimalNumberPathway))
                go_id_listPathway <- resPathway[which(lapply(resPathway, function(x) length(x)) >= optimalNumberPathway)]

                aka2Unique <- obtainDFmetadata(Object@Data[[1]], Object@DataPathways.RR)
                #order matrix
                aka2Unique <- aka2Unique[rownames(mat_corPathway),]
                aka2Unique$Cluster <- 0
                for (i in 1:length(go_id_listPathway))
                {
                  aka2Unique[go_id_listPathway[[i]], "Cluster"] <- i
                }

                clusterOrder <- sort(unique(aka2Unique[,"Cluster"]))
                clusterOrder <- c(1:(length(clusterOrder)-1), 0)
                palUnique <- pal.c[1:length(clusterOrder)]
                names(palUnique) <- clusterOrder
                palUnique["0"] <- "white"

                aka3Unique <- list(Group = groups.col,
                                   Cluster = palUnique,
                                   optimalNumberPathway = optimalNumberPathway)

                Object@plot$aka2IndependentUnique <- aka2Unique
                Object@plot$aka3IndependentUnique <- aka3Unique

                message("\nCalculating ORA and Semantic enrichment per cluster for Unique Pathways...")
                #functional_annotation <- PerformORAGOperCluster(Object, uniquePathways = FALSE, clusterIndependent = TRUE)

              }else{
                #scaling the matrix
                mat_cor <- scaleCorMatrix(Object@Data.RR)[Object@cIndependentMethod[[1]][[1]],
                                                          Object@cIndependentMethod[[1]][[1]]]
                res <- obtainDefCluster(mat_cor)
                optimalNumber <- getOptimalNumber(res)
                message(paste0("The optimal number of minimum Gene-sets for Unique Gene-set approach is ", optimalNumber))
                go_id_list <- res[which(lapply(res, function(x) length(x)) >= optimalNumber)]

                aka2 <- matrix(data = NA, nrow = nrow(Object@Data[[1]]), ncol = 2)
                colnames(aka2) <- c("Group", "Cluster")
                rownames(aka2) <- Object@Data[[1]]$RR_name

                aka2[,"Group"] <- as.character(Object@Data[[1]]$Groups)
                #order matrix
                aka2 <- aka2[rownames(mat_cor),]

                for (i in 1:length(go_id_list)) {
                  aka2[which(rownames(aka2) %in% go_id_list[[i]]), "Cluster"] <- i
                }

                aka2 <- as.data.frame(aka2)
                aka2$Cluster[is.na(aka2$Cluster)] <- 0

                if(length(unique(aka2$Cluster)) > 32)
                {
                  message("Warning: number of cluster is larger than colours supported")
                }

                #pal.c <-  getPal()

                pal <- pal.c[1:length(unique(aka2[,"Cluster"]))]
                names(pal) <- unique(aka2[,"Cluster"])
                pal["0"] <- "white"

                aka3 <- list(Group = groups.col,
                             Cluster= pal,
                             optimalNumberPathway = optimalNumber)
                #names(aka3) <-  c("Group", "Cluster")

                Object@plot$aka2Independent <- aka2
                Object@plot$aka3Independent <- aka3

                message("\nCalculating ORA and Semantic enrichment per cluster for Unique Gene-sets...")
                #functional_annotationUnique <- PerformORAGOperCluster(Object, uniquePathways = TRUE, clusterIndependent = TRUE)
              }

            } else {
              optimalNumber <- nPathways
              optimalNumberPathway <- nPathways
              message(paste0("Setting ", nPathways, " as minimum number of elements to consider a cluster...\n"))
            }


            # Object@functionalAnnotIndependent <- list(Geneset=functional_annotation,
            #                                           Pathway=functional_annotationUnique)


            message("--------------------------------------------------------------")
            message("[<<<<               SetPathway ends                >>>>>]")
            message("[============================================================]")

            return(Object)
}
)
