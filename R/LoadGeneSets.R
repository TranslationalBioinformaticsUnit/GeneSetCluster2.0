#' @title LoadGeneSets
#'
#' @description
#' Automatic loader for gene sets from the GSEA, GREAT or IPA web tools. The input can be either the file location
#' obtained from the online tools, or the GREAT or GSEA object obtained from the R tools. Also a template
#' data frame can be used to import Gene-Set analysis from other tools.
#'
#' @import readxl
#' @import stringr
#' @import rGREAT
#' @importFrom utils read.csv write.table
#' @importFrom limma strsplit2
#'
#' @param file_location A list of GREAT, gseaResult or data.frame objects or location string in a vector of files.
#' @param groupnames A vector with group names of the different gene set experiments
#' @param P.cutoff numeric Pvalue cutoff for filtering.
#' @param Mol.cutoff numeric value for minimum number of molecules.
#' @param Source Tool used to generate gene sets.
#' @param Great.Background If the Great tool was used, did the user supply a background.
#' @param type For IPA data if "Canonical_Pathways" or "Functional_annotations" were supplied.
#' @param topranks numeric with the number of ranks per group to be loaded, usefull when there is a lot of data.
#' @param structure The structure of the genes. is it SYMBOLS, ENSEMBL, NCBI etc. Used for converting when there is mutiple structure in the object.
#' @param Organism the package name for the human or mouse data, used for converting the gene structure. name of the package, currently org.Hs.eg.db and org.Mm.eg.db supported.
#' @param seperator A character used within in the string to separate genes
#'
#' @return PathwayObject
#'
#' @examples
#' # Loading from GREAT webtool results files:
#' Great.files <- c(system.file("extdata", "MM10.GREAT.KO.uGvsMac.bed.tsv",
#'                              package = "GeneSetCluster"),
#'                  system.file("extdata", "MM10.GREAT.WT.uGvsMac.bed.tsv",
#'                              package = "GeneSetCluster"))
#'
#'
#' Great.Object1 <- LoadGeneSets(file_location = Great.files,
#'                               groupnames= c("KO", "WT"),
#'                               P.cutoff = 0.05,
#'                               Mol.cutoff = 10,
#'                               Source = "Great",
#'                               Great.Background = FALSE,
#'                               type = "Canonical_Pathways",
#'                               topranks = 20,
#'                               structure = "SYMBOL",
#'                               Organism = "org.Mm.eg.db",
#'                               seperator = ",")
#'
#' # Loading from IPA web tool results files:
#' IPA.files <- c(system.file("extdata",
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
#'                          seperator = ",",
#'                          Organism = "org.Mm.eg.db")
#'
#' # Loading from GSEA results files:
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
#'
#'
#' @export
LoadGeneSets <- function(file_location,
                         groupnames,
                         P.cutoff = 0.05,
                         Mol.cutoff = 10,
                         Source,
                         Great.Background=F,
                         type = NA,
                         topranks = NA,
                         structure,
                         Organism,
                         seperator)
{

  message("[=========================================================]")
  message("[<<<<            LoadGeneSets START                  >>>>>]")
  message("-----------------------------------------------------------")


  if (!is.character(file_location))
  {
    supported_method <- c("GREAT", "GSEA", "TEMPLATE")
    if (!toupper(Source) %in% supported_method)
    {
      message(paste0(Source, " method is not supported. GREAT, GSEA and TEMPLATE are supported."))
      message("You may use the template structure.\n")
      printTemplate()
      message("\n")
      stop()
    }

    if (!length(groupnames) == length(file_location))
    {
      stop("The number of group names should be the same as inotrduced objects.")
    }

    if (toupper(Source) == "GREAT")
    {
      pathways_in <- c()
      molecules_in <- c()
      pvalues_in <- c()
      groups_in <- c()
      enirchment_in <- c()

      message("Loading GREAT objects...")
      for (i in 1:length(file_location))
      {
        df_GREAT <- obtainGREATtable(file_location[[i]])

        #filter
        df_GREAT <- df_GREAT[which(df_GREAT$p_adjust < P.cutoff),]
        df_GREAT <- df_GREAT[which(unlist(lapply(strsplit(df_GREAT$Genes, "/"), function(x) length(x)>=Mol.cutoff))),]

        pathways_in <- c(pathways_in, df_GREAT$id)
        molecules_in <- c(molecules_in, df_GREAT$Genes)
        pvalues_in <- c(pvalues_in, df_GREAT$p_adjust)
        groups_in <- c(groups_in, rep(groupnames[i], times = nrow(df_GREAT)))
        enirchment_in <- c(enirchment_in, df_GREAT$fold_enrichment_hyper)

      }

      Object <- ObjectCreator(Pathways = pathways_in,
                              Molecules = molecules_in,
                              Groups = groups_in,
                              Pvalues = pvalues_in,
                              enrichmentScore = enirchment_in,
                              structure = "ENTREZID", Type = "", sep = "/",
                              Source = "GREAT", organism = Organism)


    } else if (toupper(Source) == "GSEA")
    {

      if (!class(file_location[[1]])[1] == "gseaResult")
      {
        message("The introduced object is not gseaResult object.")
        message("You may use the template structure.\n")
        printTemplate()
        message("\n")
        stop()
      }

      pathways_in <- c()
      molecules_in <- c()
      pvalues_in <- c()
      groups_in <- c()
      enirchment_in <- c()

      message("Loading GSEA objects...")
      for (i in 1:length(file_location))
      {
        df_GSEA<- file_location[[i]]@result

        #filter
        df_GSEA <- df_GSEA[which(df_GSEA$p.adjust < P.cutoff),]
        df_GSEA <- df_GSEA[which(unlist(lapply(strsplit(df_GSEA$core_enrichment, seperator),
                                               function(x) length(x)>=Mol.cutoff))),]

        pathways_in <- c(pathways_in, df_GSEA$Groups)
        molecules_in <- c(molecules_in, df_GSEA$core_enrichment)
        pvalues_in <- c(pvalues_in, df_GSEA$p.adjust)
        groups_in <- c(groups_in, rep(groupnames[i], times = nrow(df_GSEA)))
        enirchment_in <- c(enirchment_in, df_GSEA$enrichmentScore)

      }

      Object <- ObjectCreator(Pathways = pathways_in,
                              Molecules = molecules_in,
                              Groups = groups_in,
                              Pvalues = pvalues_in,
                              enrichmentScore = enirchment_in,
                              structure = structure, Type = "", sep = seperator,
                              Source = "GSEA", organism = Organism)

    } else {
      message("TEMPLATE selected.")
      min_colnames <- c("Pathways", "Pvalues", "EnrichmentScore", "Molecules")

      if (!class(file_location[[1]]) == "data.frame")
      {
        message("The introduced object is not data.frame")
        message("You need to use the template structure.\n")
        printTemplate()
        message("\n")
        stop()
      } else if (!all(min_colnames %in% colnames(file_location[[1]])))
      {
        printTemplate()
        stop()
      }

      pathways_in <- c()
      molecules_in <- c()
      pvalues_in <- c()
      groups_in <- c()
      enirchment_in <- c()

      message("Loading TEMPLATE data frames...")
      for (i in 1:length(file_location))
      {
        df_GSEA<- file_location[[i]]

        #filter
        df_GSEA <- df_GSEA[which(df_GSEA$Pvalues < P.cutoff),]
        df_GSEA <- df_GSEA[which(unlist(lapply(strsplit(df_GSEA$Molecules, seperator), function(x) length(x)>=Mol.cutoff))),]

        pathways_in <- c(pathways_in, df_GSEA$Pathways)
        molecules_in <- c(molecules_in, df_GSEA$Molecules)
        groups_in <- c(groups_in, rep(groupnames[i], times = nrow(df_GSEA)))
        pvalues_in <- c(pvalues_in, df_GSEA$Pvalues)
        enirchment_in <- c(enirchment_in, df_GSEA$EnrichmentScore)

      }

      Object <- ObjectCreator(Pathways = pathways_in,
                              Molecules = molecules_in,
                              Groups = groups_in,
                              Pvalues = pvalues_in,
                              enrichmentScore = enirchment_in,
                              structure = structure, Type = "", sep = seperator,
                              Source = "GSEA", organism = Organism)



    }


  } else {


    for(load.i in 1:length(file_location))
    {
      message("Loading data from ", file_location[load.i])

    }
    ###########################
    ##---------IPA-----------##
    ###########################


    if(Source == "IPA")
    {
      ####################################
      ##---------Load Excel-------------##
      ####################################
      if(!length(file_location) == length(groupnames))
      {
        message("Names length dont match")
      }
      list.canonical <- list()
      for(i in 1:length(file_location))
      {

        Canonical.x <- read_excel(path = file_location[i],skip=1, sheet = 1)
        Canonical.x <- as.data.frame(Canonical.x)
        Canonical.x$MoleculesCount <- NA
        for(can.i in 1:nrow(Canonical.x))
        {
          Canonical.x[can.i,"MoleculesCount"]<- length(as.vector(limma::strsplit2(as.character(Canonical.x[can.i,"Molecules"]), split=seperator)))
        }
        list.canonical[[i]] <- as.data.frame(Canonical.x)
        names(list.canonical)[i] <- groupnames[i]
      }
      ###########################################
      ##---------filter for cutoff-------------##
      ###########################################

      for(i in 1:length(file_location))
      {
        if(P.cutoff > 1){#meaning that it is a -log10 Pvalue if smaller than 1 means it is a untransformed Pvalue
          if (nrow(list.canonical[[i]][list.canonical[[i]][,grep("P.value", colnames(list.canonical[[i]]), ignore.case = T)] > P.cutoff ,]) == 0)
          {
            message(paste0("There are not pathways that pass the filter of P value cut off for ", file_location[[i]], " file."))
            list.canonical[[i]] <- list.canonical[[i]]
          } else {
            list.canonical[[i]] <- list.canonical[[i]][list.canonical[[i]][,grep("P.value", colnames(list.canonical[[i]]), ignore.case = T)] > P.cutoff ,]
          }
        } else {
          if (nrow(list.canonical[[i]][list.canonical[[i]][,grep("P.value", colnames(list.canonical[[i]]), ignore.case = T)] > P.cutoff ,]) == 0)
          {
            message(paste0("There are not pathways that pass the filter of P value cut off for ", file_location[[i]], " file."))
            list.canonical[[i]] <- list.canonical[[i]]
          } else {
            list.canonical[[i]] <- list.canonical[[i]][list.canonical[[i]][,grep("P.value", colnames(list.canonical[[i]]), ignore.case = T)] < P.cutoff ,]
          }
        }
        if (nrow(list.canonical[[i]][list.canonical[[i]]$MoleculesCount >= Mol.cutoff,]) == 0)
        {
          message(paste0("There are not pathways that pass the filter of Mol.cutoff.\nSetting Mol.cutoff=0 for ", file_location[[i]], " file."))
          list.canonical[[i]] <- list.canonical[[i]][list.canonical[[i]]$MoleculesCount >= 0,]
        } else {
          list.canonical[[i]] <- list.canonical[[i]][list.canonical[[i]]$MoleculesCount >= Mol.cutoff,]
        }
      }
      list.canonical.f <- list()
      if(type == "Canonical_Pathways")
      {
        message("Loading IPA Canonical_Pathways")
      }
      if(type == "Functional_annotations")
      {
        message("Loading IPA Functional_annotations")
      }
      for(i in 1:length(file_location))
      {

        Data <- matrix(NA, nrow=nrow(list.canonical[[i]]), ncol = 7)
        colnames(Data) <- c("Pathways", "Molecules", "Groups","Type", "Pval", "Ratio", "MoleculesCount")
        Data <- as.data.frame(Data)

        if(type == "Canonical_Pathways")
        {
          Pathways <- list.canonical[[i]][,1]
          pval <- list.canonical[[i]]$X..log.p.value.
          ratio <- list.canonical[[i]]$Ratio
        }
        if(type == "Functional_annotations")
        {
          Pathways <- paste(list.canonical[[i]][,1], list.canonical[[i]][,2], sep="_")
          pval <- list.canonical[[i]]$p.Value
          ratio <- rep(NA, times = nrow(list.canonical[[i]]))
        }
        Data$Type <- type
        Data$Pathways <- list.canonical[[i]]$`Ingenuity Canonical Pathways`
        Data$Molecules <- list.canonical[[i]]$Molecules
        Data$Groups <- rep(groupnames[i], times =nrow(Data))
        Data$Pval <- pval
        Data$Ratio <- ratio
        Data$MoleculesCount <- list.canonical[[i]]$MoleculesCount
        list.canonical.f[[i]] <- Data
      }
      names(list.canonical.f) <- names(list.canonical)

      ####################################
      ##---------Pheno data-------------##
      ####################################
      Pdata <- matrix(NA,nrow=length(groupnames), ncol = 3)
      colnames(Pdata) <- c("Groupnames", "Length", "file_location")
      rownames(Pdata) <- paste("Experiment",1:length(file_location), sep="_" )
      Pdata <- as.data.frame(Pdata)

      Pdata[,"Groupnames"] <- groupnames
      Pdata[,"Length"] <- sapply(list.canonical.f, nrow)
      Pdata[,"file_location"] <- file_location

      ###################################
      ##---------Meta data-------------##
      ###################################
      metadata <- matrix(NA, nrow = length(groupnames), ncol = 13)
      colnames(metadata) <- c("source", "type", "structure", "organism", "Groups", "seperator", "Data",
                              "cluster.method", "highlight", "order.group", "loaded", "display", "mol.signature")
      rownames(metadata) <- paste("Experiment",1:length(groupnames), sep="_" )
      metadata <- as.data.frame(metadata)

      metadata[,"source"] <- rep(Source, times = nrow(Pdata))
      metadata[,"type"] <- rep(type, times = nrow(Pdata))
      metadata[,"structure"] <- rep(structure, times = nrow(Pdata))
      metadata[,"organism"] <- rep(Organism, times = nrow(Pdata))
      metadata[,"Groups"] <- Pdata[,"Groupnames"]
      metadata[,"loaded"] <- rep("Auto", times = nrow(Pdata))
      metadata[,"Data"] <- rep("List", times = nrow(Pdata))
      metadata[,"seperator"] <- rep(seperator, times = nrow(Pdata))
      metadata[,"display"] <- rep("Condensed", times = nrow(metadata))
      metadata[,"mol.signature"] <- rep("All", times = nrow(metadata))

      Object <-  PathwayObject(Data = list.canonical.f,
                               PData = Pdata,
                               metadata = metadata,
                               Data.RR = data.frame(),
                               plot = list(aka2 = data.frame(),
                                           aka3 = data.frame()),
                               dfTissue = data.frame())
    }

    ######################################
    ##---------GREAT--------------------##
    ######################################
    if(toupper(Source) == "GREAT")
    {
      #Split great up in to the different experiment types and give them different types
      ######################################
      ##---------Load TSV------------##
      ######################################
      if(!length(file_location) == length(groupnames))
      {
        message("Names length dont match")
      }
      list.canonical <- list()
      for(i in 1:length(file_location))
      {

        Canonical.x <- read.csv(file = file_location[i], header = T,skip=3, sep="\t")
        Canonical.x <- Canonical.x[1:(nrow(Canonical.x)-21),]
        if(Great.Background == T){
          colnames(Canonical.x)[grep("FgRegionsHit", colnames(Canonical.x))] <- "ObsRegions"
          colnames(Canonical.x)[grep("NumFgGenesHit", colnames(Canonical.x))] <- "ObsGenes"

        }
        Canonical.x <- Canonical.x[!(is.na(Canonical.x$ObsRegions)),]
        list.canonical[[i]] <- as.data.frame(Canonical.x)
        names(list.canonical)[i] <- groupnames[i]
      }

      ###########################################
      ##---------filter for cutoff-------------##
      ###########################################
      list.canonical.f <- list()
      file.locations.2 <- vector()
      message("Loading Splitting Great types passing cutoff")

      for(i in 1:length(file_location))
      {
        list.canonical[[i]] <- list.canonical[[i]][as.numeric(as.character(list.canonical[[i]]$HyperFdrQ)) <= P.cutoff,]
        list.canonical[[i]] <- list.canonical[[i]][as.numeric(as.character(list.canonical[[i]]$ObsGenes)) >= Mol.cutoff,]
        split.idx <- names(table(list.canonical[[i]]$X..Ontology))[table(list.canonical[[i]]$X..Ontology) > 0]

        if(length(split.idx) == 0)
        {
          message(paste("groupnames[i]", "has 0 pathways passing cutoff, suggest trying lower cutoff", sep=""))
          stop()
        }else{

          #split per type
          Great.Data.ls <- list()
          for(list.ii in 1:length(split.idx))
          {
            Canonical.x <- list.canonical[[i]][list.canonical[[i]]$X..Ontology == split.idx[list.ii],]

            if (!is.character(topranks) == "" & !is.na(topranks) & is.numeric(topranks))
            {
              if(nrow(Canonical.x) > topranks){
                Canonical.x <- Canonical.x[1:topranks,]
              }
            }

            Great.Data <- matrix(NA, nrow=nrow(Canonical.x), ncol = 7)
            colnames(Great.Data) <- c("Pathways", "Molecules", "Groups","Type", "Pval", "Ratio", "MoleculesCount")
            Great.Data <- as.data.frame(Great.Data)

            Great.Data$Pathways <- as.character(Canonical.x$ID)
            if(Great.Background == T)
            {
              Great.Data$Molecules <- as.character(Canonical.x$FgGeneNames)
            }else{
              Great.Data$Molecules <- as.character(Canonical.x$Genes)

            }
            Great.Data$Groups <- groupnames[i]
            Great.Data$Type <- split.idx[list.ii]
            Great.Data$Pval <- as.character(Canonical.x$HyperFdrQ)
            Great.Data$Ratio <- as.numeric(as.character(Canonical.x$ObsGenes))/as.numeric(as.character(Canonical.x$TotalGenes))
            Great.Data$MoleculesCount <- as.character(Canonical.x$ObsGenes)

            Great.Data.ls[[list.ii]] <- Great.Data
            names(Great.Data.ls)[list.ii] <- paste(groupnames[i], split.idx[list.ii],sep="_")

          }
          file.locations.1 <- rep(file_location[i], times = length(split.idx))
          file.locations.2 <- c(file.locations.2, file.locations.1)
          list.canonical.f <- do.call(c, list(list.canonical.f, Great.Data.ls))
        }
      }
      ####################################
      ##---------Pheno data-------------##
      ####################################
      Pdata <- matrix(NA,nrow=length(list.canonical.f), ncol = 3)
      colnames(Pdata) <- c("Groupnames", "Length", "file_location")
      rownames(Pdata) <- paste("Experiment",1:length(list.canonical.f), sep="_" )
      Pdata <- as.data.frame(Pdata)

      Pdata[,"Groupnames"] <- limma::strsplit2(x = names(list.canonical.f), split = "_")[,1]
      Pdata[,"Length"] <- sapply(list.canonical.f, nrow)
      Pdata[,"file_location"] <- file.locations.2

      ###################################
      ##---------Meta data-------------##
      ###################################
      metadata <- matrix(NA, nrow = length(list.canonical.f), ncol = 13)
      colnames(metadata) <- c("source", "type", "structure", "organism", "Groups", "seperator", "Data",
                              "cluster.method", "highlight", "order.group", "loaded", "display", "mol.signature")
      rownames(metadata) <- paste("Experiment",1:length(list.canonical.f), sep="_" )
      metadata <- as.data.frame(metadata)

      metadata[,"source"] <- rep(Source, times = nrow(Pdata))
      metadata[,"type"] <-  limma::strsplit2(x = names(list.canonical.f), split = "_")[,2]
      metadata[,"structure"] <- rep(structure, times = nrow(Pdata))
      metadata[,"organism"] <- rep(Organism, times = nrow(Pdata))
      metadata[,"Groups"] <-  limma::strsplit2(x = names(list.canonical.f), split = "_")[,1]
      metadata[,"loaded"] <- rep("Auto", times = nrow(Pdata))
      metadata[,"Data"] <- rep("List", times = nrow(Pdata))
      metadata[,"seperator"] <- rep(seperator, times = nrow(Pdata))
      metadata[,"display"] <- rep("Condensed", times = nrow(Pdata))
      metadata[,"mol.signature"] <- rep("All", times = nrow(Pdata))

      Object <-  PathwayObject(Data = list.canonical.f,
                               PData = Pdata,
                               metadata = metadata,
                               Data.RR = data.frame(),
                               plot = list(aka2 = data.frame(),
                                           aka3 = data.frame()),
                               dfTissue = data.frame())
    }
    ######################################
    ##---------Combine data-------------##
    ######################################
    #message("Loading IPA sheets")

    if (Source == "GSEA")
    {
      columnNames <- c("ID", "Description", "setSize", "enrichmentScore", "NES", "pvalue", "p.adjust", "qvalue", "rank", "leading_edge", "core_enrichment")
      if (str_detect(file_location[1], ".csv"))
      {

        # slots
        Pathways <- c()
        Molecules <- c()
        Groups <- c()
        Pvalues <- c()
        enrichmentScore <- c()

        for (i in 1:length(file_location))
        {
          data <- read.csv(file_location[i], sep = guessDelim(file_location[i]))
          if (!all(c("ID", "core_enrichment", "p.adjust", "NES") %in% colnames(data)))
          {
            message("The column names of the file are not correct. Please make sure to have the followings column names:")
            message(paste0(columnNames, sep="  "))
            stop()
          }
          # filtering
          data <- data[which(data[,"p.adjust"] < P.cutoff),] # padjust
          keep <- apply(data, 1, function(x) length(limma::strsplit2(x["core_enrichment"], split=seperator))) # number of molecules/genes
          data <- data[which(keep >= Mol.cutoff),]

          #Object slots
          Pathways <- c(Pathways, data$ID)
          Molecules <- c(Molecules, data$core_enrichment)
          Groups <- c(Groups, rep(groupnames[i], times=nrow(data)))
          Pvalues <- c(Pvalues, data$p.adjust)
          enrichmentScore <- c(enrichmentScore, data$NES)
        }
      }

    Object <- ObjectCreator(Pathways = Pathways,
                            Molecules = Molecules,
                            Groups = Groups,
                            Pvalues = Pvalues,
                            enrichmentScore = enrichmentScore,
                            structure = structure, Type = "", sep = seperator,
                            Source = "GSEA", organism = Organism)
    }
  }


  message("-----------------------------------------------------------")
  message("[<<<<<               LoadGeneSets END               >>>>>>]")
  message("[=========================================================]")
  message("[You may want to process CombineGeneSets next.            ]")
  message("[or merge objects using MergeObjects                     ]")
  message("[or select certain types from objects using ManageGeneSets]")

  return(Object)
}
