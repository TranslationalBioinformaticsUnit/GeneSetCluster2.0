LoadGeneSets_shiny <- function(file_location = canonical.files,
                         groupnames= c("MClust_1", "MClust_26", "MClust_3457"),
                         P.cutoff = 1.3,
                         Mol.cutoff = 5,
                         Source = "IPA",
                         Great.Background=F,
                         type = "Canonical_Pathways",
                         topranks = "",
                         structure = "SYMBOL",
                         Organism = "org.Mm.eg.db",#Human genes for converting genes to the right structure
                         seperator = ",")
{

  message("[=========================================================]")
  message("[<<<<            LoadGeneSets_shiny START                  >>>>>]")
  message("-----------------------------------------------------------")

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
        Canonical.x[can.i,"MoleculesCount"]<- length(as.vector(strsplit2(as.character(Canonical.x[can.i,"Molecules"]), split=",")))
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
        list.canonical[[i]] <- list.canonical[[i]][list.canonical[[i]][,grep("P.value", colnames(list.canonical[[i]]), ignore.case = T)] > P.cutoff ,]

      }else{
        list.canonical[[i]] <- list.canonical[[i]][list.canonical[[i]][,grep("P.value", colnames(list.canonical[[i]]), ignore.case = T)] < P.cutoff ,]

      }
      list.canonical[[i]] <- list.canonical[[i]][list.canonical[[i]]$MoleculesCount >= Mol.cutoff,]

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
      Data$Pathways <- Pathways
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

  }

  ######################################
  ##---------GREAT--------------------##
  ######################################
  if(Source == "Great")
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
      #check if it is from R or the website
      if(is.null(Canonical.x$X..Ontology)){
        #is from R
        message("File input from GREAT R package. Please, all GO supposed to be Biological Process")
        Canonical.x <- read.csv(file = file_location[i], header = T, sep="\t")
        #rename colnames
        colnames(Canonical.x)<-c("ID","description","genome_fraction","ObsRegions","fold_enrichment","p_value","HyperFdrQ","mean_tss_dist","ObsGenes","TotalGenes","fold_enrichment_hyper",	"p_value_hyper","p_adjust_hyper","Genes")
        Canonical.x$X..Ontology<-rep("GO Biological Process",nrow(Canonical.x))
      }else{
        Canonical.x <- Canonical.x[1:(nrow(Canonical.x)-21),]
      }
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
        message(paste(groupnames[i], " has 0 pathways passing cutoff, suggest trying lower cutoff", sep=""))
        return(groupnames[i])
        stop()
      }else{

        #split per type
        Great.Data.ls <- list()
        for(list.ii in 1:length(split.idx))
        {
          Canonical.x <- list.canonical[[i]][list.canonical[[i]]$X..Ontology == split.idx[list.ii],]

          if(!topranks == "")
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

    Pdata[,"Groupnames"] <- strsplit2(x = names(list.canonical.f), split = "_")[,1]
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
    metadata[,"type"] <-  strsplit2(x = names(list.canonical.f), split = "_")[,2]
    metadata[,"structure"] <- rep(structure, times = nrow(Pdata))
    metadata[,"organism"] <- rep(Organism, times = nrow(Pdata))
    metadata[,"Groups"] <-  strsplit2(x = names(list.canonical.f), split = "_")[,1]
    metadata[,"loaded"] <- rep("Auto", times = nrow(Pdata))
    metadata[,"Data"] <- rep("List", times = nrow(Pdata))
    metadata[,"seperator"] <- rep(seperator, times = nrow(Pdata))
    metadata[,"display"] <- rep("Condensed", times = nrow(Pdata))
    metadata[,"mol.signature"] <- rep("All", times = nrow(Pdata))


  }

  ######################################
  ##----------GSEA--------------------##
  ######################################
  if(Source == "GSEA")
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

      message("File input from GSEA. Please, all GO supposed to be Biological Process")
      if((length(grep("xls", file_location[i]))>0)||(length(grep("xlsx", file_location[i]))>0)){
        Canonical.x <-read_excel(file_location[i], col_names=TRUE)
      }else{
        Canonical.x <- read.csv(file = file_location[i], header = T, sep="\t")
        if (length(colnames(Canonical.x))==1){
          Canonical.x <- read.csv(file = file_location[i], header = T, sep=",")
        }
        if (length(colnames(Canonical.x))==1){
          Canonical.x <- read.csv(file = file_location[i], header = T, sep=";")
        }
      }
      Canonical.x$X..Ontology<-rep("GO Biological Process",nrow(Canonical.x))
      list.canonical[[i]] <- as.data.frame(Canonical.x)
      names(list.canonical)[i] <- groupnames[i]
    }

    ###########################################
    ##---------filter for cutoff-------------##
    ###########################################
    list.canonical.f <- list()
    file.locations.2 <- vector()
    message("Loading Splitting GSEA types passing cutoff")

    possibleSeparators<-c(",","/",";")

    for(i in 1:length(file_location))
    {
      if(is.null(list.canonical[[i]]$enrichmentScore)){
        counts<-sapply(possibleSeparators, str_count, string=list.canonical[[i]]$geneID[1])
        separator<-possibleSeparators[which(counts==max(counts))]
        #enrichGO output data
        #list.canonical[[i]] <- list.canonical[[i]][as.numeric(as.character(list.canonical[[i]]$Count)) >= Mol.cutoff,]
        list.canonical[[i]] <- list.canonical[[i]][as.numeric(str_count(as.character(list.canonical[[i]]$geneID), separator)+1) >= Mol.cutoff,]
      }else{
        counts<-sapply(possibleSeparators, str_count, string=list.canonical[[i]]$core_enrichment[1])
        separator<-possibleSeparators[which(counts==max(counts))]
        #gseGO output data
        #list.canonical[[i]] <- list.canonical[[i]][as.numeric(as.character(list.canonical[[i]]$setSize)) >= Mol.cutoff,]
        list.canonical[[i]] <- list.canonical[[i]][as.numeric(str_count(as.character(list.canonical[[i]]$core_enrichment), separator)+1) >= Mol.cutoff,]
      }
      list.canonical[[i]] <- list.canonical[[i]][as.numeric(as.character(list.canonical[[i]]$p.adjust)) <= P.cutoff,]
      split.idx <- names(table(list.canonical[[i]]$X..Ontology))[table(list.canonical[[i]]$X..Ontology) > 0]

      if(length(split.idx) == 0)
      {
        message(paste(groupnames[i], " has 0 pathways passing cutoff, suggest trying lower cutoff", sep=""))
        return(groupnames[i])
        stop()
      }else{

        #split per type
        Great.Data.ls <- list()
        for(list.ii in 1:length(split.idx))
        {
          Canonical.x <- list.canonical[[i]][list.canonical[[i]]$X..Ontology == split.idx[list.ii],]

          if(!topranks == "")
          {
            if(nrow(Canonical.x) > topranks){
              Canonical.x <- Canonical.x[1:topranks,]
            }
          }

          Great.Data <- matrix(NA, nrow=nrow(Canonical.x), ncol = 7)
          colnames(Great.Data) <- c("Pathways", "Molecules", "Groups","Type", "Pval", "Ratio", "MoleculesCount")
          Great.Data <- as.data.frame(Great.Data)


          if(is.null(list.canonical[[i]]$enrichmentScore)){
            #enrichGO output data
            Great.Data$Molecules <- as.character(Canonical.x$geneID)
            Great.Data$Ratio <- as.character(Canonical.x$GeneRatio)
            #Great.Data$MoleculesCount <- as.character(Canonical.x$Count)
            Great.Data$MoleculesCount <- str_count(as.character(Canonical.x$geneID), separator)+1

          }else{
            #gseGO output data
            Great.Data$Molecules <- as.character(Canonical.x$core_enrichment)
            #Great.Data$Ratio <- as.numeric(as.character(Canonical.x$rank))/as.numeric(as.character(Canonical.x$setSize))
            #Great.Data$MoleculesCount <- as.character(Canonical.x$setSize)
            Great.Data$MoleculesCount <- str_count(as.character(Canonical.x$core_enrichment), separator)+1
            Great.Data$Ratio <- as.numeric(as.character(Canonical.x$rank))/Great.Data$MoleculesCount
          }

          Great.Data$Pathways <- as.character(Canonical.x$ID)
          Great.Data$Groups <- groupnames[i]
          Great.Data$Type <- split.idx[list.ii]
          Great.Data$Pval <- as.character(Canonical.x$p.adjust)

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

    Pdata[,"Groupnames"] <- strsplit2(x = names(list.canonical.f), split = "_")[,1]
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
    metadata[,"type"] <-  strsplit2(x = names(list.canonical.f), split = "_")[,2]
    metadata[,"structure"] <- rep(structure, times = nrow(Pdata))
    metadata[,"organism"] <- rep(Organism, times = nrow(Pdata))
    metadata[,"Groups"] <-  strsplit2(x = names(list.canonical.f), split = "_")[,1]
    metadata[,"loaded"] <- rep("Auto", times = nrow(Pdata))
    metadata[,"Data"] <- rep("List", times = nrow(Pdata))
    metadata[,"seperator"] <- rep(seperator, times = nrow(Pdata))
    metadata[,"display"] <- rep("Condensed", times = nrow(Pdata))
    metadata[,"mol.signature"] <- rep("All", times = nrow(Pdata))

  }
  ######################################
  ##---------Combine data-------------##
  ######################################
  #message("Loading IPA sheets")

  Object <-  PathwayObject(Data = list.canonical.f,
                           PData = Pdata,
                           metadata = metadata,
                           Data.RR = data.frame(),
                           plot = list(aka2 = data.frame(),
                                       aka3 = data.frame()))


  message("-----------------------------------------------------------")
  message("[<<<<<               LoadGeneSets_shiny END               >>>>>>]")
  message("[=========================================================]")
  message("[You may want to process CombineGeneSets next.            ]")
  message("[or merge objects using MergeObjects                     ]")
  message("[or select certain types from objects using ManageGeneSets]")

  return(Object)
}
