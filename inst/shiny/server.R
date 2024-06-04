#Load R scripts
source("source/PlotTree.R")
source("source/PlotGeneSets.R")
source("source/PathwayObject.R")
source("source/LoadGeneSets_shiny.R")
source("source/HighlightGeneSets_shiny.R")
source("source/TissueExpressionPerGeneSet_shiny.R")
source("source/PlotTissueExpression_shiny.R")
source("source/GenesPerGeneSet_shiny.R")
source("source/OptimalGeneSets_shiny.R")
source("source/internalFunctions.R")
source("source/CombineGeneSets.R")
source("source/ManageGeneSets.R")
source("source/ClusterGeneSets.R")
source("source/BreakUpCluster_shiny.R")
source("source/ClusterIndependentGeneSet_shiny.R")
source("source/PlotPathwayCluster_shiny.R")
source("source/ShowPathwayCluster_shiny.R")
source("source/SetPathway_shiny.R")

###########
#VERSION 28
###########

load("databases/hpoDatabase.RData")
load("databases/mpDatabase.RData")


underlineCellsCols <- function(rows, cols){
  stopifnot(length(rows) == length(cols))
  c(
    "function(row, data, num, index){",
    sprintf("  var rows = [%s];", paste0(rows-1, collapse = ",")),
    sprintf("  var cols = [%s];", paste0(cols, collapse = ",")),
    "  for(var i = 0; i < rows.length; ++i){",
    "      $('td:eq(' + cols[i] + ')', row)",
    "        .css({'text-decoration': 'underline'});",
    "  }",
    "}"
  )
}


getcharacter<-function(list, seperator){
  nseperator<-str_count(list,seperator)
  return(nseperator)
}


calculateORA<-function(object, cluster, uniquePathways){
  if(uniquePathways){
    clustersvector<-object@plot$aka2Unique$Cluster
  }else{
    clustersvector<-object@plot$aka2$Cluster
  }
  #genes<-as.vector(unlist(sapply(strsplit(object@Data[[1]][,"Molecules"][object@plot$aka2$Cluster==cluster],object@metadata$seperator[1]),unique)))
  genes<-as.vector(unlist(sapply(strsplit(object@Data[[1]][,"Molecules"][clustersvector==cluster],object@metadata$seperator[1]),unique)))
  genes<-unique(gsub(" ","", genes))
  ora<-enrichGO(gene = genes, OrgDb = object@metadata[1,"organism"], keyType = object@metadata[1,"structure"], ont="BP", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable=TRUE)
  if (dim(ora)[1]>20) {
    totalrow=20
  }else{
    totalrow=dim(ora)[1]
  }
  oratop5<-data.frame(ora[1:totalrow,-c(8)], rep(paste("Cluster_",cluster, sep=""), nrow(ora[1:totalrow,])))
  rm(ora, genes)
  oratop5[,c(5,6,7)]<-round(oratop5[,c(5,6,7)],3)
  colnames(oratop5)[9]<-"Cluster"
  return(oratop5)
}

calculateKeywordsORA<-function(ora){
  first_rows <- ora %>%
    group_by(Cluster) %>%
    slice(1) %>%
    ungroup() %>% select(2)

  return(as.list(first_rows$Description))
}

topORA<-function(ora, cluster, top){
  return(head(ora[ora$Cluster==paste("Cluster_",cluster,sep=""),],top))
}

calculateORAindependent<-function(object, cluster, independent_info, uniquePathways){
  if(uniquePathways){
    if (checkGO(object)){
      pathway_cluster<-object@Data[[1]][,"Pathways"] %in% independent_info$ID[independent_info$Cluster==cluster]
    }else{
      pathway_cluster<-object@Data[[1]][,"Pathways"] %in% independent_info$Pathways[independent_info$Cluster==cluster]
    }
  }else{
    if (checkGO(object)){
      pathway_cluster<-object@Data[[1]][,"RR_name"] %in% paste(independent_info$Group[independent_info$Cluster==cluster],independent_info$ID[independent_info$Cluster==cluster],sep="_")
    }else{
      pathway_cluster<-object@Data[[1]][,"RR_name"] %in% paste(independent_info$Group[independent_info$Cluster==cluster],independent_info$Pathways[independent_info$Cluster==cluster],sep="_")
    }
  }
  genes<-as.vector(unlist(sapply(strsplit(object@Data[[1]][,"Molecules"][pathway_cluster],object@metadata$seperator[1]),unique)))
  genes<-unique(gsub(" ","", genes))
  ora<-enrichGO(gene = genes, OrgDb = object@metadata[1,"organism"], keyType = object@metadata[1,"structure"], ont="BP", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable=TRUE)
  if (dim(ora)[1]>10) {
    totalrow=10
  }else{
    totalrow=dim(ora)[1]
  }
  oratop10<-data.frame(ora[1:totalrow,-c(8)], rep(paste(cluster, sep=""), nrow(ora[1:totalrow,])))
  rm(ora, genes)
  oratop10[,c(5,6,7)]<-round(oratop10[,c(5,6,7)],3)
  colnames(oratop10)[9]<-"Cluster"
  return(oratop10)
}

getterm<- function(goid){
  termid<-GOTERM[[goid]]
  if(is.null(termid)) {
    return(paste(goid,"NA","NA", sep="__"))
  }else{
    return(paste(goid,Term(termid), Definition(termid), sep="__"))
  }
}


createBarPlot <- function(object, clustername, uniquePathways) {
  # Extracting categories and colors from the dataframe
  category_counts <- table(object@Data[[1]]$Groups[object@Data[[1]]$cluster==clustername])
  if(uniquePathways){
    colors <- object@plot$aka3Unique$Group[names(object@plot$aka3Unique$Group) %in% names(category_counts)]
    total<-length(unique(object@Data[[1]]$Pathways[object@Data[[1]]$cluster==clustername]))
  }else{
    colors <- object@plot$aka3$Group[names(object@plot$aka3$Group) %in% names(category_counts)]
    total<-length(object@Data[[1]]$Pathways[object@Data[[1]]$cluster==clustername])
  }
  colors <- colors[names(category_counts)]
  #create dataframe
  df<-as.data.frame(object@Data[[1]]$Groups[object@Data[[1]]$cluster==clustername])
  names(df)<-"Groups"
  df$Groups<-factor(df$Groups, levels=names(category_counts))
  #ylim_value <- max(category_counts) + 40
  title<-paste(clustername, " (Terms=", total, ")", sep = "")
  # Create a ggplot barplot
  ggplot(df, aes(x = Groups, fill = Groups)) +
    geom_bar(stat="count", width = 0.6, show.legend=FALSE) +
    scale_fill_manual(values = colors) +
    labs(title = title, y = "Count", x = "Groups") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size=16, face="bold"), panel.background=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=0.5), axis.ticks = element_line(color = "black"), axis.title.y = element_text(size=14), axis.title.x = element_blank(), axis.text.x=element_text(size=14, face="bold"), axis.text.y=element_text(size=14), legend.key= element_blank()
    ) +
    geom_text(stat = "count", aes(label = ..count..), vjust = -0.5)
}

createBarPlotIndependent <- function(object, combined_result, grouping, uniquePathways) {
  # Extracting categories and colors from the dataframe
  category_counts <- table(combined_result$Group[combined_result$Cluster==grouping])
  if(uniquePathways){
    colors <- object@plot$aka3Unique$Group[names(object@plot$aka3Unique$Group) %in% names(category_counts)]
    total<-length(unique(combined_result[combined_result$Cluster==grouping,1]))
  }else{
    colors <- object@plot$aka3$Group[names(object@plot$aka3$Group) %in% names(category_counts)]
    total<-length(combined_result[combined_result$Cluster==grouping,1])
  }
  colors <- colors[names(category_counts)]
  #create dataframe
  df<-as.data.frame(combined_result$Group[combined_result$Cluster==grouping])
  names(df)<-"Groups"
  df$Groups<-factor(df$Groups, levels=names(category_counts))
  #ylim_value <- max(category_counts) + 40
  title<-paste(grouping, " (Terms=", total, ")", sep = "")
  # Create a ggplot barplot
  ggplot(df, aes(x = Groups, fill = Groups)) +
    geom_bar(stat="count", width = 0.6, show.legend=FALSE) +
    scale_fill_manual(values = colors) +
    labs(title = title, y = "Count", x = "Groups") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size=16, face="bold"), panel.background=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=0.5), axis.ticks = element_line(color = "black"), axis.title.y = element_text(size=14), axis.title.x = element_blank(), axis.text.x=element_text(size=14, face="bold"), axis.text.y=element_text(size=14), legend.key= element_blank()
    ) +
    geom_text(stat = "count", aes(label = ..count..), vjust = -0.5)
}


transferinfo<-function(Object) {
  canonical.df<-Object@Data[[1]]
  #add in index column
  canonical.df$index<-seq_len(nrow(canonical.df))
  canonical.df_2<-Object@plot$aka2Unique
  canonical.df_2<-cbind(rownames(canonical.df_2),canonical.df_2)
  colnames(canonical.df_2)[1]<-"id"
  # merging the dataframes
  result<-merge(canonical.df, canonical.df_2, by.x="Pathways", by.y="id", all.x=TRUE)
  # Order the result by the original index
  result <- result[order(result$index), ]
  all(Object@Data[[1]]$Pathways==result$Pathways)
  #Object@Data[[1]]$cluster<-result$Cluster
  Object@Data[[1]]$cluster<-paste("Cluster_",result$Cluster, sep="")

  #return updated object
  return(Object)
}

calculategenestable<-function(Object, uniquePathways) {
  if(uniquePathways){
    clusters<-unique(Object@plot$aka2Unique$Cluster)
    clusters<-clusters[order(clusters)]
    clustersvector<-Object@plot$aka2Unique$Cluster
  }else{
    clusters<-unique(Object@plot$aka2$Cluster)
    clustersvector<-Object@plot$aka2$Cluster
  }
  j<-1
  for (i in clusters){
    genes<-as.vector(unlist(sapply(strsplit(Object@Data[[1]][,"Molecules"][clustersvector==i],Object@metadata$seperator[1]),unique)))
    genes<-gsub(" ","", genes)
    #count duplication per gene
    genesTimes<-as.data.frame(genes) %>% group_by_all() %>% count
    totalPathways<-sum(clustersvector==i)
    genesTimes$n<-round(genesTimes$n/totalPathways,6)
    colnames(genesTimes)[1]<-unique(as.character(Object@metadata[,"structure"]))
    if (j==1){
      genesFinal<-genesTimes
      colnames(genesFinal)[j+1]<- paste("Cluster_",i, " (Terms=", totalPathways,")", sep="")
    }else{
      genesFinal<-power_full_join(genesFinal, genesTimes, by=colnames(genesTimes)[1])
      colnames(genesFinal)[j+1]<- paste("Cluster_",i, " (Terms=", totalPathways,")", sep="")
    }
    j=j+1
  }
  genesFinal<-na_replace(genesFinal,0)
  return(as.data.frame(genesFinal))
}


calculateSymbol<-function(Object, genesFinal) {
  genes.df <- bitr(genesFinal[,1], fromType = unique(as.character(Object@metadata[,"structure"])), toType = c("SYMBOL"), OrgDb = unique(as.character(Object@metadata[,"organism"])))
  genes.df.merge<-merge(genesFinal, genes.df, by.x=colnames(genesFinal)[1], by.y=colnames(genes.df)[1], all.x=TRUE)
  genes.df.unique <- genes.df.merge[!duplicated(genes.df.merge[,1]), ]
  genes.df.unique<-genes.df.unique[,c(1, ncol(genes.df.unique))]
  return(genes.df.unique)
}

combineSymbolInfo<-function(genesFinal, symbolinfo) {
  genestable<-merge(genesFinal, symbolinfo, by.x=colnames(genesFinal)[1], by.y=colnames(symbolinfo)[1], all.x=TRUE)
  genestable<-genestable[,c(ncol(genestable),1:ncol(genestable)-1)]
  return(genestable)
}

calculategenestableinde<-function(Object, independent_info) {
  clusters<-unique(independent_info$Cluster)
  clusters<-clusters[order(clusters)]
  independent_info<-cbind(independent_info, paste(independent_info$Group,independent_info[,1], sep="_"))
  colnames(independent_info)[5]<-"RR_name"
  independent_info_molecules<-merge(independent_info, Object@Data[[1]][,c("RR_name", "Molecules")], by="RR_name", all.x=TRUE)
  j<-1
  for (i in clusters){
    #genes<-as.vector(unlist(sapply(strsplit(Object@Data[[1]][,"Molecules"][clustersvector==i],Object@metadata$seperator[1]),unique)))
    genes<-as.vector(unlist(sapply(strsplit(independent_info_molecules$Molecules[independent_info_molecules$Cluster==i],Object@metadata$seperator[1]),unique)))
    genes<-gsub(" ","", genes)
    #count duplication per gene
    genesTimes<-as.data.frame(genes) %>% group_by_all() %>% count
    totalPathways<-sum(independent_info_molecules$Cluster==i)
    genesTimes$n<-round(genesTimes$n/totalPathways,6)
    colnames(genesTimes)[1]<-unique(as.character(Object@metadata[,"structure"]))
    if (j==1){
      genesFinal<-genesTimes
      colnames(genesFinal)[j+1]<- paste(i, " (Terms=", totalPathways,")", sep="")
    }else{
      genesFinal<-power_full_join(genesFinal, genesTimes, by=colnames(genesTimes)[1])
      colnames(genesFinal)[j+1]<- paste(i, " (Terms=", totalPathways,")", sep="")
    }
    j=j+1
  }
  genesFinal<-na_replace(genesFinal,0)
  return(as.data.frame(genesFinal))
}

# Define server function
server <- function(input, output, session) {

  hideTab(inputId = "tabs3", target = "tabpanelresult")
  shinyjs::hide("summary.independentgroup")

  data.object<-reactive({})
  independent.info<-reactive({})
  independent.infoORA<-reactive({})
  genesfinal.info<-reactive({})

  output$tissueIntro<-renderUI({tags$p("Choose specific tissues for performing enrichment analysis (the tissues were selected from", HTML("<a href='https://gtexportal.org/home/' target='_blank'>here</a>"), "):")})

  output$independentIntro<-renderUI({tags$p("Press the buttom below to run the seriation-based analysis:")})

  ###########################
  ###########################
  # Store uploaded files
  files <- reactiveVal(list())
  files_data_group <- reactiveValues(group=NULL)

  observeEvent(input$files, {
    new_files <- input$files
    current_files <- files()
    if (length(unlist(files())) > 0){
      files(Map(c,current_files, new_files))
    }else{
      files(new_files)
    }
  })

  # Create DataTable of uploaded files
  output$fileTable <- renderDT({
    if (length(unlist(files())) > 0){
      if (is.null(files_data_group$group)){
        files_data_group$group<-as.character(seq_len(length(unlist(files()$name))))
      }else{
        if(all(!is.na(as.numeric(files_data_group$group)))){
          files_data_group$group<-as.character(seq_len(length(unlist(files()$name))))
        }else{
          if((length(files_data_group$group))==(length(unlist(files()$name)))){
          }else{
            news<-length(unlist(files()$name))-length(files_data_group$group)
            if(max(as.numeric(files_data_group$group), na.rm=TRUE)>0){
              maxvalue<-max(as.numeric(files_data_group$group), na.rm=TRUE)
              values<-(maxvalue+1):(maxvalue+news)
            }else{
              values<-1:news
            }
            currentgroup<-files_data_group$group
            files_data_group$group<-c(currentgroup, values)
          }
        }
      }
      files_data <- data.frame(
        File = unlist(files()$name),
        Group = files_data_group$group,
        Path = unlist(files()$datapath),
        Delete = '<button class="btn btn-danger btn-sm">Delete</button>',
        stringsAsFactors = FALSE
      )

      datatable(
        files_data,
        escape = FALSE,
        editable = TRUE,
        rownames = FALSE,
        options = list(
          dom = 't',
          columnDefs = list(list(targets = c(2), visible = FALSE, orderable = FALSE),
                            list(targets = c(0,1), className = "dt-center")
          )
        )
      )
    }else{
      NULL
    }
  })

  # Delete selected files
  observeEvent(input$fileTable_cell_clicked, {
    info <- input$fileTable_cell_clicked
    if (!is.null(info$value) && info$col == 3) {
      current_files <- files()
      files(lapply(current_files, function(lst) lst[-info$row]))
      files_data_group$group<-files_data_group$group[-info$row]
    }
  })

  # Observe changes in the DataTable and update the reactive value
  observeEvent(input$fileTable_cell_edit, {
    info <- input$fileTable_cell_edit
    if (!is.null(info$value) && info$col == 1) {
      # Check if the new value is different from the old value
      if (files_data_group$group[info$row] != info$value) {
        files_data_group$group[info$row] <- info$value
      }
    }
  })


  ###########################
  ###########################
  #ANNOTATIONS GO
  annotations<-reactive({})

  dataInfo_previousvalue<-reactiveVal("none")
  oraInfo_previousvalue<-reactiveVal("none")
  independentOutput_previousvalue<-reactiveVal("none")
  independentORAOutput_previousvalue<-reactiveVal("none")
  genesInfo_previousvalue<-reactiveVal("none")
  genesInfoIndependent_previousvalue<-reactiveVal("none")
  #links to quickGO
  observe({
    #Data table
    if (!is.null(input$dataInfo_cell_clicked$col)){
      if (input$dataInfo_cell_clicked$col==0){
        #if(input$dataInfo_cell_clicked$value!=""){
        if((input$dataInfo_cell_clicked$value!="")&&(input$dataInfo_cell_clicked$value!=dataInfo_previousvalue())){
          dataInfo_previousvalue(input$dataInfo_cell_clicked$value)
          js$browseURL(paste("https://www.ebi.ac.uk/QuickGO/term/",input$dataInfo_cell_clicked$value, sep=""))
        }
      }
    }
    #ORA table
    if (!is.null(input$oraInfo_cell_clicked$col)){
      if (input$oraInfo_cell_clicked$col==0){
        #if(input$oraInfo_cell_clicked$value!=""){
        if((input$oraInfo_cell_clicked$value!="")&&(input$oraInfo_cell_clicked$value!=oraInfo_previousvalue())){
          oraInfo_previousvalue(input$oraInfo_cell_clicked$value)
          js$browseURL(paste("https://www.ebi.ac.uk/QuickGO/term/",input$oraInfo_cell_clicked$value, sep=""))
        }
      }
    }
    #Independent table
    if (!is.null(input$independentOutput_cell_clicked$col)){
      if (input$independentOutput_cell_clicked$col==0){
        #if(input$oraInfo_cell_clicked$value!=""){
        if((input$independentOutput_cell_clicked$value!="")&&(input$independentOutput_cell_clicked$value!=independentOutput_previousvalue())){
          independentOutput_previousvalue(input$independentOutput_cell_clicked$value)
          js$browseURL(paste("https://www.ebi.ac.uk/QuickGO/term/",input$independentOutput_cell_clicked$value, sep=""))
        }
      }
    }
    #Independent ORA table
    if (!is.null(input$independentORAOutput_cell_clicked$col)){
      if (input$independentORAOutput_cell_clicked$col==0){
        #if(input$oraInfo_cell_clicked$value!=""){
        if((input$independentORAOutput_cell_clicked$value!="")&&(input$independentORAOutput_cell_clicked$value!=oraInfo_previousvalue())){
          oraInfo_previousvalue(input$independentORAOutput_cell_clicked$value)
          js$browseURL(paste("https://www.ebi.ac.uk/QuickGO/term/",input$independentORAOutput_cell_clicked$value, sep=""))
        }
      }
    }
    #Gene table
    if (!is.null(input$genesInfo_cell_clicked$col)){
      if (input$genesInfo_cell_clicked$col==0){
        if((input$genesInfo_cell_clicked$value!="")&&(input$genesInfo_cell_clicked$value!=genesInfo_previousvalue())){
          genesInfo_previousvalue(input$genesInfo_cell_clicked$value)
          if(unique(data.object()@metadata$organism)=="org.Mm.eg.db"){
            js$browseURL(paste("https://www.informatics.jax.org/quicksearch/summary?queryType=exactPhrase&query=", input$genesInfo_cell_clicked$value,"&submit=Quick+Search", sep=""))
          }else{
            js$browseURL(paste("https://www.genecards.org/cgi-bin/carddisp.pl?gene=",input$genesInfo_cell_clicked$value, sep=""))
          }
        }
      }
    }
    #Gene table independent
    if (!is.null(input$genesInfoIndependent_cell_clicked$col)){
      if (input$genesInfoIndependent_cell_clicked$col==0){
        if((input$genesInfoIndependent_cell_clicked$value!="")&&(input$genesInfoIndependent_cell_clicked$value!=genesInfoIndependent_previousvalue())){
          genesInfoIndependent_previousvalue(input$genesInfoIndependent_cell_clicked$value)
          if(unique(data.object()@metadata$organism)=="org.Mm.eg.db"){
            js$browseURL(paste("https://www.informatics.jax.org/quicksearch/summary?queryType=exactPhrase&query=", input$genesInfoIndependent_cell_clicked$value,"&submit=Quick+Search", sep=""))
          }else{
            js$browseURL(paste("https://www.genecards.org/cgi-bin/carddisp.pl?gene=",input$genesInfoIndependent_cell_clicked$value, sep=""))
          }
        }
      }
    }
  })

  ###########################
  ###########################
  #CLEAR
  observeEvent(input$reset, {
    runexample(FALSE)
    shinyjs::runjs("location.reload(true);")
  })


  ###########################
  ###########################
  #RUN EXAMPLE
  runexample<- reactiveVal(FALSE)

  observeEvent(input$setExample, {
    runexample(TRUE)
    if (input$exampledataset=="mousedataset"){
      updateRadioButtons(session,"source", selected = "Great")
      updateRadioButtons(session,"gene.id", selected = "SYMBOL")
      updateRadioButtons(session,"organism", selected = "org.Mm.eg.db")
      updateRadioButtons(session,"filters", selected = "none")

      tablegroup<-cbind(c("MM10.GREAT.KO.uGvsMac.bed.tsv", "MM10.GREAT.WT.uGvsMac.bed.tsv"), c("KO", "WT"))
      colnames(tablegroup)<-c("File","Group")
    }
    if (input$exampledataset=="coviddataset"){
      updateRadioButtons(session,"source", selected = "GSEA")
      updateRadioButtons(session,"gene.id", selected = "ENTREZID")
      updateRadioButtons(session,"organism", selected = "org.Hs.eg.db")
      updateRadioButtons(session,"filters", selected = "uniquepathways")

      tablegroup<-cbind(c("Covid19AI_Healthy.csv", "Covid193Mo_Healthy.csv", "Covid196Mo_Healthy.csv"), c("AI", "Mo3", "Mo6"))
      colnames(tablegroup)<-c("File","Group")
    }

    output$fileTable<-renderDataTable(datatable(tablegroup, rownames=FALSE, editable=TRUE, options=list(pageLength=25, dom='t')))

    click("printResults")
  })

  wordcloudperform<-reactiveVal(FALSE)
  wordcloudinfo<-reactive({})
  uniquepathwaysperform<-reactiveVal(FALSE)
  observe({
    if (!is.null(input$filters)&&(!runexample())){
      #if((input$filters=="uniquepathways")&&(!runexample())){
      if((input$filters=="uniquepathways")){
        uniquepathwaysperform(TRUE)
      }else{
        uniquepathwaysperform(FALSE)
      }
    }
  })

  ###########################
  ###########################
  #CALCULATE RESULTS (including duplicates, clustering=hierarchival and number of cluster=5)
  treeplotoutput<-reactive({})
  heatmapoutput<-reactive({})
  tissueoutput<-reactive({})
  tissueoutputindependent<-reactive({})
  independentoutput<-reactive({})

  ora.object<-reactive({})
  ora.objectALL<-reactive({})

  genes.symbolinfo<-reactive({})

  observeEvent(input$printResults, {


    tocontinue=TRUE
    if(((is.null(input$source))||(is.null(input$files))||(is.null(input$gene.id))||(is.null(input$organism)))&&(!runexample())){
      shinyalert(title = "Missing inputs", text = "Please, all the inputs are mandatory.", type = "error")
    }else{

      if(runexample()){
        #LOAD GREAT MOUSE EXAMPLE
        show_modal_progress_line(text="Running example")
        if (input$exampledataset=="mousedataset"){
          load("example/GSC_examplemouse.RData")
        }else{
          load("example/GSC_examplecovid.RData")
        }
        data.object<<-reactive({data_to_save})
        annotations<<-reactive({annotations_to_save})
        ora.objectALL<<-reactive({ora_to_save})
        uniquepathwaysperform<<-reactive({unique_pathways})
        update_modal_progress(value=0.3)
        if(uniquepathwaysperform()){
          clusterchoices<-paste("Cluster_", names(data.object()@plot$aka3Unique$Cluster), sep="")
          clusterchoices<-clusterchoices[order(clusterchoices)]
          clusters<-unique(data.object()@plot$aka2Unique$Cluster)
          clusters<-clusters[order(clusters)]
        }else{
          clusterchoices<-paste("Cluster_", names(data.object()@plot$aka3$Cluster), sep="")
          clusters<-unique(data.object()@plot$aka2$Cluster)
        }
        output$checkboxCluster <- renderUI({
          clusterchoices_all<-c("Uncheck all", clusterchoices)
          checkboxGroupInput("optionCluster", label = "", choices=clusterchoices_all, selected = clusterchoices, inline=TRUE)
        })

        #genes info
        genesFinal<-calculategenestable(data.object(), uniquepathwaysperform())
        if(unique(as.character(data.object()@metadata[,"structure"]))!="SYMBOL"){
          symbolinfo<-calculateSymbol(data.object(), genesFinal)
          genes.symbolinfo<<-reactive({symbolinfo})
          genesFinal<-combineSymbolInfo(genesFinal, genes.symbolinfo())
        }
        genesfinal.info<<-reactive({genesFinal})
        output$genesInfo <-renderDataTable(datatable(genesFinal, rownames = FALSE, filter="top", options = list(dom='lrtip', pageLength=10, rowCallback = JS(underlineCellsCols(c(1),c(0)))))
                                           %>% formatStyle(colnames(genesFinal)[1], color = "blue"))
        oraResultsTop<-do.call("rbind", lapply(clusters, topORA, ora=ora.objectALL(), top=as.integer(input$top)))
        ora.object<<-reactive({oraResultsTop})


        ##check if tissue is perfomed for dependent analysis
        if ((!is.null(data.object()@dfTissue))&&(length(data.object()@dfTissue>0))){

          tissueData<-data.frame(rownames(data.object()@dfTissue), data.object()@dfTissue)
          colnames(tissueData)[1]<-"Tissue"
          output$tissueOutput<-renderDataTable(datatable(tissueData, rownames = FALSE, options = list(searching=TRUE, pageLength=10)))
          #shinyjs::hide("calculateTissueEnrichment")

          output$tissueIntro<-renderUI({tags$p("Tissue enrichment results per cluster (Cluster_number..number of terms):")})
          showTab(inputId = "tabs", target = "Tissue")

          #plot tissue plots
          tissueoutput<<-reactive({PlotTissueExpression(data.object(), all = FALSE, uniquePathways=uniquepathwaysperform(), clusterIndependent=FALSE)})
          output$tissuePlot<-renderPlot(tissueoutput())
        }else{
          #hide tissue results
          hideTab(inputId = "tabs", target = "Tissue")
          hideTab(inputId = "tabs2", target = "Tissue enrichment")
        }

        independent.info<<-reactive({independent_to_save})
        independent.infoORA<<-reactive({independent_ora_to_save})
        shinyjs::hide("runIndependent")
        output$independentIntro<-renderUI({tags$p("")})
        #independent
        groupchoices<-unique(independent.info()$Cluster)
        if(checkGO(data.object())){
          output$independentOutput <-renderDataTable(datatable(independent.info(), rownames = FALSE, options = list(autoWidth=TRUE, scrollX=TRUE, scrollY=400, searching=TRUE, pageLength=25, rowCallback = JS(underlineCellsCols(c(1),c(0)))))
                                                     %>% formatStyle(colnames(independent.info())[4], backgroundColor = styleEqual(c(paste("Cluster_", names(data.object()@plot[[6]]$Cluster)[names(data.object()@plot[[6]]$Cluster)!=0], sep="")), c(data.object()@plot[[6]]$Cluster[data.object()@plot[[6]]$Cluster!="white"])))
                                                     %>% formatStyle(colnames(independent.info())[1], color = "blue"))
        }else{
          output$independentOutput <-renderDataTable(datatable(independent.info(), rownames = FALSE, options = list(autoWidth=TRUE, scrollX=TRUE, scrollY=400, searching=TRUE, pageLength=25))
                                                     %>% formatStyle(colnames(independent.info())[4], backgroundColor = styleEqual(c(paste("Cluster_", names(data.object()@plot[[6]]$Cluster)[names(data.object()@plot[[6]]$Cluster)!=0], sep="")), c(data.object()@plot[[6]]$Cluster[data.object()@plot[[6]]$Cluster!="white"]))))
        }
        #independent ORA
        output$independentORAOutput<-renderDataTable(datatable(independent.infoORA(), rownames = FALSE, options = list(autoWidth=TRUE, scrollX=TRUE, scrollY=400, searching=TRUE, pageLength=10, rowCallback = JS(underlineCellsCols(c(1),c(0)))))
                                                     %>% formatStyle(colnames(independent.infoORA())[9], backgroundColor = styleEqual(c(paste("Cluster_", names(data.object()@plot[[6]]$Cluster)[names(data.object()@plot[[6]]$Cluster)!=0], sep="")), c(data.object()@plot[[6]]$Cluster[data.object()@plot[[6]]$Cluster!="white"])))
                                                     %>% formatStyle(colnames(independent.infoORA())[1], color = "blue"))

        genesFinalIndependent<-calculategenestableinde(data.object(),independent.info())
        if(unique(as.character(data.object()@metadata[,"structure"]))!="SYMBOL"){
          genesFinalIndependent<-combineSymbolInfo(genesFinalIndependent, genes.symbolinfo())
        }
        output$genesInfoIndependent<-renderDataTable(datatable(genesFinalIndependent, rownames = FALSE, filter="top", options = list(autoWidth=TRUE, scrollX=TRUE, scrollY=400, searching=TRUE, dom='lrtip', pageLength=10, rowCallback = JS(underlineCellsCols(c(1),c(0)))))
                                                     %>% formatStyle(colnames(genesFinalIndependent)[1], color = "blue"))
        updateSelectInput(session, "summary.independentgroup", choices = groupchoices)
        shinyjs::show("summary.independentgroup")
        #create plot
        keywords_ora_independent<-calculateKeywordsORA(independent.infoORA())
        independentoutput<<-reactive({PlotPathwayCluster(data.object(), doORA = TRUE, wordcloud = FALSE, uniquePathways = uniquepathwaysperform(), keywords_ora_inde=keywords_ora_independent)})
        output$independentPlot<-renderPlot(independentoutput())
        if((data.object()@metadata[1,"organism"]=="org.Mm.eg.db")){
          #for mouse data hide tissue enrichment analysis
          hideTab(inputId = "tabsindependent", target = "Tissue")
        }
        #tissue independent
        if(!is.null(data.object()@dfTissueIndependent)){
          tissueDataIndependent<-data.frame(rownames(data.object()@dfTissueIndependent), data.object()@dfTissueIndependent)
          colnames(tissueDataIndependent)[1]<-"Tissue"
          output$tissueOutputIndependent<-renderDataTable(datatable(tissueDataIndependent, rownames = FALSE, options = list(autoWidth=TRUE, scrollX=TRUE, scrollY=400, searching=TRUE, pageLength=10)))
          tissueoutputindependent<<-reactive({PlotTissueExpression(data.object(), all = FALSE, uniquePathways=uniquepathwaysperform(), clusterIndependent=TRUE)})
          output$tissuePlotIndependent<-renderPlot(tissueoutputindependent())
        }else{
          hideTab(inputId = "tabs", target = "Tissue_S")
        }

        update_modal_progress(value=0.5)

        #download data and plots
        output$saveObject <- renderUI({
          downloadButton("save.object","Results")
        })
        #download data and plots
        output$downloadData <- renderUI({
          downloadButton("download.data","Data")
        })
        #download plots using different format
        output$downloadDataFormats <- renderUI({
          dropdownButton(inputId="download.plot.format",label="Images", circle=FALSE, radioButtons("format", label="", choices = c("jpg","png","pdf"), selected = NULL), downloadButton("download.plots","Download plots"))
        })

        nclusters<-length(clusters)
        updateNumericInput(session, inputId="cluster", value=nclusters)

        if(checkGO(data.object())){
          wordcloudperform(TRUE)
          wordcloudinfo<<-reactive({PerformORAGOperCluster(data.object(), uniquePathways=uniquepathwaysperform(), clusterIndependent=FALSE)})
          output$wordcloudCheckbox <- renderUI({checkboxInput("wordcloud", "Enable Wordcloud", value=FALSE)})
        }else{
          wordcloudperform(FALSE)
          output$wordcloudCheckbox <- renderUI({tagList()})
          #create plots
          keywords_ora<-calculateKeywordsORA(ora.objectALL())
          heatmapoutput<<-reactive({PlotGeneSets(Object = data.object(), doORA = T, wordcloud = wordcloudperform(), uniquePathways=uniquepathwaysperform(), keywords_ora=keywords_ora)})
          output$heatmap<-renderPlot(heatmapoutput())
        }
        update_modal_progress(value=0.7)

        output$barplot<-renderPlot(createBarPlot(data.object(), clusterchoices[1], uniquePathways=uniquepathwaysperform()))
        update_modal_progress(value=0.9)

        updateSelectInput(session, "breakup.cluster", choices = clusterchoices)
        updateSelectInput(session, "summary.cluster", choices = clusterchoices)

        showTab(inputId = "tabs3", target="tabpanelresult", select=TRUE)
        #updateTabsetPanel(session, "tabs3", selected = "tabpanelresult", args = list(tabpanelresult = "Results"))
        runjs('
            // Change the tab label when the button is clicked
            $("#tabs3 a[data-value=\'tabpanelresult\']").text("Results");
            ')
        shinyjs::toggle("main1")
        shinyjs::toggle("main2")
        remove_modal_progress()

      }else{
        show_modal_progress_line(text="")
        update_modal_progress(value=0.2, text="Load gene sets (1/5)") # update progress bar value
        #showModal(modalDialog("Load gene sets (1/5)", footer=NULL))
        #1-if input=great load gene set with specifc parameters and perform manage.object
        if ((input$source=="Great") || (input$source=="GSEA") || (input$source=="GSEA2")){
          if(input$source=="GSEA2"){
            load.object <- LoadGeneSets_shiny(file_location = files()$datapath,
                                              groupnames= files_data_group$group,
                                              P.cutoff = 0.05,
                                              Mol.cutoff = 20,
                                              Source = "GSEA",
                                              Great.Background = FALSE, #TRUE
                                              type = "Canonical_Pathways",
                                              topranks = "",
                                              structure = input$gene.id,
                                              Organism = input$organism,
                                              seperator = ",") #,
          }else{
            load.object <- LoadGeneSets_shiny(file_location = files()$datapath,
                                              groupnames= files_data_group$group,
                                              P.cutoff = 0.05,
                                              Mol.cutoff = 20,
                                              Source = input$source,
                                              Great.Background = FALSE, #TRUE
                                              type = "Canonical_Pathways",
                                              topranks = "",
                                              structure = input$gene.id,
                                              Organism = input$organism,
                                              seperator = ",") #,
          }

          if(is.character(load.object)){
            remove_modal_progress()
            #one of the group has 0 pathways
            shinyalert("Error", paste("Group ", load.object, " has 0 pathways passing cutoff (p<0.05)"), type = "error")
            tocontinue=FALSE
          }else{
            #showModal(modalDialog("Manage gene sets (2/5).", footer=NULL))
            update_modal_progress(value=0.4, text="Manage gene sets (2/5)")
            manage.object <- ManageGeneSets(Object = load.object,
                                            keep.type =c("Disease Ontology",
                                                         "GO Biological Process" ),
                                            exclude.type="")
          }

        }
        #1-if input=ipa load gene set with specifc parameters and don't perform manage.object
        else if ((input$source=="IPA")){

          load.object <- LoadGeneSets_shiny(file_location = input$files[,4],
                                            groupnames= tablegroup$data[,2],
                                            P.cutoff = 1.3,
                                            Mol.cutoff = 5,
                                            Source = input$source,
                                            type = "Canonical_Pathways",
                                            structure = input$gene.id,
                                            Organism = input$organism,
                                            seperator = ",")

          if(is.character(load.object)){
            remove_modal_progress()
            #one of the group has 0 pathways
            shinyalert("Error", paste("Group ", load.object, " has 0 pathways passing cutoff (p<0.05)"), type = "error")
            tocontinue=FALSE
          }else{
            manage.object<-load.object
          }
          #annotationperform=FALSE

        }
        if(tocontinue){

          rm(load.object)

          #check seperator
          possibleSeparators<-c(",","/",";")
          nelements<-length(manage.object@Data)
          for (i in (1:nelements)){
            counts<-unlist(lapply(possibleSeparators,getcharacter,list=manage.object@Data[[i]]["Molecules"]))
            manage.object@metadata$seperator[i]<-possibleSeparators[which(counts==max(counts))]
          }
          if (length(unique(manage.object@metadata$seperator))!=1){
            remove_modal_progress()
            shinyalert("Error", "Please, use the same gene-seperator (, ; /) in all data ", type = "error")

          }else{

            update_modal_progress(value=0.6, text="Combine gene sets (3/5)")
            #2-perform combine gene set
            combine.object <- CombineGeneSets(Object = manage.object, threads = 16)
            rm(manage.object)

            #2.1-calculate optimal number of cluster
            optimal.object<-OptimalGeneSets_2(object = combine.object, method = "silhouette", max_cluster= 10, cluster_method = "kmeans", main= "", uniquePathways=uniquepathwaysperform())
            nclust<-optimal.object$data

            optimalCluster<-as.numeric(nclust$clusters[which.max(nclust$y)])
            updateNumericInput(session, inputId="cluster", value=optimalCluster)

            update_modal_progress(value=0.8, text=paste("Cluster gene sets, optimal k=",optimalCluster," (4/5)", sep=""))
            #3-perform cluster gene set
            cluster.object <- ClusterGeneSets(Object = combine.object,
                                              clusters = optimalCluster,
                                              method = "kmeans", order="cluster")
            rm(combine.object)



            if(uniquepathwaysperform()){
              cluster.object <- transferinfo(cluster.object)
              clusterchoices<-paste("Cluster_", names(cluster.object@plot$aka3Unique$Cluster), sep="")
              clusters<-unique(cluster.object@plot$aka2Unique$Cluster)
              clusters<-clusters[order(clusters)]
            }else{
              cluster.object@Data[[1]]$cluster<-paste(rep("Cluster_",dim(cluster.object@Data[[1]]["cluster"])[1]), cluster.object@Data[[1]]$cluster, sep = "")
              clusters<-unique(cluster.object@plot$aka2$Cluster)
              clusterchoices<-paste("Cluster_", names(cluster.object@plot$aka3$Cluster), sep="")
            }
            clusterchoices<-clusterchoices[order(clusterchoices)]
            data.object<<-reactive({cluster.object})



            #download data and plots
            output$saveObject <- renderUI({
              downloadButton("save.object","Results")
            })
            #download data and plots
            output$downloadData <- renderUI({
              downloadButton("download.data","Data")
            })
            #download plots using different format
            output$downloadDataFormats <- renderUI({
              dropdownButton(inputId="download.plot.format",label="Images", circle=FALSE, radioButtons("format", label="", choices = c("jpg","png","pdf"), selected = NULL), downloadButton("download.plots","Download plots"))
            })

            #put check box depending the number of cluster

            output$checkboxCluster <- renderUI({
              clusterchoices_all<-c("Uncheck all", clusterchoices)
              checkboxGroupInput("optionCluster", label = "", choices=clusterchoices_all, selected = clusterchoices, inline=TRUE)
            })
            output$barplot<-renderPlot(createBarPlot(cluster.object, clusterchoices[1], uniquePathways=uniquepathwaysperform()))

            if(!is.character(cluster.object@Data[[1]]$Ratio[1])){
              cluster.object@Data[[1]]$Ratio<-round(cluster.object@Data[[1]]$Ratio,3)
            }
            #calculate annotations
            if((is.null(annotations())) && (checkGO(cluster.object)==TRUE)){
              resultAnnotations<-as.data.frame(str_split_fixed(as.data.frame(unlist(lapply(cluster.object@Data[[1]][,"Pathways"],getterm)))[,1],"__",3))
              colnames(resultAnnotations)<-c("ID", "Term", "Definition")
              resultAnnotationsFinal<-data.frame(resultAnnotations[,c(1,2)], cluster.object@Data[[1]][c("Type","Groups","cluster","Pval","Ratio")], resultAnnotations[,3])
              colnames(resultAnnotationsFinal)[c(4,5,8)]<-c("Group", "Cluster", "Definition")
              annotations<<-reactive({resultAnnotationsFinal})
            }else{
              resultAnnotationsFinal<-data.frame(cluster.object@Data[[1]][c("Pathways","Type","Groups","cluster","Ratio")])
              colnames(resultAnnotationsFinal)[c(3,4)]<-c("Group", "Cluster")
              annotations<<-reactive({resultAnnotationsFinal})
            }

            #4-ORA per cluster using enrichGO of clusterprofiler
            update_modal_progress(value=0.9, text="ORA per cluster (5/5)")

            genesFinal<-calculategenestable(cluster.object, uniquepathwaysperform())
            if(unique(as.character(cluster.object@metadata[,"structure"]))!="SYMBOL"){
              symbolinfo<-calculateSymbol(cluster.object, genesFinal)
              genes.symbolinfo<<-reactive({symbolinfo})
              genesFinal<-combineSymbolInfo(genesFinal, genes.symbolinfo())
            }

            #calulate ORA in parallel
            oraResults<-do.call("rbind", lapply(clusters, calculateORA, object=cluster.object, uniquePathways=uniquepathwaysperform()))
            ora.objectALL<<-reactive({oraResults})
            oraResultsTop<-do.call("rbind", lapply(clusters, topORA, ora=oraResults, top=as.integer(input$top)))
            ora.object<<-reactive({oraResultsTop})

            output$oraInfo<-renderDataTable(datatable(as.data.frame(oraResultsTop), rownames = FALSE, options = list(searching=TRUE, pageLength=10, rowCallback = JS(underlineCellsCols(c(1),c(0)))))
                                            %>% formatStyle(colnames(oraResultsTop)[1], color = "blue"))


            genesfinal.info<<-reactive({genesFinal})
            output$genesInfo <-renderDataTable(datatable(genesFinal, rownames = FALSE, filter="top", options = list(dom='lrtip', pageLength=10, rowCallback = JS(underlineCellsCols(c(1),c(0)))))
                                               %>% formatStyle(colnames(genesFinal)[1], color = "blue"))



            if(checkGO(cluster.object)){
              wordcloudperform(TRUE)
              wordcloudinfo<<-reactive({PerformORAGOperCluster(cluster.object, uniquePathways=uniquepathwaysperform(), clusterIndependent=FALSE)})
              output$wordcloudCheckbox <- renderUI({checkboxInput("wordcloud", "Enable Wordcloud", value=FALSE)})
            }else{
              wordcloudperform(FALSE)
              output$wordcloudCheckbox <- renderUI({tagList()})
              #plot heatmap
              keywords_ora<-calculateKeywordsORA(ora.objectALL())
              heatmapoutput<<-reactive({PlotGeneSets(cluster.object, doORA = T, wordcloud = wordcloudperform(), uniquePathways=uniquepathwaysperform(), keywords_ora=keywords_ora)})
              output$heatmap<-renderPlot(heatmapoutput())
            }


            updateSelectInput(session, "breakup.cluster", choices = clusterchoices)
            updateSelectInput(session, "summary.cluster", choices = clusterchoices)


            #hide tissue and independent results
            hideTab(inputId = "tabs", target = "Tissue")
            hideTab(inputId = "tabs", target = "Heatmap_S")
            hideTab(inputId = "tabs", target = "Tissue_S")
            hideTab(inputId = "tabsindependent", target = "ORA")
            hideTab(inputId = "tabsindependent", target = "Genes")
            hideTab(inputId = "tabsindependent", target = "Tissue")
            if((input$organism=="org.Mm.eg.db")||(is.null(input$organism))){
              #for mouse data hide tissue enrichment analysis
              hideTab(inputId = "tabs2", target = "Tissue enrichment")
            }

            showTab(inputId = "tabs3", target="tabpanelresult", select=TRUE)
            #updateTabsetPanel(session, "tabs3", selected = "tabpanelresult", args = list(tabpanelresult = "Results"))
            runjs('
            // Change the tab label when the button is clicked
            $("#tabs3 a[data-value=\'tabpanelresult\']").text("Results");
            ')
            shinyjs::toggle("main1")
            shinyjs::toggle("main2")
            remove_modal_progress()
          }
        }
      }

    }


  }) #end of printResults

  ###########################
  ###########################
  #WORD CLOUD CHECK
  observeEvent(input$wordcloud, {
    # Check if the checkbox is unchecked
    if (!input$wordcloud) {
      wordcloudperform(FALSE)
    }else{
      wordcloudperform(TRUE)
    }
    keywords_ora<-calculateKeywordsORA(ora.objectALL())
    heatmapoutput<<-reactive({PlotGeneSets(Object = data.object(), doORA = T, wordcloud = wordcloudperform(), wordclouds=wordcloudinfo(), uniquePathways=uniquepathwaysperform(), keywords_ora=keywords_ora)})
    output$heatmap<-renderPlot(heatmapoutput())
  })

  ###########################
  ###########################
  #UPLOAD FILTERS
  observe({
    if(!is.null(input$uploadfrom)){
      if(input$uploadfrom=="shiny"){
        shinyjs::hide("uploadfilters")
        shinyjs::hide("uploadnameobject")
        shinyjs::show("fileRdata")
        shinyjs::show("uploadResult")
      }
      if(input$uploadfrom=="package"){
        shinyjs::show("uploadfilters")
        shinyjs::show("uploadnameobject")
        shinyjs::show("fileRdata")
        shinyjs::show("uploadResult")
      }
    }else{
      shinyjs::hide("uploadfilters")
      shinyjs::hide("uploadnameobject")
      shinyjs::hide("fileRdata")
      shinyjs::hide("uploadResult")
    }
  })

  ###########################
  ###########################
  #DATABASE OPTIONS
  output$databases <- renderUI({
    if(!is.null(data.object())){
      if (data.object()@metadata[1,"organism"]=="org.Hs.eg.db") {
        updateSelectizeInput(session = session, inputId = "database", choices=c("","Human Phenotype Ontology (HPO)","Customize"), server=TRUE)
        selectizeInput("database", "Human databases:", choices=NULL, selected = NULL, multiple=FALSE)
      }else if (data.object()@metadata[1,"organism"]=="org.Mm.eg.db") {
        updateSelectizeInput(session = session, inputId = "database", choices=c("","Mammalian Phenotype (MP)", "Customize"), server=TRUE)
        selectizeInput("database", "Mouse databases:", choices=NULL, selected = NULL, multiple=FALSE)
      }
    }
  })
  output$databaseOptions <- renderUI({
    if(!is.null(input$database)){
      if (input$database=="Human Phenotype Ontology (HPO)") {
        updateSelectizeInput(session = session, inputId = "hpo", choices=c("",unique(hpoDatabase$hpo_id)), server=TRUE)
        selectizeInput("hpo", HTML("Choose HPO <a href='https://hpo.jax.org/app/browse/term/HP:0000118' target='_blank'>(info)</a>:"), choices=NULL, selected = NULL, multiple=FALSE)
      }else if (input$database=="Mammalian Phenotype (MP)") {
        updateSelectizeInput(session = session, inputId = "mp", choices=c("",unique(mpDatabase$mp_id)), server=TRUE)
        selectizeInput("mp", HTML("Choose MP <a href='https://www.informatics.jax.org/vocab/mp_ontology' target='_blank'>(info)</a>:"), choices=NULL, selected = NULL, multiple=FALSE)
      }else if (input$database=="Customize") {
        fileInput("fileCustomize", "One column with gene symbols", accept = c("text/csv","text/comma-separated-values,text/plain",".csv", ".txt", ".xls", ".xlsx"))
      }
    }
  })
  genescustomhighlight<-reactive({})
  observeEvent(input$fileCustomize,{
    #determine file extension
    file_ext <- file_ext(input$fileCustomize$name)
    if (file_ext %in% c("csv", "txt")) {
      dfcustomize <- read.table(input$fileCustomize$datapath, header = FALSE)
    } else {
      dfcustomize <- read_xlsx(input$fileCustomize$datapath, sheet = 1, col_names = FALSE)
    }
    colnames(dfcustomize)<-"genes"
    genescustomhighlight<<-reactive({dfcustomize})
    output$databaseOptionTitle<-renderText({paste("Number of genes = ", dim(dfcustomize)[1])})
    output$databaseGenes<-renderText({dfcustomize$genes})
  })

  observe({
    if((input$hpo!="")&&(!is.null(input$hpo))){
      output$databaseOptionTitle<-renderText({hpoDatabase$hpo_id[hpoDatabase$hpo_id==input$hpo][1]})
      output$databaseGenes<-renderText({hpoDatabase$gene_symbol[hpoDatabase$hpo_id==input$hpo]})
    }
    if((input$mp!="")&&(!is.null(input$mp))){
      output$databaseOptionTitle<-renderText({mpDatabase$mp_id[mpDatabase$mp_id==input$mp][1]})
      output$databaseGenes<-renderText({mpDatabase$gene_symbol[mpDatabase$mp_id==input$mp]})
    }
    if(!is.null(input$summary.cluster)&&(!is.null(data.object()))){
      output$barplot<-renderPlot(createBarPlot(data.object(), input$summary.cluster, uniquePathways=uniquepathwaysperform()))
    }
    if(!is.null(input$summary.independentgroup)&&(!is.null(data.object()))&&(!is.null(independent.info()))){
      output$barplotIndependent<-renderPlot(createBarPlotIndependent(data.object(), independent.info(), input$summary.independentgroup, uniquePathways=uniquepathwaysperform()))
    }
  })
  observeEvent(input$clearHighlight, {
    if(!is.null(input$database)){
      if (input$database=="Human Phenotype Ontology (HPO)") {
        shinyjs::hide("hpo")
      }else if (input$database=="Mammalian Phenotype (MP)") {
        shinyjs::hide("mp")
      }else if (input$database=="Customize") {
        shinyjs::hide("fileCustomize")
      }
      if (data.object()@metadata[1,"organism"]=="org.Hs.eg.db") {
        updateSelectizeInput(session = session, inputId = "database", choices=c("","Human Phenotype Ontology (HPO)","Customize"),selected=NULL, server=TRUE)
      }else{
        updateSelectizeInput(session = session, inputId = "database", choices=c("","Mammalian Phenotype (MP)", "Customize"), server=TRUE)
      }
      output$highlightInfo<-renderDataTable(NULL)
      output$databaseOptionTitle<-renderText({NULL})
      output$databaseGenes<-renderText({NULL})
    }
  })

  observeEvent(input$calculateHighlight, {
    if((input$database=="")|((input$database=="Mammalian Phenotype (MP)")&&(input$mp==""))|((input$database=="Human Phenotype Ontology (HPO)")&&(input$hpo==""))){
      shinyalert(title = "Wrong database", text = "Please, select a database and a gene set.", type = "error")
    }else{

      show_modal_progress_line(value=0.2, text="Calculate highlight genes")
      update_modal_progress(value=0.3)
      if((input$hpo!="")&&(!is.null(input$hpo))){
        genes<-hpoDatabase$gene_symbol[hpoDatabase$hpo_id==input$hpo]
      }
      if((input$mp!="")&&(!is.null(input$mp))){
        genes<-mpDatabase$gene_symbol[mpDatabase$mp_id==input$mp]
      }
      if((!is.null(genescustomhighlight))&&(input$database=="Customize")){
        genes<-genescustomhighlight()$genes
      }
      update_modal_progress(value=0.6)
      if(unique(as.character(data.object()@metadata[,"structure"]))!="SYMBOL"){
        highlight.object<-HighlightGeneSets_shiny(data.object(), highligt.genes = genes, name="", genesinfo=genes.symbolinfo())
      }else{
        highlight.object<-HighlightGeneSets_shiny(data.object(), highligt.genes = genes, name="", genesinfo=NULL)
      }
      data.object<<-reactive({highlight.object})
      resultHighlight<-(t(as.data.frame(str_split(unique(paste(data.object()@Data[[1]]$cluster, round(data.object()@Data[[1]]$Highlight.mean, 3)))," "))))
      update_modal_progress(value=0.8)
      colnames(resultHighlight)<-c("Cluster", "Mean")
      resultHighlight<-resultHighlight[order(resultHighlight[,1]),]
      output$highlightInfo<-renderDataTable(datatable(resultHighlight, rownames = FALSE, options = list(scrollY=400, dom='t')))
      remove_modal_progress()
    }
  })

  ###########################
  ###########################
  #FILTER TOP ORA
  observeEvent(input$top, {
    if(!is.null(data.object())){
      if(uniquepathwaysperform()){
        clusters<-unique(data.object()@plot$aka2Unique$Cluster)
        clusters<-clusters[order(clusters)]
      }else{
        clusters<-unique(data.object()@plot$aka2$Cluster)
      }

      oraResultsTop<-do.call("rbind", lapply(clusters, topORA, ora=ora.objectALL(), top=as.integer(input$top)))
      ora.object<<-reactive({oraResultsTop})
      oraFiltered<-as.data.frame(ora.object())
      oraFiltered<-oraFiltered[oraFiltered$Cluster %in% input$optionCluster,]
      #output$oraInfo<-renderDataTable(datatable(oraFiltered, rownames = FALSE, options = list(searching=TRUE, pageLength=10, rowCallback = JS(underlineCellsCols(c(1),c(0))), columnDefs = list(list(targets = c(7), visible = FALSE))))
      #                                %>% formatStyle(colnames(oraFiltered)[1], color = "blue"))
      output$oraInfo<-renderDataTable(datatable(oraFiltered, rownames = FALSE, options = list(searching=TRUE, pageLength=10, rowCallback = JS(underlineCellsCols(c(1),c(0)))))
                                      %>% formatStyle(colnames(oraFiltered)[1], color = "blue"))

    }

  })

  ###########################
  ###########################
  #FILTER BY CLUSTER -> ORA and table of GO
  observe({

    if((is.null(input$optionCluster)) && (!is.null(data.object()))){
      if ("Uncheck all" %in% input$optionCluster) {
        updateCheckboxGroupInput(session, "optionCluster", selected = character(0))
      }

      #if any option is selected print NULL
      tableFiltered<-annotations()
      tableFiltered<-tableFiltered[tableFiltered$Cluster %in% c("Cluster_X"),]

      if("Term" %in% colnames(tableFiltered)){
        output$dataInfo <-renderDataTable(datatable(tableFiltered, rownames = FALSE, options = list(autoWidth=TRUE, scrollX=TRUE, scrollY=400, searching=TRUE, pageLength=25, rowCallback = JS(underlineCellsCols(c(1),c(0)))))
                                          %>% formatStyle(colnames(tableFiltered)[1], color = "blue"))
      }else{
        output$dataInfo <-renderDataTable(datatable(tableFiltered, rownames = FALSE, options = list(autoWidth=TRUE, scrollX=TRUE, scrollY=400, searching=TRUE, pageLength=25)))
      }

      oraFiltered<-as.data.frame(ora.object())
      oraFiltered<-oraFiltered[oraFiltered$Cluster %in% c("Cluster_X"),]
      output$oraInfo<-renderDataTable(datatable(oraFiltered, rownames = FALSE, options = list(searching=TRUE, pageLength=10, rowCallback = JS(underlineCellsCols(c(1),c(0)))))
                                      %>% formatStyle(colnames(oraFiltered)[1], color = "blue"))

    }
  })

  observeEvent(input$optionCluster,{

    if ("Uncheck all" %in% input$optionCluster) {
      updateCheckboxGroupInput(session, "optionCluster", selected = character(0))
    }

    tableFiltered<-annotations()
    tableFiltered<-tableFiltered[tableFiltered$Cluster %in% input$optionCluster,]

    if("Term" %in% colnames(tableFiltered)){
      if(uniquepathwaysperform()){
        output$dataInfo <-renderDataTable(datatable(tableFiltered, rownames = FALSE, options = list(autoWidth=TRUE, scrollX=TRUE, scrollY=400, searching=TRUE, pageLength=25, rowCallback = JS(underlineCellsCols(c(1),c(0)))))
                                          %>% formatStyle(colnames(tableFiltered)[5], backgroundColor = styleEqual(c(paste("Cluster_", names(data.object()@plot$aka3Unique$Cluster), sep="")), c(data.object()@plot$aka3Unique$Cluster)))
                                          %>% formatStyle(colnames(tableFiltered)[1], color = "blue"))
      }else{
        output$dataInfo <-renderDataTable(datatable(tableFiltered, rownames = FALSE, options = list(autoWidth=TRUE, scrollX=TRUE, scrollY=400, searching=TRUE, pageLength=25, rowCallback = JS(underlineCellsCols(c(1),c(0)))))
                                          %>% formatStyle(colnames(tableFiltered)[5], backgroundColor = styleEqual(c(paste("Cluster_", names(data.object()@plot$aka3$Cluster), sep="")), c(data.object()@plot$aka3$Cluster)))
                                          %>% formatStyle(colnames(tableFiltered)[1], color = "blue"))
      }
    }else{
      output$dataInfo <-renderDataTable(datatable(tableFiltered, rownames = FALSE, options = list(autoWidth=TRUE, scrollX=TRUE, scrollY=400, searching=TRUE, pageLength=25)))
    }

    oraFiltered<-as.data.frame(ora.object())
    oraFiltered<-oraFiltered[oraFiltered$Cluster %in% input$optionCluster,]
    if(uniquepathwaysperform()){
      output$oraInfo<-renderDataTable(datatable(oraFiltered, rownames = FALSE, options = list(searching=TRUE, pageLength=10, rowCallback = JS(underlineCellsCols(c(1),c(0)))))
                                      %>% formatStyle(colnames(oraFiltered)[9], backgroundColor = styleEqual(c(paste("Cluster_", names(data.object()@plot$aka3Unique$Cluster), sep="")), c(data.object()@plot$aka3Unique$Cluster)))
                                      %>% formatStyle(colnames(oraFiltered)[1], color = "blue"))
    }else{
      output$oraInfo<-renderDataTable(datatable(oraFiltered, rownames = FALSE, options = list(searching=TRUE, pageLength=10, rowCallback = JS(underlineCellsCols(c(1),c(0)))))
                                      %>% formatStyle(colnames(oraFiltered)[9], backgroundColor = styleEqual(c(paste("Cluster_", names(data.object()@plot$aka3$Cluster), sep="")), c(data.object()@plot$aka3$Cluster)))
                                      %>% formatStyle(colnames(oraFiltered)[1], color = "blue"))
    }

  })


  ###########################
  ###########################
  #RECALCULATING CLUSTERING
  observeEvent(input$recalculateClustering, {
    show_modal_progress_line(value=0.3, text="Reclustering (1/2)")
    update_modal_progress(value=0.4, text="Reclustering (1/2)") # update progress bar value
    new.object<-ClusterGeneSets(Object = data.object(),
                                clusters = input$cluster,
                                method = "kmeans", order="cluster")
    if(uniquepathwaysperform()){
      new.object <- transferinfo(new.object)
      clusterchoices<-paste("Cluster_", names(new.object@plot$aka3Unique$Cluster), sep="")
      clusters<-unique(new.object@plot$aka2Unique$Cluster)
      clusters<-clusters[order(clusters)]
    }else{
      #change names of the cluster with Cluster_1
      new.object@Data[[1]]$cluster<-paste(rep("Cluster_",dim(new.object@Data[[1]]["cluster"])[1]), new.object@Data[[1]]$cluster, sep = "")
      clusterchoices<-paste("Cluster_", names(new.object@plot$aka3$Cluster), sep="")
      clusters<-unique(new.object@plot$aka2$Cluster)
    }
    clusterchoices<-clusterchoices[order(clusterchoices)]
    data.object<<-reactive({new.object})

    update_modal_progress(value=0.5, text="Reclustering (1/2)")


    #show new data results
    annotationsUpdate<-annotations()
    annotationsUpdate$Cluster<-data.object()@Data[[1]]$cluster
    annotations<<-reactive({annotationsUpdate})

    update_modal_progress(value=0.8, text="ORA per cluster (2/2)") # update progress bar value

    genesFinal<-calculategenestable(data.object(), uniquepathwaysperform())
    if(unique(as.character(data.object()@metadata[,"structure"]))!="SYMBOL"){
      genesFinal<-combineSymbolInfo(genesFinal, genes.symbolinfo())
    }
    #recalculate ORA per cluster
    oraResults<-do.call("rbind", lapply(clusters, calculateORA, object=data.object(), uniquePathways=uniquepathwaysperform()))
    ora.objectALL<<-reactive({oraResults})

    oraResultsTop<-do.call("rbind", lapply(clusters, topORA, ora=oraResults, top=as.integer(input$top)))
    ora.object<<-reactive({oraResultsTop})

    genesfinal.info<<-reactive({genesFinal})
    output$genesInfo <-renderDataTable(datatable(genesFinal, rownames = FALSE, filter="top", options = list(pageLength=10, rowCallback = JS(underlineCellsCols(c(1),c(0)))))
                                       %>% formatStyle(colnames(genesFinal)[1], color = "blue"))


    #plot new heatmap
    keywords_ora<-calculateKeywordsORA(ora.objectALL())
    if(wordcloudperform()){
      wordcloudinfo<<-reactive({PerformORAGOperCluster(cluster.object, uniquePathways=uniquepathwaysperform(), clusterIndependent=FALSE)})
    }
    heatmapoutput<<-reactive({PlotGeneSets(Object = data.object(), doORA = T, wordcloud = wordcloudperform(), wordclouds=wordcloudinfo(), uniquePathways=uniquepathwaysperform(), keywords_ora=keywords_ora)})
    output$heatmap<-renderPlot(heatmapoutput())

    update_modal_progress(value=0.9, text="ORA per cluster (2/2)") # update progress bar value

    #put check box of cluster
    output$checkboxCluster <- renderUI({
      clusterchoices_all<-c("Uncheck all", clusterchoices)
      checkboxGroupInput("optionCluster", label = "", choices=clusterchoices_all, selected = clusterchoices, inline=TRUE)
    })
    output$barplot<-renderPlot(createBarPlot(data.object(), clusterchoices[1], uniquePathways=uniquepathwaysperform()))

    updateSelectInput(session, "breakup.cluster", choices = clusterchoices)
    updateSelectInput(session, "summary.cluster", choices = clusterchoices)

    #reset tissue results
    shinyjs::show("calculateTissueEnrichment")
    hideTab(inputId = "tabs", target = "Tissue")
    output$tissueOutput<-renderDataTable(NULL)
    output$tissueIntro<-renderUI({tags$p("Choose specific tissues for performing enrichment analysis (the tissues were selected from", HTML("<a href='https://gtexportal.org/home/' target='_blank'>here</a>"), "):")})
    remove_modal_progress()

  })
  ###########################
  ###########################
  #PERFORM ORA PER GENES
  #shinyjs::hide("downloadORAgenes")
  oragenesfiltered<-reactive({})
  observeEvent(input$performORAgenes, {
    numbergenes<-dim(genesfinal.info()[input$genesInfo_rows_all,])[1]
    if(numbergenes>0){
      show_modal_progress_line(value=0.3, text="Performing ORA")
      update_modal_progress(value=0.5, text="Performing ORA")
      oragenestable<-enrichGO(gene = genesfinal.info()[input$genesInfo_rows_all,1], OrgDb = data.object()@metadata[1,"organism"], keyType = "SYMBOL", ont="BP", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable=TRUE)
      update_modal_progress(value=0.8, text="Performing ORA")
      oragenesfiltered<<-reactive({oragenestable})
      runjs('$("#downloadORAgenes")[0].click();')
      remove_modal_progress()
    }else{
      shinyalert("Error", "No genes slected", type = "error")
    }
  })

  output$downloadORAgenes<-downloadHandler(
    filename = paste("ORAgenes_",Sys.time(),".csv", sep=""),
    content = function(file) {
      write.csv(oragenesfiltered(), file, row.names=FALSE)
    }
  )

  ###########################
  ###########################
  #CALCULATE TISSUE ENRICHMENT CLUSTERING
  observeEvent(input$calculateTissueEnrichment, {
    if (is.null(input$tissuesselected)){
      shinyalert(title = "Missing tissue", text = "Please, select 1 or more tissues.", type = "error")
    }else{
      show_modal_progress_line(value=0.3, text="Calculate tissue enrichment")
      update_modal_progress(value=0.5, text="Calculate tissue enrichment")
      load("databases/TissueLocalDatabase15.RData")
      load("databases/dic.rda")
      #shinyjs::hide("calculateTissueEnrichment")

      localDatabase15<-localDatabase15[input$tissuesselected]
      tissue.object <-TissueExpressionPerGeneSet(data.object(), localDatabase = localDatabase15, dic=dic, uniquePathways=uniquepathwaysperform(), clusterIndependent=FALSE)
      tissue.object@dfTissue<-round(tissue.object@dfTissue,3)
      tissue.object@dfTissue<-na_replace(tissue.object@dfTissue,0)
      data.object<<-reactive({tissue.object})

      update_modal_progress(value=0.7, text="Calculate tissue enrichment")
      tissueData<-data.frame(rownames(tissue.object@dfTissue), tissue.object@dfTissue)
      colnames(tissueData)[1]<-"Tissue"
      output$tissueOutput<-renderDataTable(datatable(tissueData, rownames = FALSE, options = list(searching=TRUE, pageLength=10)))
      #shinyjs::hide("calculateTissueEnrichment")

      output$tissueIntro<-renderUI({tags$p("Tissue enrichment results per cluster (Cluster_number..number of genes):")})
      showTab(inputId = "tabs", target = "Tissue", select=TRUE)

      update_modal_progress(value=0.9)
      #plot tissue plots
      tissueoutput<<-reactive({PlotTissueExpression(tissue.object, all = FALSE, uniquePathways=uniquepathwaysperform(), clusterIndependent=FALSE)})
      output$tissuePlot<-renderPlot(tissueoutput())



      rm(localDatabase15)
      rm(dic)
      remove_modal_progress()
    }

  })

  ###########################
  #RUN INDEPENDENT ANALYSIS
  observeEvent(input$runIndependent, {
    shinyjs::hide("runIndependent")
    show_modal_progress_line(value=0.3, text="Run seriation-based analysis")
    update_modal_progress(value=0.5, text="Run seriation-based analysis")

    output$independentIntro<-renderUI({tags$p("")})
    independent.object<-ClusterIndependentGeneSet(data.object(), uniquePathways = uniquepathwaysperform(),  nPathways = "optimal")
    data.object<<-reactive({independent.object})

    update_modal_progress(value=0.6, text="Run seriation-based analysis")
    if(checkGO(data.object())){
      wordcloudindependent<-TRUE
    }else{
      wordcloudindependent<-FALSE
    }
    #pathway information per grouping
    resultindependent<-ShowPathwayCluster(independent.object,  uniquePathways = uniquepathwaysperform())
    annotated_resultindependent <- lapply(seq_along(resultindependent), function(i,uniquePathways = uniquepathwaysperform(), goids=wordcloudindependent) {
      go_list <- resultindependent[[i]]
      if(uniquePathways){
        if(goids){
          subset_df <- annotations()[annotations()$ID %in% go_list, c("ID","Term","Group")]
          annotated_terms <- merge(subset_df, data.frame(ID = go_list, Cluster = paste("Cluster_",i, sep="")), by = "ID", all.y = TRUE)
        }else{
          subset_df <- annotations()[annotations()$Pathways %in% go_list, c("Pathways","Group")]
          annotated_terms <- merge(subset_df, data.frame(Pathways = go_list, Cluster = paste("Cluster_",i, sep="")), by = "Pathways", all.y = TRUE)
        }
      }else{
        if(goids){
          subset_df <- annotations()[paste(annotations()$Group, annotations()$ID, sep="_") %in% go_list, c("ID","Term","Group")]
          subset_df<-cbind(paste(subset_df$Group, subset_df$ID, sep="_"), subset_df)
          colnames(subset_df)[1]<-"rowname"
          annotated_terms <- merge(subset_df, data.frame(ID = go_list, Cluster = paste("Cluster_",i, sep="")), by.x="rowname", by.y = "ID", all.y = TRUE)
          annotated_terms<-annotated_terms[,-1]
        }else{
          subset_df <- annotations()[paste(annotations()$Group, annotations()$Pathways, sep="_") %in% go_list, c("ID","Term","Group")]
          subset_df<-cbind(paste(subset_df$Group, subset_df$ID, sep="_"), subset_df)
          colnames(subset_df)[1]<-"rowname"
          annotated_terms <- merge(subset_df, data.frame(Pathways = go_list, Cluster = paste("Cluster_",i, sep="")), by.x="rowname", by.y = "Pathways", all.y = TRUE)
          annotated_terms<-annotated_terms[,-1]
        }
      }

      return(annotated_terms)
    })
    combined_resultindependent<-do.call(rbind, annotated_resultindependent)
    groupchoices<-unique(combined_resultindependent$Cluster)
    independent.info<<-reactive({combined_resultindependent})
    updateSelectInput(session, "summary.independentgroup", choices = groupchoices)
    shinyjs::show("summary.independentgroup")

    if(checkGO(data.object())){
      output$independentOutput <-renderDataTable(datatable(combined_resultindependent, rownames = FALSE, options = list(autoWidth=TRUE, scrollX=TRUE, scrollY=400, searching=TRUE, pageLength=25, rowCallback = JS(underlineCellsCols(c(1),c(0)))))
                                                 %>% formatStyle(colnames(combined_resultindependent)[4], backgroundColor = styleEqual(c(paste("Cluster_", names(data.object()@plot[[6]]$Cluster)[names(data.object()@plot[[6]]$Cluster)!=0], sep="")), c(data.object()@plot[[6]]$Cluster[data.object()@plot[[6]]$Cluster!="white"])))
                                                 %>% formatStyle(colnames(combined_resultindependent)[1], color = "blue"))
    }else{
      output$independentOutput <-renderDataTable(datatable(combined_resultindependent, rownames = FALSE, options = list(autoWidth=TRUE, scrollX=TRUE, scrollY=400, searching=TRUE, pageLength=25))
                                                 %>% formatStyle(colnames(combined_resultindependent)[4], backgroundColor = styleEqual(c(paste("Cluster_", names(data.object()@plot[[6]]$Cluster)[names(data.object()@plot[[6]]$Cluster)!=0], sep="")), c(data.object()@plot[[6]]$Cluster[data.object()@plot[[6]]$Cluster!="white"]))))
    }

    genesFinalIndependent<-calculategenestableinde(data.object(),combined_resultindependent)
    if(unique(as.character(data.object()@metadata[,"structure"]))!="SYMBOL"){
      genesFinalIndependent<-combineSymbolInfo(genesFinalIndependent, genes.symbolinfo())
    }
    output$genesInfoIndependent<-renderDataTable(datatable(genesFinalIndependent, rownames = FALSE, filter="top", options = list(autoWidth=TRUE, scrollX=TRUE, scrollY=400, searching=TRUE, dom='lrtip', pageLength=10, rowCallback = JS(underlineCellsCols(c(1),c(0)))))
                                                 %>% formatStyle(colnames(genesFinalIndependent)[1], color = "blue"))

    #output$independentOutput<-renderDataTable(datatable(combined_resultindependent, rownames = FALSE, options = list(searching=TRUE, pageLength=10)))
    output$barplotIndependent<-renderPlot(createBarPlotIndependent(data.object(), combined_resultindependent, groupchoices[1], uniquePathways=uniquepathwaysperform()))

    #Calculate ORA in independent results
    oraResultsIndependent<-do.call("rbind", lapply(groupchoices, calculateORAindependent, object=data.object(), independent_info=combined_resultindependent, uniquePathways=uniquepathwaysperform()))
    independent.infoORA<<-reactive({oraResultsIndependent})
    output$independentORAOutput<-renderDataTable(datatable(oraResultsIndependent, rownames = FALSE, options = list(autoWidth=TRUE, scrollX=TRUE, scrollY=400, searching=TRUE, pageLength=10, rowCallback = JS(underlineCellsCols(c(1),c(0)))))
                                                 %>% formatStyle(colnames(oraResultsIndependent)[9], backgroundColor = styleEqual(c(paste("Cluster_", names(data.object()@plot[[6]]$Cluster)[names(data.object()@plot[[6]]$Cluster)!=0], sep="")), c(data.object()@plot[[6]]$Cluster[data.object()@plot[[6]]$Cluster!="white"])))
                                                 %>% formatStyle(colnames(oraResultsIndependent)[1], color = "blue"))

    update_modal_progress(value=0.8, text="Run seriation-based analysis")
    #create plot
    keywords_ora_independent<-calculateKeywordsORA(independent.infoORA())
    independentoutput<<-reactive({PlotPathwayCluster(independent.object, doORA = TRUE, wordcloud = FALSE, uniquePathways = uniquepathwaysperform(), keywords_ora_inde=keywords_ora_independent)})
    output$independentPlot<-renderPlot(independentoutput())
    update_modal_progress(value=0.9, text="Run seriation-based analysis")
    remove_modal_progress()
    showTab(inputId = "tabs", target = "Heatmap_S", select=TRUE)
    if((data.object()@metadata[1,"organism"]!="org.Mm.eg.db")){
      #for mouse data hide tissue enrichment analysis
      showTab(inputId = "tabsindependent", target = "Tissue")
    }
    showTab(inputId = "tabsindependent", target = "ORA")
    showTab(inputId = "tabsindependent", target = "Genes")
  })

  ###########################
  ###########################
  #CALCULATE TISSUE ENRICHMENT INDEPENDENT
  observeEvent(input$calculateTissueEnrichmentIndependent, {
    if (is.null(input$tissuesselectedindependent)){
      shinyalert(title = "Missing tissue", text = "Please, select 1 or more tissues.", type = "error")
    }else{
      show_modal_progress_line(value=0.3, text="Calculate tissue enrichment")
      update_modal_progress(value=0.5, text="Calculate tissue enrichment")
      load("databases/TissueLocalDatabase15.RData")
      load("databases/dic.rda")

      localDatabase15<-localDatabase15[input$tissuesselectedindependent]
      tissue.objectindependent <- TissueExpressionPerGeneSet(data.object(), localDatabase = localDatabase15, dic=dic, uniquePathways=uniquepathwaysperform(), clusterIndependent=TRUE)
      tissue.objectindependent@dfTissueIndependent<-round(tissue.objectindependent@dfTissueIndependent,3)
      tissue.objectindependent@dfTissueIndependent<-na_replace(tissue.objectindependent@dfTissueIndependent,0)
      data.object<<-reactive({tissue.objectindependent})

      update_modal_progress(value=0.7, text="Calculate tissue enrichment")
      tissueDataIndependent<-data.frame(rownames(tissue.objectindependent@dfTissueIndependent), tissue.objectindependent@dfTissueIndependent)
      colnames(tissueDataIndependent)[1]<-"Tissue"
      output$tissueOutputIndependent<-renderDataTable(datatable(tissueDataIndependent, rownames = FALSE, options = list(autoWidth=TRUE, scrollX=TRUE, scrollY=400, searching=TRUE, pageLength=10)))

      #output$tissueIntro<-renderUI({tags$p("Tissue enrichment results per cluster (Cluster_number..number of genes):")})

      showTab(inputId = "tabs", target = "Tissue_S", select=TRUE)

      update_modal_progress(value=0.9)
      #plot tissue plots
      tissueoutputindependent<<-reactive({PlotTissueExpression(tissue.objectindependent, all = FALSE, uniquePathways=uniquepathwaysperform(), clusterIndependent=TRUE)})
      output$tissuePlotIndependent<-renderPlot(tissueoutputindependent())

      rm(localDatabase15)
      rm(dic)
      remove_modal_progress()
    }

  })


  ###########################
  #BREAK UP CLUSTER
  observeEvent(input$breakup, {

    show_modal_progress_line(text="Calculate subclusters (1/2)")
    update_modal_progress(value=0.2)
    if (input$nbreakup.cluster=="Automatic"){
      new.object<-data.object()

      if(uniquepathwaysperform()){
        RR<-new.object@DataPathways.RR
        breakup.cluster<-sub("Cluster_","",input$breakup.cluster)
        new.object@DataPathways.RR <- new.object@DataPathways.RR[new.object@plot$aka2Unique$Cluster==breakup.cluster,new.object@plot$aka2Unique$Cluster==breakup.cluster]
      }else{
        new.object@Data.RR <- new.object@Data.RR[new.object@Data[[1]]$cluster==input$breakup.cluster,new.object@Data[[1]]$cluster==input$breakup.cluster]
      }
      #calculate optimal cluster
      optimal.object<-OptimalGeneSets_2(object = new.object, method = "silhouette", max_cluster= 10, cluster_method = "kmeans", main= "", uniquePathways=uniquepathwaysperform())
      nclust<-optimal.object$data
      optimalCluster<-as.numeric(nclust$clusters[which.max(nclust$y)])
    }else{
      optimalCluster<-input$nbreakup.cluster
    }
    update_modal_progress(value=0.5)
    data_break_up <- BreakUpCluster(Object = data.object(), breakup.cluster = input$breakup.cluster, sub.cluster=optimalCluster, uniquePathways=uniquepathwaysperform())
    if(uniquepathwaysperform()){
      data_break_up <- transferinfo(data_break_up)
      clusterchoices<-paste("Cluster_", names(data_break_up@plot$aka3Unique$Cluster), sep="")
      clusterchoices<-clusterchoices[order(clusterchoices)]
      clusters<-unique(data_break_up@plot$aka2Unique$Cluster)
      clusters<-clusters[order(clusters)]
    }else{
      clusterchoices<-paste("Cluster_", names(data_break_up@plot$aka3$Cluster), sep="")
      clusters<-unique(data_break_up@plot$aka2$Cluster)
    }

    data.object<<-reactive({data_break_up})

    annotationsUpdate<-annotations()
    if("Term" %in% colnames(annotations())){
      #wordcloudperform=T
      annotationsUpdate$Cluster<-data_break_up@Data[[1]]$cluster[which(data_break_up@Data[[1]]$Pathways %in% annotationsUpdate$ID)]
    }else{
      #wordcloudperform=F
      annotationsUpdate$Cluster<-data_break_up@Data[[1]]$cluster[which(data_break_up@Data[[1]]$Pathways %in% annotationsUpdate$Pathways)]
    }

    annotations<<-reactive({annotationsUpdate})

    genesFinal<-calculategenestable(data.object(), uniquepathwaysperform())
    if(unique(as.character(data.object()@metadata[,"structure"]))!="SYMBOL"){
      genesFinal<-combineSymbolInfo(genesFinal, genes.symbolinfo())
    }
    #genes info
    genesfinal.info<<-reactive({genesFinal})
    output$genesInfo <-renderDataTable(datatable(genesFinal, rownames = FALSE, filter="top", options = list(dom='lrtip', pageLength=10, rowCallback = JS(underlineCellsCols(c(1),c(0)))))
                                       %>% formatStyle(colnames(genesFinal)[1], color = "blue"))

    #put check box of cluster
    output$checkboxCluster <- renderUI({
      clusterchoices_all<-c("Uncheck all", clusterchoices)
      checkboxGroupInput("optionCluster", label = "", choices=clusterchoices_all, selected = clusterchoices, inline=TRUE)
    })

    update_modal_progress(value=0.7)

    output$barplot<-renderPlot(createBarPlot(data.object(), clusterchoices[1], uniquePathways=uniquepathwaysperform()))

    updateSelectInput(session, "breakup.cluster", choices = clusterchoices)
    updateSelectInput(session, "summary.cluster", choices = clusterchoices)

    #showModal(modalDialog("ORA per newcluster (2/2)", footer=NULL))
    update_modal_progress(text="ORA per new subcluster (2/2)", value=0.8)
    #calculate ORA
    clusterpattern<-paste(sub("Cluster_","",input$breakup.cluster), ".", sep="")
    oraTemp<-ora.objectALL()
    oraTemp<-oraTemp[oraTemp$Cluster!=input$breakup.cluster,]
    oraResults<-do.call("rbind", lapply(grep(clusterpattern, unique(clusters), value = TRUE), calculateORA, object=data.object(), uniquePathways=uniquepathwaysperform()))
    oraResults<-rbind(oraResults, oraTemp)
    oraResults<-oraResults[order(oraResults$Cluster),]
    ora.objectALL<<-reactive({oraResults})

    oraTopTemp<-ora.object()
    oraTopTemp<-oraTopTemp[oraTopTemp$Cluster!=input$breakup.cluster,]
    oraResultsTop<-do.call("rbind", lapply(grep(clusterpattern,unique(clusters), value = TRUE), topORA, ora=oraResults, top=as.integer(input$top)))
    oraResultsTop<-rbind(oraResultsTop, oraTopTemp)
    oraResultsTop<-oraResultsTop[order(oraResultsTop$Cluster),]
    ora.object<<-reactive({oraResultsTop})

    #create plots
    if(wordcloudperform()){
      wordcloudinfo<<-reactive({PerformORAGOperCluster(data.object(), uniquePathways=uniquepathwaysperform(), clusterIndependent=FALSE)})
    }
    keywords_ora<-calculateKeywordsORA(ora.objectALL())
    heatmapoutput<<-reactive({PlotGeneSets(Object = data.object(), doORA = T, wordcloud = wordcloudperform(), wordclouds=wordcloudinfo(), uniquePathways=uniquepathwaysperform(), keywords_ora=keywords_ora)})
    output$heatmap<-renderPlot(heatmapoutput())

    #reset tissue results
    shinyjs::show("calculateTissueEnrichment")
    hideTab(inputId = "tabs", target = "Tissue")
    output$tissueOutput<-renderDataTable(NULL)
    output$tissueIntro<-renderUI({tags$p("Choose specific tissues for performing enrichment analysis (the tissues were selected from", HTML("<a href='https://gtexportal.org/home/' target='_blank'>here</a>"), "):")})
    remove_modal_progress()

  })

  ###########################
  #UPLOAD RESULTS
  observeEvent(input$uploadResult, {

    if(is.null(input$fileRdata$datapath)){
      shinyalert(title = "Missing file", text = "Please, insert Rdata file.", type = "error")
    }else{

      if((input$uploadfrom=="package")&&((is.null(input$uploadfilters))||(is.null(input$uploadnameobject)))){
        shinyalert(title = "Missing input data", text = "Please, specify object name and pathays included.", type = "error")
      }else{
        show_modal_progress_line(text="Uploading")
        load(input$fileRdata$datapath)

        #load data
        if(input$uploadfrom=="shiny"){
          data.object<<-reactive({data_to_save})
          annotations<<-reactive({annotations_to_save})
          ora.objectALL<<-reactive({ora_to_save})
          uniquepathwaysperform<<-reactive({unique_pathways})

        }else{
          data.object<<-reactive({get(input$uploadnameobject)})
          if(input$uploadfilters=="uniquepathways"){
            uniquepathwaysperform(TRUE)
          }else{
            uniquepathwaysperform(FALSE)
          }
          if(uniquepathwaysperform()){
            updateobject <- transferinfo(data.object())
          }else{
            updateobject<-data.object()
            updateobject@Data[[1]]$cluster<-paste(rep("Cluster_",dim(updateobject@Data[[1]]["cluster"])[1]), updateobject@Data[[1]]$cluster, sep = "")
          }
          if(!is.character(updateobject@Data[[1]]$Ratio[1])){
            updateobject@Data[[1]]$Ratio<-round(updateobject@Data[[1]]$Ratio,3)
          }
          data.object<<-reactive({updateobject})
        }

        update_modal_progress(value=0.2)
        #########
        if(uniquepathwaysperform()){
          clusterchoices<-paste("Cluster_", names(data.object()@plot$aka3Unique$Cluster), sep="")
          clusterchoices<-clusterchoices[order(clusterchoices)]
          clusters<-unique(data.object()@plot$aka2Unique$Cluster)
          clusters<-clusters[order(clusters)]
        }else{
          clusterchoices<-paste("Cluster_", names(data.object()@plot$aka3$Cluster), sep="")
          clusters<-unique(data.object()@plot$aka2$Cluster)
        }

        #put check box of cluster
        output$checkboxCluster <- renderUI({
          clusterchoices_all<-c("Uncheck all", clusterchoices)
          checkboxGroupInput("optionCluster", label = "", choices=clusterchoices_all, selected = clusterchoices, inline=TRUE)
        })

        #genes info
        genesFinal<-calculategenestable(data.object(), uniquepathwaysperform())
        if(unique(as.character(data.object()@metadata[,"structure"]))!="SYMBOL"){
          symbolinfo<-calculateSymbol(data.object(), genesFinal)
          genes.symbolinfo<<-reactive({symbolinfo})
          genesFinal<-combineSymbolInfo(genesFinal, genes.symbolinfo())
        }
        genesfinal.info<<-reactive({genesFinal})
        output$genesInfo <-renderDataTable(datatable(genesFinal, rownames = FALSE, filter="top", options = list(dom='lrtip', pageLength=10, rowCallback = JS(underlineCellsCols(c(1),c(0)))))
                                           %>% formatStyle(colnames(genesFinal)[1], color = "blue"))

        if(input$uploadfrom=="shiny"){
          #ora results
          oraResultsTop<-do.call("rbind", lapply(clusters, topORA, ora=ora.objectALL(), top=as.integer(input$top)))
          ora.object<<-reactive({oraResultsTop})
          if(!is.null(independent_to_save)){
            independent.info<<-reactive({independent_to_save})
            independent.infoORA<<-reactive({independent_ora_to_save})
            shinyjs::hide("runIndependent")
            output$independentIntro<-renderUI({tags$p("")})
            #independent
            groupchoices<-unique(independent.info()$Cluster)
            if(checkGO(data.object())){
              output$independentOutput <-renderDataTable(datatable(independent.info(), rownames = FALSE, options = list(autoWidth=TRUE, scrollX=TRUE, scrollY=400, searching=TRUE, pageLength=25, rowCallback = JS(underlineCellsCols(c(1),c(0)))))
                                                         %>% formatStyle(colnames(independent.info())[4], backgroundColor = styleEqual(c(paste("Cluster_", names(data.object()@plot[[6]]$Cluster)[names(data.object()@plot[[6]]$Cluster)!=0], sep="")), c(data.object()@plot[[6]]$Cluster[data.object()@plot[[6]]$Cluster!="white"])))
                                                         %>% formatStyle(colnames(independent.info())[1], color = "blue"))

            }else{
              output$independentOutput <-renderDataTable(datatable(independent.info(), rownames = FALSE, options = list(autoWidth=TRUE, scrollX=TRUE, scrollY=400, searching=TRUE, pageLength=25))
                                                         %>% formatStyle(colnames(independent.info())[4], backgroundColor = styleEqual(c(paste("Cluster_", names(data.object()@plot[[6]]$Cluster)[names(data.object()@plot[[6]]$Cluster)!=0], sep="")), c(data.object()@plot[[6]]$Cluster[data.object()@plot[[6]]$Cluster!="white"]))))
            }
            #independent ORA
            output$independentORAOutput<-renderDataTable(datatable(independent.infoORA(), rownames = FALSE, options = list(autoWidth=TRUE, scrollX=TRUE, scrollY=400, searching=TRUE, pageLength=10, rowCallback = JS(underlineCellsCols(c(1),c(0)))))
                                                         %>% formatStyle(colnames(independent.infoORA())[9], backgroundColor = styleEqual(c(paste("Cluster_", names(data.object()@plot[[6]]$Cluster)[names(data.object()@plot[[6]]$Cluster)!=0], sep="")), c(data.object()@plot[[6]]$Cluster[data.object()@plot[[6]]$Cluster!="white"])))
                                                         %>% formatStyle(colnames(independent.infoORA())[1], color = "blue"))

            genesFinalIndependent<-calculategenestableinde(data.object(),independent.info())
            if(unique(as.character(data.object()@metadata[,"structure"]))!="SYMBOL"){
              genesFinalIndependent<-combineSymbolInfo(genesFinalIndependent, genes.symbolinfo())
            }
            output$genesInfoIndependent<-renderDataTable(datatable(genesFinalIndependent, rownames = FALSE, filter="top", options = list(autoWidth=TRUE, scrollX=TRUE, scrollY=400, searching=TRUE, dom='lrtip', pageLength=10, rowCallback = JS(underlineCellsCols(c(1),c(0)))))
                                                         %>% formatStyle(colnames(genesFinalIndependent)[1], color = "blue"))
            updateSelectInput(session, "summary.independentgroup", choices = groupchoices)
            shinyjs::show("summary.independentgroup")
            #create plot
            keywords_ora_independent<-calculateKeywordsORA(independent.infoORA())
            independentoutput<<-reactive({PlotPathwayCluster(data.object(), doORA = TRUE, wordcloud = FALSE, uniquePathways = uniquepathwaysperform(), keywords_ora_inde=keywords_ora_independent)})
            output$independentPlot<-renderPlot(independentoutput())
            if((data.object()@metadata[1,"organism"]=="org.Mm.eg.db")){
              #for mouse data hide tissue enrichment analysis
              hideTab(inputId = "tabsindependent", target = "Tissue")
            }

          }else{
            hideTab(inputId = "tabs", target = "Heatmap_S")
            hideTab(inputId = "tabs", target = "Tissue_S")
            hideTab(inputId = "tabsindependent", target = "ORA")
            hideTab(inputId = "tabsindependent", target = "Genes")
            hideTab(inputId = "tabsindependent", target = "Tissue")
          }
        }else{
          #R packageeee
          #calculate annotations
          if((is.null(annotations())) && (checkGO(data.object())==TRUE)){
            resultAnnotations<-as.data.frame(str_split_fixed(as.data.frame(unlist(lapply(data.object()@Data[[1]][,"Pathways"],getterm)))[,1],"__",3))
            colnames(resultAnnotations)<-c("ID", "Term", "Definition")
            resultAnnotationsFinal<-data.frame(resultAnnotations[,c(1,2)], data.object()@Data[[1]][c("Type","Groups","cluster","Pval","Ratio")], resultAnnotations[,3])
            colnames(resultAnnotationsFinal)[c(4,5,8)]<-c("Group", "Cluster", "Definition")
            annotations<<-reactive({resultAnnotationsFinal})
          }else{
            resultAnnotationsFinal<-data.frame(data.object()@Data[[1]][c("Pathways","Type","Groups","cluster","Ratio")])
            colnames(resultAnnotationsFinal)[c(3,4)]<-c("Group", "Cluster")
            annotations<<-reactive({resultAnnotationsFinal})
          }
          #calulate ORA in parallel
          oraResults<-do.call("rbind", lapply(clusters, calculateORA, object=data.object(), uniquePathways=uniquepathwaysperform()))
          ora.objectALL<<-reactive({oraResults})
          oraResultsTop<-do.call("rbind", lapply(clusters, topORA, ora=oraResults, top=as.integer(input$top)))
          ora.object<<-reactive({oraResultsTop})
          #add independent info
          if (!is.null(data.object()@cIndependentMethod)){
            shinyjs::hide("runIndependent")
            output$independentIntro<-renderUI({tags$p("")})
            if(uniquepathwaysperform()){
              temp.object<-data.object()
              temp.object@cIndependentMethod[[1]][[1]]<-temp.object@cIndependentMethod[[2]][[1]]
              data.object<<-reactive({temp.object})
            }
            if(checkGO(data.object())){
              wordcloudindependent<-TRUE
            }else{
              wordcloudindependent<-FALSE
            }
            resultindependent<-ShowPathwayCluster(data.object(), uniquePathways = uniquepathwaysperform())
            annotated_resultindependent <- lapply(seq_along(resultindependent), function(i,uniquePathways = uniquepathwaysperform(), goids=wordcloudindependent) {
              go_list <- resultindependent[[i]]
              if(uniquePathways){
                if(goids){
                  subset_df <- annotations()[annotations()$ID %in% go_list, c("ID","Term","Group")]
                  annotated_terms <- merge(subset_df, data.frame(ID = go_list, Cluster = paste("Cluster_",i, sep="")), by = "ID", all.y = TRUE)
                }else{
                  subset_df <- annotations()[annotations()$Pathways %in% go_list, c("Pathways","Group")]
                  annotated_terms <- merge(subset_df, data.frame(Pathways = go_list, Cluster = paste("Cluster_",i, sep="")), by = "Pathways", all.y = TRUE)
                }
              }else{
                if(goids){
                  subset_df <- annotations()[paste(annotations()$Group, annotations()$ID, sep="_") %in% go_list, c("ID","Term","Group")]
                  subset_df<-cbind(paste(subset_df$Group, subset_df$ID, sep="_"), subset_df)
                  colnames(subset_df)[1]<-"rowname"
                  annotated_terms <- merge(subset_df, data.frame(ID = go_list, Cluster = paste("Cluster_",i, sep="")), by.x="rowname", by.y = "ID", all.y = TRUE)
                  annotated_terms<-annotated_terms[,-1]
                }else{
                  subset_df <- annotations()[paste(annotations()$Group, annotations()$Pathways, sep="_") %in% go_list, c("ID","Term","Group")]
                  subset_df<-cbind(paste(subset_df$Group, subset_df$ID, sep="_"), subset_df)
                  colnames(subset_df)[1]<-"rowname"
                  annotated_terms <- merge(subset_df, data.frame(Pathways = go_list, Cluster = paste("Cluster_",i, sep="")), by.x="rowname", by.y = "Pathways", all.y = TRUE)
                  annotated_terms<-annotated_terms[,-1]
                }
              }

              return(annotated_terms)
            })
            combined_resultindependent<-do.call(rbind, annotated_resultindependent)
            groupchoices<-unique(combined_resultindependent$Cluster)

            independent.info<<-reactive({combined_resultindependent})
            output$independentOutput<-renderDataTable((datatable(combined_resultindependent, rownames = FALSE, options = list(searching=TRUE, pageLength=10)))
                                                      %>% formatStyle(colnames(independent.info())[4], backgroundColor = styleEqual(c(paste("Cluster_", names(data.object()@plot[[6]]$Cluster)[names(data.object()@plot[[6]]$Cluster)!=0], sep="")), c(data.object()@plot[[6]]$Cluster[data.object()@plot[[6]]$Cluster!="white"]))))

            genesFinalIndependent<-calculategenestableinde(data.object(),independent.info())
            if(unique(as.character(data.object()@metadata[,"structure"]))!="SYMBOL"){
              genesFinalIndependent<-combineSymbolInfo(genesFinalIndependent, genes.symbolinfo())
            }
            output$genesInfoIndependent<-renderDataTable(datatable(genesFinalIndependent, rownames = FALSE, filter="top", options = list(autoWidth=TRUE, scrollX=TRUE, scrollY=400, searching=TRUE, dom='lrtip', pageLength=10, rowCallback = JS(underlineCellsCols(c(1),c(0)))))
                                                         %>% formatStyle(colnames(genesFinalIndependent)[1], color = "blue"))
            oraResultsIndependent<-do.call("rbind", lapply(groupchoices, calculateORAindependent, object=data.object(), independent_info=combined_resultindependent, uniquePathways=uniquepathwaysperform()))
            independent.infoORA<<-reactive({oraResultsIndependent})
            output$independentORAOutput<-renderDataTable(datatable(oraResultsIndependent, rownames = FALSE, options = list(autoWidth=TRUE, scrollX=TRUE, scrollY=400, searching=TRUE, pageLength=10, rowCallback = JS(underlineCellsCols(c(1),c(0)))))
                                                         %>% formatStyle(colnames(oraResultsIndependent)[9], backgroundColor = styleEqual(c(paste("Cluster_", names(data.object()@plot[[6]]$Cluster)[names(data.object()@plot[[6]]$Cluster)!=0], sep="")), c(data.object()@plot[[6]]$Cluster[data.object()@plot[[6]]$Cluster!="white"])))
                                                         %>% formatStyle(colnames(oraResultsIndependent)[1], color = "blue"))


            updateSelectInput(session, "summary.independentgroup", choices = groupchoices)
            shinyjs::show("summary.independentgroup")

            #create plot
            keywords_ora_independent<-calculateKeywordsORA(independent.infoORA())
            independentoutput<<-reactive({PlotPathwayCluster(data.object(), doORA = TRUE, wordcloud = FALSE, uniquePathways = uniquepathwaysperform(), keywords_ora_inde=keywords_ora_independent)})
            output$independentPlot<-renderPlot(independentoutput())
          }else{
            hideTab(inputId = "tabs", target = "Heatmap_S")
          }
        }

        update_modal_progress(value=0.5)
        #tissue independent
        if(!is.null(data.object()@dfTissueIndependent)){
          temp.object<-data.object()
          if(input$uploadfrom=="package"){
            if(uniquepathwaysperform()){
              temp.object@dfTissueIndependent<- temp.object@dfTissueIndependent$Pathway
            }else{
              temp.object@dfTissueIndependent<- temp.object@dfTissueIndependent$Geneset
            }
          }
          temp.object@dfTissueIndependent<-round(temp.object@dfTissueIndependent,3)
          temp.object@dfTissueIndependent<-na_replace(temp.object@dfTissueIndependent,0)
          data.object<<-reactive({temp.object})
          tissueDataIndependent<-data.frame(rownames(data.object()@dfTissueIndependent), data.object()@dfTissueIndependent)
          colnames(tissueDataIndependent)[1]<-"Tissue"
          output$tissueOutputIndependent<-renderDataTable(datatable(tissueDataIndependent, rownames = FALSE, options = list(autoWidth=TRUE, scrollX=TRUE, scrollY=400, searching=TRUE, pageLength=10)))
          tissueoutputindependent<<-reactive({PlotTissueExpression(data.object(), all = FALSE, uniquePathways=uniquepathwaysperform(), clusterIndependent=TRUE)})
          output$tissuePlotIndependent<-renderPlot(tissueoutputindependent())
        }else{
          hideTab(inputId = "tabs", target = "Tissue_S")
        }

        #download data and plots
        output$saveObject <- renderUI({
          downloadButton("save.object","Results")
        })
        #download data and plots
        output$downloadData <- renderUI({
          downloadButton("download.data","Data")
        })
        #download plots using different format
        output$downloadDataFormats <- renderUI({
          dropdownButton(inputId="download.plot.format",label="Images", circle=FALSE, radioButtons("format", label="", choices = c("jpg","png","pdf"), selected = NULL), downloadButton("download.plots","Download plots"))
        })

        nclusters<-length(clusters)
        updateNumericInput(session, inputId="cluster", value=nclusters)

        if(checkGO(data.object())){
          wordcloudperform(TRUE)
          wordcloudinfo<<-reactive({PerformORAGOperCluster(data.object(), uniquePathways=uniquepathwaysperform(), clusterIndependent=FALSE)})
          output$wordcloudCheckbox <- renderUI({checkboxInput("wordcloud", "Enable Wordcloud", value=FALSE)})
        }else{
          wordcloudperform(FALSE)
          output$wordcloudCheckbox <- renderUI({tagList()})
          #create plots
          keywords_ora<-calculateKeywordsORA(ora.objectALL())
          heatmapoutput<<-reactive({PlotGeneSets(Object = data.object(), doORA = T, wordcloud = wordcloudperform(), uniquePathways=uniquepathwaysperform(), keywords_ora=keywords_ora)})
          output$heatmap<-renderPlot(heatmapoutput())
        }
        update_modal_progress(value=0.7)

        output$barplot<-renderPlot(createBarPlot(data.object(), clusterchoices[1], uniquePathways=uniquepathwaysperform()))
        update_modal_progress(value=0.9)

        updateSelectInput(session, "breakup.cluster", choices = clusterchoices)
        updateSelectInput(session, "summary.cluster", choices = clusterchoices)

        ##check if tissue is perfomed for dependent analysis
        if ((!is.null(data.object()@dfTissue))&&(length(data.object()@dfTissue)>0)){
          temp.object<-data.object()
          if(input$uploadfrom=="package"){
            if(uniquepathwaysperform()){
              temp.object@dfTissue<- temp.object@dfTissue$Pathway
            }else{
              temp.object@dfTissue<- temp.object@dfTissue$Geneset
            }
          }
          temp.object@dfTissue<-round(temp.object@dfTissue,3)
          temp.object@dfTissue<-na_replace(temp.object@dfTissue,0)
          data.object<<-reactive({temp.object})
          tissueData<-data.frame(rownames(data.object()@dfTissue), data.object()@dfTissue)
          colnames(tissueData)[1]<-"Tissue"
          output$tissueOutput<-renderDataTable(datatable(tissueData, rownames = FALSE, options = list(searching=TRUE, pageLength=10)))
          #shinyjs::hide("calculateTissueEnrichment")

          output$tissueIntro<-renderUI({tags$p("Tissue enrichment results per cluster (Cluster_number..number of terms):")})
          showTab(inputId = "tabs", target = "Tissue")

          #plot tissue plots
          tissueoutput<<-reactive({PlotTissueExpression(data.object(), all = FALSE, uniquePathways=uniquepathwaysperform(), clusterIndependent=FALSE)})
          output$tissuePlot<-renderPlot(tissueoutput())
        }else{
          #hide tissue results
          hideTab(inputId = "tabs", target = "Tissue")
        }

        showTab(inputId = "tabs3", target="tabpanelresult", select=TRUE)
        #updateTabsetPanel(session, "tabs3", selected = "tabpanelresult", args = list(tabpanelresult = "Results"))
        runjs('
              // Change the tab label when the button is clicked
              $("#tabs3 a[data-value=\'tabpanelresult\']").text("Results");
              ')

        if(unique(data.object()@metadata$organism)=="org.Mm.eg.db"){
          #for mouse data hide tissue enrichment analysis
          hideTab(inputId = "tabs2", target = "Tissue enrichment")
        }
        shinyjs::toggle("main1")
        shinyjs::toggle("main2")
        remove_modal_progress()

        #reset upload inputs
        updateRadioButtons(session,"uploadfrom", selected = character(0))
        updateRadioButtons(session,"uploadfilters", selected = character(0))
        updateTextInput(session,"uploadnameobject", value = "")
        reset('fileRdata')
      }
    }
  })

  ###########################
  ###########################
  #DOWNLOAD
  output$save.object <- downloadHandler(
    filename = function() {
      paste("GSC_", format(Sys.Date(), "%Y%m%d"), ".RData", sep = "")
    },
    content = function(file) {
      # Access the value of the reactive expression within a reactive context
      data_to_save <- data.object()
      annotations_to_save <- annotations()
      ora_to_save <- ora.objectALL()
      unique_pathways <- uniquepathwaysperform()
      independent_to_save<-independent.info()
      independent_ora_to_save<-independent.infoORA()
      # Check the structure of data_to_save
      #str(data_to_save)
      # Save the data to the file
      save(data_to_save, annotations_to_save, ora_to_save, unique_pathways, independent_to_save, independent_ora_to_save, file = file)
    }
  )

  output$downloadTemplate <- downloadHandler(
    filename = "template.xls",
    content = function(file) {
      file.copy("source/template.xls",file)
    }
  )

  output$downloadUserGuide <- downloadHandler(
    filename = "UserGuide.pdf",
    content = function(file) {
      file.copy("source/UserGuide.pdf",file)
    }
  )

  output$download.plots<-downloadHandler(
    file = function() {
      paste('plots-', Sys.Date(), '.zip', sep = '')
    },
    content = function(file) {

      nofile<-c()
      filenames<-c(paste("heatmap_plot.",input$format, sep=""),paste("tissue_plot.",input$format, sep=""), paste("heatmapS_plot.",input$format, sep=""), paste("tissueS_plot.",input$format, sep=""))
      # Save the heatmap plot
      #ggsave(plot = as.ggplot(heatmapoutput()), filename = filenames[1], width = 30, unit = "cm", device = input$format)
      ggsave(plot = as.ggplot(heatmapoutput()@ht_list[[1]]), filename = filenames[1], width = 30, height= 20, unit = "cm", device = input$format)
      if(!is.null(tissueoutput())){
        # Save the tissue plot
        ggsave(plot = tissueoutput(), filename = filenames[2], width = 30, height= 20, unit = "cm", device = input$format)
      }else{
        nofile<-c(nofile,2)
      }
      if(!is.null(independentoutput())){
        # Save the independent
        ggsave(plot = as.ggplot(independentoutput()@ht_list[[1]]), filename = filenames[3], width = 30, height= 20, unit = "cm", device = input$format)
      }else{

        nofile<-c(nofile,3)
      }
      if(!is.null(tissueoutputindependent())){
        # Save the tissue plot
        ggsave(plot = tissueoutputindependent(), filename = filenames[4], width = 30, height= 20, unit = "cm", device = input$format)
      }else{
        nofile<-c(nofile,4)
      }
      if(!is.null(nofile)){
        filenames<-filenames[-nofile]
      }
      # Create the zip file
      zip::zipr(zipfile = file, files = filenames)
    }
  )

  output$download.data<-downloadHandler(
    filename = function() {
      paste0("data_", Sys.Date(), ".zip")
    },
    content = function(file) {
      # Create a temporary directory
      temp_dir <- tempdir()

      if(!is.null(independent.info())){
        filenames <- c("dataClassic.csv", "dataSeriation.csv")
        #paths <- file.path(temp_dir, filenames)
        write.csv(independent.info(), filenames[2], row.names=FALSE)
      }else{
        # Save each data frame as a CSV file in the temporary directory
        filenames <- c("dataClassic.csv")
      }
      if (unique(data.object()@metadata$source)=="IPA"){
        write.csv(data.object()@Data[[1]][,c(1:9)], filenames[1], row.names=FALSE) # add parentheses to data arg if reactive
      }else{
        write.csv(data.object()@Data[[1]][,c(1:7,10)], filenames[1], row.names=FALSE) # add parentheses to data arg if reactive
      }

      # Zip the files
      zip(file, filenames)

    }
  )

  output$downloadGSE111385 <- downloadHandler(
    filename = "GSE111385_files.zip",
    content = function(file) {
      file.copy("example/GSE111385_files.zip",file)
    }
  )

  output$downloadGSE198256 <- downloadHandler(
    filename = "GSE198256_files.zip",
    content = function(file) {
      file.copy("example/GSE198256_files.zip",file)
    }
  )

} # server



