library(shiny)
library(shinythemes)
library(DT)
library(shinydashboard)
library(shinyWidgets)
library(shinyBS)
library(shinyjs)
library(ggplot2)
library(gridExtra)
library(htmltools)
library(shinyFiles)

#data
options(shiny.maxRequestSize=50*1024^3)

js <- "$(document).ready(function(){
    var $treeplot = $('#hello li > a[data-value=treeplot_val]').parent();
    $treeplot.removeClass('active').addClass('hide');
    $('#printPlot').on('click', function(){
      $treeplot.removeClass('hide');
    });
  });
"
js_code <- "
shinyjs.browseURL = function(url) {
  window.open(url,'_blank');
}
"
tissues<-c("Adipose Tissue (Subcutaneous)","Artery (Aorta)","Brain (Cortex)","Colon (Sigmoid)","Kidney (Cortex)","Liver","Lung","Pancreas","Skeletal Muscle","Skin (Lower Leg)","Small Intestine (Terminal Ileum)","Stomach","Thyroid","Vagina","Whole Blood")

# Define UI
ui <- fluidPage(theme = shinytheme("flatly"),
                tags$head(
                  tags$script(HTML(js)),
                  tags$link(rel = "stylesheet", href = "https://fonts.googleapis.com/css2?family=Pacifico&display=swap"),
                  tags$style(
                    HTML("
                      .covered-element {
                        margin-top: 70px; /* Adjust the value as needed */
                      }
                      .covered-element2 {
                        margin-top: 0px; /* Adjust the value as needed */
                      }
                    ")
                  ),
                ),
                useShinyjs(),
                useShinydashboard(),
                extendShinyjs(text = js_code, functions = 'browseURL'),
                navbarPage(
                  #"GeneSetCluster",
                  title = div(HTML("GeneSet<br>Cluster"), style = "font-weight: bold; font-family: 'Consolas', monospace; font-size: 17px; color: white;"),
                  windowTitle = "GeneSetCluster",
                  id="navbar",
                  position = c("fixed-top"),
                  header = tagList(useShinydashboard()),
                  tabPanel("Main",
                           #fluidRow(
                             div(class="covered-element",
                             tabsetPanel(id="tabs3", type="pills",
                                tabPanel("New",
                                   br(),
                                   column(width=3,
                                         div(style="display:inline-block;",downloadButton("downloadTemplate", "Download template", style = "font-size: 12px; padding: 6px; background-color: #7EBC72; border-color: #7EBC72; color: white;")),
                                         div(style="display:inline-block;", dropdownButton(status="runExamples",label="Examples", circle=FALSE, radioButtons("exampledataset", label="", choices = c("GSE111385 (GREAT)"="mousedataset","GSE198256 (GSEA)"="coviddataset"), selected = NULL), actionButton(inputId="setExample", label="Run example", icon=icon("play", lib = "glyphicon"), style = "font-size: 12px; padding: 6px; background-color: #70B7D7; border-color: #70B7D7; color: white;"), width="200px")),
                                         tags$head(tags$style(
                                           ".btn-runExamples{
                                              background-color: #70B7D7 !important;
                                              border: solid; border-radius: 4px; border-width: 1px; border-color: #70B7D7;
                                              color: white; text-align: center; font-weight:bold; font-size: 12px;
                                            }
                                           .btn-runExamples .shiny-input-container{
                                              max-width: 13px;
                                           }
                                           "
                                         )),
                                         br(),
                                         radioButtons("source", label = h4("Source", bsButton("infosource", label="", icon=icon("info"), style="info", size="extra-small")), choices=c("GREAT"="Great","IPA","GSEA"="GSEA", "Template"="GSEA2"), selected = character(0), inline=TRUE),
                                         bsPopover(id="infosource", title="Source", content="GREAT, IPA, GSEA or template inputs.", placement = "right", trigger="hover"),
                                         radioButtons("gene.id", label = h4("Gene ID", bsButton("infogeneid", label="", icon=icon("info"), style="info", size="extra-small")), choices=c("Ensembl ID"="ENSEMBL", "Symbol"="SYMBOL", "Entrez ID"="ENTREZID"), selected = character(0), inline=TRUE),
                                         bsPopover(id="infogeneid", title="Gene ID", content="Gene identification used.", placement = "right", trigger="hover"),
                                         radioButtons("organism", label = h4("Organism", bsButton("infoorganism", label="", icon=icon("info"), style="info", size="extra-small")), choices=c("Homo sapiens"="org.Hs.eg.db", "Mus musculus"="org.Mm.eg.db"), selected = character(0), inline=TRUE),
                                         bsPopover(id="infoorganism", title="Organism", content="Organism of your data.", placement = "right", trigger="hover"),
                                         radioButtons("filters", label = h4("Approaches", bsButton("infofilters", label="", icon=icon("info"), style="info", size="extra-small")), choices=c("Unique GeneSet"="none","Unique Pathways"="uniquepathways"), selected = character(0), inline=TRUE),
                                         bsPopover(id="infofilters", title="Approaches", content="Unique_GeneSet:include all genesets independently. Unique_Pathways: combine genesets with same label.", placement = "right", trigger="hover"),
                                         fileInput("files", h4("Upload files",bsButton("infofiles", label="", icon=icon("info"), style="info", size="extra-small")), multiple=TRUE, accept = c("text/csv","text/comma-separated-values,text/plain",".csv", ".txt", ".xls", ".xlsx")),
                                         bsPopover(id="infofiles", title="Files input", content=".txt, .csv or .xls files.", placement = "right", trigger="hover"),
                                         DTOutput("fileTable"),
                                         br(),
                                         actionButton(inputId="printResults", label="Run analysis", icon=icon("play", lib = "glyphicon")),
                                   ),
                                ),
                                tabPanel("Upload",
                                    br(),
                                    column(width=3,
                                        radioButtons("uploadfrom", label = h4("From:", bsButton("infouploadfrom", label="", icon=icon("info"), style="info", size="extra-small")), choices=c("Shiny"="shiny","R package"="package"), selected = character(0), inline=TRUE),
                                        bsPopover(id="infouploadfrom", title="From", content="Extracted from R package or from the shiny application.", placement = "right", trigger="hover"),
                                        radioButtons("uploadfilters", label = h4("Approaches", bsButton("infouploadfilters", label="", icon=icon("info"), style="info", size="extra-small")), choices=c("Unique GeneSet"="none","Unique Pathways"="uniquepathways"), selected = character(0), inline=TRUE),
                                        bsPopover(id="infouploadfilters", title="Approaches", content="Unique_GeneSet:include all geneset independently. Unique_Pathways: combine geneset with same label.", placement = "right", trigger="hover"),
                                        textInput("uploadnameobject", label = h4("Object name", bsButton("infonameobject", label="", icon=icon("info"), style="info", size="extra-small"))),
                                        bsPopover(id="infonameobject", title="Name", content="Name of the GeneSetCluster object.", placement = "right", trigger="hover"),
                                        fileInput("fileRdata", "Choose a file", accept = ".RData"),
                                        actionButton(inputId="uploadResult", label="Load", icon=icon("play", lib = "glyphicon"))
                                    ),
                                ),
                                tabPanel("", value="tabpanelresult",
                                     br(),
                                     conditionalPanel(
                                       condition = "input.tabs != 'Heatmap_S' && input.tabs != 'Tissue_S'",
                                       column(width=3,
                                          selectInput("summary.cluster","Select cluster (classic):", choices = NULL),
                                          plotOutput("barplot", width="80%"),
                                       ),
                                     ),
                                     conditionalPanel(
                                       condition = "input.tabs == 'Heatmap_S' || input.tabs == 'Tissue_S'",
                                       column(width=3,
                                              selectInput("summary.independentgroup","Select cluster (seriation):", choices = NULL),
                                              plotOutput("barplotIndependent", width="80%"),
                                       )
                                     ),
                                ),
                             ),
                             ),
                           #),
                           conditionalPanel(
                             condition = "(input.tabs3 == 'New' || input.tabs3 == 'Upload') && !(input.tabs == 'Heatmap_S' || input.tabs == 'Tissue_S' || input.tabs == 'Tissue') && !(output.heatmap)",
                             div(style = "text-align: left; margin-left: 40px;",
                                 img(src = "GeneSetClusterLogo_small.png", height = "150px", style = "float: left; margin-right: 20px;"),
                                 h1(style = "font-family: 'Consolas', monospace; margin-top: 20px; margin-bottom: 10px; font-size: 70px", "Gene Set Cluster"),
                                 h4(style = "font-family: 'Consolas', monospace; padding-left: 20px; margin-bottom: 40px; font-size: 25px", "Summarizing and integrating genesets results")
                             ),
                           ),
                             div(id="main1", class="covered-element2",
                               column(width=9,
                                     fluidRow(
                                       column(10,
                                         actionButton(inputId="reset", label="Reset", icon=icon("refresh", lib = "glyphicon")),
                                         p("plots", style="color:white"),
                                         tabsetPanel(id="tabs", type="tabs",
                                                     tabPanel("Heatmap",
                                                          uiOutput("wordcloudCheckbox"),
                                                          plotOutput("heatmap", height = "500px"),
                                                     ),
                                                     tabPanel("Tissue",
                                                          plotOutput("tissuePlot", height = "500px"),
                                                     ),
                                                     tabPanel("Heatmap_S",
                                                          uiOutput("wordcloudCheckboxIndependent"),
                                                          plotOutput("independentPlot", height = "500px")
                                                     ),
                                                     tabPanel("Tissue_S",
                                                          plotOutput("tissuePlotIndependent", height = "500px")
                                                     ),
                                         ),
                                       ),
                                       column(2,
                                         fluidRow(div(
                                           div(style="display:inline-block;", p("ai", style="color:white"), width=6),
                                           div(style="display:inline-block;", uiOutput("saveObject"), width=6),
                                           div(style="display:inline-block;", uiOutput("downloadData"), width=6),
                                         )),
                                         br(),
                                         div(style="display:inline-block",uiOutput("downloadDataFormats"), br(), width=6),
                                         conditionalPanel(
                                           condition = "input.tabs != 'Heatmap_S' && input.tabs != 'Tissue_S'",
                                           box(title=h4("ReClustering", style="margin-top:-10px"), style = "font-size:13px", width=15, status="primary", collapsible = TRUE, collapsed=TRUE,
                                            numericInput("cluster","Set number of clusters:",value=NULL,min=2, width="200px"),
                                            actionButton(inputId="recalculateClustering", label="Recalculate", icon=icon("play", lib = "glyphicon")),
                                           ),
                                           box(title=h4("SubClustering", style="margin-top:-10px"), style = "font-size:13px", width=15, status="primary", collapsible = TRUE, collapsed=TRUE,
                                             selectInput("breakup.cluster","Select cluster:", choices = NULL),
                                             selectInput("nbreakup.cluster","Select number cluster:", choices = c("Automatic", 2:20), selected = c("Automatic")),
                                             actionButton(inputId="breakup", label="Split", icon=icon("play", lib = "glyphicon")),
                                           )
                                         ),
                                       ),
                                     ),
                                  ),#column main panel
                              ) %>% shinyjs::hidden(),
                            #),
                            div(id="main2",
                              fluidRow(
                                 column(width=12,
                                     hr(style = "border-top: 1px solid #000000;"),
                                     conditionalPanel(
                                       condition = "input.tabs2 != 'Seriation'",
                                        uiOutput("checkboxCluster"),
                                     ),
                                     tabsetPanel(id="tabs2", type="pills",
                                                 tabPanel("Data",
                                                      br(),
                                                      column(9,
                                                        DTOutput("dataInfo"),
                                                      ),
                                                      column(3,
                                                        uiOutput("databases"),
                                                        uiOutput("databaseOptionInformation"),
                                                        uiOutput("databaseOptions"),
                                                        textOutput("databaseOptionTitle"),
                                                        verbatimTextOutput("databaseGenes"),
                                                        fluidRow(div(
                                                          div(style="display:inline-block;", p("hola", style="color:white"), width=6),
                                                          div(style="display:inline-block;", actionButton(inputId="calculateHighlight", label="Calculate", icon=icon("play", lib = "glyphicon")), width=6),
                                                          div(style="display:inline-block;", actionButton(inputId="clearHighlight", label="Clear", icon=icon("refresh", lib = "glyphicon")), width=6),
                                                        )),
                                                        br(),
                                                        DTOutput("highlightInfo"),
                                                      )
                                                 ),
                                                 tabPanel("ORA",
                                                      div(style="display:inline-block;vertical-align:center;horizontal-align:center;", class="row-fluid", strong("Top"),),
                                                      div(style="display:inline-block; padding-left:10px; width:auto", class="row-fluid", selectInput("top", "", choices=c(5,10,15,20), selected = c(5), width="60px")),
                                                      div(style="display:inline-block;vertical-align:center;horizontal-align:center;padding-left:10px", class="row-fluid", strong("per cluster:")),
                                                      br(),
                                                      br(),
                                                      DTOutput("oraInfo"),
                                                 ),
                                                 tabPanel("Genes",
                                                          fluidRow(column(12,
                                                            div(style="display:inline-block; float:right",actionButton(inputId="performORAgenes", label="ORA", icon=icon("play", lib = "glyphicon"))),
                                                          )),
                                                          br(),
                                                          p("Table of gene frequency distribution per cluster:"),
                                                          DTOutput("genesInfo"),
                                                          downloadButton("downloadORAgenes", "",style = "font-size: 2px; background-color: white; border-color: white; color: white;"),
                                                 ),
                                                 tabPanel("Tissue enrichment",
                                                      br(),
                                                      uiOutput("tissueIntro"),
                                                      fluidRow(div(
                                                        div(style="display:inline-block;", pickerInput(inputId = "tissuesselected", label="Tissues:", choices=tissues, selected=NULL, multiple=TRUE)),
                                                        div(style="display:inline-block;", actionButton(inputId="calculateTissueEnrichment", label="Run", icon=icon("play", lib = "glyphicon"))),
                                                      )),
                                                      DTOutput("tissueOutput"),
                                                 ),
                                                 tabPanel("Seriation",
                                                          br(),
                                                          uiOutput("independentIntro"),
                                                          actionButton(inputId="runIndependent", label="Run", icon=icon("play", lib = "glyphicon")),
                                                          fluidRow(
                                                            column(5,
                                                                   DTOutput("independentOutput"),
                                                            ),
                                                            column(7,
                                                              tabsetPanel(id="tabsindependent", type="pills",
                                                                tabPanel("ORA",
                                                                    DTOutput("independentORAOutput"),
                                                                  ),
                                                                tabPanel("Genes",
                                                                    DTOutput("genesInfoIndependent"),
                                                                ),
                                                                tabPanel("Tissue",
                                                                     fluidRow(div(
                                                                       div(style="display:inline-block;", pickerInput(inputId = "tissuesselectedindependent", label="Tissues:", choices=tissues, selected=NULL, multiple=TRUE)),
                                                                       div(style="display:inline-block;", actionButton(inputId="calculateTissueEnrichmentIndependent", label="Run", icon=icon("play", lib = "glyphicon"))),
                                                                     )),
                                                                     DTOutput("tissueOutputIndependent"),
                                                                  ),
                                                                ),
                                                            ),
                                                          )
                                                 ),
                                     ),
                                ),
                              ),
                            ) %>% shinyjs::hidden(),
                  ), # Navbar 1, tabPanel
                  tabPanel("About",  br(), br(), br(),
                             h2("Documentation"),
                               fluidRow(
                                 div(style="display:inline-block;", p("gene", style="color:white")),
                                 div(style="display:inline-block;",h4("Shiny user guide:")),
                                 div(style="display:inline-block;",downloadButton("downloadUserGuide", "User guide", style = "font-size: 12px; padding: 6px; background-color: #E4B5BB; border-color: #E4B5BB; color: black;")),
                               ),
                               fluidRow(
                                 div(style="display:inline-block;", p("gene", style="color:white")),
                                 div(style="display:inline-block;",h4("Examples:")),
                                 div(style="display:inline-block;",downloadButton("downloadGSE111385", "GSE111385", style = "font-size: 12px; padding: 6px; background-color: #E4D2B5; border-color: #E4D2B5; color: black;")),
                                 div(style="display:inline-block;",downloadButton("downloadGSE198256", "GSE198256", style = "font-size: 12px; padding: 6px; background-color: #E4D2B5; border-color: #E4D2B5; color: black;")),
                               ),
                               fluidRow(
                                 div(style="display:inline-block;", p("gene", style="color:white")),
                                 div(style="display:inline-block;",h4("R package vignette:")),
                                 div(style="display:inline-block;",tags$a(href = "https://htmlpreview.github.io/?https://github.com/TranslationalBioinformaticsUnit/GeneSetCluster2.0/blob/main/inst/extdata/GeneSetCluster_vignette.html", "link", target = "_blank"))
                               ),
                               fluidRow(
                                 div(style="display:inline-block;", p("gene", style="color:white")),
                                 div(style="display:inline-block;",h4("GitHub:")),
                                 div(style="display:inline-block;",tags$a(href = "https://github.com/TranslationalBioinformaticsUnit/GeneSetCluster2.0", "link", target = "_blank"))
                               ),
                           br(),
                           h2("References"),
                             p("If you use GeneSetCluster, please cite the following paper"),
                             p("[1] Ewing, E., Planell-Picola, N., Jagodic, M. et al. GeneSetCluster: a tool for summarizing and integrating gene-set analysis results. BMC Bioinformatics 21, 443 (2020).", tags$a(href = "https://doi.org/10.1186/s12859-020-03784-z", "https://doi.org/10.1186/s12859-020-03784-z",target = "_blank")),
                           br(),
                           h2("Contact"),
                            p(tags$b("Ewoud Ewing: "), tags$a(href = "mailto:ewoud.ewing@ki.se", "ewoud.ewing@ki.se")),
                            p("Department of Clinical Neuroscience"),
                            p("Center for Molecular Medicine"),
                            p("Karolinska Institutet"),
                           br(),
                           p(tags$b("Asier Ortega-Legarreta: "), tags$a(href = "mailto:aortegal@navarra.es", "aortegal@navarra.es")),
                           p("Translational Bioinformatic Unit"),
                           p("Navarrabiomed"),
                          ),
                ) # navbarPage
) # fluidPage

