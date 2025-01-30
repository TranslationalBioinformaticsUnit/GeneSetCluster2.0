# Shiny applciation
# library(shiny)
# library(shinythemes)
# library(shinydashboard)
# library(shinyWidgets)
# library(shinyBS)
# library(shinyjs)
# library(shinyalert)
# library(shinyFiles)
# library(shinybusy)
# library(DT)
# library(ggplot2)
# library(gridExtra)
# library(htmltools)
# library(ggplotify)
# library(tools)
# library(powerjoin)
# library(imputeTS)
# library(htmltools)
# library(reshape2)
# library(clusterProfiler)
# library(factoextra)
# library(GGally)
# library(org.Hs.eg.db)
# library(org.Mm.eg.db)
# library(dplyr)
# library(stringr)
# library(jsonlite)
# library(doParallel)
# library(parallel)
# library(httr)
# library(utils)
# library(readxl)
# library(pbapply)
# library(RColorBrewer)
# library(patchwork)
# library(grid)
# library(simplifyEnrichment)
# library(ggwordcloud)
# library(ComplexHeatmap)
# library(colorRamp2)
# library(bigstatsr)
# library(seriation)
# library(GO.db)
# library(limma)

dependencies <- c(
  "shiny", "shinythemes", "shinydashboard", "shinyWidgets", "shinyBS", "shinyjs",
  "shinyalert", "shinyFiles", "shinybusy", "DT", "ggplot2", "gridExtra",
  "htmltools", "ggplotify", "tools", "powerjoin", "imputeTS", "reshape2",
  "clusterProfiler", "factoextra", "GGally", "org.Hs.eg.db", "org.Mm.eg.db",
  "dplyr", "stringr", "jsonlite", "doParallel", "parallel", "httr", "utils",
  "readxl", "pbapply", "RColorBrewer", "patchwork", "grid", "simplifyEnrichment",
  "ggwordcloud", "ComplexHeatmap", "colorRamp2", "bigstatsr", "seriation",
  "GO.db", "limma"
)

unavailable_packages <- suppressMessages(lapply(dependencies, require,
                                                character.only = TRUE))

if (!all(unlist(unavailable_packages)))
{
  not_installed <- dependencies[which(unlist(unavailable_packages)==FALSE)]
  message("The following package are not intalled and are necessary for the shiny application:")
  message(paste0(not_installed, sep="\n"))
  stop("Please install the requiered packages for the shiny application.")
}

suppressMessages(source("server.R"))
suppressMessages(source("ui.R"))

shinyApp(ui = ui, server = server)
