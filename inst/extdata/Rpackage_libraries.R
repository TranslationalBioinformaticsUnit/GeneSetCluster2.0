#############################################
####### R PACKAGE (MANDATORY) ###############
#############################################

packages <- c("limma", "stats", "methods", "RColorBrewer", "clusterProfiler",
               "readxl", "cluster", "factoextra", "stringr", "AnnotationDbi", "ComplexHeatmap",
              "GetoptLong", "bigstatsr", "colorRamp2", "doParallel", "foreach", "ggplot2",
              "grid", "parallel", "patchwork", "pbapply", "reshape2", "seriation", "graphics",
              "simplifyEnrichment", "slam", "utils", "grDevices", "rGREAT",
              "org.Hs.eg.db", "org.Mm.eg.db", "STRINGdb", "GGally", "network", "httr", "jsonlite",
              "dplyr", "sna", "shiny")
#check if packages are already installed
packagecheck <- match( packages, utils::installed.packages()[,1] )
packagestoinstall <- packages[ is.na( packagecheck ) ]

#install CRAN packages
if( length( packagestoinstall ) > 0L ) {
  utils::install.packages( packagestoinstall,
                           repos = "http://cran.us.r-project.org"
  )
} else {
  print( "All requested packages from CRAN already installed" )
}
packagecheck <- match( packages, utils::installed.packages()[,1] )
packagestoinstall <- packages[ is.na( packagecheck ) ]
#install BIOCONDUCTOR packages
if( length( packagestoinstall) > 0L ) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(packagestoinstall)	
} else {
  print( "All requested packages from BIOCONDUCTOR already installed" )
}
for( package in packages ) {
  suppressPackageStartupMessages(
    library( package, character.only = TRUE, quietly = TRUE )
  )
}


#############################################
######### R SHINY (OPTIONAL) ###############
#############################################
packages <- c("shinythemes", "shinydashboard", "shinyWidgets", "shinyBS", "shinyjs",
              "shinyalert", "shinyFiles", "shinybusy", "DT", "GO.db", "htmltools", "ggplotify",
              "gridExtra", "ggwordcloud", "tools", "powerjoin", "imputeTS")
#check if packages are already installed
packagecheck <- match( packages, utils::installed.packages()[,1] )
packagestoinstall <- packages[ is.na( packagecheck ) ]

#install CRAN packages
if( length( packagestoinstall ) > 0L ) {
  utils::install.packages( packagestoinstall,
                           repos = "http://cran.us.r-project.org"
  )
} else {
  print( "All requested packages from CRAN already installed" )
}
packagecheck <- match( packages, utils::installed.packages()[,1] )
packagestoinstall <- packages[ is.na( packagecheck ) ]
#install BIOCONDUCTOR packages
if( length( packagestoinstall) > 0L ) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(packagestoinstall)	
} else {
  print( "All requested packages from BIOCONDUCTOR already installed" )
}
for( package in packages ) {
  suppressPackageStartupMessages(
    library( package, character.only = TRUE, quietly = TRUE )
  )
}

