#' WriteGeneSets
#'
#' Exports the pathways into a csv file
#'
#' @param Object is a PathwayObject.
#' @param file_location where to write the files to
#' @param name name to be added to the files, what experiments are these
#' @param write what to write, either "Data", "RR" or "Both"
#'
#' @return written tables in the folder designated.
#'
#'
setGeneric(name="WriteGeneSets",
           def=function(Object, file_location = "~/Project9/test/", name = "IPA_20181123", write = "Both")
           {
             standardGeneric("WriteGeneSets")
           }
)
#' WriteGeneSets
#'
#' Exports the pathways into a csv file
#'
#' @import utils
#' @param Object is a PathwayObject.
#' @param PathwayObject  a PathwayObject
#' @param file_location where to write the files to
#' @param name name to be added to the files, what experiments are these
#' @param write what to write, either "Data", "RR" or "Both"
#'
#' @return written tables in the folder designated.
#'
setMethod(f="WriteGeneSets",
          signature="PathwayObject",
          definition=function(Object, file_location = "~/Project9/test/", name = "IPA_20181123", write = "Both")
          {
            if(write=="Data")
            {
              write.table(Object@Data[[1]],
                          file= paste(file_location,name,"_Pathway",".csv",sep=""),
                          sep=";",
                          row.names = F)
            }
            if(write=="RR")
            {
              write.table(Object@Data.RR,
                          file= paste(file_location,name,"_RR",".csv",sep=""),
                          sep=";",
                          row.names = F)
            }
            if(write=="Both")
            {
              write.table(Object@Data[[1]],
                          file= paste(file_location,name,"_Pathway",".csv",sep=""),
                          sep=";",
                          row.names = F)
              write.table(Object@Data.RR,
                          file= paste(file_location,name,"_RR",".csv",sep=""),
                          sep=";",
                          row.names = F)
            }
          }
)
