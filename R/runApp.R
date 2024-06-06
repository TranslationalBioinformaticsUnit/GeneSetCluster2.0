#' @title run_app
#'
#' @description
#' Runs the Shiny application
#'
#' @importFrom shiny runApp
#'
#'
#' @return NULL
#'
#' @export
#' @examples
#' \dontrun{
#' run_app()
#' }
#'

run_app <- function()
{
  path <- system.file("shiny", "app.R", package = "GeneSetCluster")
  shiny::runApp(path)
}
