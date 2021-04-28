
#' Shiny App for FRmatch Demo
#'
#' This function launches the Shiny App for FRmatch Demo.
#'
#' @import shiny
#' @export


runShiny <- function() {
  appDir <- system.file("shinyapp", package = "FRmatch")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `FRmatch`.", call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal")
}
