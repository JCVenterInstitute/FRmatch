
#' Shiny App for FRmatch Demo
#'
#' This function lauches the Shiny App for FRmatch Demo.
#'
#' @export


runShiny <- function() {
  appDir <- system.file("shinyapp", package = "FRmatch")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `FRmatch`.", call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal")
}
