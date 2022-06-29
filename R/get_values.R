
#' Utility function to get the plotted values
#'
#' Get the plotted values from an output matrix if \code{return.value = T}.
#'
#' @param output The output from \code{plot_FRmatch()} or \code{plot_FRmatch()} if \code{return.value=T}.
#'
#' @return A data frame of row names, column names, and values.
#'
#' @seealso \code{\link[FRmatch]{plot_FRmatch}, \link[FRmatch]{plot_bi_FRmatch}}.
#'
#' @export


get_values <- function(output){
  df <- output %>% as_tibble(rownames = NA) %>% rownames_to_column(var = "rowname") %>%
    pivot_longer(!rowname, names_to = "colname", values_to = "value") %>% filter(value>0)
  return(df)
}
