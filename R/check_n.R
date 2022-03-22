##' Checks if enough data is available
##' @description Returns FALSE for less than 7 observations
##' @param trend_variable Vector with values
##' @param min_n Minimum amount of observations for trenda
##' @importFrom stats na.omit
##' @return logical
check_n <- function(trend_variable, min_n = 7) {
  length(na.omit(trend_variable)) >= min_n
}
