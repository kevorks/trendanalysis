##' Extract phi from the gls model
##' @param mod gls model
##' @importFrom stats coef
##' @return phi parameter of correlation
extract_phi <- function(mod) {
  # the function is specified for gls objects
  if (inherits(mod, "gls")) {
    x <- summary(mod$modelStruct)$corStruct
    class(x) <- attr(x, "oClass")
    phi <- round(coef(x, unconstrained = FALSE), 2)
    # phi is 0 for linear model
  } else {
    phi <- 0
  }
  if (length(phi) == 1) {phi[2] <- 0}
  return(as.character(phi))
}
