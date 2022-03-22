##' Extrahiere phi aus gls - Modell
##' @param mod gls-Modell
##' @importFrom stats coef
##' @return phi Parameter der Korrelation
extract_phi <- function(mod) {
  if (inherits(mod, "gls")) { # Die Funktion ist fuer gls-Objekt spezifiziert
    x <- summary(mod$modelStruct)$corStruct
    class(x) <- attr(x, "oClass")
    phi <- round(coef(x, unconstrained = FALSE), 2)
  } else { # fuer lineares Modell betrdgt phi eben 0
    phi <- 0
  }
  if (length(phi) == 1) {phi[2] <- 0}
  return(as.character(phi))
}
