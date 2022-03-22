#' Copy of predictmeans::CookD, but with ability to include the used formula
#' and correlation structure so that the update function works. The model needs
#' to be fitted with ML.
#'
#' @inheritParams predictmeans::CookD
#' @param used_formula formula used when fitting the model
#' @param cor_structure correlation structure used when fitting the
#' @param data data used for the fit
#'
#' @details Rather quick and dirty method; not exported function from package
#' predictmeans is used. Maybe better implementation if scoping is adapted.
#' @importFrom grDevices dev.new
#' @importFrom graphics points text
#' @importFrom stats update vcov
#' @return vector with Cook's distance
CookD_gls <- function (model, group = NULL, plot = TRUE, idn = 3, newwd = TRUE,
                       used_formula, cor_structure, data)
{
  if (class(model)[1] != "gls") stop("The model has to be of class 'gls'")
  if (model$method != "ML") stop("The model has to be fitted with method 'ML'")

  mdf <- data

  mp <- mymodelparm_default_imp(model)
  beta0 <- mp$coef
  vcovb <- mp$vcov
  vb.inv <- solve(vcovb)
  if (is.null(group) || group %in% c("NULL", "")) {
    rn <- rownames(mdf)
    LOOmp <- lapply(rn, function(x) mymodelparm_default_imp(gls(model = used_formula,
                                                                data = mdf[rn != x, ],
                                                                correlation = cor_structure,
                                                                method = "ML")))
  }
  else {
    rn <- unique(mdf[, group])
    LOOmp <- lapply(rn, function(x) {
      rind <- mdf[, group] != x
      mymodelparm_default_imp(update(model, data = mdf[rind, ]))
    })
  }
  LOObeta <- sapply(LOOmp, function(x) x$coef)
  rK <- t(LOObeta - beta0)
  CookD <- diag(rK %*% tcrossprod(vb.inv, rK)/length(beta0))
  names(CookD) <- rn
  if (plot) {
    if (newwd)
      dev.new()
    outD <- CookD >= sort(CookD, decreasing = TRUE)[idn]
    labid <- names(CookD)
    plot(CookD, xlab = "Obs. number", col = "blue", ylim = c(0,
                                                             max(CookD) + 0.005), main = "Cook's Distance", ylab = "Cook's distance",
         type = "h")
    text((1:length(CookD))[outD], CookD[outD], labid[outD],
         pos = 3)
    points((1:length(CookD))[outD], CookD[outD], pch = 16,
           col = "blue")
  }
  return(invisible(CookD))
}
