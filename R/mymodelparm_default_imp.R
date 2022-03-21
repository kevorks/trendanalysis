mymodelparm_default_imp <- function (model, coef. = coef, vcov. = vcov, df = NULL, ...)
{
  beta <- try(coef.(model))
  if (inherits(beta, "try-error"))
    stop("no ", sQuote("coef"), " method for ", sQuote("model"),
         " found!")
  sigma <- try(vcov.(model))
  if (inherits(sigma, "try-error"))
    stop("no ", sQuote("vcov"), " method for ", sQuote("model"),
         " found!")
  sigma <- as.matrix(sigma)
  if (any(length(beta) != dim(sigma)))
    beta = na.omit(beta)
  if (is.null(df)) {
    df <- 0
    if (inherits(model, "aov") || inherits(model, "lm") ||
        inherits(model, "glm")) {
      class(model) <- "lm"
      df <- summary(model)$df[2]
    }
    if (inherits(model, "gls")) {
      dd <- model$dims
      df <- dd[["N"]] - dd[["p"]]
    }
    if (inherits(model, "parm"))
      df <- model$df
  }
  else {
    if (df < 0)
      stop(sQuote("df"), " is not positive")
  }
  ocoef <- coef.(model)
  if (inherits(model, "aov"))
    ocoef <- model$coefficients
  estimable <- rep(TRUE, length(ocoef))
  if (any(is.na(ocoef))) {
    estimable[is.na(ocoef)] <- FALSE
    beta <- ocoef[estimable]
    if (dim(sigma)[1] == length(estimable))
      sigma <- sigma[estimable, estimable]
  }
  if (length(beta) != ncol(sigma) || nrow(sigma) != sum(estimable))
    stop("could not extract coefficients and covariance matrix from ",
         sQuote("model"))
  RET <- list(coef = beta, vcov = sigma, df = df, estimable = estimable)
  class(RET) <- "mymodelparm"
  RET
}
