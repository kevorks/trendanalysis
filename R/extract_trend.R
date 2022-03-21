##' Method for extract trends from model
##' @param model gls or lm model
##' @param alpha Level of significance
##' @param ... zur Zeit nicht benutzt
##' @return String with trendresults
extract_trend <- function(model, alpha, ...) {
  UseMethod("extract_trend")
}

##' Extract trend from lm-objects
##' @param model lm model
##' @param alpha Level of significance
##' @return String with trendresults from lm
extract_trend.lm <- function(model, alpha = 0.05) {
  ## extract p-wert for linear term
  smod <- summary(model)
  beta0 <- smod$coefficients["(Intercept)", "Estimate"]
  beta1 <- smod$coefficients["Zeit", "Estimate"]
  beta2 <- as.character(0)
  p_beta1 <- smod$coefficients["Zeit", "Pr(>|t|)"]

  ## if the model includes a quadratic term a quadratic trend is present
  quadratic <- "I(Zeit^2)" %in% rownames(smod$coefficients)

  if (quadratic){
    # since the function fit_trend was already checked for significance we can assume
    # a quadratic trend
    beta2 <- smod$coefficients["I(Zeit^2)", "Estimate"]
    extrem <- (-beta1 / 2 / beta2)
    maxi <- max(as.numeric(model$model$Zeit))
    mini <- min(as.numeric(model$model$Zeit))
    if (extrem >= mini && extrem <= maxi) {
      trend <- ifelse(beta2 > 0, "Quadratischer Trend(u)", "Quadratischer Trend(n)")
    } else {
      trend <- ifelse(model$fitted.values[1] < model$fitted.values[length(model$fitted.values)],
                      "Quadratischer Trend(+)", "Quadratischer Trend(-)")
    }
  } else {
    # Check if p-values could be calculated
    if (!is.na(p_beta1)) {
      if (p_beta1 < alpha) {
        trend <- ifelse(beta1 < 0, "Fallender Trend", "Steigender Trend")
      } else {
        trend <- "Kein signifikanter Trend"
      }
    } else {
      # if a p-value cant be cant be calculated
      trend <- "es konnte kein p-Wert berechnet werden"
    }
  }
  return(list(beta0, beta1, beta2, trend))
}

##' Extract trend from gls-objects
##' @param model gls model
##' @param alpha level of significance
##' @return String with trendresults from gls
extract_trend.gls <- function(model, alpha = 0.05){

  ## extract p-value for linear term
  smod <- summary(model)
  beta0 <- smod$tTable["(Intercept)", "Value"]
  beta1 <- smod$tTable["Zeit", "Value"]
  beta2 <- as.character(0)
  p_beta1 <- smod$tTable["Zeit", "p-value"]

  ## Quadratic trend if a quadratic term is present
  quadratic <- "I(Zeit^2)" %in% rownames(smod$tTable)

  if (quadratic) {
    # since the function fit_trend was already checked for significance we can assume
    # a quadratic trend
    beta2 <- smod$tTable["I(Zeit^2)", "Value"]
    extrem <- (-beta1 / 2 / beta2)
    maxi <- max(as.numeric(names(model$fitted)))
    mini <- min(as.numeric(names(model$fitted)))
    if (extrem >= mini && extrem <= maxi) {
      trend <- ifelse(beta2 > 0, "Quadratischer Trend(u)", "Quadratischer Trend(n)")
    } else {
      trend <- ifelse(model$fitted[1] < model$fitted[length(model$fitted)],
                      "Quadratischer Trend(+)", "Quadratischer Trend(-)")
    }
  } else {
    if (p_beta1 < alpha) {
      trend <- ifelse(beta1 < 0, "Fallender Trend", "Steigender Trend")
    } else {
      trend <- "Kein signifikanter Trend"
    }
  }
  return(list(beta0, beta1, beta2, trend))
}
