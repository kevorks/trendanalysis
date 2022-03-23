##' Method for extract trends from model
##' @param model gls or lm model
##' @param alpha Level of significance
##' @param ... placeholder
##' @return String with trend results
extract_trend <- function(model, alpha, ...) {
  UseMethod("extract_trend")
}

##' Extract trend from lm-objects
##' @param model lm model
##' @param alpha Level of significance
##' @return String with trend results from lm
extract_trend.lm <- function(model, alpha = 0.05) {
  ## extract p-value for linear term
  smod <- summary(model)
  beta0 <- smod$coefficients["(Intercept)", "Estimate"]
  beta1 <- smod$coefficients["Time", "Estimate"]
  beta2 <- as.character(0)
  p_beta1 <- smod$coefficients["Time", "Pr(>|t|)"]

  ## if the model includes a quadratic term a quadratic trend is present
  quadratic <- "I(Time^2)" %in% rownames(smod$coefficients)

  if (quadratic){
    # since the function fit_trend was already checked for significance we can assume
    # a quadratic trend
    beta2 <- smod$coefficients["I(Time^2)", "Estimate"]
    extrem <- (-beta1 / 2 / beta2)
    maxi <- max(as.numeric(model$model$Time))
    mini <- min(as.numeric(model$model$Time))
    if (extrem >= mini && extrem <= maxi) {
      trend <- ifelse(beta2 > 0, "Quadratic trend(u)", "Quadratic trend(n)")
    } else {
      trend <- ifelse(model$fitted.values[1] < model$fitted.values[length(model$fitted.values)],
                      "Quadratic trend(+)", "Quadratic trend(-)")
    }
  } else {
    # Check if p-values could be calculated
    if (!is.na(p_beta1)) {
      if (p_beta1 < alpha) {
        trend <- ifelse(beta1 < 0, "Falling trend", "Rising trend")
      } else {
        trend <- "No significant trend"
      }
    } else {
      # if a p-value cant be cant be calculated
      trend <- "no p-value could be calculated"
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
  beta1 <- smod$tTable["Time", "Value"]
  beta2 <- as.character(0)
  p_beta1 <- smod$tTable["Time", "p-value"]

  ## Quadratic trend if a quadratic term is present
  quadratic <- "I(Time^2)" %in% rownames(smod$tTable)

  if (quadratic) {
    # since the function fit_trend was already checked for significance we can assume
    # a quadratic trend
    beta2 <- smod$tTable["I(Time^2)", "Value"]
    extrem <- (-beta1 / 2 / beta2)
    maxi <- max(as.numeric(names(model$fitted)))
    mini <- min(as.numeric(names(model$fitted)))
    if (extrem >= mini && extrem <= maxi) {
      trend <- ifelse(beta2 > 0, "Quadratic trend(u)", "Quadratic trend(n)")
    } else {
      trend <- ifelse(model$fitted[1] < model$fitted[length(model$fitted)],
                      "Quadratic trend(+)", "Quadratic trend(-)")
    }
  } else {
    if (p_beta1 < alpha) {
      trend <- ifelse(beta1 < 0, "Falling trend", "Rising trend")
    } else {
      trend <- "No significant trend"
    }
  }
  return(list(beta0, beta1, beta2, trend))
}
