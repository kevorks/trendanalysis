#' Define influential observations according to Cook's distance
#'
#' @description All observations whose Cook's distance is greater than the threshold
#' are specified#'
#' @param mod method lm or gls model
#' @param threshold Threshold at which observations are considered influential
#' @param ... placeholder
#' @return List with indices and observations
infl_observations <- function(mod, threshold = 0.5, ...) {
  UseMethod("infl_observations")
}

##' Used in method infl_observations as lm model
##' @param mod method used is lm
##' @param threshold level of significance
##' @param data data which are used for to fit the model
##' @param ... placeholder
##' @importFrom stats cooks.distance
infl_observations.lm <- function(mod, threshold = 0.5, data, ...) {
  # Cook's distance berechnen
  cooks_distance <- cooks.distance(mod)
  # Index der auffälligen Beobachtungen bestimmen
  index <- which(cooks_distance > threshold)
  attributes(index) <- NULL
  cooks_distance_selected <- cooks_distance[index]
  # die Werte der auffälligen Beobachtungen bestimmen
  values <- mod$model[index, 1]
  # das Jahr der auffälligen Beobachtung bestimmen
  jahr <- data[index, "Jahr"]

  list(infl_obs_index = index, infl_obs_value = values,
       infl_obs_cookd = cooks_distance_selected,
       infl_obs_jahr = jahr)
}

##' Used in method infl_observations as gls model
##' @param mod method used is gls
##' @param used_formula placeholder
##' @param cor_structure placeholder
##' @param threshold placeholder
##' @param data placeholder
##' @param varname placeholder
infl_observations.gls <- function(mod, used_formula, cor_structure,
                                  threshold = 0.5, data, varname) {
  # Cook's distance calculated
  cooks_distance <- CookD_gls(mod, used_formula = used_formula,
                              cor_structure = cor_structure, data = data,
                              plot = FALSE)
  # Index of noticable observations
  index <- which(cooks_distance > threshold)
  attributes(index) <- NULL
  cooks_distance_selected <- cooks_distance[index]
  # get data from used model
  model_matrix <- data
  # determin values of noticable observations
  values <- model_matrix[index, varname]
  # determin Jahr of noticable observations
  jahr <- data[index, "Jahr"]

  list(infl_obs_index = index, infl_obs_value = values,
       infl_obs_cookd = cooks_distance_selected,
       infl_obs_jahr = jahr)
}
