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
  # calculate Cook's distance
  cooks_distance <- cooks.distance(mod)
  # determin index of noticeable observations
  index <- which(cooks_distance > threshold)
  attributes(index) <- NULL
  cooks_distance_selected <- cooks_distance[index]
  # determin the value of noticeable observations
  values <- mod$model[index, 1]
  # determin the year of noticeable observations
  year <- data[index, "Year"]

  list(infl_obs_index = index, infl_obs_value = values,
       infl_obs_cookd = cooks_distance_selected,
       infl_obs_year = year)
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
  # calculate Cook's distance
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
  # determin year of noticable observations
  year <- data[index, "Year"]

  list(infl_obs_index = index, infl_obs_value = values,
       infl_obs_cookd = cooks_distance_selected,
       infl_obs_year = year)
}
