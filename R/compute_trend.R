<<<<<<< HEAD
##' Calculates the trend for data and variables
##' @description Checks whether enough data is available and runs the model
##' @param dat Dataset
##' @param varname list with variablenames from the dataset
##' @param log_trans logical indicates if data is log-transformed
##' @param calc_infl_obs logical includes influential observations
##' @importFrom stats na.omit
##' @importFrom ggplot2 ggplot geom_point aes_string
##' @return A list with trendgraphs, variablenames, model and information about the results
compute_trend <- function(dat, varname, log_trans = FALSE, calc_infl_obs = TRUE){

  ## Code variable as numeric
  dat[, varname] <- as.numeric(dat[ ,varname])
  dat <- na.omit(dat[, c("Time", "Jahr", "ID", varname)])

  ## Check if enough datapoints are available
  enough_data <- check_n(dat[, varname, drop = TRUE])

  ## Results if not enough data is available
  if (!enough_data) {
    return(list(plot = ggplot(dat) + geom_point(aes_string(x = "Jahr", y = varname)),
                varname = varname,
                mod = NULL, trend = "No trend analysis possible",
                phi = c(NA, NA),
                beta = list(beta0 = NA, beta1 = NA, beta2 = NA),
                rsq = NA,
                infl_obs_index = NA,
                infl_obs_value = NA,
                infl_obs_cookd = NA,
                infl_obs_year = NA))
  }
  ## Fiting the model
  model_data <- fit_trend(dat, varname = varname)
  mod <- model_data[["mod"]]
  rsq <- round(1 - (sum(mod$residuals ^ 2) / sum(((mod$residuals + mod$fitted) - mean(mod$residuals + mod$fitted)) ^ 2)), 4)
  if (calc_infl_obs) {
    infl_obs <- infl_observations(mod = mod,
                                  used_formula =
                                    model_data[["used_formula"]],
                                  cor_struct =
                                    model_data[["correlation_structure"]],
                                  data = dat,
                                  varname = varname)
  } else {
    infl_obs <- list(infl_obs_index = NA, infl_obs_value = NA,
                     infl_obs_cookd = NA,
                     infl_obs_year = NA)
  }

  list(plot = plot_trend(mod, df = dat, log_trans = log_trans), varname = varname,
       mod = mod, trend = extract_trend(mod)[4], phi = extract_phi(mod),
       beta = extract_trend(mod)[1:3], rsq = max(0,rsq),
       infl_obs_index = paste0(infl_obs[[1]], collapse = ", "),
       infl_obs_value = paste(infl_obs[[2]], collapse = ", "),
       infl_obs_cookd = paste(infl_obs[[3]], collapse = ", "),
       infl_obs_year = paste(infl_obs[[4]], collapse = ", "))
}
=======
##' Calculates the trend for data and variables
##' @description Checks whether enough data is available and runs the model
##' @param dat Dataset
##' @param varname list with variablenames from the dataset
##' @param log_trans logical indicates if data is log-transformed
##' @param calc_infl_obs logical includes influential observations
##' @importFrom stats na.omit
##' @importFrom ggplot2 ggplot geom_point aes_string
##' @return A list with trendgraphs, variablenames, model and information about the results
compute_trend <- function(dat, varname, log_trans = FALSE, calc_infl_obs = TRUE){

  ## Code variable as numeric
  dat[, varname] <- as.numeric(dat[ ,varname])
  dat <- na.omit(dat[, c("Time", "Jahr", "ID", varname)])

  ## Check if enough datapoints are available
  enough_data <- check_n(dat[, varname, drop = TRUE])

  ## Results if not enough data is available
  if (!enough_data) {
    return(list(plot = ggplot(dat) + geom_point(aes_string(x = "Jahr", y = varname)),
                varname = varname,
                mod = NULL, trend = "No trend analysis possible",
                phi = c(NA, NA),
                beta = list(beta0 = NA, beta1 = NA, beta2 = NA),
                rsq = NA,
                infl_obs_index = NA,
                infl_obs_value = NA,
                infl_obs_cookd = NA,
                infl_obs_year = NA))
  }
  ## Fiting the model
  model_data <- fit_trend(dat, varname = varname)
  mod <- model_data[["mod"]]
  rsq <- round(1 - (sum(mod$residuals ^ 2) / sum(((mod$residuals + mod$fitted) - mean(mod$residuals + mod$fitted)) ^ 2)), 4)
  if (calc_infl_obs) {
    infl_obs <- infl_observations(mod = mod,
                                  used_formula =
                                    model_data[["used_formula"]],
                                  cor_struct =
                                    model_data[["correlation_structure"]],
                                  data = dat,
                                  varname = varname)
  } else {
    infl_obs <- list(infl_obs_index = NA, infl_obs_value = NA,
                     infl_obs_cookd = NA,
                     infl_obs_year = NA)
  }

  list(plot = plot_trend(mod, df = dat, log_trans = log_trans), varname = varname,
       mod = mod, trend = extract_trend(mod)[4], phi = extract_phi(mod),
       beta = extract_trend(mod)[1:3], rsq = max(0,rsq),
       infl_obs_index = paste0(infl_obs[[1]], collapse = ", "),
       infl_obs_value = paste(infl_obs[[2]], collapse = ", "),
       infl_obs_cookd = paste(infl_obs[[3]], collapse = ", "),
       infl_obs_year = paste(infl_obs[[4]], collapse = ", "))
}
>>>>>>> 3009302fe227bbffce700fdebf0b01c92f77b4a7
