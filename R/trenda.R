#' Function to perform trend analyses in time series and identify outliers using
#' Cook's Distance
#' @description performs a series of analysis to find trends
#' in time-series related data by first checking the number of observations,
#' performing the Generalized Durbin-Watson-Test to test for autocorrelation,
#' performing AR(1) or ARMA to get the correlation structure before performing
#' a GLS and further deciding on the trends of an observation by looking at significance
#' and beta values
#' @param df file where data is stored. The function requires data to
#' be stored in a semicolon delimited csv file format and therefore only creates a
#' list of csv files from the specified folder.
#' @param log_trans logical. Set to TRUE if data is already log10-transformed.
#' If create_dir = TRUE a new folder named “_results_standard_log” will be created.
#' Observables will be transformed back to their original form and plots will be
#' adjusted for normal values.
#' @param plot_graphs logical. Set to TRUE by default. Will plot to the environment.
#' The plots summarizes the information about the data the trend
#' @param calc_infl_obs logical: If set to TRUE the function will calculate the
#' influential observations using Cook's Distance. The data is stored in a list
#' with information
#' @details The first column of the data frame is
#' expected to be a time series of type int followed by any number of
#' columns with numeric values.
#' Based on the number of observations (n<7, 7 >= n < 13, n >= 13)
#' a linear equation or linear & quadratic equation is is formed.
#' A Durbin-Watson-Test with max.lag set to 2 is performed to find significant autocorrelation.
#' If p-value of first-order is smaller than  0.05 we reject H0 and assume significant autocorrelation
#' of first-order and perform an autoregressive model AR(1) to get the correlation structure for the GLS model.
#' If p >= 0.05 the second p-value of the DWT is inspected and checked again for significance.
#' If p >= 0.05 a linear regression model is fitted. For p < 0.05 we assume
#' autocorrelated  errors of second-order and perform an autoregressive-moving-average model
#' to get the correlation structure for GLS.
#' @return If calc_infl_obs is set to true a list is returned containing two data.frames.
#'  The first data frame informs about the findings of the first run and suggested
#'  which values are influential using Cook's Distance. In a second run the influential
#'  observations are removed from the data frame and only those variables are reevaluated
#'  and returned in the second data frame.
#'
#' @examples
#' ## Generate a data frame with to variables. One value of height_m is set
#' ## to be extremely out of bound to illustrate how trenda() handles outliers.
#'\dontrun{ rnd_preci <- data.frame(Year = c(1991:2020), precip_mm = rnorm(30, 770, 50),
#'                                height_m = c(rnorm(10, 10, 1), 100,
#'                                rnorm(9, 15, 1),
#'                                rnorm(10, 20, 1)))
#'
#' trenda(rnd_preci, plot_graphs = TRUE, log_trans = FALSE, calc_infl_obs = TRUE)
#' }
#' @importFrom utils str write.table read.table read.csv2
#' @importFrom grDevices jpeg dev.off
#' @export



trenda <- function(df, plot_graphs = FALSE, log_trans = FALSE, calc_infl_obs = FALSE) {
  trend_files <- deparse(substitute(df))
  # Create empty data frame
  ResTab <- data.frame(matrix(ncol = 13, nrow = 0))
  colnames(ResTab) <- c("File", "Index", "Beta0", "Beta1", "Beta2", "Phi1", "Phi2",
                        "Trend", "rSquared", "observations_to_remove_index",
                        "observations_to_remove_value",
                        "cooks_distance_of_observation",
                        "observations_to_remove_year"
  )
  # Create empty data frame for calc_infl_obs = TRUE
  ResTab2 <- data.frame(matrix(ncol = 13, nrow = 0))
  colnames(ResTab2) <- c("File", "Index", "Beta0", "Beta1", "Beta2", "Phi1", "Phi2",
                        "Trend", "rSquared", "observations_to_remove_index",
                        "observations_to_remove_value",
                        "cooks_distance_of_observation",
                        "observations_to_remove_year"
  )
  list <- list()

  #set.seed(1)

  f_r <- names(df[1])
  z <- 1
  output <- lapply(trend_files, function(trend_file) {

  # some data cause wrong calculations if they are separated by ","
    dat <- data.frame(apply(apply(df, 2, gsub, patt=",", replace="."), 2, as.numeric))
    print(str(dat))
    print(all(apply(dat, 2, is.numeric)))
    # Dots are replaced for the plot_trend() function
    names(dat) <- gsub("\\.", "\\_", names(dat))
    varnames <- setdiff(names(dat), c("ID", "Time", f_r))
    # ID needed for GLS
    dat["ID"] <- 1
    # time beginns at zero
    dat["Time"] <- dat[1] - min(dat[1])

    for(varname in varnames) {

    print(varname)
    #calculate trend
    res <- compute_trend(dat, varname, log_trans = log_trans, calc_infl_obs = calc_infl_obs)

    ResTab[z, "File"] <<- trend_file
    ResTab[z, "Index"] <<- varname
    ResTab[z, "Beta0"] <<- round(as.numeric(res$beta[1]),4)
    ResTab[z, "Beta1"] <<- round(as.numeric(res$beta[2]),4)
    ResTab[z, "Beta2"] <<- round(as.numeric(res$beta[3]),4)
    ResTab[z, "Phi1"] <<- res$phi[1]
    ResTab[z, "Phi2"] <<- res$phi[2]
    ResTab[z, "Trend"] <<- res$trend
    ResTab[z, "rSquared"] <<- res$rsq
    ResTab[z, "observations_to_remove_index"] <<- res$infl_obs_index
    ResTab[z, "observations_to_remove_value"] <<- res$infl_obs_value
    ResTab[z, "cooks_distance_of_observation"] <<- res$infl_obs_cookd
    ResTab[z, "observations_to_remove_year"] <<- res$infl_obs_year
    if(plot_graphs) {
    plot(res$plot)

  }
  z <<- z + 1

}


})
result_infl <- ResTab
m <- 1
if(calc_infl_obs) {
  output <- lapply(trend_files, function(trend_file) {

    dat <- data.frame(apply(apply(df, 2, gsub, patt=",", replace="."), 2, as.numeric))
        print(str(dat))
    f_r <- names(dat[1])

    names(dat) <- gsub("\\.", "\\_", names(dat))

    varnames <- setdiff(names(dat), c("ID", "Time", f_r))

    dat$ID <- 1

    dat["Time"] <- dat[1] - min(dat[1])

    for(varname in varnames){
      listr2 <- list()
      print(varname)
      index <- result_infl[result_infl$File == trend_file &
                             result_infl$Index == varname,
                           "observations_to_remove_index"]
      if (!is.na(index)) {
        index <- as.character(index)
        index <- strsplit(index, ",")
        index <- as.vector(sapply(index, as.numeric))
      }

      if (!is.list(index) && !is.na(index)) {
        temp_dat <- na.omit(dat[, c(f_r, "Time", "ID", varname)])
        temp_dat <- temp_dat[-index, ]

        print(str(temp_dat))

        res2 <- compute_trend(temp_dat, varname, log_trans = log_trans, calc_infl_obs = FALSE)

        ResTab2[m,"File"] <<- trend_file
        ResTab2[m,"Index"] <<- varname
        ResTab2[m,"Beta0"] <<- round(as.numeric(res2$beta[1]),4)
        ResTab2[m,"Beta1"] <<- round(as.numeric(res2$beta[2]),4)
        ResTab2[m,"Beta2"] <<- round(as.numeric(res2$beta[3]),4)
        ResTab2[m,"Phi1"] <<- res2$phi[1]
        ResTab2[m,"Phi2"] <<- res2$phi[2]
        ResTab2[m,"Trend"] <<- res2$trend
        ResTab2[m,"rSquared"] <<- res2$rsq
        ResTab2[m, "observations_to_remove_index"] <<- res2$infl_obs_index
        ResTab2[m, "observations_to_remove_value"] <<- res2$infl_obs_value
        ResTab2[m, "cooks_distance_of_observation"] <<- res2$infl_obs_cookd
        ResTab2[m, "observations_to_remove_year"] <<- res2$infl_obs_year

        if (log_trans) {
          name_file <- substr(trend_file, start = 6, stop = nchar(trend_file) - 4)
        } else {
          name_file <- suppressWarnings(abbreviate(trend_file, 8))
        }
        if (plot_graphs) {
          plot(res2$plot)
        }

        m <<- m + 1
      }
    }
  })
}

if (calc_infl_obs){
  list(ResTab,
       ResTab2)
} else {
  ResTab
}
}
