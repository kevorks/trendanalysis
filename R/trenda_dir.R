#' Perform trend analyses in time series related climate related data.
#'
#' @description trenda_dir() is a use-case specific addition to the function trenda().
#' It was created for bulk-processing time series related data stored as
#' comma delimited csv files in a folder.
#' Each file is checked for
#' It performs a series of analysis to find trends
#' in time-series related data by first checking the number of observations,
#' performing the Generalized Durbin-Watson-Test to test for autocorrelation errors,
#' performing AR(1) or ARMA to get the correlation strucutre before peforming
#' a GLS
#' @param data_dir Directory where data is stored. The function requires data to
#' be stored in a comma delimited csv file format and therefore only creates a
#' list of csv files from the specified folder.
#' @param log_trans logical. Set to TRUE if data is already log10-transformed.
#' If create_dir = TRUE a new folder named “_results_standard_log” will be created.
#' Observables will be transformed back to their original form and plots will be
#' adjusted for normal values.
#' @param calc_infl_obs logical. TRUE removes the influential observations
#' generated and stored in the file “res_tab_standard_”, Also creates a new
#' folder named “_results_standard_obs”
#' @param create_dir logical. Creates a subfolder in the data_dir specified path
#' if set to TRUE. Print to console if set to FALSE If log_trans and calc_infl_obs
#' is set to TRUE the following folders are created:
#' “_results_standard”,
#' “_results_standard_obs”,
#' “_resulst_log”,
#' “_results_log_obs”
#' @details Provided data_dir all semicolon delimited csv files are analysed
#' and summarized in a csv file created in a subfolder. The first column of each
#' file is expected to be a time series of type int followed by any number of
#' columns with numeric values.
#'
#' Based on the number of observations (n<7, 7 >= n < 13, n >= 13)
#' a linear equation or linear & quadratic equation is is formed.
#' A Durbin-Watson-Test with max.lag set to 2 is performed to find significant autocorrelation.
#' If p-value of first-order is smaller than  0.05 we reject H0 and assume significant autocorrelation
#' of first-order and perform an autoregressive model AR(1) to get the correlation structure for the GLS model.
#' If p >= 0.05 the second p-value of the DWT is inspected and checked again for significance.
#' If p >= 0.05 a linear regression model is fitted. For p < 0.05 we assume
#' autocorrelated  errors of second-order and perform an autoregressive-moving-average model
#' to get the correlation structure for GLS.
#' @examples # Assuming the folder with the data is placed in the working directory.
#' \dontrun{
#' list.files(path = "./files_folder")
#' }
#' # The analysis can be started by entering the path of the folder:
#' \dontrun{
#' trenda_dir("./files_folder/", log_trans = FALSE, create_dir = TRUE, calc_infl_obs = TRUE)
#' }
#' @return the function does not return a value but stores the plot and table in
#' the assigned directories
#' @importFrom utils str write.table read.table read.csv2
#' @importFrom grDevices jpeg dev.off
#' @export

trenda_dir <- function(data_dir, log_trans = FALSE, create_dir = TRUE, calc_infl_obs = TRUE) {
  if(create_dir) {
    if (!log_trans) {
     dir.create(paste0(data_dir, Sys.Date(), "_results_standard"))
     plot_dir <- paste0(data_dir, Sys.Date(), "_results_standard/")
     result_name <- paste0(plot_dir, "res_tab_standard_")
   } else {
     dir.create(paste0(data_dir, Sys.Date(), "_results_log"))
     plot_dir <- paste0(data_dir, Sys.Date(), "_results_log/")
     result_name <- paste0(plot_dir, "res_tab_log_")
  }
  }

#  plot_dir <- dir.create(paste(data_dir, ))
  trend_files <- list.files(data_dir, pattern = "*.csv")

  ResTab <- data.frame(matrix(ncol = 13, nrow = 0))
  colnames(ResTab) <- c("File", "Index", "Beta0", "Beta1", "Beta2", "Phi1", "Phi2",
                        "Trend", "rSquared", "observations_to_remove_index",
                        "observations_to_remove_value",
                        "cooks_distance_of_observation",
                        "observations_to_remove_year"
  )

  # since generalized Durbin-Watson-Test is random -> set seed
  set.seed(1)

  z <- 1

  output <- lapply(trend_files, function(trend_file) {
    # Exceltitle
    print(sprintf("Data:----------%s", trend_file))

    # read dataset
    dat <- read.csv2(sprintf("%s%s", data_dir, trend_file), header = TRUE, sep = ";", fileEncoding = "WINDOWS-1252")
    f_r <- names(dat[1])
    # show structure of data to be analysed
    print(str(dat))
    print(all(apply(dat, 2, is.numeric)))

    # dots have to be replaced in the variable names for the plot function
    names(dat) <- gsub("\\.", "\\_", names(dat))

    # all variables except for ID, Time and Jahr are environmental indicators
    varnames <- setdiff(names(dat), c("ID", "Time", f_r))

    # ID is needed for the correlation structure in the GLS-function
    dat$ID <- 1

    # Time begins at 0
    dat["Time"] <- dat[1] - min(dat[1])

    for(varname in varnames){
      print(varname)

      # calculate trend
      res <- compute_trend(dat, varname, log_trans = log_trans)

      ## save results
      ResTab[z,"File"] <<- trend_file
      ResTab[z,"Index"] <<- varname
      ResTab[z,"Beta0"] <<- round(as.numeric(res$beta[1]),4)
      ResTab[z,"Beta1"] <<- round(as.numeric(res$beta[2]),4)
      ResTab[z,"Beta2"] <<- round(as.numeric(res$beta[3]),4)
      ResTab[z,"Phi1"] <<- res$phi[1]
      ResTab[z,"Phi2"] <<- res$phi[2]
      ResTab[z,"Trend"] <<- res$trend
      ResTab[z,"rSquared"] <<- res$rsq
      ResTab[z, "observations_to_remove_index"] <<- res$infl_obs_index
      ResTab[z, "observations_to_remove_value"] <<- res$infl_obs_value
      ResTab[z, "cooks_distance_of_observation"] <<- res$infl_obs_cookd
      ResTab[z, "observations_to_remove_year"] <<- res$infl_obs_year

      # save trend graphs
      if (log_trans) {
        name_file <- substr(trend_file, start = 6, stop = nchar(trend_file) - 4)
      } else {
        name_file <- suppressWarnings(abbreviate(trend_file, 8))
      }
      if (create_dir) {
      jpeg(filename = sprintf("%s%s_%s.jpeg", plot_dir, name_file, varname),
           width = 800, height = 800, quality = 640000)
      plot(res$plot)
      dev.off()
      }
      z <<- z + 1

    }

  })
  #summary(ResTab)
  ResTab = data.frame(apply(ResTab, 2, gsub, patt="\\.", replace=","))
  if (create_dir) {
  write.table(ResTab,
              paste0(result_name, format(Sys.Date(), "%Y-%m-%d"), ".csv"),
              sep = ";", dec = ",", row.names = FALSE)
  } else {
    #ResTab <- ResTab
    return(ResTab)
  }
  # check if all files were used
  trend_files == unique(ResTab$File)
  #perform trenda_obs
  if (calc_infl_obs) {
  trenda_obs(data_dir, log_trans, calc_infl_obs = FALSE, create_dir = create_dir, res_tab_file = paste0(result_name, format(Sys.Date(), "%Y-%m-%d"), ".csv"))
  }
}
