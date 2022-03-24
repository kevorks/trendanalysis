#' Function to perform trend analyses in time series
#'
#' @description Given a folder with csv files
#' @param data_dir directory where the data is stored
#' @param log_trans logical if log-transformed values are used (TRUE leads to
#' the adjustments for plot function)
#' @return the function does not return a value but stores the plot and table in
#' the assigned directories
#' @importFrom utils str write.table read.table read.csv2
#' @importFrom grDevices jpeg dev.off
#' @export
#'

trenda <- function(data_dir, log_trans = FALSE, create_dir = TRUE) {
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
    dat <- read.csv2(sprintf("%s%s", data_dir, trend_file), header = TRUE, sep = ";")

    # show structure of data to be analysed
    print(str(dat))
    print(all(apply(dat, 2, is.numeric)))

    # dots have to be replaced in the variablenames for the plot function
    names(dat) <- gsub("\\.", "\\_", names(dat))

    # all variables except for ID, Time and Jahr are environmental indicators
    varnames <- setdiff(names(dat), c("ID", "Time", "Jahr"))

    # ID is needed for the correlationstructure in the gls-function
    dat$ID <- 1

    # Time begins at 0
    dat$Time <- dat$Jahr - min(dat$Jahr)

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
        name_file <- abbreviate(trend_file, 8)
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

  summary(ResTab)
  if (create_dir) {
  write.table(ResTab,
              paste0(result_name, format(Sys.Date(), "%Y-%m-%d"), ".csv"),
              sep = ";", dec = ",", row.names = FALSE)
  } else {
    return(ResTab)
  }
  # check if all files were used
  trend_files == unique(ResTab$File)
}
